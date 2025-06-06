# Python Wrapper for phobius.py with Asynchronous Concurrency Control
#
# This script processes a multi-sequence FASTA file by calling the
# phobius.py script for each sequence. It uses Python's asyncio library
# to manage a constant pool of concurrent jobs, ensuring that no more than
# a specified number of jobs are running at any one time. This is the
# idiomatic and most efficient way to handle many I/O-bound tasks in Python.

import asyncio
import os
import sys
import argparse
from pathlib import Path
import pandas as pd
import re
from typing import List, Tuple, Optional

# --- Configuration ---
# You can adjust these default values if needed.
DEFAULT_MAX_CONCURRENT_JOBS = 30
DEFAULT_PYTHON_EXECUTABLE = "python"


def parse_fasta(file_path: Path) -> List[Tuple[str, str]]:
    """
    Parses a FASTA file and returns a list of (header, sequence) tuples.
    This simple parser handles multi-line sequences.
    """
    sequences = []
    current_header = None
    current_sequence = []
    with file_path.open("r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header:
                    sequences.append((current_header, "".join(current_sequence)))
                current_header = line[1:].split()[
                    0
                ]  # Use the first part of the header as ID
                current_sequence = []
            else:
                current_sequence.append(line)
    if current_header:
        sequences.append((current_header, "".join(current_sequence)))
    return sequences


async def parse_phobius_output(
    output_file: Path, sequence_id: str
) -> Optional[pd.DataFrame]:
    """
    Parses the output file from phobius.py and returns a pandas DataFrame.
    Returns None if the file doesn't exist or no features are found.
    """
    if not output_file.exists():
        print(
            f"Warning: Output file not found for sequence {sequence_id}",
            file=sys.stderr,
        )
        return None

    try:
        content = output_file.read_text()
    except Exception as e:
        print(
            f"Warning: Could not read output file for {sequence_id}: {e}",
            file=sys.stderr,
        )
        return None

    # Regex to find all feature lines (FT) in the output
    # FT <TYPE> <START> <END> <DESCRIPTION>
    feature_pattern = re.compile(r"^FT\s+(\S+)\s+(\d+)\s+(\d+)\s+(.*)$", re.MULTILINE)
    matches = feature_pattern.findall(content)

    if not matches:
        return None

    # Create a DataFrame from the found features
    df = pd.DataFrame(matches, columns=["feature_type", "start", "end", "description"])
    df["sequence_id"] = sequence_id
    df["start"] = pd.to_numeric(df["start"])
    df["end"] = pd.to_numeric(df["end"])
    df["description"] = df["description"].str.strip()

    # Reorder columns for a clean output
    return df[["sequence_id", "feature_type", "description", "start", "end"]]


async def run_single_phobius_job(
    semaphore: asyncio.Semaphore,
    sequence_id: str,
    sequence: str,
    temp_dir: Path,
    python_script_path: Path,
    user_email: str,
    python_executable: str,
    pbar: str = "tqdm",
) -> Optional[pd.DataFrame]:
    """
    Runs a single phobius.py job for one sequence with a retry mechanism.

    This coroutine acquires a slot from the semaphore, runs the external
    process, retries on failure with exponential backoff, parses the result
    on success, and finally releases the semaphore slot.
    """
    MAX_RETRIES = 3  # Maximum number of retries for each job
    async with semaphore:
        safe_seq_id = re.sub(r"[^a-zA-Z0-9_-]", "_", sequence_id)
        job_name = f"phobius_{safe_seq_id}_{os.getpid()}_{id(sequence)}"
        temp_fasta = temp_dir / f"{job_name}.fasta"
        temp_outfile_base = temp_dir / job_name

        with temp_fasta.open("w") as f:
            f.write(f">{sequence_id}\n{sequence}\n")

        args = [
            python_executable,
            str(python_script_path),
            "--email",
            user_email,
            "--outfile",
            str(temp_outfile_base),
            str(temp_fasta),
        ]

        process = None
        for attempt in range(MAX_RETRIES + 1):
            if attempt > 0:
                # Exponential backoff: 2, 4, 8 seconds
                wait_time = 2**attempt
                print(
                    f"Retrying Phobius for {sequence_id} in {wait_time}s (attempt {attempt}/{MAX_RETRIES})",
                    file=sys.stderr,
                )
                await asyncio.sleep(wait_time)

            process = await asyncio.create_subprocess_exec(
                *args, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE
            )
            await (
                process.wait()
            )  # Use wait() instead of communicate() to check returncode first

            if process.returncode == 0:
                if attempt > 0:
                    print(f"Phobius for {sequence_id} succeeded on attempt {attempt}.")
                break  # Success, exit the retry loop
            else:
                stderr_output = await process.stderr.read()  # type: ignore
                print(
                    f"Error running Phobius for {sequence_id} (attempt {attempt}):\n{stderr_output.decode()}",
                    file=sys.stderr,
                )

        result_df = None
        # After the loop, check the final status of the process
        if process and process.returncode == 0:
            expected_result_file = temp_outfile_base.with_suffix(".out.txt")
            result_df = await parse_phobius_output(expected_result_file, sequence_id)
        else:
            print(
                f"Phobius job for sequence {sequence_id} failed after {MAX_RETRIES} retries.",
                file=sys.stderr,
            )

        # Clean up temporary files
        if temp_fasta.exists():
            temp_fasta.unlink()
        if "expected_result_file" in locals() and expected_result_file.exists():
            expected_result_file.unlink()

        pbar.update(1)  # type: ignore
        return result_df


async def main(args: argparse.Namespace):
    """
    The main asynchronous function to orchestrate the entire process.
    """
    # Validate input paths
    input_fasta = Path(args.input_fasta_file)
    python_script = Path(args.python_script_path)
    output_csv = Path(args.output_csv_file)

    if not input_fasta.is_file():
        print(f"Error: Input FASTA file not found at {input_fasta}", file=sys.stderr)
        sys.exit(1)
    if not python_script.is_file():
        print(
            f"Error: Phobius python script not found at {python_script}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Create a temporary directory for intermediate files
    temp_dir = Path("./phobius_temp")
    temp_dir.mkdir(exist_ok=True)
    print(f"Using temporary directory: {temp_dir.resolve()}")

    # Parse all sequences from the input file
    sequences = parse_fasta(input_fasta)
    if not sequences:
        print("No sequences found in the input file. Exiting.", file=sys.stderr)
        sys.exit(0)

    print(f"Found {len(sequences)} sequences to process.")

    # A semaphore will limit the number of concurrent jobs
    semaphore = asyncio.Semaphore(args.max_concurrent_jobs)

    # Try to import tqdm for a progress bar, but don't fail if it's not there
    try:
        from tqdm import tqdm

        pbar = tqdm(total=len(sequences), desc="Processing sequences")
    except ImportError:
        # Create a dummy object if tqdm is not installed
        class DummyPbar:
            def update(self, n):
                pass

            def __enter__(self):
                return self

            def __exit__(self, *args):
                pass

        pbar = DummyPbar()

    # Create a list of asynchronous tasks, one for each sequence
    with pbar:
        tasks = [
            run_single_phobius_job(
                semaphore,
                seq_id,
                seq,
                temp_dir,
                python_script,
                args.user_email,
                args.python_executable,
                pbar,  # type: ignore
            )
            for seq_id, seq in sequences
        ]

        # asyncio.gather runs all tasks concurrently and collects their results
        results = await asyncio.gather(*tasks)

    # Filter out any None results from failed jobs
    successful_results = [df for df in results if df is not None and not df.empty]

    if not successful_results:
        print("\nProcessing complete, but no features were found or parsed.")
    else:
        # Concatenate all result DataFrames into a single one
        final_df = pd.concat(successful_results, ignore_index=True)
        final_df.to_csv(output_csv, index=False)
        print(f"\nResults successfully saved to {output_csv}")

    # Clean up the temporary directory
    try:
        for f in temp_dir.iterdir():
            f.unlink()
        temp_dir.rmdir()
        print("Temporary directory cleaned up.")
    except Exception as e:
        print(
            f"Warning: Could not fully clean up temporary directory {temp_dir}: {e}",
            file=sys.stderr,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run Phobius in parallel for a multi-sequence FASTA file using asyncio.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "input_fasta_file", help="Path to the input multi-sequence FASTA file."
    )
    parser.add_argument(
        "output_csv_file", help="Path to write the final, collated CSV results."
    )
    parser.add_argument(
        "user_email", help="Your email address, required by the EBI API."
    )
    parser.add_argument("python_script_path", help="Path to the phobius.py script.")
    parser.add_argument(
        "--max_concurrent_jobs",
        type=int,
        default=DEFAULT_MAX_CONCURRENT_JOBS,
        help="The maximum number of phobius.py jobs to run in parallel.",
    )
    parser.add_argument(
        "--python_executable",
        default=DEFAULT_PYTHON_EXECUTABLE,
        help="The command to run Python (e.g., 'python' or 'python3').",
    )

    # Note: For this to run, Python 3.7+ is required for asyncio.
    # The 'pandas' and 'tqdm' libraries are also recommended.
    # You can install them with:
    # pip install pandas tqdm

    args = parser.parse_args()
    asyncio.run(main(args))
