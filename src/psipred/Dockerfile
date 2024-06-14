# Use the official Ubuntu base image
FROM amd64/ubuntu:20.04

# Set environment variables to prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
    build-essential \
    wget \
    tcsh \
    csh \
    gcc \
    perl \
    unzip \
    gzip \
    vim \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Download and install classic BLAST
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-x64-linux.tar.gz && \
    tar -xzvf blast-2.2.26-x64-linux.tar.gz && \
    rm blast-2.2.26-x64-linux.tar.gz && \
    cp blast-2.2.26/bin/* /usr/local/bin && \
    rm -rf blast-2.2.26

# Download and install PSIPRED
RUN wget http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/psipred.4.02.tar.gz && \
    tar -xzvf psipred.4.02.tar.gz && \
    rm psipred.4.02.tar.gz

# Set the working directory to PSIPRED source, compile PSIPRED
RUN cd /psipred/src && make && make install
ENV PATH="/psipred:/psipred/bin:${PATH}"

# Download UniRef90 data bank and set up BLAST database
RUN wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz && \
    gunzip -v uniref90.fasta.gz && \
    rm -r uniref90.fasta.gz

RUN formatdb -p T -o T -t uniref90 -i uniref90.fasta 


# Create a mount point
VOLUME [ "/shared" ]

# Set the default command to run a shell
CMD ["/bin/bash"]
