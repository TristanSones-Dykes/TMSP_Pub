#!/bin/bash
################################################################################
##
## This script demonstrates the following workflow:
##      
##   1. Parse CLAs and read Gene IDs from a file
##   2. Search for genes based on Gene IDs
##   3. Perform a GO Enrichment Analysis on the search result
##   4. Download the analysis result as a tabular file
##
################################################################################

#---------------------------------------------------------------
# parse command line arguments
#---------------------------------------------------------------

while getopts t:i:o: flag
do
    case "${flag}" in
        t) taxonomy=${OPTARG};;
        i) inputFile=${OPTARG};;
        o) outputDir=${OPTARG};;
    esac
done

# validate input file
if [ ! -r "$inputFile" ]; then
    echo "Input file not found: $inputFile"
    exit 1
fi

# create output directory if it doesn't exist
mkdir -p "$outputDir"

#---------------------------------------------------------------
# script config
#---------------------------------------------------------------

webappUrl="https://fungidb.org/fungidb"
cookieJar="$outputDir/cookie-jar"
outputFile="$outputDir/goEnrichmentResult.tsv"

#---------------------------------------------------------------
# search config
#---------------------------------------------------------------

geneIds=$(jq -R -s -c 'split("\n")[:-1]' < "$inputFile")

#---------------------------------------------------------------
# analysis config
#---------------------------------------------------------------

ontologies="Cellular Component"
evidenceCodes="[\\\"Computed\\\",\\\"Curated\\\"]" # JSON array (double escapes this time)
goSubset="Yes"
pValueCutoff="0.05"


################################################################################
##
##  Succession of CURL commands to download GO Enrichment result
##
################################################################################

#---------------------------------------------------------------
# establish a guest user and cookies needed to act as that user
#---------------------------------------------------------------

echo "Establishing guest user (creates cookie jar)"
curl -s -S -c $cookieJar \
     $webappUrl/service/users/current > /dev/null

#---------------------------------------------------------------
# create a dataset param ID containing specified gene IDs
#---------------------------------------------------------------

echo "Creating input dataset"
datasetData="{\"sourceType\":\"idList\",\"sourceContent\":{\"ids\":$geneIds}}"
datasetId=`curl -s -S -b $cookieJar \
     -X POST \
     -H "Content-Type:application/json" \
     -d "$datasetData" \
     $webappUrl/service/users/current/datasets | \
     jq .id`

#---------------------------------------------------------------
# create a step for this search's results
#---------------------------------------------------------------

echo "Creating search step"
stepData="{\"searchName\":\"GeneByLocusTag\",\"searchConfig\":{\"parameters\":{\"ds_gene_ids\":\"${datasetId}\"},\"wdkWeight\":10},\"customName\":\"IDs List\"}"
stepId=`curl -s -S -b $cookieJar \
     -X POST \
     -H "Content-Type:application/json" \
     -d "$stepData" \
     $webappUrl/service/users/current/steps | \
     jq .id`

#---------------------------------------------------------------
# create an analysis on this step
#---------------------------------------------------------------

echo "Creating analysis"
analysisData="{\"analysisName\":\"go-enrichment\",\"displayName\":\"Gene Ontology Enrichment\",\"parameters\":{\"organism\":\"${taxonomy}\",\"goAssociationsOntologies\":\"${ontologies}\",\"goEvidenceCodes\":\"${evidenceCodes}\",\"goSubset\":\"${goSubset}\",\"pValueCutoff\":\"${pValueCutoff}\"}}"
analysisId=`curl -s -S -b $cookieJar \
     -X POST \
     -H "Content-Type:application/json" \
     -d "$analysisData" \
     $webappUrl/service/users/current/steps/$stepId/analyses | \
     jq .analysisId`


#---------------------------------------------------------------
# start the analysis job
#---------------------------------------------------------------

echo "\nStarting the analysis job..."
start_analysis_url="$webappUrl/service/users/current/steps/$stepId/analyses/$analysisId/result"
start_analysis_response=$(curl -s -S -b $cookieJar \
     -X POST \
     -H "Content-Type: application/json" \
     $start_analysis_url)

# Check the HTTP status code
http_status=$(echo "$start_analysis_response" | jq -r '.status')

if [ "$http_status" == "202" ]; then
    echo "Analysis job started successfully"
elif [ "$http_status" == "400" ] || [ "$http_status" == "406" ]; then
    echo "Invalid request or unacceptable analysis status"
    exit 1
else
    echo "Analysis job started with status: $http_status"
fi

#---------------------------------------------------------------
# wait for analysis to complete
#---------------------------------------------------------------

while :; do
  echo "Checking analysis job status"
  jobStatus=`curl -s -S -b $cookieJar \
       $webappUrl/service/users/current/steps/$stepId/analyses/$analysisId/result/status | \
       jq .status`
  if [[ "$jobStatus" == '"RUNNING"' ]]; then
    echo "Job still running; checking again in 2 seconds"
    sleep 2
    continue
  fi
  if [[ "$jobStatus" == '"COMPLETE"' ]]; then
    echo "Job complete. Ready to download results."
    break
  fi
  echo "Job failed: $jobStatus"
  exit 1
done

#---------------------------------------------------------------
# download result
#---------------------------------------------------------------

echo "Results downloaded to $outputFile"
curl -s -S -b $cookieJar \
     -o $outputFile \
     $webappUrl/service/users/current/steps/$stepId/analyses/$analysisId/resources?path=hiddenGoEnrichmentResult.tsv

#---------------------------------------------------------------
# clean up cookie jar
#---------------------------------------------------------------

rm $cookieJar

echo "Done"

