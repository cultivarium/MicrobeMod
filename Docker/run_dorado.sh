#!/bin/bash

OPTIND=1
use_aws_credentials=0
use_gpu=0

while getopts "h?agl:" OPT; do
  case "$OPT" in
  h|\?)
    echo "$ ./run_dorado.sh -g -l <library_name>"
    echo "Options:"
    echo "-a: pass AWS credentials specified as environment variables to the aws command"
    echo "-g: use GPU"
    echo "-l: library name"
    echo "Environment variables:"
    echo "BUCKET_LOCATION: S3 path with the input pod5 files (required)"
    echo "DORADO_MODEL: dorado model complex (default: sup)"
    echo "DORADO_ARGS: dorado basecaller flags after the reads input (default: --recursive --modified-bases 4mC_5mC 6mA)"
    exit 0
    ;;
  a)
    use_aws_credentials=1
    ;;
  g)
    use_gpu=1
    ;;
  l)
    library_name=$OPTARG
    ;;
  esac
done

: "${DORADO_MODEL:=sup}"
: "${DORADO_ARGS:=--recursive --modified-bases 4mC_5mC 6mA}"

if [[ -z $BUCKET_LOCATION ]]; then
  echo "Please provide a bucket location"
  exit 1
fi

if [[ -z $library_name ]]; then
  echo "Please provide a library name"
  exit 1
fi

if [[ $use_aws_credentials == 1 ]]; then
  export AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID
  export AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY
  export AWS_SESSION_TOKEN=$AWS_SESSION_TOKEN
fi

mkdir tmp$library_name
cd tmp$library_name

aws s3 cp --recursive $BUCKET_LOCATION $library_name
 
if [[ $use_gpu == 1 ]]; then
  echo "Using GPU"
  /home/ubuntu/dorado-2.0.0-linux-x64/bin/dorado basecaller $DORADO_MODEL $library_name $DORADO_ARGS > $library_name.d3.bam
else
  echo "Using CPU"
  /home/ubuntu/dorado-2.0.0-linux-x64/bin/dorado basecaller -x cpu $DORADO_MODEL $library_name $DORADO_ARGS > $library_name.d3.bam
fi

aws s3 cp $library_name.d3.bam $BUCKET_LOCATION

cd ..
rm -rf tmp$library_name
