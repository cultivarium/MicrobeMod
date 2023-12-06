OPTIND=1
bam_file=""
file_dir="./"
output_prefix="microbemod_output"
reference_genome=""
subcommand="call_methylation"
threads=1

while getopts "h?b:d:o:r:s:t:" OPT; do
  case "$OPT" in
  h|\?)
    echo "$ ./microbemod.sh -s <step>"
    echo "Options:"
    echo "-b: BAM file"
    echo "-d: dir for input and output files, default ./"
    echo "-o: output file prefix, default microbemod_output"
    echo "-r: reference genome"
    echo "-s: subcommand, either call_methylation or annotate_rm"
    echo "-t: threads, default 1"
    exit 0
    ;;
  b)
    bam_file=$OPTARG
    ;;
  d)
    file_dir=$OPTARG
    ;;
  o)
    output_prefix=$OPTARG
    ;;
  r)
    reference_genome=$OPTARG
    ;;
  s)
    subcommand=$OPTARG
    ;;
  t)
    threads=$OPTARG
    ;;
  esac
done

if [[ $subcommand != "call_methylation" ]] && [[ $subcommand != "annotate_rm" ]]; then
  echo "Unknown subcommand: $subcommand"
  exit 1
fi

if [[ -z $reference_genome ]]; then
  echo "Please provide a reference genome with the -r option"
  exit 1
fi

if [[ $subcommand == "call_methylation" ]] && [[ -z $bam_file ]]; then
  echo "Please provide a BAM file with the -b option"
  exit 1
fi

if [[ $subcommand == "call_methylation" ]]; then
  python ./bin/MicrobeMod call_methylation -b /files/$bam_file -r /files/$reference_genome -o $output_prefix -t $threads
fi

if [[ $subcommand == "annotate_rm" ]]; then
  python ./bin/MicrobeMod annotate_rm -f /files/$reference_genome -o $output_prefix -t $threads
fi

cp $output_prefix* $file_dir
