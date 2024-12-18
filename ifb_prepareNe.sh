#! /bin/bash
#SBATCH --job-name=PREPARENE
#SBATCH -o /shared/projects/plant_lewontin_paradox/vcf_processing/genouest_log/prepare_Ne/%x.%j.out
#SBATCH -e /shared/projects/plant_lewontin_paradox/vcf_processing/genouest_log/prepare_Ne/%x.%j.err
#SBATCH -p long
#SBATCH --time=2-00
#SBATCH --mem=10G

. /shared/projects/plant_lewontin_paradox/paths.sh

echo "Command :
$snp_scripts/go_snparcher.sh $@" >&2

unlock=false
rerun=false

while getopts "i:o:ur" opt; do
  case $opt in
    i)
      config="$OPTARG" # path/to/config.yaml
      ;;
    o)
      dir="$OPTARG" # path/to/analysis
      ;;
    u)
      unlock=true
      ;;
    r)
      rerun=true
      ;;
    \?)
      echo "Usage : ne/scripts/prepare_smcpp.sh [-i config_path] [-o dir] [-u]"
      exit 1
      ;;
  esac
done

module load snakemake/8.9.0

if [ $unlock == true ]; then
  snakemake -s $vcf_scripts/snakemake_prepareNe/workflow/Snakefile \
    -d $dir \
    --unlock
  exit 1
fi

mkdir -p $dir/config
cp $config $dir/config/config.yml

if [ $rerun == true ]; then
  snakemake -s $vcf_scripts/snakemake_prepareNe/workflow/Snakefile \
  -d $dir \
  --workflow-profile $vcf_scripts/snakemake_prepareNe/profiles/slurm \
  --use-conda \
  --conda-frontend conda \
  --rerun-incomplete
else
  snakemake -s $vcf_scripts/snakemake_prepareNe/workflow/Snakefile \
    -d $dir \
    --workflow-profile $vcf_scripts/snakemake_prepareNe/profiles/slurm \
    --use-conda \
    --conda-frontend conda
fi