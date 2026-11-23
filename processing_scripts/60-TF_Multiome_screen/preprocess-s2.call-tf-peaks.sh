######################################
# STEP 2:
######################################
# Call peaks on gRNA cells for each
# individual TF
######################################

# Path to the project directory (must contain data/targets.numbered.tsv
# and the s1 output data/atac-rna.preprocessed-s1.signac.rds)
work_dir="<path to project directory>"
scripts_dir="<path to processing_scripts/60-TF_Multiome_screen>"

cd "${work_dir}"

# Activate the conda environment with R 4.3.2, Seurat, Signac and MACS2
source "<path to conda>/etc/profile.d/conda.sh"
conda activate "<environment name>"

LOOKUP_FILE="${work_dir}/data/targets.numbered.tsv"

# Check if the lookup file exists and is readable
if [[ ! -r "${LOOKUP_FILE}" ]]; then
    echo "Cannot find ${LOOKUP_FILE}" >&2
    exit 1
fi

# Call peaks for each TF (column 1 = index, column 2 = TF name).
# Each TF is independent, so these can also be run in parallel,
# e.g. as a job array over the 60 rows (about 25 GB memory each).
while read -r idx TF_NAME; do
    echo "Calling peaks for ${TF_NAME}"
    mkdir -p "${work_dir}/data/peaks/${TF_NAME}"
    Rscript "${scripts_dir}/preprocess-s2.call-tf-peaks.R" --tf "${TF_NAME}"
done < "${LOOKUP_FILE}"

echo "Complete."
