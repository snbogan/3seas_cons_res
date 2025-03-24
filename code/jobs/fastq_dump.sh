#!/bin/bash
#SBATCH --account=pi-jkoc
#SBATCH --partition=lab-colibri
#SBATCH --qos=pi-jkoc
#SBATCH --job-name=download_sra
#SBATCH --output=out/fastq-dump/fastq-dump%A_%a.out
#SBATCH --error=error/fastq-dump/fastq-dump%A_%a.err
#SBATCH --array=1-240
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=snbogan@ucsc.edu
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=ALL
#SBATCH --mem=10G
#SBATCH --cpus-per-task=2

# Load required module
module load sratoolkit

# Set working directory
cd /hb/home/snbogan/3seas/data

# Define output directory
OUTDIR="/hb/home/snbogan/3seas/data/fastqs"
mkdir -p "$OUTDIR"

# List of all accessions embedded in script
ACCESSIONS=(
SRR21814094 SRR21814095 SRR21814096 SRR21814097 SRR21814098 SRR21814099 SRR21814100
SRR21814101 SRR21814102 SRR21814103 SRR21814104 SRR21814105 SRR21814106 SRR21814107
SRR21814108 SRR21814109 SRR21814110 SRR21814111 SRR21814112 SRR21814113 SRR21814114
SRR21814115 SRR21814116 SRR21814117 SRR21814118 SRR21814119 SRR21814120 SRR21814121
SRR21814122 SRR21814123 SRR21814124 SRR21814125 SRR21814126 SRR21814127 SRR21814128
SRR21814129 SRR21814130 SRR21814131 SRR21814132 SRR21814133 SRR21814134 SRR21814135
SRR21814136 SRR21814137 SRR21814138 SRR21814139 SRR21814140 SRR21814141 SRR21814142
SRR21814143 SRR21814144 SRR21814145 SRR21814146 SRR21814147 SRR21814148 SRR21814149
SRR21814150 SRR21814151 SRR21814152 SRR21814153 SRR21814154 SRR21814155 SRR21814156
SRR21814157 SRR21814158 SRR21814159 SRR21814160 SRR21814161 SRR21814162 SRR21814163
SRR21814164 SRR21814165 SRR21814166 SRR21814167 SRR21814168 SRR21814169 SRR21814170
SRR21814171 SRR21814172 SRR21814173 SRR21814174 SRR21814175 SRR21814176 SRR21814177
SRR21814178 SRR21814179 SRR21814180 SRR21814181 SRR21814182 SRR21814183 SRR21814184
SRR21814185 SRR21814186 SRR21814187 SRR21814188 SRR21814189 SRR21814190 SRR21814191
SRR21814192 SRR21814193 SRR21814194 SRR21814195 SRR21814196 SRR21814197 SRR21814198
SRR21814199 SRR21814200 SRR21814201 SRR21814202 SRR21814203 SRR21814204 SRR21814205
SRR21814206 SRR21814207 SRR21814208 SRR21814209 SRR21814210 SRR21814211 SRR21814212
SRR21814213 SRR21814214 SRR21814215 SRR21814216 SRR21814217 SRR21814218 SRR21814219
SRR21814220 SRR21814221 SRR21814222 SRR21814223 SRR21814224 SRR21814225 SRR21814226
SRR21814227 SRR21814228 SRR21814229 SRR21814230 SRR21814231 SRR21814232 SRR21814233
SRR21814234 SRR21814235 SRR21814236 SRR21814237 SRR21814238 SRR21814239 SRR21814240
SRR21814241 SRR21814242 SRR21814243 SRR21814244 SRR21814245 SRR21814246 SRR21814247
SRR21814248 SRR21814249 SRR21814250 SRR21814251 SRR21814252 SRR21814253 SRR21814254
SRR21814255 SRR21814256 SRR21814257 SRR21814258 SRR21814259 SRR21814260 SRR21814261
SRR21814262 SRR21814263 SRR21814264 SRR21814265 SRR21814266 SRR21814267 SRR21814268
SRR21814269 SRR21814270 SRR21814271 SRR21814272 SRR21814273 SRR21814274 SRR21814275
SRR21814276 SRR21814277 SRR21814278 SRR21814279 SRR21814280 SRR21814281 SRR21814282
SRR21814283 SRR21814284 SRR21814285 SRR21814286 SRR21814287 SRR21814288 SRR21814289
SRR21814290 SRR21814291 SRR21814292 SRR21814293 SRR21814294 SRR21814295 SRR21814296
SRR21814297 SRR21814298 SRR21814299 SRR21814300 SRR21814301 SRR21814302 SRR21814303
SRR21814304 SRR21814305 SRR21814306 SRR21814307 SRR21814308 SRR21814309 SRR21814310
SRR21814311 SRR21814312 SRR21814313 SRR21814314 SRR21814315 SRR21814316 SRR21814317
SRR21814318 SRR21814319 SRR21814320 SRR21814321 SRR21814322 SRR21814323 SRR21814324
SRR21814325 SRR21814326 SRR21814327 SRR21814328 SRR21814329 SRR21814330 SRR21814331
SRR21814332 SRR21814333
)

# Get accession for this array task
ACCESSION="${ACCESSIONS[$SLURM_ARRAY_TASK_ID-1]}"

# Verify we got a valid accession
if [[ -z "$ACCESSION" ]]; then
    echo "ERROR: No accession found for array index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "Processing accession: $ACCESSION (Task $SLURM_ARRAY_TASK_ID of ${#ACCESSIONS[@]})"

# Step 1: Prefetch the SRA file
echo "Prefetching $ACCESSION..."
prefetch --max-size 100G "$ACCESSION" --output-directory "$OUTDIR"

# Check if prefetch succeeded
if [[ $? -ne 0 ]]; then
    echo "ERROR: prefetch failed for $ACCESSION"
    exit 1
fi

# Step 2: Convert to FASTQ
echo "Converting $ACCESSION to FASTQ..."
fastq-dump --gzip --split-files "$OUTDIR/$ACCESSION/$ACCESSION.sra" -O "$OUTDIR"

# Check if conversion succeeded
if [[ $? -ne 0 ]]; then
    echo "ERROR: fastq-dump failed for $ACCESSION"
    exit 1
fi

# Cleanup: Remove the .sra file and directory
echo "Cleaning up..."
rm -rf "$OUTDIR/$ACCESSION"

echo "Successfully processed $ACCESSION"
