#!/bin/bash
##########################################################
# QIIME2 16S Amplicon Analysis Pipeline (Template)
# Applicable: 16S paired-end sequencing data
# Environments: qc_preprocess (QC/trimming) + qiime2-2025.7 (core analysis)
# Database: SILVA138.2_SSURef_NR99_uniform_classifier_full-length.qza
##########################################################
# First use: Run qiime2-16s-pipeline_install_en.sh to install environments,
#            then run this script for data analysis
##########################################################

## 0. Parameter settings (modify according to your data)
wd=~/amplicon_analysis
metadata=metadata.txt
# classifiers/SILVA138.2_SSURef_NR99_uniform_classifier_full-length.qza must be pre-downloaded to classifiers/
# Download: https://pan.baidu.com/s/1yFOvSJXofc6H7AahHSJq2A?pwd=5gdz
#      Or obtain from the SILVA website: https://www.arb-silva.de/
# Note: SILVA138.2 classifier is only compatible with qiime2 ≥ 2023.x (this pipeline uses qiime2-2025.7)
#       If using qiime2-2023.2, please use silva-138-99-nb-classifier.qza instead

## 1. Directory initialization and file placement
# Before running this script, place the following files in the specified locations:
#   - Raw sequencing data: *.fq.gz into seq/ directory (naming: sampleName_1.fq.gz / sampleName_2.fq.gz)
#   - Metadata file: metadata.txt in the working directory
#     Format: TSV (tab-separated), first row is header, first column is sample ID, must include Group column
#     Example (see: examples/metadata.txt):
#       SampleID  Group  ...
#       WT1       WT     ...
#       KO1       KO     ...
mkdir -p ${wd}/{seq,trimmed,qiime2,logs,results/{fastqc_raw,fastqc_trimmed,cutadapt_logs,export}}
cd ${wd}

# ln /path/to/raw/seq/*.fq.gz seq/
# ln /path/to/metadata.txt ./

## 2. Environment verification
conda activate qc_preprocess
fastqc --version
cutadapt --version
multiqc --version

conda activate qiime2-2025.7
qiime --version

## 3. Raw data quality control
conda activate qc_preprocess
cd ${wd}

fastqc seq/*_1.fq.gz seq/*_2.fq.gz -t 8 -o results/fastqc_raw
multiqc results/fastqc_raw/ -o results/fastqc_raw/ -n multiqc_report_raw.html
# Check: results/fastqc_raw/multiqc_report_raw.html

## 4. Primer trimming
conda activate qc_preprocess
cd ${wd}

for f in seq/*_1.fq.gz; do
    base=$(basename "${f}" _1.fq.gz)
    r="seq/${base}_2.fq.gz"
    [ ! -f "${r}" ] && echo "Warning: reverse file missing ${r}" && continue

    cutadapt \
        -g GTGCCAGCMGCCGCGG \
        -G CCGTCAATTCMTTTRAGTTT \
        --pair-filter=any \
        --minimum-length 150 \
        --quality-cutoff 20,20 \
        --max-n 0 \
        -j 8 \
        -o trimmed/"${base}_1.fq.gz" \
        -p trimmed/"${base}_2.fq.gz" \
        "${f}" "${r}" \
        > results/cutadapt_logs/"${base}.log" 2>&1
done
# Check: results/cutadapt_logs/ for individual sample logs, focus on Pairs written percentage

## 4a. Primer trimming results summary
conda activate qc_preprocess
cd ${wd}

echo "=== Sequence Pass Rate Summary ===" > results/cutadapt_logs/summary_report.txt

for f in seq/*_1.fq.gz; do
    base=$(basename "${f}" _1.fq.gz)
    log="results/cutadapt_logs/${base}.log"

    [ ! -f "${log}" ] && echo "Warning: log file missing ${log}" | tee -a results/cutadapt_logs/summary_report.txt && continue

    echo "Sample: ${base}" >> results/cutadapt_logs/summary_report.txt
    grep "Total read pairs processed" "${log}" >> results/cutadapt_logs/summary_report.txt
    grep "Pairs written" "${log}" >> results/cutadapt_logs/summary_report.txt
    echo "" >> results/cutadapt_logs/summary_report.txt
done
# Check: results/cutadapt_logs/summary_report.txt for Pairs written percentage per sample

## 5. Post-trimming quality control
conda activate qc_preprocess
cd ${wd}

fastqc trimmed/*_1.fq.gz trimmed/*_2.fq.gz -t 8 -o results/fastqc_trimmed
multiqc results/fastqc_trimmed/ -o results/fastqc_trimmed/ -n multiqc_report_trimmed.html
# Check: results/fastqc_trimmed/multiqc_report_trimmed.html, mean quality should be ≥Q25

## 6. Generate QIIME2 manifest file
conda activate qiime2-2025.7
cd ${wd}

awk -v pwd="$PWD" 'BEGIN {FS=OFS="\t"}
    NR==1 {print "sample-id\tforward-absolute-filepath\treverse-absolute-filepath"}
    NR>1 {print $1, pwd "/trimmed/" $1 "_1.fq.gz", pwd "/trimmed/" $1 "_2.fq.gz"}' \
    ${metadata} > manifest
# Check: head manifest, verify paths and sample IDs

## 7. Import into QIIME2
conda activate qiime2-2025.7
cd ${wd}

time qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path manifest \
    --output-path qiime2/demux.qza \
    --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
    --i-data qiime2/demux.qza \
    --o-visualization qiime2/demux.qzv
# Check: qiime2/demux.qzv, view sequencing depth distribution across samples

## 8. DADA2 denoising (adjust parameters based on data quality)
conda activate qiime2-2025.7
cd ${wd}

nohup qiime dada2 denoise-paired \
    --i-demultiplexed-seqs qiime2/demux.qza \
    --p-trim-left-f 0 --p-trim-left-r 0 \
    --p-trunc-len-f 230 --p-trunc-len-r 217 \
    --p-max-ee-f 5.0 --p-max-ee-r 5.0 \
    --p-trunc-q 10 \
    --p-n-threads 8 \
    --p-chimera-method consensus \
    --o-table qiime2/table.qza \
    --o-representative-sequences qiime2/rep-seqs.qza \
    --o-denoising-stats qiime2/denoising-stats.qza \
    > logs/dada2.log 2>&1 &
echo "DADA2 running in background, check progress: tail -f logs/dada2.log"
echo "Continue with downstream commands after completion"
# Note: 'pooled' was removed in QIIME2 2025.4; 2025.7 supports 'consensus' and 'none' only

qiime metadata tabulate \
    --m-input-file qiime2/denoising-stats.qza \
    --o-visualization qiime2/denoising-stats.qzv

qiime feature-table summarize \
    --i-table qiime2/table.qza \
    --o-visualization qiime2/table.qzv \
    --m-sample-metadata-file ${metadata}

qiime feature-table tabulate-seqs \
    --i-data qiime2/rep-seqs.qza \
    --o-visualization qiime2/rep-seqs.qzv
# Check: qiime2/table.qzv for sample sequencing depth, qiime2/denoising-stats.qzv for denoising efficiency

## 9. Taxonomic classification
conda activate qiime2-2025.7
cd ${wd}

nohup qiime feature-classifier classify-sklearn \
    --i-classifier classifiers/SILVA138.2_SSURef_NR99_uniform_classifier_full-length.qza \
    --i-reads qiime2/rep-seqs.qza \
    --o-classification qiime2/taxonomy.qza \
    --p-n-jobs 8 \
    > logs/classify.log 2>&1 &
echo "Taxonomic classification running in background, check progress: tail -f logs/classify.log"

qiime metadata tabulate \
    --m-input-file qiime2/taxonomy.qza \
    --o-visualization qiime2/taxonomy.qzv

## 10. Phylogenetic tree construction
conda activate qiime2-2025.7
cd ${wd}

nohup qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences qiime2/rep-seqs.qza \
    --o-alignment qiime2/aligned-rep-seqs.qza \
    --o-masked-alignment qiime2/masked-aligned-rep-seqs.qza \
    --o-tree qiime2/unrooted-tree.qza \
    --o-rooted-tree qiime2/rooted-tree.qza \
    --p-n-threads 8 \
    > logs/tree.log 2>&1 &
echo "Phylogenetic tree running in background, check progress: tail -f logs/tree.log"

## 11. Alpha rarefaction curve (select sampling depth)
conda activate qiime2-2025.7
cd ${wd}

# --p-max-depth refers to the maximum sample sequence count in table.qzv, modify as needed
qiime diversity alpha-rarefaction \
    --i-table qiime2/table.qza \
    --i-phylogeny qiime2/rooted-tree.qza \
    --p-max-depth 11576 \
    --m-metadata-file ${metadata} \
    --o-visualization results/alpha-rarefaction.qzv
# Check: results/alpha-rarefaction.qzv, observe curve plateau
# The plateau depth is the appropriate sampling depth, also refer to sample sequence distribution in table.qzv
# 【Replace sampling_depth below with the selected value】

sampling_depth=11576

## 12. Core diversity analysis (Alpha + Beta)
conda activate qiime2-2025.7
cd ${wd}

rm -rf results/core-metrics-results
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny qiime2/rooted-tree.qza \
    --i-table qiime2/table.qza \
    --p-sampling-depth ${sampling_depth} \
    --m-metadata-file ${metadata} \
    --output-dir results/core-metrics-results
# Outputs: faith_pd / shannon / observed_features / evenness (alpha)
#          unweighted_unifrac / weighted_unifrac / bray_curtis / jaccard (beta)

## 13. Taxonomic composition bar plot
conda activate qiime2-2025.7
cd ${wd}

qiime taxa barplot \
    --i-table qiime2/table.qza \
    --i-taxonomy qiime2/taxonomy.qza \
    --m-metadata-file ${metadata} \
    --o-visualization results/taxa-bar-plots.qzv

## 14. Export QIIME2 results to plain text (for downstream R/Python analysis)
conda activate qiime2-2025.7
cd ${wd}

mkdir -p results/export

# 14a. Rarefied ASV abundance table (export first, rename to avoid overwrite by subsequent feature-table.biom)
qiime tools export \
    --input-path results/core-metrics-results/rarefied_table.qza \
    --output-path results/export
mv results/export/feature-table.biom results/export/rarefied_table.biom
biom convert -i results/export/rarefied_table.biom \
    -o results/export/rarefied_table.tsv --to-tsv
# Output: results/export/rarefied_table.tsv (rarefied ASV abundance table)

# 14b. Raw ASV/OTU feature table
qiime tools export --input-path qiime2/table.qza \
    --output-path results/export
biom convert -i results/export/feature-table.biom \
    -o results/export/feature-table.tsv --to-tsv
# View total sequence count (total_frequency)
biom summarize-table -i results/export/feature-table.biom
# Output: results/export/feature-table.tsv (rows = ASVs, columns = samples, values = sequence counts)

# 14c. Taxonomic classification
qiime tools export --input-path qiime2/taxonomy.qza \
    --output-path results/export
# Output: results/export/taxonomy.tsv (ASV ID → taxonomy + confidence)

# 14d. Representative sequences
qiime tools export --input-path qiime2/rep-seqs.qza \
    --output-path results/export
# Output: results/export/dna-sequences.fasta

# 14e. Denoising statistics
qiime tools export --input-path qiime2/denoising-stats.qza \
    --output-path results/export

# 14f. Alpha diversity (Faith PD / Shannon / Observed Features / Evenness)
qiime tools export \
    --input-path results/core-metrics-results/faith_pd_vector.qza \
    --output-path results/export
qiime tools export \
    --input-path results/core-metrics-results/shannon_vector.qza \
    --output-path results/export
qiime tools export \
    --input-path results/core-metrics-results/observed_features_vector.qza \
    --output-path results/export
qiime tools export \
    --input-path results/core-metrics-results/evenness_vector.qza \
    --output-path results/export

# 14g. Beta diversity distance matrices
qiime tools export \
    --input-path results/core-metrics-results/unweighted_unifrac_distance_matrix.qza \
    --output-path results/export
qiime tools export \
    --input-path results/core-metrics-results/weighted_unifrac_distance_matrix.qza \
    --output-path results/export
qiime tools export \
    --input-path results/core-metrics-results/bray_curtis_distance_matrix.qza \
    --output-path results/export
qiime tools export \
    --input-path results/core-metrics-results/jaccard_distance_matrix.qza \
    --output-path results/export

# 14h. PCoA coordinates
qiime tools export \
    --input-path results/core-metrics-results/unweighted_unifrac_pcoa_results.qza \
    --output-path results/export
qiime tools export \
    --input-path results/core-metrics-results/weighted_unifrac_pcoa_results.qza \
    --output-path results/export
qiime tools export \
    --input-path results/core-metrics-results/bray_curtis_pcoa_results.qza \
    --output-path results/export
qiime tools export \
    --input-path results/core-metrics-results/jaccard_pcoa_results.qza \
    --output-path results/export

# 14i. Phylogenetic tree (Newick format)
qiime tools export \
    --input-path qiime2/rooted-tree.qza \
    --output-path results/export

## 15. One-step R visualization (copy R script to project root first)
conda activate qiime2-2025.7
cd ${wd}

# Copy examples/QIIME2_16S_visualization_EN.R to $(pwd), modify setwd(),
# then the script reads from results/export/ and outputs to subdirectories,
# also prepares FAPROTAX/PICRUSt2 input files
if [ -f "QIIME2_16S_visualization_EN.R" ]; then
    Rscript QIIME2_16S_visualization_EN.R
elif [ -f "QIIME2_16S_visualization.R" ]; then
    Rscript QIIME2_16S_visualization.R
else
    echo "Please copy examples/QIIME2_16S_visualization_EN.R or examples/QIIME2_16S_visualization.R to $(pwd)"
fi

## 16. FAPROTAX functional prediction
conda activate qiime2-2025.7
cd ${wd}

cd ${wd}/results/export/faprotax

# FAPROTAX script path (modify as needed)
# Download: http://www.loucalab.com/archive/FAPROTAX/lib/php/index.php?section=Download
sd=~/db/EasyMicrobiome/script/FAPROTAX_1.2.12

# Prepare BIOM input with taxonomy metadata
biom add-metadata \
    -i rarefied_table.biom \
    --observation-metadata-fp taxonomy.tsv \
    -o rarefied_tax.biom \
    --sc-separated taxonomy \
    --observation-header OTUID,taxonomy

# FAPROTAX collapse
python ${sd}/collapse_table.py \
    -i rarefied_tax.biom \
    -g ${sd}/FAPROTAX.txt \
    --collapse_by_metadata 'taxonomy' \
    -v --force \
    -o faprotax.txt \
    -r faprotax_report.txt

# Generate presence/absence matrix
grep '*' -B 1 faprotax_report.txt | grep -v -P '^--$' > faprotax_report.clean
perl ${sd}/../faprotax_report_sum.pl \
    -i faprotax_report.clean \
    -o faprotax_report

echo "FAPROTAX done: results/export/faprotax/faprotax_report.txt"

## 17. PICRUSt2 functional prediction
conda activate picrust2
cd ${wd}

cd ${wd}/results/export/picrust2

# Run PICRUSt2 pipeline (background, wait for completion before annotations)
nohup picrust2_pipeline.py -s dna-sequences.fasta -i feature-table.tsv \
    -o ./out -p 16 > picrust2.log 2>&1 &
echo "PICRUSt2 running in background, check: tail -f results/export/picrust2/picrust2.log"

# Add EC/KO/Pathway annotations (run after PICRUSt2 completes)
cd out
add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
  -o METACYC.tsv
add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
  -o EC.tsv
add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
  -o KO.tsv

# KEGG hierarchy merge (requires EasyMicrobiome database)
db=~/db/EasyMicrobiome/
python3 ${db}/script/summarizeAbundance.py \
    -i KO.tsv \
    -m ${db}/kegg/KO1-4.txt \
    -c 2,3,4 -s ',+,+,' -n raw \
    -o KEGG
wc -l KEGG*

# Export and functional prediction complete: results/export/
# Key files:
#   rarefied_table.tsv       — Rarefied ASV abundance table
#   feature-table.tsv        — ASV/OTU abundance table
#   taxonomy.tsv             — Taxonomic classification
#   dna-sequences.fasta      — Representative sequences
#   *.tsv (alpha)            — Alpha diversity (Faith PD, Shannon, etc.)
#   *_distance_matrix.tsv    — Beta diversity distance matrices
#   *_pcoa_results.tsv       — PCoA ordination coordinates
#   faprotax/                — FAPROTAX functional prediction results
#   picrust2/                — PICRUSt2 functional prediction results
#
# One-step visualization: copy examples/QIIME2_16S_visualization_EN.R to
# project root, modify setwd(), run Rscript, reads results/export/ automatically
