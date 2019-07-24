# Supplementary Note 3: DADA2, QIIME2 and QIIME1 pipelines

Following workflow recommendations for each pipeline, DADA2 was used with SILVA v132 NR database, QIIME2 with Greengenes database clustered at 99% identity and Deblur denoiser, and QIIME1 with Greengenes database clustered at 97% identity both using Open Reference and De Novo OTU picking. The main text presents the (recommended) workflow using Open Reference OTUs.

## Mock: DADA2 SILVA pipeline

```R
#--------------- DADA2 pipeline for Zheng 2015 et al. -------------------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dada2))

main_path = '~/data1/'
set_path = file.path(main_path, 'zheng-2015')
setwd(set_path)

# Follow the DADA2 pipeline:
# https://benjjneb.github.io/dada2/tutorial.html

# Define some helper functions
fasta_writer = function (meta, seqs, output, ncpu=detectCores()-1) {
  # Write FASTA sequences from meta data and sequences
  writeChar(paste(paste0('>', meta),
                  seqs,
                  sep='\n', collapse='\n'
            ),
            output, eos=NULL)
}

#--------------- Filter and trim files ----------------------------------
fq_path = file.path(set_path,'fastq_v3v4')
fq_fwd = sort(list.files(fq_path, pattern='.*R1.*gz', full.names=T))
fq_rev = sort(list.files(fq_path, pattern='.*R2.*gz', full.names=T))
sample_ids = sapply(strsplit(basename(fq_fwd), '_', fixed=T), `[`, 1)
# Check quality profiles
plotQualityProfile(fq_fwd[1:2])
plotQualityProfile(fq_rev[1:2])
# Trim last 100 nts for reverse reads. Generate output filenames:
fq_fwd_fil = file.path(set_path, 'dada2_analysis/filtered',
                       paste0(sample_ids, '_R1_filtered.fastq'))
fq_rev_fil = file.path(set_path, 'dada2_analysis/filtered',
                       paste0(sample_ids, '_R2_filtered.fastq'))
# Filter and trim
ft_out = filterAndTrim(fq_fwd, fq_fwd_fil, fq_rev, fq_rev_fil,
                       trimLeft=c(22,22), truncLen=c(300,200), maxN=0, maxEE=c(2,2),
                       truncQ=2, rm.phix=T, compress=T, multithread=T)
# Check the number of retained reads after paired filter
head(ft_out)
# Most reads retained.

#--------------------- DADA2 denoising ----------------------------------
# Learn errors for fwd and rev reads separately
# This step takes few hours on 2016 Macbook Pro. In my tests with other
# data sets, it works equally well to use much less than 1e6 reads for
# learning error rates. 10-100k read sample was still fine.
err_fwd = learnErrors(fq_fwd_fil, multithread=T)
err_rev = learnErrors(fq_rev_fil, multithread=T)

# Save intermediate R files, so we can resume later if needed.
saveRDS(err_fwd, 'dada2_analysis/err_fwd')
saveRDS(err_rev, 'dada2_analysis/err_rev')

# Plot errors to check how well learned errors fit
plotErrors(err_fwd, nominalQ=T)
plotErrors(err_rev, nominalQ=T)

# Looks good. Run dada2.
derepFs = derepFastq(fq_fwd_fil, verbose=T)
derepRs = derepFastq(fq_rev_fil, verbose=T)
# Name the derep-class objects by the sample names
names(derepFs) = sample_ids
names(derepRs) = sample_ids
dadaFs = dada(derepFs, err=err_fwd, multithread=T)
dadaRs = dada(derepRs, err=err_rev, multithread=T)
saveRDS(dadaFs, 'dada2_analysis/dadaFs')
saveRDS(dadaRs, 'dada2_analysis/dadaRs')

# Merge reads
mergers = mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=T)
saveRDS(mergers, 'dada2_analysis/mergers')

# Generate table with sequence counts
seqtab = makeSequenceTable(mergers)
saveRDS(seqtab, 'dada2_analysis/seqtab')

# Chimera removal
seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus",
                                   multithread=T, verbose=T)
saveRDS(seqtab.nochim, 'dada2_analysis/seqtab.nochim')
# For multiple samples collapse:
# seqtab.nochim.coll = collapseNoMismatch(seqtab.nochim)

# Fraction of chimeric reads
cat('Fraction of non-chimeric reads: ', sum(seqtab.nochim)/sum(seqtab), '\n')
# 94%. That looks good.

# Assign taxonomy with SILVA
taxa = assignTaxonomy(seqtab.nochim,
                      'dada2_analysis/db/silva_nr_v132_train_set.fa.gz',
                      multithread=T, verbose=T)
saveRDS(taxa, 'dada2_analysis/taxa_before_addSpecies')
# Add species assignment
taxa = addSpecies(taxa, 'dada2_analysis/db/silva_species_assignment_v132.fa.gz',
                  allowMultiple=T, verbose=T)
saveRDS(taxa, 'dada2_analysis/taxa')

# Print assignments
taxa.print = taxa
rownames(taxa.print) = NULL

# Save taxonomy abundance and FASTA sequences
tax.dt = as.data.table(taxa, keep.rownames=T)
write.table(tax.dt, 'dada2_analysis/taxonomy.txt',
            sep='\t', row.names=F, quote=F)
ab.dt = as.data.table(t(seqtab.nochim), keep.rownames=T)
write.table(ab.dt,
            'dada2_analysis/abundances.txt',
            sep='\t', row.names=F, quote=F)

# Write FASTA file with sequences
seqs.dt = ab.dt[, .(rn, id=1:.N)]
fasta_writer(seqs.dt[, id], seqs.dt[, rn],
             'dada2_analysis/sequences.fasta')

# We manually BLAST these 20 sequences against NCBI database to identify
# the 20 mock species/strains.
# blastn -query sequences.fasta \
#        -db /data1/igor/himap/db/V3-V4_341F-805R_hang22_uniq \
#        -num_threads 4 -outfmt 6 > sequences_blast.txt
# Then we used these assignments to check which sequence was found or missing.
```


## Mock: QIIME2 GG99 pipeline

In the Zheng et al. 2015 mock community dataset PCR primers were still present (unlike the DIABIMMUNE data, where this step is skipped!). We found that trimming first 22 nt instead of using cutadapt (can used through QIIME2) in this case removes PCR primers more accurately. This is an important step, because we noticed that a PCR primer with a mismatch to the exact sequence can be more common than the one with an exact match. This single nt mismatch might make taxonomy prediction more difficult later down the pipeline.

```sh
# Trim PCR primers manually (tried cutadapt within QIIME2 and it doesn't remove all of them for some reason)

mkdir fastq_trimmed
vsearch --fastx_filter ../fastq_v3v4/V3V4Rep1_S5_L001_R1_001.fastq.gz --fastq_stripleft 22 --fastqout fastq_trimmed/V3V4Rep1_S5_L001_R1_001.fastq
vsearch --fastx_filter ../fastq_v3v4/V3V4Rep1_S5_L001_R2_001.fastq.gz --fastq_stripleft 22 --fastqout fastq_trimmed/V3V4Rep1_S5_L001_R2_001.fastq
vsearch --fastx_filter ../fastq_v3v4/V3V4Rep2_S6_L001_R1_001.fastq.gz --fastq_stripleft 22 --fastqout fastq_trimmed/V3V4Rep2_S6_L001_R1_001.fastq
vsearch --fastx_filter ../fastq_v3v4/V3V4Rep2_S6_L001_R2_001.fastq.gz --fastq_stripleft 22 --fastqout fastq_trimmed/V3V4Rep2_S6_L001_R2_001.fastq
vsearch --fastx_filter ../fastq_v3v4/V3V4Rep3_S7_L001_R1_001.fastq.gz --fastq_stripleft 22 --fastqout fastq_trimmed/V3V4Rep3_S7_L001_R1_001.fastq
vsearch --fastx_filter ../fastq_v3v4/V3V4Rep3_S7_L001_R2_001.fastq.gz --fastq_stripleft 22 --fastqout fastq_trimmed/V3V4Rep3_S7_L001_R2_001.fastq
vsearch --fastx_filter ../fastq_v3v4/V3V4Rep4_S8_L001_R1_001.fastq.gz --fastq_stripleft 22 --fastqout fastq_trimmed/V3V4Rep4_S7_L001_R1_001.fastq
vsearch --fastx_filter ../fastq_v3v4/V3V4Rep4_S8_L001_R2_001.fastq.gz --fastq_stripleft 22 --fastqout fastq_trimmed/V3V4Rep4_S7_L001_R2_001.fastq

cd fastq_trimmed
gzip *.fastq
cd ..
```

Manually create a manifest.txt files:

```
sample-id,absolute-filepath,direction
V3V4Rep1,/Users/igor/cloud/research/microbiome/zheng-2015/qiime2_analysis/fastq_trimmed/V3V4Rep1_S5_L001_R1_001.fastq.gz,forward
V3V4Rep1,/Users/igor/cloud/research/microbiome/zheng-2015/qiime2_analysis/fastq_trimmed/V3V4Rep1_S5_L001_R2_001.fastq.gz,reverse
V3V4Rep2,/Users/igor/cloud/research/microbiome/zheng-2015/qiime2_analysis/fastq_trimmed/V3V4Rep2_S6_L001_R1_001.fastq.gz,forward
V3V4Rep2,/Users/igor/cloud/research/microbiome/zheng-2015/qiime2_analysis/fastq_trimmed/V3V4Rep2_S6_L001_R2_001.fastq.gz,reverse
V3V4Rep3,/Users/igor/cloud/research/microbiome/zheng-2015/qiime2_analysis/fastq_trimmed/V3V4Rep3_S7_L001_R1_001.fastq.gz,forward
V3V4Rep3,/Users/igor/cloud/research/microbiome/zheng-2015/qiime2_analysis/fastq_trimmed/V3V4Rep3_S7_L001_R2_001.fastq.gz,reverse
V3V4Rep4,/Users/igor/cloud/research/microbiome/zheng-2015/qiime2_analysis/fastq_trimmed/V3V4Rep4_S8_L001_R1_001.fastq.gz,forward
V3V4Rep4,/Users/igor/cloud/research/microbiome/zheng-2015/qiime2_analysis/fastq_trimmed/V3V4Rep4_S8_L001_R2_001.fastq.gz,reverse
```

```sh
# Import paired-end FASTQ files into artifact
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path manifest_trimmed.txt \
    --output-path qiime_import_trimmed \
    --source-format PairedEndFastqManifestPhred33

# Merge reads
qiime vsearch join-pairs \
    --i-demultiplexed-seqs qiime_import_trimmed.qza \
    --o-joined-sequences qiime_import_trimmed_merged.qza

# Quality filter merged reads
qiime quality-filter q-score-joined \
    --i-demux qiime_import_trimmed_merged.qza \
    --o-filtered-sequences qiime_import_trimmed_merged_filtered.qza \
    --o-filter-stats qiime_import_trimmed_merged_filtered_stats.qza

# Run Deblur denoising
qiime deblur denoise-16S \
    --i-demultiplexed-seqs qiime_import_trimmed_merged_filtered.qza \
    --p-trim-length 393 \
    --o-representative-sequences qiime_import_trimmed_merged_filtered_rep-seqs.qza \
    --o-table qiime_import_trimmed_merged_filtered_table.qza \
    --p-sample-stats \
    --o-stats qiime_import_trimmed_merged_filtered_deblur_stats.qza

# Export data from artifacts into normal files
qiime tools export qiime_import_trimmed_merged_filtered_table.qza \
    --output-dir qiime_import_trimmed_merged_filtered_feature_table
qiime tools export qiime_import_trimmed_merged_filtered_rep-seqs.qza \
    --output-dir qiime_import_trimmed_merged_filtered_feature_table
biom convert -i qiime_import_trimmed_merged_filtered_feature_table/feature-table.biom \
    -o qiime_import_trimmed_merged_filtered_feature_table/feature-table.txt --to-tsv
```

For taxonomic classification, we have to do a bit more work since we couldn’t find precomputed artifact file for v3-v4 region on QIIME2 website. There’s only precomputed files for either full-length 16S sequences or for weirdly short v4 region (120 nt?) at https://docs.qiime2.org/2018.2/data-resources/ (scroll a bit down).
Let’s follow the guide from here https://docs.qiime2.org/2018.2/tutorials/feature-classifier/  to train Naive Bayes classifier on V3-V4 region of 99% OTU Greengenes 13.8 sequences. First download Greengenes 13.8 99% OTU files:

```sh
# Download fasta and taxonomy from Greengenes
mkdir greengenes
cd greengenes
wget ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
tar xvfz gg_13_8_otus.tar.gz
rm gg_13_8_otus.tar.gz
mv gg_13_8_otus/rep_set/99_otus.fasta .
mv gg_13_8_otus/taxonomy/99_otu_taxonomy.txt .
rm -rf gg_13_8_otus
cd ..

# Import FASTA and txt files as QIIME2 artifacts
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path greengenes/99_otus.fasta \
  --output-path 99_otus.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --source-format HeaderlessTSVTaxonomyFormat \
  --input-path greengenes/99_otu_taxonomy.txt \
  --output-path ref-taxonomy.qza
```

We will also need primer sequences for V3-V4 region. Forward primer 341F 5’-3’ sequence: `GACAGCCTACGGGNGGCWGCAG`. Reverse primer 805R 3’-5’ sequence: `GACTACCAGGGTATCTAATC`. Now we can extract reference reads and truncate to 419 nt:

```sh
# Extract V3-V4 regions
qiime feature-classifier extract-reads \
    --i-sequences 99_otus.qza \
    --p-f-primer GACAGCCTACGGGNGGCWGCAG \
    --p-r-primer GACTACCAGGGTATCTAATC \
    --p-trunc-len 419 \
    --o-reads 99_otus_refseqs.qza

# Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads 99_otus_refseqs.qza \
    --i-reference-taxonomy ref-taxonomy.qza \
    --o-classifier 99_otus_v3-v4_341f-805r_classifier.qza 
```

Now run this feature classifier to generate taxonomy for the sequences in our data:

```sh
qiime feature-classifier classify-sklearn \
  --i-classifier 99_otus_v3-v4_341f-805r_classifier.qza \
  --i-reads qiime_import_trimmed_merged_filtered_rep-seqs.qza \
  --o-classification qiime_import_trimmed_merged_filtered_taxonomy.qza

# Export taxonomy to tab-delimited file
qiime tools export qiime_import_trimmed_merged_filtered_taxonomy.qza \
  --output-dir qiime_import_trimmed_merged_filtered_taxonomy
```


## Mock: QIIME1 GG97 pipeline
Start macOS session using `macqiime`. This uses QIIME 1.9. Run everything from zheng-2015 folder. USEARCH 6.1 needs to be manually downloaded from http://www.drive5.com and be present in PATH as usearch61 executable for the chimera removal part.

### Merge reads
Start from PCR primer trimmed reads that we prepared from QIIME 2 analysis. The merging is done using (the default) fastq-join algorithm:

```sh
join_paired_ends.py -f qiime2_analysis/fastq_trimmed/V3V4Rep1_S5_L001_R1_001.fastq.gz \
                    -r qiime2_analysis/fastq_trimmed/V3V4Rep1_S5_L001_R2_001.fastq.gz \
                    -o qiime1_analysis/trimmed_merged/V3V4Rep1_merged.fastq

join_paired_ends.py -f qiime2_analysis/fastq_trimmed/V3V4Rep2_S6_L001_R1_001.fastq.gz \
                    -r qiime2_analysis/fastq_trimmed/V3V4Rep2_S6_L001_R2_001.fastq.gz \
                    -o qiime1_analysis/trimmed_merged/V3V4Rep2_merged.fastq

join_paired_ends.py -f qiime2_analysis/fastq_trimmed/V3V4Rep3_S7_L001_R1_001.fastq.gz \
                    -r qiime2_analysis/fastq_trimmed/V3V4Rep3_S7_L001_R2_001.fastq.gz \
                    -o qiime1_analysis/trimmed_merged/V3V4Rep3_merged.fastq

join_paired_ends.py -f qiime2_analysis/fastq_trimmed/V3V4Rep4_S8_L001_R1_001.fastq.gz \
                    -r qiime2_analysis/fastq_trimmed/V3V4Rep4_S8_L001_R2_001.fastq.gz \
                    -o qiime1_analysis/trimmed_merged/V3V4Rep4_merged.fastq
```

### Quality control and filtering

```sh
split_libraries_fastq.py -i qiime1_analysis/trimmed_merged/V3V4Rep1_merged.fastq/fastqjoin.join.fastq \
                         --sample_ids V3V4Rep1 \
                         -o qiime1_analysis/V3V4Rep1_quality_filtered_q20/ -q 19 \
                         --barcode_type 'not-barcoded' --phred_offset=33

split_libraries_fastq.py -i qiime1_analysis/trimmed_merged/V3V4Rep2_merged.fastq/fastqjoin.join.fastq \
                         --sample_ids V3V4Rep2 \
                         -o qiime1_analysis/V3V4Rep2_quality_filtered_q20/ -q 19 \
                         --barcode_type 'not-barcoded' --phred_offset=33

split_libraries_fastq.py -i qiime1_analysis/trimmed_merged/V3V4Rep3_merged.fastq/fastqjoin.join.fastq \
                         --sample_ids V3V4Rep3 \
                         -o qiime1_analysis/V3V4Rep3_quality_filtered_q20/ -q 19 \
                         --barcode_type 'not-barcoded' --phred_offset=33

split_libraries_fastq.py -i qiime1_analysis/trimmed_merged/V3V4Rep4_merged.fastq/fastqjoin.join.fastq \
                         --sample_ids V3V4Rep4 \
                         -o qiime1_analysis/V3V4Rep4_quality_filtered_q20/ -q 19 \
                         --barcode_type 'not-barcoded' --phred_offset=33
```

### Chimera removal

```sh
wget https://drive5.com/uchime/gold.fa --no-check-certificate

# Make sure usearch61 (this needs to be the name of the executable) is in %PATH%
identify_chimeric_seqs.py -m usearch61 -i V3V4Rep1_quality_filtered_q20/seqs.fna -r gold.fa -o qiime_chimeras1/
identify_chimeric_seqs.py -m usearch61 -i V3V4Rep2_quality_filtered_q20/seqs.fna -r gold.fa -o qiime_chimeras2/
identify_chimeric_seqs.py -m usearch61 -i V3V4Rep3_quality_filtered_q20/seqs.fna -r gold.fa -o qiime_chimeras3/
identify_chimeric_seqs.py -m usearch61 -i V3V4Rep4_quality_filtered_q20/seqs.fna -r gold.fa -o qiime_chimeras4/

mkdir qiime_nochim
filter_fasta.py -f V3V4Rep1_quality_filtered_q20/seqs.fna \
                -o qiime_nochim/V3V4Rep1_nochim.fasta \
                -s qiime_chimeras1/chimeras.txt -n

filter_fasta.py -f V3V4Rep2_quality_filtered_q20/seqs.fna \
                -o qiime_nochim/V3V4Rep2_nochim.fasta \
                -s qiime_chimeras2/chimeras.txt -n

filter_fasta.py -f V3V4Rep3_quality_filtered_q20/seqs.fna \
                -o qiime_nochim/V3V4Rep3_nochim.fasta \
                -s qiime_chimeras3/chimeras.txt -n

filter_fasta.py -f V3V4Rep4_quality_filtered_q20/seqs.fna \
                -o qiime_nochim/V3V4Rep4_nochim.fasta \
                -s qiime_chimeras4/chimeras.txt -n
```

### Pick open reference OTUs

```sh
pick_open_reference_otus.py -i qiime_nochim/V3V4Rep1_nochim.fasta -o V3V4Rep1_openref_otus
pick_open_reference_otus.py -i qiime_nochim/V3V4Rep2_nochim.fasta -o V3V4Rep2_openref_otus
pick_open_reference_otus.py -i qiime_nochim/V3V4Rep3_nochim.fasta -o V3V4Rep3_openref_otus
pick_open_reference_otus.py -i qiime_nochim/V3V4Rep4_nochim.fasta -o V3V4Rep4_openref_otus
```

### Pick de novo OTUs
We show Open reference OTU picking in the main text, but de novo OTU picking results look very similar.

```sh
pick_de_novo_otus.py -i qiime_nochim/V3V4Rep1_nochim.fasta -o V3V4Rep1_denovo_otus
pick_de_novo_otus.py -i qiime_nochim/V3V4Rep2_nochim.fasta -o V3V4Rep2_denovo_otus
pick_de_novo_otus.py -i qiime_nochim/V3V4Rep3_nochim.fasta -o V3V4Rep3_denovo_otus
pick_de_novo_otus.py -i qiime_nochim/V3V4Rep4_nochim.fasta -o V3V4Rep4_denovo_otus
```

### Extract counts for each sample id

For open ref OTUs:

```sh
awk '{ print $1, NF-1 }' V3V4Rep1_openref_otus/final_otu_map_mc2.txt > V3V4Rep1_otu_counts_openref.txt
awk '{ print $1, NF-1 }' V3V4Rep2_openref_otus/final_otu_map_mc2.txt > V3V4Rep2_otu_counts_openref.txt
awk '{ print $1, NF-1 }' V3V4Rep3_openref_otus/final_otu_map_mc2.txt > V3V4Rep3_otu_counts_openref.txt
awk '{ print $1, NF-1 }' V3V4Rep4_openref_otus/final_otu_map_mc2.txt > V3V4Rep4_otu_counts_openref.txt
```

For de novo OTUs:

```sh
awk '{ print $1, NF-1 }' V3V4Rep1_denovo_otus/uclust_picked_otus/V3V4Rep1_nochim_otus.txt > V3V4Rep1_denovo_otus/otu_counts.txt
awk '{ print $1, NF-1 }' V3V4Rep2_denovo_otus/uclust_picked_otus/V3V4Rep2_nochim_otus.txt > V3V4Rep2_denovo_otus/otu_counts.txt
awk '{ print $1, NF-1 }' V3V4Rep3_denovo_otus/uclust_picked_otus/V3V4Rep3_nochim_otus.txt > V3V4Rep3_denovo_otus/otu_counts.txt
awk '{ print $1, NF-1 }' V3V4Rep4_denovo_otus/uclust_picked_otus/V3V4Rep4_nochim_otus.txt > V3V4Rep4_denovo_otus/otu_counts.txt
```

## DIABIMMUNE: DADA2 SILVA pipeline

Used latest (at the time of writing) DADA2 v1.8, analysis tutorial followed at https://benjjneb.github.io/dada2/tutorial.html (accessed May 14th, 2018). As per instructions, we downloaded `silva_nr_v132_train_set.fa.gz` and `silva_species_assignment_v132.fa.gz` files for taxonomic classification. The following is the R code used to produce results in the main text:

```R
library(dada2)
library(data.table)

main_path = '/data1/igor/diabimmune'
fq_fwd = sort(dir(file.path(main_path, '16s_fastq_all'), 'R1', full.names=T))
fq_rev = sort(dir(file.path(main_path, '16s_fastq_all'), 'R2', full.names=T))

sample.names = sapply(strsplit(basename(fq_fwd), "_"), `[`, 1)

# Use this to guesstimte truncLen
plotQualityProfile(fq_fwd[1:2])
plotQualityProfile(fq_rev[1:2])

# Filter and trim
filtFs = file.path(main_path, 'dada2', 'filtered', paste0(sample.names, "_F_filt.fastq.gz"))
filtRs = file.path(main_path, 'dada2', 'filtered', paste0(sample.names, "_R_filt.fastq.gz"))
filtOut = filterAndTrim(fq_fwd, filtFs, fq_rev, filtRs, truncLen=c(150,150),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

# Learn errors
errF = learnErrors(filtFs, multithread=TRUE)
errR = learnErrors(filtRs, multithread=TRUE)

# Iterate 1 by 1 sample since there are too many to load all into memory
dadaFs = list()
dadaRs = list()
mergers = list()
for (i in 1:length(sample.names)) {
  # Forward read
  cat('Processing sample ', i , ' out of ', length(sample.names), '\n')
  derepFs = list(derepFastq(filtFs[i], verbose=F))
  names(derepFs) = sample.names[i]
  dadaFs[[i]] = dada(derepFs[[1]], err=errF, multithread=TRUE)
  # Reverse read
  derepRs = list(derepFastq(filtRs[i], verbose=F))
  names(derepRs) = sample.names[i]
  dadaRs[[i]] = dada(derepRs[[1]], err=errR, multithread=TRUE)
  # Merge
  mergers[[i]] = mergePairs(dadaFs[[i]], derepFs[[1]], dadaRs[[i]], derepRs[[1]], verbose=TRUE)
}

# Save DADA results
saveRDS(dadaFs, file.path(main_path, 'dada2', 'dadaFs'))
saveRDS(dadaRs, file.path(main_path, 'dada2', 'dadaRs'))
saveRDS(mergers, file.path(main_path, 'dada2', 'mergers'))

# Generate sequence count table
seqtab = makeSequenceTable(mergers)
saveRDS(seqtab, file.path(main_path, 'dada2', 'seqtab'))
dim(seqtab)

# Remove chimeras
seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
rownames(seqtab.nochim) = sample.names
saveRDS(seqtab.nochim, file.path(main_path, 'dada2', 'seqtab.nochim'))

# Assign taxonomy
taxa = assignTaxonomy(seqtab.nochim, '/data1/igor/databases/silva_nr_v132_train_set.fa.gz',
                      multithread=TRUE)
taxa = addSpecies(taxa, '/data1/igor/databases/silva_species_assignment_v132.fa.gz',
                  allowMultiple=T)
saveRDS(taxa, file.path(main_path, 'dada2', 'taxa'))

# Convert matrices to data tables
seqtab.dt = melt(as.data.table(seqtab.nochim, keep.rownames=T),
                 variable.name='sequence', value.name='count', id.vars='rn')
seqtab.dt[, seq_id := .GRP, by=sequence]
seqtab.dt = seqtab.dt[count > 0]
seqtab.dt[, sequence := NULL]

taxa.dt = as.data.table(taxa, keep.rownames=T)
taxa.dt = merge(taxa.dt, unique(seqtab.dt[, .(seq_id, sequence)]),
                by.x='rn', by.y='sequence')
taxa.dt[, rn := NULL]
setorderv(taxa.dt, c('seq_id'))

# Write tables
write.table(seqtab.dt, file.path(main_path, 'dada2', 'seqtab.dt.txt'),
            sep='\t', quote=F, row.names=F)

write.table(taxa.dt, file.path(main_path, 'dada2', 'taxa.dt.txt'),
            sep='\t', quote=F, row.names=F)
```


## DIABIMMUNE: QIIME2 GG99 pipeline

### Manifest
The analysis was done using QIIME version 2018.4. First we need to generate a manifest file (metadata) for QIIME2 to use. We generate it from the raw FASTQ files (from the folder 16_fastq_all) using this bash script `generate_manifest`:

```sh
#!/bin/bash
#
# Arguments: $1 absolute filepath to scan for fastq files, e.g. "/data1/igor/diabimmune/16s_fastq_all/"
#            $2 separator for sample id, e.g. "_"
#            $3 forward read substring, e.g. "_R1"
#            $4 reverse reads substring, e.g. "_R2"

echo "sample-id,absolute-filepath,direction"
for f in $(find $1 -name "*.fastq*"); do
        base=${f##*/}
        if [[ $base = *"$3"* ]]; then
                direction="forward"
        fi
        if [[ $base = *"$4"* ]]; then
                direction="reverse"
        fi
        echo "${base%%_*},$f,$direction"
done
```

Now run the script:

```sh
./generate_manifest /data1/igor/diabimmune/16s_fastq_all/ _ _R1 _R2 > manifest.txt
```

Run the entire pipeline:

```sh
# Import FASTQ files into artifact:
# This part takes A LONG time (about a day and a half on our server)
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path manifest.txt \
    --output-path diab_import \
    --source-format PairedEndFastqManifestPhred33

# This is a large dataset and QIIME2 copies stuff over and over into and out
# of the temp folder, which will run out of space on our 50GB boot partition.
# Change the temp folder to the large storage volume.
export TMPDIR="/data1/igor/diabimmune/qiime2/tmp/"

# Merge reads
qiime vsearch join-pairs \
    --i-demultiplexed-seqs diab_import.qza \
    --o-joined-sequences diab_merged.qza

# Quality filter merged reads
qiime quality-filter q-score-joined \
    --i-demux diab_merged.qza \
    --o-filtered-sequences diab_filtered.qza \
    --o-filter-stats diab_filtered_stats.qza

# Run Deblur denoising
qiime deblur denoise-16S \
    --i-demultiplexed-seqs diab_filtered.qza \
    --p-trim-length 251 \
    --o-representative-sequences diab_repseqs.qza \
    --o-table diab_table.qza \
    --p-sample-stats \
    --o-stats diab_deblur_stats.qza

# Export data from artifacts into normal files
qiime tools export diab_table.qza \
    --output-dir export_feature_table
qiime tools export diab_repseqs.qza \
    --output-dir export_feature_table
biom convert -i export_feature_table/feature-table.biom \
        -o diab_feature-table.txt --to-tsv
```

For taxonomic classification, we have to do a bit more work since I couldn’t find precomputed artifact file for the normal V4 region on QIIME2 website. There’s only precomputed files for either full-length 16S sequences or for an unusually short "V4" region (120 nt) at https://docs.qiime2.org/2018.2/data-resources/ (scroll a bit down).
Let’s follow the guide from here https://docs.qiime2.org/2018.2/tutorials/feature-classifier/ to train Naive Bayes classifier on V4 region of 99% OTU Greengenes 13.8 sequences. Download the Greengenes 13.8 99% OTU files (FASTA from rep_set and taxonomy) into greengenes folder and then run this:

```sh
# Import FASTA and txt files as QIIME2 artifacts
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path greengenes/99_otus.fasta \
  --output-path 99_otus.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --source-format HeaderlessTSVTaxonomyFormat \
  --input-path greengenes/99_otu_taxonomy.txt \
  --output-path ref-taxonomy.qza
```


We will also need primer sequences for V4 region. Forward primer 515F 5’-3’ sequence: `GTGCCAGCMGCCGCGGTAA`. Reverse primer 805R 3’-5’ sequence: `GACTACCAGGGTATCTAATC`. Now we can extract reference reads and truncate to 251 nt because that’s what we used in the HiMAP pipeline:

```sh
# Extract V3-V4 regions
qiime feature-classifier extract-reads \
    --i-sequences 99_otus.qza \
    --p-f-primer GTGCCAGCMGCCGCGGTAA \
    --p-r-primer GACTACCAGGGTATCTAATC \
    --p-trunc-len 251 \
    --o-reads 99_otus_refseqs.qza

# Train the classifier
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads 99_otus_refseqs.qza \
    --i-reference-taxonomy ref-taxonomy.qza \
    --o-classifier 99_otus_v4_515f-805r_classifier.qza 
```

Now run this feature classifier to generate taxonomy for the sequences in our data:

```sh
qiime feature-classifier classify-sklearn \
  --i-classifier 99_otus_v4_515f-805r_classifier.qza \
  --i-reads diab_repseqs.qza \
  --o-classification diab_taxonomy.qza

# Export taxonomy to tab-delimited file
qiime tools export diab_taxonomy.qza \
  --output-dir export_taxonomy
```

## DIABIMMUNE: QIIME1 GG97 pipeline

All raw paired-end FASTQ files (*R1* and *R2*) are placed in “16s_fastq_all” folder. Then the following commands were ran. The analysis was done using QIIME version 1.9 (as part of the Anaconda installation), ran natively under Ubuntu Linux 16.04.4 LTS. 

### Merge reads

```sh
mkdir merged
# This will generate a separate folder for each read pair
for f in $(find ../16s_fastq_all -name "*_R1*.fastq*"); do
    base=${f##*/}
    base=${base%%_*}
    join_paired_ends.py -f $f \
                    -r ${f/_R1/_R2} \
                    -o merged/${base}
done

# Extract files from all folders into merged, cleanup unmerged
for f in $(find merged -name "*.join.fastq"); do
    base=${f#*/}
    base=${base%%/*}
    mv $f merged/${base}.fastq
    rm -rf merged/${base}
done
```

### Quality filtering

```sh
mkdir filtered_q20
for f in $(find merged -maxdepth 1 -name "*.fastq"); do
    base=${f%.*} # remove extension
    base=${base#*/}
    echo $f
    echo $base
    split_libraries_fastq.py -i $f \
                                 --sample_ids $base \
                                 -o filtered_q20/${base} -q 19 \
                                 --barcode_type 'not-barcoded' --phred_offset=33
done
```

### Chimera removal

```sh
wget https://drive5.com/uchime/gold.fa --no-check-certificate
mkdir chimeras

for f in $(find filtered_q20 -name "seqs.fna"); do
    base=${f#*/}
    base=${base%/*}
    identify_chimeric_seqs.py -m usearch61 -i $f -r gold.fa -o chimeras/$base
done

mkdir nonchim
for f in $(find filtered_q20 -name "seqs.fna"); do
    base=${f##*/}
    base=${base%/*}
    filter_fasta.py -f $f \
                    -o nonchim/${base}.fasta \
                    -s chimeras/${base}/chimeras.txt -n
done
```


### Pick De Novo and Open Ref OTUs

```sh
mkdir otus_openref
mkdir otus_denovo

for f in $(find nonchim -name "*.fasta"); do
    base=${f##*/}
    base=${base%%.*}
    # Pick OTUs
    pick_open_reference_otus.py -i $f -o otus_openref/$base
    pick_de_novo_otus.py -i $f -o otus_denovo/$base
    # Export tables for R
    awk '{ print $1, NF-1 }' otus_openref/${base}/final_otu_map_mc2.txt > otus_openref/${base}_otus_openref_counts.txt
    awk '{ print $1, NF-1 }' otus_denovo/${base}/uclust_picked_otus/${base}_otus.txt > otus_denovo/${base}_otus_denovo_counts.txt
done
```

### Export tables for R

Run this from qiime1 folder to make a single OTU and single taxonomy table.

For De Novo OTU picking:

```sh
out="qiime1_diabimmune_otu_table.txt"
echo -e "sample_id\totu_id\tcount" > $out
for f in $(find otus_denovo -name "*_otus_denovo_counts.txt"); do
    sampleid=${f##*/}
    sampleid=${sampleid%%_*}
    awk -F" " -v sid=$sampleid '{ print sid"\t"sid"-"$1"\t"$2 }' $f >> $out
done

out="qiime1_diabimmune_tax_table.txt"
echo -e "sample_id\totu_id\ttaxonomy" > $out
for f in $(find otus_denovo -name "*_rep_set_tax_assignments.txt"); do
    sampleid=${f##*/}
    sampleid=${sampleid%%_*}
    awk -F "\t" -v sid=$sampleid '{ print sid"\t"sid"-"$1"\t"$2 }' $f >> $out
done
```

For Open Reference OTU picking:

```sh
out="qiime1_diabimmune_otu_table_openref.txt"
echo -e "sample_id\totu_id\tcount" > $out
for f in $(find otus_openref -name "*_otus_openref_counts.txt"); do
    sampleid=${f##*/}
    sampleid=${sampleid%%_*}
    awk -F" " -v sid=$sampleid '{
        # if (substr($1,1,3) == "New") { print sid"\t"sid"-"$1"\t"$2 }
        # else { print sid"\t"$1"\t"$2 }
        print sid"\t"sid"-"$1"\t"$2
    }' $f >> $out
done

out="qiime1_diabimmune_tax_table_openref.txt"
echo -e "sample_id\totu_id\ttaxonomy" > $out
for f in $(find otus_openref -name "rep_set_tax_assignments.txt"); do
    sampleid=${f#*/}
    sampleid=${sampleid%%/*};
        awk -F "\t" -v sid=$sampleid '{ 
        # if (substr($1,1,3) == "New") { print sid"\t"sid"-"$1"\t"$2 }
        # else { print sid"\t"$1"\t"$2 }
        print sid"\t"sid"-"$1"\t"$2
    }' $f >> $out
done
```