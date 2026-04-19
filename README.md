# Practical-RNA-seq-analysis
RNA-seq is a high-throughput sequencing technique used to quantify and compare transcript abundance, enabling the analysis of gene expression and transcriptomic changes across biological conditions.

# What is RNA-seq?

**RNA-seq (RNA sequencing)** is a method used to study **RNA molecules** in a sample.

It tells us:
* which genes are being expressed
* how much they are expressed
* which genes are upregulated or downregulated
* whether different conditions affect gene expression

It can be said that **DNA is the blueprint**, but **RNA shows which parts of the blueprint are active**.
Inside a cell:
**DNA → RNA → Protein**

By RNA-seq analysis we can measure the **RNA level**, so it helps us understand gene activity.
For example:
* healthy tissue vs diseased tissue
* treated vs untreated cells
* mutant vs wild type
* one time point vs another

# Types of RNA studied

RNA-seq can target different RNA types:

* **mRNA**: protein-coding genes
* **rRNA**: ribosomal RNA, usually removed because it is too abundant
* **tRNA**
* **miRNA / small RNA** in special protocols
* total RNA in broader experiments

Most common gene expression studies focus on **mRNA**.

# Main goal of RNA-seq

* quantify gene expression
* Compare the expression between groups
* identify differentially expressed genes
* discover new transcripts or splice variants
* study pathway changes

# General experimental concept

The experiment usually follows this logic:

1. collect biological sample
2. extract total RNA
3. check RNA quality
4. prepare sequencing library
5. sequence the cDNA
6. analyze reads with bioinformatics
7. interpret biological meaning

# Why RNA is converted to cDNA

RNA is less stable.So during library preparation, RNA is converted into **complementary DNA (cDNA)**, because DNA is more stable and suitable for sequencing.

# Important concepts in RNA-seq

## 1. Gene expression

Gene expression means how actively a gene is producing RNA.

If a gene has many reads mapped to it, it usually means that gene is highly expressed.

## 2. Reads

The sequencer produces many short sequences called **reads**.

These reads come from RNA fragments in the sample.

## 3. Read counts

After mapping reads to a genome or transcriptome, we count how many reads belong to each gene.

This gives a **count matrix**, where:

* rows = genes
* columns = samples
* values = counts

## 4. Transcriptome

The **transcriptome** is the complete set of RNA molecules present in a cell or sample at a given time.

It changes depending on condition, tissue, and environment.

# RNA-seq workflow theory

## A. RNA extraction

RNA is extracted from cells or tissues. This must be high quality because degraded RNA gives poor results.

Important quality terms:

* **concentration**
* **purity**
* **integrity**

## B. RNA quality

RNA quality is very important.

Common measurements:

* **A260/A280** for purity
* **RIN (RNA Integrity Number)** for integrity

High-quality RNA is preferred for reliable sequencing.

## C. Library preparation

The RNA is prepared for sequencing.

This often includes:

* rRNA depletion or poly-A selection
* fragmentation
* cDNA synthesis
* adapter ligation
* PCR amplification

### Poly-A selection

Used mainly for eukaryotic mRNA.
It enriches mRNA because mRNA has a poly-A tail.

### rRNA depletion

Removes ribosomal RNA and keeps more RNA types.
Useful when studying non-coding RNA or degraded samples.

# Sequencing strategy

## Single-end vs paired-end

### Single-end

Each fragment is read from one side only.

* cheaper
* simpler
* less information

### Paired-end

Each fragment is read from both ends.

* better alignment
* better transcript discovery
* better splice junction detection

## Read length

Common read lengths are 50 bp, 75 bp, 100 bp, 150 bp.

Longer reads often improve mapping and transcript analysis.

# Raw data format

The raw sequencing data usually comes in **FASTQ** format.

FASTQ contains:

* sequence ID
* nucleotide sequence
* quality score

The quality score indicates confidence in each base.

# Quality control theory

Before analysis, raw reads are checked for:

* low-quality bases
* adapter contamination
* GC bias
* overrepresented sequences
* sequence duplication

This step ensures the data are reliable.

# Trimming

Sometimes reads need trimming to remove:

* adapter sequences
* poor-quality ends
* technical artifacts

Not all datasets need aggressive trimming, but it depends on quality.

# Alignment or pseudoalignment

After QC, reads are assigned to genes or transcripts.

## Alignment

Reads are mapped to a reference genome or transcriptome.

# Reference genome vs transcriptome

## Reference genome

Full genome sequence of the organism

Used when:

* mapping reads to genomic coordinates
* studying splice junctions
* discovering novel transcripts

## Reference transcriptome

Known transcript sequences

Used when:

* estimating transcript abundance
* faster expression analysis

# Counting reads

Once reads are assigned, gene or transcript abundance is measured.

There are two main levels:

## Gene-level quantification

Counts all reads belonging to a gene

Useful for:

* differential gene expression

## Transcript-level quantification

Counts reads for individual isoforms

Useful for:

* alternative splicing
* isoform expression

# Normalization theory

Raw counts cannot be directly compared because samples differ in:

* sequencing depth
* RNA composition
* library size

So we normalize.

## Why normalization is needed

Example:

* Sample A has 10 million reads
* Sample B has 30 million reads

Even if the biology is same, raw counts will look larger in B.
Normalization corrects this.

## Common normalization ideas

### CPM

**Counts per million**
Adjusts for sequencing depth

### FPKM / RPKM

Adjusts for:

* sequencing depth
* gene length

Used for expression estimation, but not ideal for differential expression testing.

### TPM

**Transcripts per million**
Better for comparing transcript abundance within and between samples.

### DESeq2 size-factor normalization

Widely used for differential expression testing.

# Differential expression theory

This is one of the main goals of RNA-seq.

We compare groups to find genes that significantly change.

Examples:

* control vs treatment
* healthy vs diseased
* before vs after exposure

## Upregulated gene

Expression is higher in one group.

## Downregulated gene

Expression is lower in one group.

# Statistical idea behind differential expression

RNA-seq count data are not normally distributed.
They are discrete counts and often show biological variability.

So methods such as DESeq2 or edgeR use models based on the **negative binomial distribution**.

They estimate:

* mean count
* dispersion
* fold change
* significance

# Important statistical terms

## Fold change

How much a gene changes between conditions.

Example:

* 2-fold increase = expression doubled
* 0.5-fold = expression halved

Usually shown as **log2 fold change**

* log2FC = 1 means 2× increase
* log2FC = -1 means 2× decrease

## p-value

Probability of observing the result by chance

## Adjusted p-value / FDR

Because thousands of genes are tested at once, multiple testing correction is needed.

**FDR** helps reduce false positives.

Usually significant genes are chosen by:

* adjusted p-value < 0.05
* and a log2 fold change cutoff

# Biological replicates

These are very important.

A biological replicate means an independent sample from the same condition.

Why needed?

Because biological systems vary naturally.
Replicates help distinguish real biological change from random variation.

Without enough replicates, differential expression becomes weak or unreliable.

# Batch effect

A **batch effect** is unwanted variation caused by technical differences, such as:

* sequencing on different days
* different reagent lots
* different operators
* different machines

Batch effects can hide true biology or create false differences.

# Exploratory data analysis

Before testing differential expression, we check sample relationships.

## PCA

Principal component analysis reduces complex data into major variation patterns.

It helps answer:

* do samples cluster by treatment?
* is any sample an outlier?
* is there batch effect?

## Hierarchical clustering / heatmap

Shows similarity among samples or genes.

# Volcano plot theory

A volcano plot combines:

* x-axis = log2 fold change
* y-axis = statistical significance

It helps quickly identify genes with:

* strong fold change
* strong statistical support

# Heatmap theory

Heatmaps display expression patterns across samples.

Often used for:

* top differentially expressed genes
* pathway-related genes
* clustering samples

# Functional interpretation

After identifying differentially expressed genes, we ask:

**What do these genes mean biologically?**

This is done by enrichment analysis.

## GO analysis

Gene Ontology classifies genes into:

* biological process
* molecular function
* cellular component

## KEGG / pathway analysis

Shows which pathways are enriched, such as:

* immune response
* metabolism
* cell cycle
* apoptosis

# Alternative splicing

In eukaryotes, one gene can produce multiple transcript isoforms.

RNA-seq can detect:

* exon skipping
* intron retention
* alternative splice sites

This gives information beyond total gene expression.

# Strandedness

Some RNA-seq libraries preserve strand information. This tells which DNA strand the RNA came from. Important because overlapping genes on opposite strands can otherwise be confusing.

# Depth of sequencing

**Sequencing depth** means total number of reads per sample. More depth usually means better detection of low-abundance genes. Needed depth depends on study aim:

* basic gene expression: moderate depth
* isoform detection: higher depth
* rare transcript detection: higher depth

# Common outputs of RNA-seq analysis

Typical results include:

* QC report
* alignment summary
* count matrix
* normalized expression matrix
* PCA plot
* heatmap
* volcano plot
* list of differentially expressed genes
* enriched pathways

# Advantages of RNA-seq

Compared with older methods like microarray, RNA-seq has many advantages:

* higher sensitivity
* wider dynamic range
* no need for predesigned probes
* can detect novel transcripts
* can study splice variants
* works for many organisms

# Limitations of RNA-seq

It also has limitations:

* requires good RNA quality
* bioinformatics can be complex
* batch effects can influence results
* expensive compared with some simple methods
* transcript abundance does not always equal protein abundance

# Bulk RNA-seq vs single-cell RNA-seq

## Bulk RNA-seq

RNA from many cells is pooled together.

Gives average expression across all cells.

## Single-cell RNA-seq

RNA is measured at the level of individual cells.

Reveals cell-to-cell heterogeneity.

Bulk RNA-seq is simpler and more common for many experiments.

# RNA-seq in bacteria vs eukaryotes

There are some differences.

## Eukaryotic RNA-seq

* mRNA often has poly-A tail
* introns and splicing are important
* transcript isoforms are common

## Bacterial RNA-seq

* usually no poly-A tail enrichment
* rRNA depletion is more common
* splicing is rare
* operons may be present

# Important things to remember in interpretation

1. high read count does not always mean biological importance
2. significant genes should be interpreted in biological context
3. technical quality strongly affects results
4. replicates are essential
5. RNA change does not always mean protein change

# Very short summary

RNA-seq is used to measure gene expression by sequencing RNA-derived cDNA.
The main theory is:

* extract RNA
* sequence it
* map or quantify reads
* count reads per gene
* normalize data
* statistically compare groups
* identify biological meaning
