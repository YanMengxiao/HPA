HPA (Haplotype-based Phylogenetic Analysis)
===========================================

HPA is developed to resolve the genetic origin of polyploid. It
distinguishes the relationship between polyploid species and its
polyploid relatives based on high-throughput phylogenetic analysis of
the homologous haplotypes. The demo data Y601 is hexaploid sweet potato
and Y413 is its tetraploid relative.

Dependencies
------------

The code is implemented in R 4.0.2 and can be run from linux shell. To
run HPA the following libraries and tools should be available:

samtools 1.10

MAFFT v7.471

megacc 10.1.8

iqtree 1.6.12

R package: ape 5.4-1, phangorn 2.5.5, plyr 1.8.6, stringr 1.4.0

Data preparation
----------------

Illumina raw reads were mapped to reference genome

``` shell
bwa index -p pasi3.fa
# hexaploid sweet potato
bwa mem Y601_R1.fastq.gz Y601_R2.fastq.gz -t 30 |samtools view -@ 5 -Su - | samtools sort -@ 5 - -T 500 -o Y601.bam &
# tetraploid relative
bwa mem Y413_R1.fastq.gz Y413_R2.fastq.gz -t 30 |samtools view -@ 5 -Su - | samtools sort -@ 5 - -T 500 -o Y413.bam &
# The command lines for other accessions are similar, taking two accessions as an example.
```

PCR duplicates were removed

``` shell
# hexaploid sweet potato
samtools rmdup Y601.bam Y601.rmdup.bam &
# tetraploid relative
samtools rmdup Y413.bam Y413.rmdup.bam &
# The command lines for other accessions are similar, taking two accessions as an example.
```

Variant calling

``` shell
fasta_generate_regions.py pasi3.fa.fai 100000 > pasi3.100kb &
# hexaploid sweet potato
freebayes-parallel pasi3.100kb 48 -f pasi3.fa -p 6 -F 0.08 Y601.bam > Y601.rmdup.vcf &
# tetraploid relative
freebayes-parallel pasi3.100kb 48 -f pasi3.fa -p 4 -F 0.08 Y413.bam > Y413.rmdup.vcf &
# The command lines for other accessions are similar, taking two accessions as an example.
```

Haplotyping

``` shell
# hexaploid sweet potato
python ranbow.py hap -mode index -par /project/sweet/liubeiti/Y601/Y601.par 
for i in {0..64}; do python ranbow.py hap -mode hap -par /project/sweet/liubeiti/Y601/Y601.par -processorIndex $i >/project/sweet/liubeiti/Y601/Ranbow/$i.log & done
# tetraploid relative
python ranbow.py hap -mode index -par /project/sweet/liubeiti/Y413/Y413.par 
for i in {0..64}; do python ranbow.py hap -mode hap -par /project/sweet/liubeiti/Y413/Y413.par -processorIndex $i >/project/sweet/liubeiti/Y413/Ranbow/$i.log & done
# The command lines for other accessions are similar, taking two accessions as an example.
```

``` shell
# parameters for Y601
-ploidy 6
-noProcessor    303
-bamFile    /project/sweet/liubeiti/Y601/Y601.bwa.bam
-refFile    /project/sweet/b47/pasi/pasi3.fa
-vcfFile    /project/sweet/liubeiti/Y601/Y601.6.08.bwa.vcf
-selectedScf    /project/sweet/pasi3.100kb
-outputFolderBase   /project/sweet/liubeiti/Y601/Ranbow/
# parameters for Y413 
-ploidy 6
-noProcessor    303
-bamFile    /project/sweet/liubeiti/Y413/Y413.bwa.bam
-refFile    /project/sweet/b47/pasi/pasi3.fa
-vcfFile    /project/sweet/liubeiti/Y413/Y413.6.08.bwa.vcf
-selectedScf    /project/sweet/pasi3.100kb
-outputFolderBase   /project/sweet/liubeiti/Y413/RanbowLG1/  
```

Index bam file
--------------

``` shell
samtools index ./data/Y601.hap.bam
samtools index ./data/Y413.hap.bam
```

Simplify the bam file
---------------------

Usage:

split\_bam.sh *Reference* *Wild relative* *Params* *Chromosome*

Options:

*Reference* The reference accession name

*Wild relative* The wild relative name

*Params* Path of parameter file

*Chromosome* The chromosome name

-h Help  
<br/>

Deal with a single chromosome

``` shell
bash src/split_bam.sh Y601 Y413 ./data/params LG1
```

Batch processing for all chromosomes

``` shell
tail -15 ~/genomics/genome/genome.info | awk '{print $1}' | awk '{print "bash src/split_bam.sh Y601 Y413 ./data/params " $i " &"}'
```

    ## bash src/split_bam.sh Y601 Y413 ./data/params LG1 &
    ##  bash src/split_bam.sh Y601 Y413 ./data/params LG2 &
    ##  bash src/split_bam.sh Y601 Y413 ./data/params LG3 &
    ##  bash src/split_bam.sh Y601 Y413 ./data/params LG4 &
    ##  bash src/split_bam.sh Y601 Y413 ./data/params LG5 &
    ##  bash src/split_bam.sh Y601 Y413 ./data/params LG6 &
    ##  bash src/split_bam.sh Y601 Y413 ./data/params LG7 &
    ##  bash src/split_bam.sh Y601 Y413 ./data/params LG8 &
    ##  bash src/split_bam.sh Y601 Y413 ./data/params LG9 &
    ##  bash src/split_bam.sh Y601 Y413 ./data/params LG10 &
    ##  bash src/split_bam.sh Y601 Y413 ./data/params LG11 &
    ##  bash src/split_bam.sh Y601 Y413 ./data/params LG12 &
    ##  bash src/split_bam.sh Y601 Y413 ./data/params LG13 &
    ##  bash src/split_bam.sh Y601 Y413 ./data/params LG14 &
    ##  bash src/split_bam.sh Y601 Y413 ./data/params LG15 &

The syntenic haplotype blocks extraction
----------------------------------------

Usage:

Rscript haplotype\_extraction.R *Reference* *Wild relative* *Params*
*Chromosome*

Options:

*Reference* The reference accession name

*Wild relative* The wild relative name

*Params* Path of parameter file

*Chromosome*(optional) The chromosomes to be analysed, and seperated by
comma. If not named, all chromosomes will be analysed.

-h Help  
<br/>

Deal with a single chromosome

``` shell
Rscript src/haplotype_extraction.R Y601 Y413 ./data/params LG1
```

Batch processing for all chromosomes

``` shell
tail -15 ~/genomics/genome/genome.info | awk '{print $1}' | awk '{print "Rscript src/haplotype_extraction.R Y601 Y413 ./data/params " $i " &"}'
```

    ## Rscript src/haplotype_extraction.R Y601 Y413 ./data/params LG1 &
    ##  Rscript src/haplotype_extraction.R Y601 Y413 ./data/params LG2 &
    ##  Rscript src/haplotype_extraction.R Y601 Y413 ./data/params LG3 &
    ##  Rscript src/haplotype_extraction.R Y601 Y413 ./data/params LG4 &
    ##  Rscript src/haplotype_extraction.R Y601 Y413 ./data/params LG5 &
    ##  Rscript src/haplotype_extraction.R Y601 Y413 ./data/params LG6 &
    ##  Rscript src/haplotype_extraction.R Y601 Y413 ./data/params LG7 &
    ##  Rscript src/haplotype_extraction.R Y601 Y413 ./data/params LG8 &
    ##  Rscript src/haplotype_extraction.R Y601 Y413 ./data/params LG9 &
    ##  Rscript src/haplotype_extraction.R Y601 Y413 ./data/params LG10 &
    ##  Rscript src/haplotype_extraction.R Y601 Y413 ./data/params LG11 &
    ##  Rscript src/haplotype_extraction.R Y601 Y413 ./data/params LG12 &
    ##  Rscript src/haplotype_extraction.R Y601 Y413 ./data/params LG13 &
    ##  Rscript src/haplotype_extraction.R Y601 Y413 ./data/params LG14 &
    ##  Rscript src/haplotype_extraction.R Y601 Y413 ./data/params LG15 &

High-throughput phylogenetic reconstruction
-------------------------------------------

Usage:

tree\_building.sh *Path* *Thread*

Options:

*Path* The parameter file for MEGA-CC UPGMA tree

*Thread* The thread number

-h Help

``` shell
bash src/tree_building.sh ../../src/infer_UPGMA_nucleotide.mao 128
```

Topological analysis
--------------------

Usage:

Rscript topology\_analysis.R *Reference* *Wild relative* *Params*
*Chromosome*

Options:

*Reference* The reference accession name

*Wild relative* The wild relative name

*Params* Path of parameter file

*Chromosome*(optional) The chromosomes to be analysed, and seperated by
comma. If not named, all chromosomes will be analysed.

-h Help  
<br/>

Analyzing all chromosomes

``` shell
Rscript src/topology_analysis.R Y601 Y413 ./data/params
```

Analyzing specific chromosomes

``` shell
Rscript src/topology_analysis.R Y601 Y413 ./data/params LG1,LG2
```

Gene conversion analysis
------------------------

Usage:

Rscript gene\_conversion\_analysis.R *Reference* *Wild relative*
*Params* *Chromosome*

Options:

*Reference* The reference accession name

*Wild relative* The wild relative name

*Params* Path of parameter file

*Chromosome*(optional) The chromosomes to be analysed, and seperated by
comma. If not named, all chromosomes will be analysed.

-h Help  
<br/>

Analyzing all chromosomes

``` shell
Rscript src/gene_conversion_analysis.R Y601 Y413 ./data/params
```

Analyzing specific chromosomes

``` shell
Rscript src/gene_conversion_analysis.R Y601 Y413 ./data/params LG1,LG2
```
