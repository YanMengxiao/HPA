#!/bin/bash
###
### Syntenic haplotype blocks extraction for polyploid and its wild relative
###
### Usage:
###   split_bam.sh <Reference> <Wild relative> <Chromosome> <Params> 
###
### Options:
###   <Reference>   The reference accession name
###   <Wild relative>   The wild relative name
###   <Params>    Path of parameter file
###   <Chromosome>   The chromosome name
###   -h   Help
### 
### Example:
###   src/split_bam.sh Y601 Y413 params LG1

help() {
	awk -F'### ' '/^###/ { print $2 }' "$0"
}

if [[ $# == 0 ]] || [[ "$1" == "-h" ]]; then
	help
	exit 1
fi

RefFile=`grep Reference_path $3 | awk '{print $2}'`
WRFile=`grep Wild_relative_path $3 | awk '{print $2}'`

#samtools index $RefFile
samtools view $RefFile $4 | awk '{if($1!~"ploidy_1") print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}' > $1_$4.sam
#samtools index $WRFile
samtools view $WRFile $4 | awk '{if($1!~"ploidy_1") print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}' > $2_$4.sam
