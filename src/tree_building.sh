#!/bin/bash
###
### High-throughput phylogenetic reconstruction
###
### Usage:
###   tree_building.sh <path> <thread>
###
### Options:
###   <path>   The parameter file for MEGA-CC UPGMA tree
###   <thread>   The thread number
###   -h   Help
### Example:
###   bash tree_building.sh ../../src/infer_UPGMA_nucleotide.mao 64

help() {
	awk -F'### ' '/^###/ { print $2 }' "$0"
}

if [[ $# == 0 ]] || [[ "$1" == "-h" ]]; then
	help
	exit 1
fi

function ceil(){
  floor=`echo "scale=0;$1/1"|bc -l ` # 向下取整
  add=`awk -v num1=$floor -v num2=$1 'BEGIN{print(num1<num2)?"1":"0"}'`
  echo `expr $floor  + $add`
}

alignment
mkdir fas
cd fa    
ls | awk '{print "mafft " $i " > ../fas/" $i ".fas"}' > align.sh 
total=`wc -l align.sh | cut -f 1 -d " "`
thread=$2
num=`echo "sclae=1;$total/$thread" | bc -l`
num=`ceil $num`
split -l $num align.sh
ls x* | awk '{print "bash " $i " &"}' > para_align.sh
echo "wait" >> para_align.sh
bash para_align.sh
rm x*
rm *sh*

# iqtree
cd ../fas
mkdir ../tree
# ls | awk '{print "iqtree -s " $i " -m TEST -safe -pre ../tree/" $i}' > iqtree.sh
ls | awk '{print "mkdir ../tree/" $i "\niqtree -s " $i " -m TEST -safe -pre ../tree/" $i "/" $i}' > iqtree.sh
total=`wc -l iqtree.sh | cut -f 1 -d " "`
thread=$2
num=`echo "sclae=1;$total/$thread" | bc -l`
num=`ceil $num`
split -l $num iqtree.sh
ls x* | awk '{print "bash " $i " &"}' > para_tree.sh
echo "wait" >> para_tree.sh
bash para_tree.sh > iqtree.log

# UPGMA
rm x*
rm *sh*
mkdir ../tree_UPGMA
ls | awk -v val=$1 '{print "megacc -a " val "  -d " $i " -o ../tree_UPGMA -n"}' > megacc.sh
#head  megacc.sh
total=`wc -l megacc.sh | cut -f 1 -d " "`
thread=$2
num=`echo "sclae=1;$total/$thread" | bc -l`
num=`ceil $num`
split -l $num megacc.sh
ls x* | awk '{print "bash " $i " &"}' > para_tree.sh
echo "wait" >> para_tree.sh
bash para_tree.sh > UPGMA.log
rm x*
rm *sh*
cd ..

