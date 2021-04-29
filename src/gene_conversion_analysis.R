# Help
Args <- commandArgs(T)
if (length(Args)==0 | Args[1]=="-h") {
  cat("Gene conversion analysis

Usage:
  Rscript gene_conversion_analysis.R <Reference> <Wild relative> <Params> <Chromosome>

Options:
  <Reference>   The reference accession name
  <Wild relative>   The wild relative name
  <Params>    Path of parameter file
  <Chromosome>(optional)   The chromosomes to be analysed, and seperated by comma. If not named, all chromosomes will be analysed.
  -h   Help

Example:
  Rscript src/gene_conversion_analysis.R Y601 Y413 ./data/params LG1,LG2")
} else {
  options(warn=-1)
  library(ape)
  library( phangorn)
  options(warn=1)
  
  ##############################################################
  ################   Gene conversion analysis  #################
  ##############################################################
  params_path <- Args[3]
  params <- read.delim(params_path, header = F, as.is = T,sep=" ")
  # Gene_conversion <- params$V2[grep("Gene_conversion",params$V1)]
  # if (Gene_conversion) {
    source("~/genomics/origin/HPA/Y601/test/src/functions.R") 
    # Read commandArgs and config
    message("1. Read commandArgs and config")
    ref <- Args[1]
    wr <- Args[2]
    # GenomeInfo <- read.delim(params$V2[grep("Genome_info_path",params$V1)], header =T, as.is = T)
    # chrs <- GenomeInfo$Chr
        message("done")
    # Gene conversion analysis
    message("2. Gene conversion analysis")
  # }
  ### gene conversion UPGMA
  RefPloidy <- as.numeric(params$V2[grep("Reference_ploid",params$V1)])
  WRPloidy <- as.numeric(params$V2[grep("Wild_relative_ploid",params$V1)])
  # Extract blocks of specific ploid
  message("2.1 Extract blocks of specific ploid")
  statistic_all <- read.delim(paste(ref,wr,"tree_topology_info.txt",sep="_"), as.is = T, header=T)
  statistic_all <- statistic_all[statistic_all$RefPloidy==RefPloidy & statistic_all$WRPloidy==WRPloidy,]
  statistic_all <- statistic_all[statistic_all$Monophy == statistic_all$UMonophy,]
  if (length(Args)>3) {
    chrs <- unlist(strsplit(Args[4],","))
    statistic_all <- statistic_all[statistic_all$Chr %in% chrs,]
  }
  # Extract blocks within gene region
  message("2.2 Extract blocks within gene region")
  path <- params$V2[grep("Genome_annotation_path",params$V1)]
  rows <- extract_gene_region(path,statistic_all)
  statistic_all <- statistic_all[unique(rows),]
  
  # Gene conversion identified by tree topology
  message("2.2 Gene conversion identified by tree topology")
  topology_upgma <- gene_conversion_upgma(TreeFolder="tree_UPGMA",TreeSuffix="_.*nwk")
  head(topology_upgma)
  write.table(topology_upgma, paste0(wr,"_gene_conversion.txt"), quote = F, row.names = F)
  message(paste0("The result is saved as: ", paste(ref,wr,"topology_info_ploidy_based_summary.txt",sep="_")))
  
  message("The ratio from B2 to B1:")
  cat(length(which(topology_upgma$Type=="tetra2di"))/nrow(topology_upgma),"\n")
  message("The ratio from B1 to B2:")
  cat(length(which(topology_upgma$Type=="di2tetra"))/nrow(topology_upgma),"\n")
  message("How many trees are used for gene conversion analyse?")
  cat(nrow(topology_upgma),"\n")
  
  ### gene conversion IQTREE
  if (F) {
    topology_iqtree <- gene_conversion_upgma(TreeFolder="tree",TreeFolderPattern=".fa.fas")
    write.table(topology_iqtree, paste0(wr,"_gene_conversion_iqtree.txt"), quote = F, row.names = F)
    message("The ratio from B2 to B1:")
    length(which(topology_upgma$Type=="tetra2di"))/nrow(topology_upgma)
    message("The ratio from B1 to B2:")
    length(which(topology_upgma$Type=="di2tetra"))/nrow(topology_upgma)
    message("How many trees are used for gene conversion analyse?")
    nrow(topology_upgma)
  }
}



