# Help
Args <- commandArgs(T)
if (length(Args)==0 | Args[1]=="-h") {
  cat("Phylogenetic topology analysis

Usage:
  Rscript topology_analysis.R <Reference> <Wild relative> <Params> <Chromosome>

Options:
  <Reference>   The reference accession name
  <Wild relative>   The wild relative name
  <Params>    Path of parameter file
  <Chromosome>(optional)   The chromosomes to be analysed, and seperated by comma. If not named, all chromosomes will be analysed.
  -h   Help

Example:
  Rscript src/topology_analysis.R Y601 Y413 ./data/params LG1,LG2")
} else {
  options(warn=-1)
  library(ape)
  library(phangorn)
  # library(pegas)
  library(plyr)
  options(warn=1)
  
  ##############################################################
  ################   Tree topology analysis  ###################
  ##############################################################
  params_path <- Args[3]
  params <- read.delim(params_path, header = F, as.is = T,sep=" ")
  Topology_analysis <- params$V2[grep("Topology",params$V1)]
  script_path <- params$V2[grep("script_path",params$V1)]
  # if (Topology_analysis) {
    source(paste0(script_path,"/functions.R"))
    # Read commandArgs and config
    message("1. Read commandArgs and config")
    ref <- Args[1]
    wr <- Args[2]
    params_path <- Args[3]
    GenomeInfo <- read.delim(params$V2[grep("Genome_info_path",params$V1)], header =T, as.is = T)
    chrs <- GenomeInfo$Chr
    if (length(Args)>3) {
      chrs <- unlist(strsplit(Args[4],","))
    }
    message("Topology analysis is performed for the following chromosomes: \n")
    chrs
    # Topology analysis
    message("2. Topology analysis")
    statistic_all <- NULL
    ploid_tree_all <- NULL
    for (j in 1:length(chrs)) {
      cat(paste0(chrs[j]),"\n")
      statistic_chr <- topology_analysis_iqtree(chrs[j])
      statistic_chr <- topology_analysis_upgma(chrs[j])
      statistic_chr$Chr <- chrs[j]
      statistic_chr$WR <- wr
      statistic_chr <- statistic_chr[order(statistic_chr$Ploidy),]
      statistic_all <- rbind(statistic_chr,statistic_all)
      statistic_chr <- statistic_chr[statistic_chr$Monophy==statistic_chr$UMonophy,]
      ploid_tree <- summarise_ploid(chrs[j],statistic_chr)
      ploid_tree_all <- rbind(ploid_tree_all,ploid_tree)
    }
    write.table(statistic_all,paste(ref,wr,"tree_topology_info.txt",sep="_"),row.names = F,quote=F, sep="\t")
    write.table(as.matrix(ploid_tree_all),paste(ref,wr,"topology_info_ploidy_based_summary.txt",sep="_"),row.names = F,quote=F, sep="\t")
    
    message("Topological result of each syntenic haplotype block")
    head(statistic_chr[,c(1,14,18:21)])
    message(paste0("The result is saved as: ", paste(ref,wr,"tree_topology_info.txt",sep="_")))
    message("Statistic summary of LG1")
    ploid_tree
    message(paste0("The result is saved as: ", paste(ref,wr,"topology_info_ploidy_based_summary.txt",sep="_")))
    message("done")
  }
# }

