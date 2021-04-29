# Help
Args <- commandArgs(T)
if (length(Args)==0 | Args[1]=="-h") {
  cat("Phylogenetic topology

Usage:
  Rscript haplotype_extraction.R <Reference> <Wild relative> <Params> <Chromosome>

Options:
  <Reference>   The reference accession name
  <Wild relative>   The wild relative name
  <Params>    Path of parameter file
  <Chromosome>   The chromosomes to be analysed
  -h   Help

Example:
  Rscript src/haplotype_extraction.R Y601 Y413 ./data/params LG1")
} else {
  library(stringr)
  options(warn=-1)
  options(stringsAsFactors=F)
  options(scipen=200)
  options(warn=1)
  
  params_path <- Args[3]
  params <- read.delim(params_path, header = F, as.is = T,sep=" ")
  Topology_analysis <- params$V2[grep("Topology",params$V1)]
  script_path <- params$V2[grep("script_path",params$V1)]
  source(paste0(script_path,"/functions.R"))
  
  #########################################################################################
  ############ Haplotype blocks extraction of reference and wild relative  ################
  #########################################################################################
  ref <- Args[1]
  wr <- Args[2]
  chr <- Args[4]
  message("1.1 Haplotype blocks extraction of reference")
  # read in sam file for reference
  message(" Read in sam file for reference")
  sam_ref <- read_sam(ref,chr)
  head(sam_ref[,c(1,2,3,8,6)])
  message(paste0("The number of haplotype blocks for reference is: ",nrow(block_ref)))
  
  # extraction haplotype blocks for reference
  message("Haplotype blocks for reference")
  block_ref <- extract_block(sam_ref)
  head(block_ref)
  write.table(block_ref,paste0(ref,"_",chr,"_WR_blocks.txt"),row.names = F,quote=F, sep="\t")
  message(paste0("The result is saved as: ", paste0(ref,"_",chr,"_ref_blocks.txt")))
  
  message("1.2 Haplotype blocks extraction of wild relative")
  # read in sam file for wild relative
  message(" Read in sam file for wild relative")
  sam_wr <- read_sam(wr,chr)
  head(sam_wr[,c(1,2,3,8,6)])
  message(paste0("The number of haplotype blocks for wild relative is: ",nrow(block_wr)))
  # extraction haplotype blocks for wild relative
  message("Haplotype blocks for wild relative")
  block_wr <- extract_block(sam_wr)
  head(block_wr)
  write.table(block_wr,paste0(wr,"_",chr,"_WR_blocks.txt"),row.names = F,quote=F, sep="\t")
  message(paste0("The result is saved as: ", paste0(wr,"_",chr,"_WR_blocks.txt")))
  message("Done")
  
  #########################################################################################
  ##############  Homologous blocks between reference and wild relative  ##################
  #########################################################################################
  
  message("2. Homologous blocks between reference and wild relative")
  
  # extract the syntenic haplotype blocks between reference and wild relative
  message("2.1 Extract the syntenic haplotype blocks between reference and wild relative")
  SyntenicBlock <- extract_homologous_block(block_ref,block_wr)
  head(SyntenicBlock)
  
  # Remove blocks with shorter sequence than length threshhold
  message("2.2 Remove blocks with shorter sequence than length threshhold")
  LenThreshhold <- 3
  SyntenicBlock <- SyntenicBlock[SyntenicBlock$Length >= LenThreshhold,]
  
  # save as table
  write.table(SyntenicBlock,paste(ref,wr,chr,"ref_WR_all_blocks.txt",sep="_"),row.names = F,quote=F, sep="\t")
  message(paste0("The result is saved as: ", paste(ref,wr,chr,"ref_WR_all_blocks.txt",sep="_")))
  
  # reomove duplicated syntenic haplotype blocks
  message("2.2 Reomove duplicated syntenic haplotype blocks")
  ## Step1: deal with duplicated region with different ploidy
  ## sort by RefPloid & WRPloid, retain the first row of each duplicated region
  SyntenicBlock <- rm_duplicated(SyntenicBlock)
  nrow(SyntenicBlock)
  
  ## Step2: deal with total SyntenicBlock
  ## sort by RefPloid & WRPloid & Length, retain the first row of each total SyntenicBlock region
  SyntenicBlock <- cbind(SyntenicBlock,matrix(NA,nrow(SyntenicBlock),2))
  colnames(SyntenicBlock)[(ncol(SyntenicBlock)-1):ncol(SyntenicBlock)] <- c("DeltaSE","DeltaE")
  SyntenicBlock$DeltaSE[1:(nrow(SyntenicBlock)-1)] <- SyntenicBlock$Start[2:nrow(SyntenicBlock)] - SyntenicBlock$End[1:(nrow(SyntenicBlock)-1)]
  SyntenicBlock$DeltaE[1:(nrow(SyntenicBlock)-1)] <- SyntenicBlock$End[2:nrow(SyntenicBlock)] - SyntenicBlock$End[1:(nrow(SyntenicBlock)-1)]
  SyntenicBlock <- rm_total_overlap(SyntenicBlock)
  nrow(SyntenicBlock)

  ## Step3: deal with partial SyntenicBlock
  ## delete the overlap region for one of duplicated blocks
  SyntenicBlock <- sub_partial_overlap(SyntenicBlock)
  nrow(SyntenicBlock)
  SyntenicBlock$Length <- SyntenicBlock$End - SyntenicBlock$Start + 1
  
  # Remove blocks with shorter sequence than length threshhold
  LenThreshhold <- 3
  SyntenicBlock <- SyntenicBlock[SyntenicBlock$Length >= LenThreshhold,]
  nrow(SyntenicBlock)
  
  # save as table
  write.table(SyntenicBlock,paste(ref,wr,chr,"ref_WR_filtered_blocks.txt",sep="_"),row.names = F,quote=F, sep="\t")
  message(paste0("The result is saved as: ", paste(ref,wr,chr,"ref_WR_filtered_blocks.txt",sep="_")))
  
  # update the sam table for reference and wild relative
  sam_ref <- sam_ref[sam_ref$Block %in% SyntenicBlock$RefBlock,]
  sam_wr <- sam_wr[sam_wr$Block %in% SyntenicBlock$WRBlock,]
  
  message("Done")
  
  ##############################################################
  ################   Calculate cigar position  #################
  ##############################################################
  
  message("3. Calculate cigar position")
  ## Calculate cigar position for reference species
  message("3.1 Calculate cigar position for reference species")
  cigar_ref <- extract_cigar(sam_ref)
  
  ## Calculate cigar position for wild relative
  message("3.2 Calculate cigar position for wild relative")
  cigar_wr <- extract_cigar(sam_wr)
  
  message("Done")
  
  ##############################################################
  ####################   Extract sequence  #####################
  ##############################################################
  message("4. Extract sequence")
  
  # Extract sequence of reference species
  message("4.1 Extract sequence of reference species")
  seqs_ref <- extract_sequence_ref(SyntenicBlock,sam_ref,cigar_ref,ref)
  
  # Extract sequence of wild relative
  message("4.2 Extract sequence of  wild relative")
  seqs_wr <- extract_sequence_wr(SyntenicBlock,sam_wr,cigar_wr,wr)
  message(paste0("The number of syntenic haplotype block: "),length(seqs_ref)) 
  # save sequences of syntenic haplotypes as fasta file
  
  SyntenicBlock$SeqName <- NA
  seqs <- list()
  # create directory fa to store fasta file
  if (!file.exists("./fa")) {system("mkdir fa")}
  for (i in 1:nrow(SyntenicBlock)) {
    seqs[[i]] <- rbind(seqs_wr[[i]],seqs_ref[[i]])
    # remove blocks with all identical haplotype sequences
    if (length(unique(seqs[[i]]$Seq)) != 1) {
      SyntenicBlock$SeqName[i] <- paste0(paste(ref,wr,chr,SyntenicBlock$Start[i],SyntenicBlock$End[i],i,sep = "_"),".fa")
      write.table(seqs[[i]][,c("Name","Seq")],file=paste0("./fa/",SyntenicBlock$SeqName[i]),sep="\n",
                  quote = F, row.names = F,col.names = F)
    }
  }
  SyntenicBlock <- SyntenicBlock[!is.na(SyntenicBlock$SeqName),]
  write.table(SyntenicBlock, paste(ref,wr,chr,"ref_WR_final_blocks.txt",sep="_"), quote = F, row.names = F)
  message(paste0("The final number of syntenic haplotype block: "),nrow(SyntenicBlock))
  message(paste0("The result is saved as: ", paste(ref,wr,chr,"ref_WR_final_blocks.txt",sep="_")))
  message("Done\n")
}

