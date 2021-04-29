
read_sam <- function(sample,chr) {
  sam_wr <- read.delim(paste0(paste(sample,chr,sep = "_"),".sam"),header = F,as.is= T)
  colnames(sam_wr) <- c("Hap","Ref","StartPos","Cigar","Seq")
  sam_wr$Block <- gsub("_Hap_[0-9]","",sam_wr$Hap)
  # calculate the end position
  sam_wr <- cbind(sam_wr,matrix(NA,nrow(sam_wr),1))
  colnames(sam_wr)[ncol(sam_wr)] <- c("M")
  for (i in 1:nrow(sam_wr)) {
    sam_wr$M[i] <- sum(as.numeric(gsub("M","",unlist(str_extract_all(sam_wr$Cigar[i],"[0-9]*M")))))
  }
  sam_wr$EndPos <- sam_wr$StartPos + sam_wr$M - 1
  sam_wr
}

extract_block <- function(sam) {
  blocks <- unique(sam$Block)
  SampleBlock <- data.frame(matrix(NA,length(blocks),4),as.is=T)
  colnames(SampleBlock) <- c("Ploidy","Block","Ref","Start","End")
  SampleBlock$Block <- blocks
  SampleBlock$Ploidy <- as.numeric(gsub("ploidy_","",str_extract(SampleBlock$Block,"ploidy_[1-6]")))
  for (i in 1:length(blocks)) {
    rows <- which(sam$Block==blocks[i])
    # remove the blocks which all haplotypes not overlap
    if (length(rows) == SampleBlock$Ploidy[i]) {
      SampleBlock$Start[i] <- max(sam$StartPos[rows])
      SampleBlock$End[i] <- min(sam$EndPos[rows])
      SampleBlock$Ref[i] <- sam$Ref[rows[1]]
    }
  }
  SampleBlock$Length <- SampleBlock$End - SampleBlock$Start + 1
  SampleBlock <- SampleBlock[SampleBlock$Length>0 & !is.na(SampleBlock$Length),]
  SampleBlock
}

extract_homologous_block <- function (block_ref,block_wr) {
  block_all <- cbind(block_wr[,c("Ploidy","Block","Ref","Start","End")],matrix(NA,nrow(block_wr),6))
  colnames(block_all) <- c("WRPloidy","WRBlock","Ref","WRStart","WREnd","RefPloidy","RefBlock","RefStart","RefEnd","Start","End")
  # for more than one overlap region
  block_add <- data.frame(matrix(NA,1,11))
  colnames(block_add) <- colnames(block_all)
  num <- nrow(block_wr)
  # num stand for different ploid block of c4
  for (i in 1:num) {
    # find out which row of block_ref has overlap region with block_wr[i]
    rows <- which((block_all$WRStart[i]>=block_ref$Start & block_all$WRStart[i]<block_ref$End) | 
                    (block_all$WREnd[i]>block_ref$Start & block_all$WREnd[i]<=block_ref$End) | 
                    (block_all$WRStart[i]<=block_ref$Start & block_all$WREnd[i]>=block_ref$Start & block_all$WRStart[i]<=block_ref$End & block_all$WREnd[i]>=block_ref$End))
    # block_wr[i] has one verlap region with c6
    if (length(rows)==1) {
      # add overlaped c6 info
      # block_all$RefPloidy[i] <- block_ref$Ploidy[rows]
      block_all$RefBlock[i] <- block_ref$Block[rows]
      block_all$RefStart[i] <- block_ref$Start[rows]
      block_all$RefEnd[i] <- block_ref$End[rows]
      # calculate overlap region
      block_all$Start[i] <- max(block_all$WRStart[i],block_all$RefStart[i])
      block_all$End[i] <- min(block_all$WREnd[i],block_all$RefEnd[i])
    }
    # block_wr[i] has more than one verlap region with c6
    if (length(rows)>1) {
      # no <- c(no,i)
      # for the first overlap region add in table block_all
      # block_all$RefPloidy[i] <- block_ref$Ploidy[rows[1]]
      block_all$RefBlock[i] <- block_ref$Block[rows[1]]
      block_all$RefStart[i] <- block_ref$Start[rows[1]]
      block_all$RefEnd[i] <- block_ref$End[rows[1]]
      block_all$Start[i] <- max(block_all$WRStart[i],block_all$RefStart[i])
      block_all$End[i] <- min(block_all$WREnd[i],block_all$RefEnd[i])
      ## for the second and more overlap region add in table block_add
      for (j in 2:length(rows)) {
        block_add$WRBlock[1] <- block_all$WRBlock[i]
        block_add$Ref[1] <- block_all$Ref[i]
        block_add$WRStart[1] <- block_all$WRStart[i]
        block_add$WREnd[1] <- block_all$WREnd[i]
        # block_all$RefPloidy[i] <- block_ref$Ploidy[rows[j]]
        block_add$RefBlock[1] <- block_ref$Block[rows[j]]
        block_add$RefStart[1] <- block_ref$Start[rows[j]]
        block_add$RefEnd[1] <- block_ref$End[rows[j]]
        block_add$Start <- max(block_add$WRStart,block_add$RefStart)
        block_add$End <- min(block_add$WREnd,block_add$RefEnd)
        ## combine table block_all & block_add: add each overlap to block_all
        block_all <- rbind(block_all,block_add)
      }
    }
  }
  # Step2: remove empty record which c4 has no overlap region with c6
  block_all <- block_all[!is.na(block_all$RefBlock),]
  # Step3: rank block_all by decreasing order of RefPloidy,WRPloidy and block_all$Length
  block_all$Length <- block_all$End - block_all$Start + 1
  block_all$RefPloidy <- as.numeric(gsub("ploidy_","",str_extract(block_all$RefBlock,"ploidy_[1-6]")))
  block_all$WRPloidy <- as.numeric(gsub("ploidy_","",str_extract(block_all$WRBlock,"ploidy_[1-6]")))
  block_all <- block_all[order(block_all$RefPloidy,block_all$WRPloidy,block_all$Length, decreasing = T),]
  block_all
}

rm_duplicated <- function( block_all) {
  block_all <-  block_all[order( block_all$Start,- block_all$RefPloid,- block_all$WRPloid,- block_all$Length),]
  block_all$Ploidy <- paste0( block_all$RefPloidy,"_", block_all$WRPloidy)
  rownames( block_all) <- 1:nrow( block_all)
  block_all$Pos <- paste0( block_all$Start,"_", block_all$End)
  dulppos <- unique( block_all$Pos[duplicated( block_all$Pos)])
  nrow( block_all)
  rmrow <- NULL
  for (i in dulppos) {
    sub <-  block_all[ block_all$Pos == i,]
    rmrow <- c(rmrow,rownames(sub)[-1])
  }
  block_all <-  block_all[!rownames( block_all) %in% rmrow,]
  rownames( block_all) <- 1:nrow( block_all)
  block_all
}


rm_total_overlap <- function( block_all) {
  rmrow <- NULL
  # count <- 1
  rows <- which( block_all$DeltaSE[1:(nrow( block_all)-1)] < 0 &  block_all$DeltaE[1:(nrow( block_all)-1)] <= 0)
  if (length(rows) != 0) {
    for(i in rows) {
      sub <-  block_all[i:(i+1),]
      sub <- sub[order(-sub$RefPloidy,-sub$WRPloidy,-sub$Length),]
      rmrow <- c(rmrow, rownames(sub)[2])
    }
    block_all <-  block_all[!rownames( block_all) %in% rmrow,]
    rownames( block_all) <- 1:nrow( block_all)
    block_all$DeltaSE[1:(nrow( block_all)-1)] <-  block_all$Start[2:nrow( block_all)] -  block_all$End[1:(nrow( block_all)-1)]
    block_all$DeltaE[1:(nrow( block_all)-1)] <-  block_all$End[2:nrow( block_all)] -  block_all$End[1:(nrow( block_all)-1)]
    # count <- count + 1
    block_all <- rm_total_overlap( block_all)
  }
  block_all 
}

sub_partial_overlap <- function( block_all) {
  rows <- which( block_all$DeltaSE[1:(nrow( block_all)-1)] < 0 &  block_all$DeltaE[1:(nrow( block_all)-1)] >= 0)
  for (i in rows) {
    block_all$Start[i+1] <-  block_all$Start[i+1] -  block_all$DeltaSE[i] +1
  }
  block_all$DeltaSE[1:(nrow( block_all)-1)] <-  block_all$Start[2:nrow( block_all)] -  block_all$End[1:(nrow( block_all)-1)]
  block_all$DeltaE[1:(nrow( block_all)-1)] <-  block_all$End[2:nrow( block_all)] -  block_all$End[1:(nrow( block_all)-1)]
  block_all
}

extract_cigar <- function(sam) {
  options(warn=-1)
  # Step1: replace D/I by DA/IA, using A as the marker to seperate cigar event
  sam$Cigar <- gsub("D","DA",sam$Cigar)
  sam$Cigar <- gsub("I","IA",sam$Cigar)
  # step2: for haps with D/I, calculate the extraction position
  dein_ref <- which(str_detect(sam$Cigar,"D|I"))
  cigar <- list()
  if (length(dein_ref) > 0) {
    for (i in dein_ref) {
      # split cigar, get all M/I/D events
      cigars <- strsplit(sam$Cigar[i],"A")[[1]]
      # calculate D/I events, minus the 1st M event
      rows <- length(cigars)-1
      # store D/I events for every haplotype; if no event, table should be NULL
      cigar[[i]] <- data.frame(matrix(NA,rows,0))
      # first row is M, not INS/DEL, remove 1st row
      cigar[[i]]$Cigar <- cigars[-(rows+1)]
      cigar[[i]]$IDpos <- as.numeric(gsub("M[0-9]*[D,I]","",cigar[[i]]$Cigar)) 
      # if cigar pattern without M (e.g.,100M10D11I), NA is introduced and get warning message. 
      # Without M, no M length needed to add. Therefore, Spos for NA is 0. The warning is fixed.
      cigar[[i]]$IDpos[is.na(cigar[[i]]$IDpos)] <- 0
      for (x in 1:rows) {
        # for pattern like 166M4I, extract event length
        if (str_detect(cigar[[i]]$Cigar[x],"M")) {
          if (str_detect(cigar[[i]]$Cigar[x],"I"))
            cigar[[i]]$Length[x] <- as.numeric(gsub("[0-9]*M([0-9]*)[D,I]","\\1",cigar[[i]]$Cigar[x]))
          else 
            cigar[[i]]$Length[x] <- -as.numeric(gsub("[0-9]*M([0-9]*)[D,I]","\\1",cigar[[i]]$Cigar[x]))
        }
        # for pattern like 4I, extract event length
        else {
          if (str_detect(cigar[[i]]$Cigar[x],"I"))
            cigar[[i]]$Length[x] <- as.numeric(gsub("[D,I]","",cigar[[i]]$Cigar[x]))
          else 
            cigar[[i]]$Length[x] <- -as.numeric(gsub("[D,I]","",cigar[[i]]$Cigar[x]))
        } 
      }
      cigar[[i]]$CLength <- cigar[[i]]$Length
      if (rows>1) {
        for (j in 2:rows) {
          # calculate event start position
          cigar[[i]]$IDpos[j] <- cigar[[i]]$IDpos[j-1]+cigar[[i]]$IDpos[j]
          # calculate cumulative event length
          cigar[[i]]$CLength[j] <- cigar[[i]]$CLength[j-1]+cigar[[i]]$CLength[j]
        }
      }
    }
    names(cigar) <- sam$Hap[1:dein_ref[length(dein_ref)]]
  }
  options(warn=1)
  cigar
}

extract_sequence_ref <- function(SyntenicBlock,sam_ref,cigar_ref,ref) {
  seqs_ref <- list()
  for (i in 1:nrow(SyntenicBlock)) {
    RefPloidy <- SyntenicBlock$RefPloidy[i]
    seqs_ref[[i]] <- data.frame(matrix(NA,RefPloidy,0))
    seqs_ref[[i]]$WRBlock <- rep(SyntenicBlock$WRBlock[i],RefPloidy)
    seqs_ref[[i]]$RefBlock <- rep(SyntenicBlock$RefBlock[i],RefPloidy)
    seqs_ref[[i]]$Hap <- sam_ref$Hap[sam_ref$Block == SyntenicBlock$RefBlock[i]]
    seqs_ref[[i]][,4:5] <- sam_ref[match(seqs_ref[[i]]$Hap,sam_ref$Hap),c("StartPos","Seq")]
    seqs_ref[[i]]$StartSeq <- SyntenicBlock$Start[i]-seqs_ref[[i]]$StartPos+1
    seqs_ref[[i]]$EndSeq <- SyntenicBlock$End[i]-seqs_ref[[i]]$StartPos+1
    for (j in 1:RefPloidy) {
      # if haplotype without D/I event, StartSeq & EndSeq don't change
      if (!is.null(cigar_ref[[seqs_ref[[i]]$Hap[j]]])) {
        # calculate seqs_ref[[i]]$StartSeq
        diff <- seqs_ref[[i]]$StartSeq[j]-cigar_ref[[seqs_ref[[i]]$Hap[j]]]$IDpos
        # <0 D/I events occurred after StartSeq, >=0 D/I events occurred before StartSeq
        # all(diff<0), all D/I events occurred after StartSeq, so StartSeq & EndSeq don't change
        if (!all(diff<=0)) {
          # find out the last influenced event
          rows <- max(which(diff==min(diff[diff>0])))
          # calculate StartSeq by add the length of insertion and minus the length ofdeletion
          seqs_ref[[i]]$StartSeq[j] <- seqs_ref[[i]]$StartSeq[j]+cigar_ref[[seqs_ref[[i]]$Hap[j]]]$CLength[rows]
        }
        # calculate seqs_ref[[i]]$EndSeq, same principle with seqs_ref[[i]]$StartSeq
        diff <- seqs_ref[[i]]$EndSeq[j]-cigar_ref[[seqs_ref[[i]]$Hap[j]]]$IDpos
        if (!all(diff<=0)) {
          rows <- max(which(diff==min(diff[diff>0])))
          seqs_ref[[i]]$EndSeq[j] <- seqs_ref[[i]]$EndSeq[j]+cigar_ref[[seqs_ref[[i]]$Hap[j]]]$CLength[rows]
        }
      }
    }
    seqs_ref[[i]]$Seq <- substr(seqs_ref[[i]]$Seq,seqs_ref[[i]]$StartSeq,seqs_ref[[i]]$EndSeq)
    seqs_ref[[i]]$Name <- paste0(">",ref,"_",seqs_ref[[i]]$Hap,"_",seqs_ref[[i]]$WRBlock)
  }
  names(seqs_ref) <- paste(SyntenicBlock$RefBlock,SyntenicBlock$WRBlock,sep="_")
  seqs_ref
}

extract_sequence_wr <- function(SyntenicBlock,sam_wr,cigar_wr,wr) {
  seqs_wr <- list()
  for (i in 1:nrow(SyntenicBlock)) {
    WRPloidy <- SyntenicBlock$WRPloidy[i]
    seqs_wr[[i]] <- data.frame(matrix(NA,WRPloidy,0))
    seqs_wr[[i]]$WRBlock <- rep(SyntenicBlock$WRBlock[i],WRPloidy)
    seqs_wr[[i]]$RefBlock <- rep(SyntenicBlock$RefBlock[i],WRPloidy)
    seqs_wr[[i]]$Hap <- sam_wr$Hap[sam_wr$Block == SyntenicBlock$WRBlock[i]]
    seqs_wr[[i]][,4:5] <- sam_wr[match(seqs_wr[[i]]$Hap,sam_wr$Hap),c("StartPos","Seq")]
    seqs_wr[[i]]$StartSeq <- SyntenicBlock$Start[i]-seqs_wr[[i]]$StartPos+1
    seqs_wr[[i]]$EndSeq <- SyntenicBlock$End[i]-seqs_wr[[i]]$StartPos+1
    for (j in 1:WRPloidy) {
      if (!is.null(cigar_wr[[seqs_wr[[i]]$Hap[j]]])) {
        # calculate seqs_wr[[i]]$StartSeq
        diff <- seqs_wr[[i]]$StartSeq[j]-cigar_wr[[seqs_wr[[i]]$Hap[j]]]$IDpos
        # <0 D/I events occurred after StartSeq, >=0 D/I events occurred before StartSeq
        # all(diff<0), all D/I events occurred after StartSeq, so StartSeq & EndSeq don't change
        if (!all(diff<=0)) {
          # find out the last influenced event
          rows <- max(which(diff==min(diff[diff>0])))
          # calculate StartSeq by add the length of insertion and minus the length ofdeletion
          seqs_wr[[i]]$StartSeq[j] <- seqs_wr[[i]]$StartSeq[j]+cigar_wr[[seqs_wr[[i]]$Hap[j]]]$CLength[rows]
        }
        # calculate seqs_wr[[i]]$EndSeq, same principle with seqs_wr[[i]]$StartSeq
        diff <- seqs_wr[[i]]$EndSeq[j]-cigar_wr[[seqs_wr[[i]]$Hap[j]]]$IDpos
        if (!all(diff<=0)) {
          rows <- max(which(diff==min(diff[diff>0])))
          seqs_wr[[i]]$EndSeq[j] <- seqs_wr[[i]]$EndSeq[j]+cigar_wr[[seqs_wr[[i]]$Hap[j]]]$CLength[rows]
        }
      }
    }
    seqs_wr[[i]]$Seq <- substr(seqs_wr[[i]]$Seq,seqs_wr[[i]]$StartSeq,seqs_wr[[i]]$EndSeq)
    seqs_wr[[i]]$Name <- paste0(">",wr,"_",seqs_wr[[i]]$Hap,"_",seqs_wr[[i]]$WRBlock)
  }
  names(seqs_wr) <- paste(SyntenicBlock$WRBlock,SyntenicBlock$RefBlock,sep="_")
  seqs_wr
}

na2zero <- function(x) {x[is.na(x)] <- 0; return(x)}
null2zero <- function(x) {x[length(x)==0] <- 0; return(x)}

topology_analysis_iqtree <- function(chr) {
  dirs <- list.files(path="tree",pattern=paste0(".*",chr))
  statistic <- data.frame(matrix(NA,length(dirs),5))
  colnames(statistic) <- c("SeqName","Monophy","BranchLenT","BranchLenA","BranchLen")
  # statistic$SeqName <- gsub("-[0-9]*[.]nwk","",trees)
  statistic$SeqName <- gsub("[.]fas","",dirs)
  for (i in 1:length(dirs)) {
    # statistic[i,2:4] <- unlist(strsplit(trees[i],"_"))[2:4]
    path <- paste0("tree/",dirs[i],"/",dirs[i],".treefile")
    if (file.exists(path)) {
      tree <- read.tree(path)
      tree$tip.label[grep(ref,tree$tip.label)] <- "T"
      # statistic$HapNo4[i] <- length(grep("T",tree$tip.label))
      tree$tip.label[grep(wr,tree$tip.label)] <- "A"
      # statistic$HapNo6[i] <- length(grep("A",tree$tip.label))
      if (is.monophyletic(tree,which(tree$tip.label=="T"))) {
        statistic$Monophy[i] <- TRUE
        statistic$BranchLenT[i] <- null2zero(na2zero(tree$edge.length[which.edge(tree,mrca.phylo(tree,which(tree$tip.label=="T")))]))
        statistic$BranchLenA[i] <- null2zero(na2zero(tree$edge.length[which.edge(tree,mrca.phylo(tree,which(tree$tip.label=="A")))]))
        statistic$BranchLen[i] <- statistic$BranchLenT[i] + statistic$BranchLenA[i]
      } else {
        statistic$Monophy[i] <- FALSE
        statistic$BranchLen[i] <- 0
      }
    }
  }
  SyntenicBlock <- read.delim(paste(ref,wr,chr,"ref_WR_final_blocks.txt",sep="_"),sep=" ",as.is=T)
  statistic <- merge(SyntenicBlock,statistic[,c(1,2,5)],by="SeqName",all=F)
  statistic <- statistic[!is.na(statistic$Monophy),]
  statistic
}

topology_analysis_upgma <- function(chr) {
  trees <- list.files(path="tree_UPGMA", pattern=paste0(chr,"_.*nwk"))
  # trees <- trees[1:1000]
  statistic <- data.frame(matrix(NA,length(trees),5))
  colnames(statistic) <- c("SeqName","Monophy","BranchLenT","BranchLenA","BranchLen")
  statistic$SeqName <- gsub("-[0-9]*[.]nwk","",trees)
  null <- NULL
  for (i in 1:length(trees)) {
    # statistic[i,2:4] <- unlist(strsplit(trees[i],"_"))[2:4]
    tree <- read.tree(paste0("tree_UPGMA/",trees[i]))
    # remove trees without branch length, which has
    if (!is.null(tree$edge.length)) {
      tree$tip.label[grep(ref,tree$tip.label)] <- "T"
      # statistic$HapNo4[i] <- length(grep("T",tree$tip.label))
      tree$tip.label[grep(wr,tree$tip.label)] <- "A"
      # statistic$HapNo6[i] <- length(grep("A",tree$tip.label))
      if (is.monophyletic(tree,which(tree$tip.label=="T"))) {
        statistic$Monophy[i] <- TRUE
        statistic$BranchLenT[i] <- null2zero(na2zero(tree$edge.length[which.edge(tree,mrca.phylo(tree,which(tree$tip.label=="T")))]))
        statistic$BranchLenA[i] <- null2zero(na2zero(tree$edge.length[which.edge(tree,mrca.phylo(tree,which(tree$tip.label=="A")))]))
        statistic$BranchLen[i] <- statistic$BranchLenT[i] + statistic$BranchLenA[i]
        # statistic$BranchLen[i] <- na2zero(tree$edge.length[which.edge(tree,mrca.phylo(tree,which(tree$tip.label=="T")))]) +
        #   na2zero(tree$edge.length[which.edge(tree,mrca.phylo(tree,which(tree$tip.label=="A")))])
      } else {
        statistic$Monophy[i] <- FALSE
        statistic$BranchLen[i] <- 0
      }
      null <- c(null,i)
    }
  }
  colnames(statistic) <- c("SeqName","UMonophy","UBranchLenT","UBranchLenA","UBranchLen")
  statistic <- merge(statistic_chr,statistic[,c(1,2,5)],by="SeqName",all=F)
  statistic <- statistic[!is.na(statistic$UMonophy),]
  statistic
}

summarise_ploid <- function(chr,statistic_chr) {
  ploid_tree <- data.frame(matrix(NA,length(unique(statistic_chr$Ploidy))+1,0))
  ploid_tree$WR <- wr
  ploid_tree$Chr <- chr
  ## IQTREE
  ploid_Monophy <- lapply(split(statistic_chr$Monophy, statistic_chr$Ploidy),function(x) {length(which(x==TRUE))/length(x)})
  nrows <- nrow(ploid_tree)-1
  ploid_tree$Ploidy <- c(names(ploid_Monophy),"All")
  ploid_tree$Num[nrows+1] <- nrow(statistic_chr)
  ploid_tree$Num[1:nrows] <- lapply(split(statistic_chr$Ploidy, statistic_chr$Ploidy),length)
  ploid_tree$Percent <- unlist(ploid_tree$Num)/nrow(statistic_chr)
  ploid_tree$MeanLength[nrows+1] <- mean(statistic_chr$Length)
  ploid_tree$MeanLength[1:nrows] <- lapply(split(statistic_chr$Length, statistic_chr$Ploidy),mean)
  ploid_tree$GenomeCoverage[nrows+1] <- sum(statistic_chr$Length)/GenomeInfo$Length[j]
  ploid_tree$GenomeCoverage[1:nrows] <- unlist(lapply(split(statistic_chr$Length, statistic_chr$Ploidy),sum))/GenomeInfo$Length[j]
  ploid_tree$IQMonoRatio[nrows+1] <- length(which(statistic_chr$Monophy==TRUE))/length(trees)
  ploid_tree$IQMonoRatio[1:nrows] <- unlist(ploid_Monophy)
  ploid_tree$IQBranchLen[nrows+1] <- mean(statistic_chr$BranchLen)
  ploid_tree$IQBranchLen[1:nrows] <- unlist(lapply(split(statistic_chr$BranchLen, statistic_chr$Ploidy),function(x) {mean(x,na.rm = T)}))
  ## UPGMA
  ploid_Monophy <- lapply(split(statistic_chr$UMonophy, statistic_chr$Ploidy),function(x) {length(which(x==TRUE))/length(x)})
  ploid_tree$UPGMAMonoRatio[nrows+1] <- length(which(statistic_chr$UMonophy==TRUE))/length(trees)
  ploid_tree$UPGMAMonoRatio[1:nrows] <- unlist(ploid_Monophy)
  ploid_tree$UPGMABranchLen[nrows+1] <- mean(statistic_chr$UBranchLen,na.rm = T)
  ploid_tree$UPGMABranchLen[1:nrows] <- unlist(lapply(split(statistic_chr$UBranchLen, statistic_chr$Ploidy),function(x) {mean(x,na.rm = T)}))
  ploid_tree
}

extract_gene_region <- function(path,statistic_all) {
  gff <- read.delim(path,header=F,as.is=T,sep="")
  gff <- gff[gff$V3=="gene",c(1,4,5,9)]
  colnames(gff) <- c("Chr","Start","End","Gene")
  chrs <- unique(statistic_all$Chr)
  rows <- NULL
  for (j in 1:length(chrs)) {
    sub <- gff[gff$Chr==chrs[j],]
    statistic_chr <- statistic_all[statistic_all$Chr==chrs[j],]
    for (i in 1:nrow(sub)) {
      num <- rownames(statistic_chr)[which(statistic_chr$Start>=sub$Start[i] & statistic_chr$Start<sub$End[i] |
                                             statistic_chr$End>sub$Start[i] & statistic_chr$End<=sub$End[i] |
                                             statistic_chr$Start<=sub$Start[i] & statistic_chr$End>=sub$End[i])]
      rows <- c(rows,num)
    }
  }
  rows
}

gene_conversion_upgma <- function(TreeFolder,TreeSuffix) {
  trees <- list.files(path=TreeFolder, pattern=TreeSuffix)
  trees <- trees[gsub("-.*","",trees) %in% statistic_all$SeqName]
  topology <- data.frame(matrix(NA,length(trees),5))
  colnames(topology) <- c("Tree","Tmono","Topo","DeRefl","Type")
  topology$Tree <- trees
  for (i in 1:length(trees)) {
    tree <- read.tree(paste0(TreeFolder,"/",trees[i]))
    tree$tip.label[grep(ref,tree$tip.label)] <- "T"
    tree$tip.label[grep(wr,tree$tip.label)] <- "A"
    # Topology
    node <- Children(tree,11)
    tips <- Descendants(tree, node)
    topology$Topo[i] <- paste(sort(unlist(lapply(tips,length))),collapse = ":")
    # Componant of small branch
    topology$DeRefl[i] <- paste(sort(tree$tip.label[Descendants(tree,node[which.min(lapply(tips,length))])[[1]]]),collapse = "")
    # Count monophy of T
    node <- node[which(unlist(lapply(Descendants(tree, node),function(x){all(tree$tip.label[x]=="T")})))]
    if(length(node)>0) {
      topology$Tmono[i] <- length(Descendants(tree, node)[[1]])
    }
  }
  topology$Tmono[is.na(topology$Tmono)] <- 0
  topology$Type[which(topology$DeRefl %in% c("TTT","ATTT","TTTT"))] <- "tetra2di"
  topology$Type[which(topology$DeRefl %in% c("T","AT","AAT","AAAT"))] <- "di2tetra"
  topology$Type[is.na(topology$Type)] <- "Others"
  topology
}

gene_conversion_iqtree <- function(TreeFolder,TreeFolderPattern) {
  trees <- list.files(TreeFolder, TreeFolderPattern)
  trees <- trees[gsub(".fas","",trees) %in% statistic_all$SeqName]
  # trees <- trees[gsub(".fas.treefile","",trees) %in% statistic_all$SeqName]
  topology <- data.frame(matrix(NA,length(trees),5))
  colnames(topology) <- c("Tree","Tmono","Topo","DeRefl","Type")
  topology$Tree <- trees
  for (i in 1:length(trees)) {
    tree <- read.tree(paste0(TreeFolder,"/",trees[i],"/",trees[i],".treefile"))
    # tree <- read.tree(paste0("tree/",trees[i]))
    tree$tip.label[grep(ref,tree$tip.label)] <- "T"
    tree$tip.label[grep(wr,tree$tip.label)] <- "A"
    # Topology
    node <- Children(tree,11)
    tips <- Descendants(tree, node)
    topology$Topo[i] <- paste(sort(unlist(lapply(tips,length))),collapse = ":")
    # Componant of small branch
    tt <- table(tree$tip.label[Descendants(tree,node[which.max(lapply(tips,length))])[[1]]])
    if(dim(tt)>1) {
      tt <- table(tree$tip.label)-tt
      topology$DeRefl[i] <- paste(c(rep(names(tt)[1],tt[1]),rep(names(tt)[2],tt[2])),collapse = "")
    } else {
      tt <- table(tree$tip.label)[names(tt)]-tt
      topology$DeRefl[i] <- paste(rep(names(tt),tt),collapse = "")
    }
    # Count monophy of T
    node <- node[which(unlist(lapply(Descendants(tree, node),function(x){all(tree$tip.label[x]=="T")})))]
    if(length(node)>0) {
      topology$Tmono[i] <- length(Descendants(tree, node)[[1]])
    }
  }
  topology$Tmono[is.na(topology$Tmono)] <- 0
  topology$Type[which(topology$DeRefl %in% c("TTT","ATTT","TTTT"))] <- "tetra2di"
  topology$Type[which(topology$DeRefl %in% c("T","AT","AAT","AAAT"))] <- "di2tetra"
  topology
}