#############################
##### SNP array package #####
#############################
# Author: Lachlan McIntosh  #
rm(list=ls())
library(ggplot2)
library(PSCBS)
library(data.table)
library(mgcv)
library(cluster)
library(gridExtra)


#### unpaired test dataset: SAMPLE 13 ####
# args=(commandArgs(trailingOnly = TRUE))

args=(commandArgs(TRUE))
##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  a = 1
  b = c(1,1,1)
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

#sample <-  as.character(args[1])
# sample <- 1
# PARENT_FOLDER = paste("/home/users/allstaff/lmcintosh/P2_LEON/ILO2.58-8359/","ascat_",as.character(sample),"/",sep="")
# filename = paste(PARENT_FOLDER,"/raw_data.txt",sep="")
# filename = commandArgs(trailingOnly=T)[1]

# filename="/export/share/prkfs2/shared/bioinf-data/Papenfuss_lab/projects/melanoma/SNP_arrays/data/raw/ILO2.58-6865/ILO2_58-6865/ILO2_58-6865_FinalReport1.txt"

# load functions




# filename="/wehisan/home/allstaff/l/lmcintosh/run_pres/ILO2.58-7453/ascat_6/raw_data.txt"
setwd("/wehisan/home/allstaff/l/lmcintosh/SAP")
source(paste(getwd(),"/utils.R",sep=""))
#source("/wehisan/home/allstaff/l/lmcintosh/SAP/utils.R")

#="./ILO2.58-6865/ascat_10/raw_data.txt"
BASE = "/home/users/allstaff/lmcintosh/run_pres"
filename <- paste(BASE,substring(filename,2),sep="")
data <- read_raw_illumina_file(filename)
data <- preprocess_raw_data(data)

result <- do_some_seg(data,0)
data <- result$data
segments <- result$segments
result <-  reformat(data,segments)
data <- result$data
segments <- result$segments

#### find het SNPs if the SNPs are not annotated ####
# not required for the illumina data
# find_all_het_SNPs(data)
#### find het SNPs if the SNPs are not annotated ####
GC <- get_GC_file()
data <- merge_data_with_GC(data,GC)


temp_data <- data
temp_seg <- segments

data <- temp_data
segments <- temp_seg
# segments$seg_name <- rownames(segments)
# data[!complete.cases(data$CT),]
# data <- data[complete.cases(data$CT),]

segr <- 1
pr <- 0
old_segments <- list(segments)
segments <- get_new_seg_estimates(data=data,segments=segments,
                                  new_snp_name = "CT",
                                  new_seg_name = paste("segmentCN_pre",as.character(pr),sep=""))
data <- put_seg_estimates_into_data_array(data=data,segments=segments,
                                          segment_name = paste("segmentCN_pre",as.character(pr),sep=""),
                                          new_snp_name = paste("segmentCN_pre",as.character(pr),sep=""))
segments[,"length"] <- segments[,"tcnNbrOfLoci"]
segments$TCN <- segments[,paste("segmentCN_pre",as.character(pr),sep="")]
# interesting <- which(!(complete.cases(segments$TCN) & complete.cases(segments$length) & complete.cases(segments$chromosome)))
# segments[sort(c(interesting+1,interesting,interesting-1)),]
# data[which(data$segment == 152),]



local_adj_thresh = min(max(0.005*pr,0.02),0.06)
global_adj_thresh = max(local_adj_thresh -0.03,0)
segments <- cluster_ACN2(dat=segments,iterations=30,local_adj_thresh=local_adj_thresh,
                         global_adj_thresh=global_adj_thresh,p_close_local = 0.5,p_close_global=0.5,
                         TCN=paste("segmentCN_pre",as.character(pr),sep=""),
                         TCN_se=paste("segmentCN_pre",as.character(pr),"_stderr",sep=""))
data <- put_seg_estimates_into_data_array(data=data,segments=segments,
                                          segment_name = "TCN",
                                          new_snp_name = paste("segmentCN_pre_clustered",as.character(pr),sep=""))
segments[,paste("segmentCN_pre_clustered",as.character(pr),sep="")] <- segments$TCN

# g2 <- get_TCN_track_precision(segments,paste("segmentCN_pre",as.character(pr),sep="")) + ylim(0,5)
# g3 <- get_TCN_track_precision(segments,paste("segmentCN_pre_clustered",as.character(pr),sep="")) + ylim(0,5)
# # g4 <- get_TCN_track_precision(segments,paste("segmentCN_pre","0",sep="")) + ylim(0,5)
# grid.arrange(g2,g3)

# interesting <- which(abs(segments[,paste("segmentCN_pre_clustered",as.character(pr),sep="")] - segments[,paste("segmentCN_pre",as.character(pr),sep="")]) >local_adj_thresh)
# segments[sort(c(interesting,interesting+1,interesting-1)),]
# after we finish the clustering algorithm we should actually move the snps as otherwise the other things will fit this arbritrary change in segment values?


pr

SEG_ROUNDS = 4
PR_ROUNDS = 4
setwd("/wehisan/home/allstaff/l/lmcintosh/SAP")
source(paste(getwd(),"/utils.R",sep=""))
try(while(segr < SEG_ROUNDS){
  while(pr < PR_ROUNDS*segr){
    pr <- pr + 1
    print(pr)
    result <-  renormalise2(data,segments,pr)
    data <- result$data
    segments <- result$segments
    print("gam started")
    gam <- gam(CNsnp/CNseg ~ 0 + te(GC50) + te(GC150) + te(GC500) + te(GC2500)+te(GC500,CNseg), data=data)
    print("gam_done")
    if(!is.null(gam$na.action)){
      data <- data[-gam$na.action,]
    }

    summary(gam)
    data$CNsnp <- gam$residuals*data$CNseg
    data$CNsnp <- data$CNsnp/median(data$CNsnp,na.rm=TRUE)*2

    data[,paste("CN_pre",as.character(pr),sep="")] <- data$CNsnp
    segments <- get_new_seg_estimates(data=data,segments=segments,
                                      new_snp_name = paste("CN_pre",as.character(pr),sep=""),
                                      new_seg_name = paste("segmentCN_pre",as.character(pr),sep=""))
    data <- put_seg_estimates_into_data_array(data=data,segments=segments,
                                              segment_name = paste("segmentCN_pre",as.character(pr),sep=""),
                                              new_snp_name = paste("segmentCN_pre",as.character(pr),sep=""))

    local_adj_thresh = min(max(0.1*pr,0.02),0.12)
    global_adj_thresh = max(local_adj_thresh -0.03,0)/3
    segments <- cluster_ACN2(dat=segments,iterations=30,local_adj_thresh=local_adj_thresh,
                            global_adj_thresh=global_adj_thresh,p_close_local = 0.5,p_close_global=0.5,
                            TCN=paste("segmentCN_pre",as.character(pr),sep=""),
                            TCN_se=paste("segmentCN_pre",as.character(pr),"_stderr",sep=""))
    data <- put_seg_estimates_into_data_array(data=data,segments=segments,
                                              segment_name = "TCN",
                                              new_snp_name = paste("segmentCN_pre_clustered",as.character(pr),sep=""))
    segments[,paste("segmentCN_pre_clustered",as.character(pr),sep="")] <- segments$TCN

    quantile(data[,paste("segmentCN_pre",as.character(pr-1),sep="")],c(0.5,0.9))
    quantile(data[,paste("segmentCN_pre",as.character(pr-1),sep="")],c(0.5,0.9))

    # g1 <- get_TCN_track_precision(segments,paste("segmentCN_pre",as.character(pr-1),sep="")) + ylim(0,5)
    # g2 <- get_TCN_track_precision(segments,paste("segmentCN_pre",as.character(pr),sep="")) + ylim(0,5)
    # g3 <- get_TCN_track_precision(segments,paste("segmentCN_pre_clustered",as.character(pr-1),sep="")) + ylim(0,5)
    # g4 <- get_TCN_track_precision(segments,paste("segmentCN_pre_nogam",as.character(pr),sep="")) + ylim(0,5)
    # # g4 <- get_TCN_track_precision(segments,paste("segmentCN_pre","0",sep="")) + ylim(0,5)
    # grid.arrange(g1,g2,g3,g4)
  }
  segr <- segr+1
  old_segments[[segr]] <- segments
  pr <- pr + 1
  result <- do_some_seg(data,pr,paste("segmentCN_pre_clustered",as.character(pr-1),sep=""))
  segments <- result$segments
  fit <- result$fit
  data <- result$data
  segments <- get_new_seg_estimates(data=data,segments=segments,
                                    new_snp_name = paste("CN_pre",as.character(pr-1),sep=""),
                                    new_seg_name = paste("segmentCN_pre",as.character(pr),sep=""))
  data <- put_seg_estimates_into_data_array(data=data,segments=segments,
                                            segment_name = paste("segmentCN_pre",as.character(pr),sep=""),
                                            new_snp_name = paste("segmentCN_pre",as.character(pr),sep=""))
  local_adj_thresh = min(max(0.1*pr,0.02),0.12)
  global_adj_thresh = max(local_adj_thresh -0.03,0)/3
  segments <- cluster_ACN2(dat=segments,iterations=30,local_adj_thresh=local_adj_thresh,
                           global_adj_thresh=global_adj_thresh,p_close_local = 0.5,p_close_global=0.5,
                           TCN=paste("segmentCN_pre",as.character(pr),sep=""),
                           TCN_se=paste("segmentCN_pre",as.character(pr),"_stderr",sep=""))
  data <- put_seg_estimates_into_data_array(data=data,segments=segments,
                                            segment_name = "TCN",
                                            new_snp_name = paste("segmentCN_pre_clustered",as.character(pr),sep=""))
  segments[,paste("segmentCN_pre_clustered",as.character(pr),sep="")] <- segments$TCN
})

for(round in 0:pr){
  name <- as.character(round)
  if(pr%%PR_ROUNDS == 0) {
    segments <- get_new_seg_estimates(data=data,segments=segments,
      new_snp_name = paste("segmentCN_pre",name,sep=""),
      new_seg_name = paste("segmentCN_pre",name,sep=""))
  }else{
    segments <- get_new_seg_estimates(data=data,segments=segments,
      new_snp_name = paste("CT_pre",name,sep=""),
      new_seg_name = paste("segmentCN_pre",name,sep=""))
  }
  print(name)
}


saveRDS(segments, file=paste(filename,".segmented",sep=""))

# to load these later use:
# bar <- readRDS(file="data.Rda")

# now we want to save old_segments




# if we put in there that data values shrink toward their segment values we might just get it perfect.
# but leave it for now.
# i can't think what the implications of that would be.


# head(GC)
#
# local_adj_thresh = 0.12
# global_adj_thresh = 0.04
# segments <- cluster_ACN2(dat=segments,iterations=30,local_adj_thresh=local_adj_thresh,
#                          global_adj_thresh=global_adj_thresh,p_close_local = 0.5,p_close_global=0.5,
#                          TCN=paste("segmentCN_pre",as.character(pr),sep=""),
#                          TCN_se=paste("segmentCN_pre",as.character(pr),"_stderr",sep=""))
# data <- put_seg_estimates_into_data_array(data=data,segments=segments,
#                                           segment_name = "TCN",
#                                           new_snp_name = paste("segmentCN_pre_clustered",as.character(pr),sep=""))
# segments[,paste("segmentCN_pre_clustered",as.character(pr),sep="")] <- segments$TCN
#
# get_TCN_track_precision(segments,paste("segmentCN_pre_clustered",as.character(pr),sep="")) + ylim(0,5)
# get_TCN_track_precision(segments,paste("segmentCN_pre",as.character(pr),sep="")) + ylim(0,5)

# i notice that the start or end of chromosomes are often up - is there a way to check if this is biological or real?







# g1 <- get_TCN_track_precision(old_segments[[1]],"segmentCN_pre0") + ylim(0,4)
# g2 <- get_TCN_track_precision(segments,"segmentCN_pre9") + ylim(0,4)
# grid.arrange(g1,g2)
#
# g1 <-get_TCN_track(segments,"segmentCN_pre0")+ylim(0,4)
# g2 <-get_TCN_track(segments,"segmentCN_pre1")+ylim(0,4)
# g3 <-get_TCN_track(segments,"segmentCN_pre_nogam1")+ylim(0,4)
# g4 <-get_TCN_track(segments,"segmentCN_pre3")+ylim(0,4)
# grid.arrange(g1,g3)
# grid.arrange(g1,g2,g3,g4)


# # do a saftey look on everything:
# g1 <-get_TCN_track_precision(old_segments[[2]],"segmentCN_pre0")+ylim(0,4)
# g2 <-get_TCN_track_precision(old_segments[[2]],"segmentCN_pre1")+ylim(0,4)
# g3 <-get_TCN_track_precision(old_segments[[2]],"segmentCN_pre2")+ylim(0,4)
# g4 <-get_TCN_track_precision(old_segments[[2]],"segmentCN_pre3")+ylim(0,4)
# g5 <-get_TCN_track_precision(old_segments[[2]],"segmentCN_pre4")+ylim(0,4)
# grid.arrange(g1,g2,g3,g4,g5)
#
# g1 <-get_TCN_track_precision(old_segments[[3]],"segmentCN_pre5")+ylim(0,4)
# g2 <-get_TCN_track_precision(old_segments[[3]],"segmentCN_pre6")+ylim(0,4)
# g3 <-get_TCN_track_precision(old_segments[[3]],"segmentCN_pre7")+ylim(0,4)
# g4 <-get_TCN_track_precision(old_segments[[3]],"segmentCN_pre8")+ylim(0,4)
# grid.arrange(g1,g2,g3,g4)
#
# g1 <-get_TCN_track_precision(old_segments[[4]],"segmentCN_pre9")+ylim(0,4)
# g2 <-get_TCN_track_precision(old_segments[[4]],"segmentCN_pre10")+ylim(0,4)
# g3 <-get_TCN_track_precision(old_segments[[4]],"segmentCN_pre11")+ylim(0,4)
# g4 <-get_TCN_track_precision(old_segments[[4]],"segmentCN_pre12")+ylim(0,4)
# grid.arrange(g1,g2,g3,g4)
#
# g1 <-get_TCN_track_precision(segments,"segmentCN_pre13")+ylim(0,4)
# g2 <-get_TCN_track_precision(segments,"segmentCN_pre14")+ylim(0,4)
# g3 <-get_TCN_track_precision(segments,"segmentCN_pre15")+ylim(0,4)
# g4 <-get_TCN_track_precision(segments,"segmentCN_pre16")+ylim(0,4)
# grid.arrange(g1,g2,g3,g4)
#
#
# # in terms of the bias inference, forget it,
# # in terms of allelic balance, forget it,
# # just do soem stuff
#
# # tidyness code
#
# g1 <-get_TCN_track(segments,"segmentCN_pre1")+ylim(0,4)
# g2 <-get_TCN_track(segments,"segmentCN_pre2")+ylim(0,4)
# g3 <-get_TCN_track(segments,"segmentCN_pre3")+ylim(0,4)
# g4 <-get_TCN_track(segments,"segmentCN_pre4")+ylim(0,4)
# grid.arrange(g1,g2,g3,g4)
#
# g1 <-get_TCN_track(segments,"segmentCN_pre5")+ylim(0,4)
# g2 <-get_TCN_track(segments,"segmentCN_pre6")+ylim(0,4)
# g3 <-get_TCN_track(segments,"segmentCN_pre7")+ylim(0,4)
# g4 <-get_TCN_track(segments,"segmentCN_pre8")+ylim(0,4)
# grid.arrange(g1,g2,g3,g4)
#
# g1 <-get_TCN_track(segments,"segmentCN_pre9")+ylim(0,4)
# g2 <-get_TCN_track(segments,"segmentCN_pre10")+ylim(0,4)
# g3 <-get_TCN_track(segments,"segmentCN_pre11")+ylim(0,4)
# g4 <-get_TCN_track(segments,"segmentCN_pre12")+ylim(0,4)
# grid.arrange(g1,g2,g3,g4)
#
# g1 <-get_TCN_track(segments,"segmentCN_pre13")+ylim(0,4)
# g2 <-get_TCN_track(segments,"segmentCN_pre14")+ylim(0,4)
# g3 <-get_TCN_track(segments,"segmentCN_pre15")+ylim(0,4)
# g4 <-get_TCN_track(segments,"segmentCN_pre16")+ylim(0,4)
# grid.arrange(g1,g2,g3,g4)
