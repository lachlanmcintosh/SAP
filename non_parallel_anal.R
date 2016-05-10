#############################
##### SNP array package #####
#############################
# Author: Lachlan McIntosh
rm(list=ls())
library(ggplot2)
library(PSCBS)
library(data.table)
library(mgcv)
library(cluster)


#### unpaired test dataset: SAMPLE 13 ####
args=(commandArgs(TRUE))
#sample <-  as.character(args[1])
sample <- 1
PARENT_FOLDER = paste("/home/users/allstaff/lmcintosh/P2_LEON/ILO2.58-8359/","ascat_",as.character(sample),"/",sep="")
filename = paste(PARENT_FOLDER,"/raw_data.txt",sep="")
# load functions
setwd("/wehisan/home/allstaff/l/lmcintosh/SAP")
source(paste(getwd(),"/utils.R",sep=""))
#source("/wehisan/home/allstaff/l/lmcintosh/SAP/utils.R")

data <- read_raw_illumina_file(filename)
data <- preprocess_raw_data(data)

result <- do_some_seg(data,0)
data <- result$data
segments <- result$segments
data <-  reformat(data,segments)

#### find het SNPs if the SNPs are not annotated ####
# not required for the illumina data
# find_all_het_SNPs(data)
#### find het SNPs if the SNPs are not annotated ####
GC <- get_GC_file()
data <- merge_data_with_GC(data,GC)

get_TCN_track <- function(mat,seg_est_name){
  g <- ggplot(mat) +
    geom_segment(aes_string(x = "tcnStart", y = seg_est_name, xend = "tcnEnd", yend = seg_est_name),size=2)+
    facet_grid(.~chromosome,scales = "free_x", space = "free")
  return(g)
}
temp_data <- data
data <- temp_data

segr <- 1
pr <- 0
old_segments <- list(segments)
segments <- get_new_seg_estimates(data=data,segments=segments,
                                  new_snp_name = "CT",
                                  new_seg_name = paste("segmentCN_pre",as.character(pr),sep=""))
data <- put_seg_estimates_into_data_array(data=data,segments=segments,
                                          segment_name = paste("segmentCN_pre",as.character(pr),sep=""),
                                          new_snp_name = paste("segmentCN_pre",as.character(pr),sep=""))
pr
while(segr < 4){ # could potentially do this until convergence.... but might just shrink all to zero.
  while(pr < 4*segr-1){
    pr <- pr + 1
    print(pr)
    result <-  renormalise(data,segments,pr)
    data <- result$data
    segments <- result$segments
    print("gam started")
    gam <- gam(CNsnp/CNseg ~ 0+s(GC50) + s(GC150) + s(GC500) + s(GC2500) , data=data) #+ ti(GC50,CNseg)+ ti(GC150,CNseg) + ti(GC500,CNseg) + ti(GC2500,CNseg), data=data)
    print("gam_done")
    if(!is.null(gam$na.action)){
      data <- data[-gam$na.action,]
    }

    data$CNsnp <- gam$residuals*data$CNseg
    data$CNsnp <- data$CNsnp/median(data$CNsnp,na.rm=TRUE)*2

    data[,paste("CN_pre",as.character(pr),sep="")] <- data$CNsnp
    segments <- get_new_seg_estimates(data=data,segments=segments,
                                      new_snp_name = paste("CN_pre",as.character(pr),sep=""),
                                      new_seg_name = paste("segmentCN_pre",as.character(pr),sep=""))
    data <- put_seg_estimates_into_data_array(data=data,segments=segments,
              segment_name = paste("segmentCN_pre",as.character(pr),sep=""),
              new_snp_name = paste("segmentCN_pre",as.character(pr),sep=""))
    adj_thresh = min(max(0.005*pr,0.02),0.04)
    segments <- cluster_ACN(dat=segments,iterations=10,local_adj_thresh=adj_thresh,
                            global_adj_thresh=adj_thresh/2,p_close_local = 0.2,p_close_global=0.2,
                            TCN=paste("segmentCN_pre",as.character(pr),sep=""),
                            TCN_se=paste("segmentCN_pre",as.character(pr),"_stderr",sep=""))

    g1 <- get_TCN_track(segments,paste("segmentCN_pre",as.character(pr-1),sep="")) + ylim(0,5)
    g2 <- get_TCN_track(segments,paste("segmentCN_pre",as.character(pr),sep="")) + ylim(0,5)
    g3 <- get_TCN_track(segments,paste("segmentCN_pre_nogam",as.character(pr),sep="")) + ylim(0,5)
    g4 <- get_TCN_track(segments,paste("segmentCN_pre","0",sep="")) + ylim(0,5)
    grid.arrange(g1,g2,g3,g4)
    stop()
  }
  segr <- segr+1
  if(segr <5){
    old_segments[[segr]] <- segments
    pr <- pr + 1
    result <- do_some_seg(data,pr,paste("CN_pre",as.character(pr-1),sep=""))
    segments <- result$segments
    fit <- result$fit
    data <- result$data
    segments <- get_new_seg_estimates(data=data,segments=segments,
                                      new_snp_name = paste("CN_pre",as.character(pr),sep=""),
                                      new_seg_name = paste("segmentCN_pre",as.character(pr),sep=""))
    data <- put_seg_estimates_into_data_array(data=data,segments=segments,
                                              segment_name = paste("segmentCN_pre",as.character(pr),sep=""),
                                              new_snp_name = paste("segmentCN_pre",as.character(pr),sep=""))
  }
}


get_TCN_track(segments,paste("segmentCN_pre",as.character(0),sep="")) + ylim(0,5)

# sample <- 1
# PARENT_FOLDER = paste("/home/users/allstaff/lmcintosh/P2_LEON/ILO2.58-8359/","ascat_",as.character(sample),"/",sep="")PARENT_FOLDER
# load(paste(PARENT_FOLDER,"data.Rda",sep=""))
# load(paste(PARENT_FOLDER,"segments.Rda",sep=""))
for(pr in 0:12){
  name <- as.character(pr)
  if(pr%%4 == 0) {
    segments <- get_new_seg_estimates(data=data,segments=segments,
    new_snp_name = paste("segmentCN_pre",name,sep=""),
    new_seg_name = paste("segmentCN_pre",name,sep=""))
  }else{
    segments <- get_new_seg_estimates(data=data,segments=segments,
    new_snp_name = paste("CN_pre",name,sep=""),
    new_seg_name = paste("segmentCN_pre",name,sep=""))
  }
  print(pr)
}

library(gridExtra)
g1 <- get_TCN_track(segments,"segmentCN_pre0") + ylim(0,4)
g2 <- get_TCN_track(segments,"segmentCN_pre4") + ylim(0,4)
grid.arrange(g1,g2)

g1 <-get_TCN_track(segments,"segmentCN_pre0")+ylim(0,4)
g2 <-get_TCN_track(segments,"segmentCN_pre1")+ylim(0,4)
g3 <-get_TCN_track(segments,"segmentCN_pre_nogam1")+ylim(0,4)
g4 <-get_TCN_track(segments,"segmentCN_pre3")+ylim(0,4)
grid.arrange(g1,g3)
grid.arrange(g1,g2,g3,g4)

g1 <-get_TCN_track(segments,"segmentCN_pre4")+ylim(0,4)
g2 <-get_TCN_track(segments,"segmentCN_pre5")+ylim(0,4)
g3 <-get_TCN_track(segments,"segmentCN_pre6")+ylim(0,4)
g4 <-get_TCN_track(segments,"segmentCN_pre7")+ylim(0,4)
grid.arrange(g1,g2,g3,g4)

g1 <-get_TCN_track(segments,"segmentCN_pre8")+ylim(0,4)
g2 <-get_TCN_track(segments,"segmentCN_pre9")+ylim(0,4)
g3 <-get_TCN_track(segments,"segmentCN_pre10")+ylim(0,4)
g4 <-get_TCN_track(segments,"segmentCN_pre12")+ylim(0,4)
grid.arrange(g1,g2,g3,g4)





> get_new_seg_estimates
function(data,segments,new_snp_name,new_seg_name){
  segments[,new_seg_name] <- sapply(1:nrow(segments),function(x) mean(data[which(data$segment == x),new_snp_name],na.rm=TRUE))
  return(segments)
}


# tidyness code
for(pr in 1:16){
  name <- as.character(pr)
  if(pr%%4 == 0) {
    segments <- get_new_seg_estimates(data=data,segments=segments,
                                      new_snp_name = paste("segmentCN_pre",name,sep=""),
                                      new_seg_name = paste("segmentCN_pre",name,sep=""))
  }else{
    segments <- get_new_seg_estimates(data=data,segments=segments,
                                      new_snp_name = paste("CN_pre",name,sep=""),
                                      new_seg_name = paste("segmentCN_pre",name,sep=""))
  }
}

> for(pr in 1:16){
  +   name <- as.character(pr)
  +   if(pr%%4 == 0) {
    +     segments <- get_new_seg_estimates(data=data,segments=segments,
                                            +                                       new_snp_name = paste("segmentCN_pre",name,sep=""),
                                            +                                       new_seg_name = paste("segmentCN_pre",name,sep=""))
    +   }else{
      +   segments <- get_new_seg_estimates(data=data,segments=segments,
                                            +                                     new_snp_name = paste("CN_pre",name,sep=""),
                                            +                                     new_seg_name = paste("segmentCN_pre",name,sep=""))
      +   }
  +   print(pr)
  + }

g1 <-get_TCN_track(segments,"segmentCN_pre1")+ylim(0,4)
g2 <-get_TCN_track(segments,"segmentCN_pre2")+ylim(0,4)
g3 <-get_TCN_track(segments,"segmentCN_pre3")+ylim(0,4)
g4 <-get_TCN_track(segments,"segmentCN_pre4")+ylim(0,4)
grid.arrange(g1,g2,g3,g4)

g1 <-get_TCN_track(segments,"segmentCN_pre5")+ylim(0,4)
g2 <-get_TCN_track(segments,"segmentCN_pre6")+ylim(0,4)
g3 <-get_TCN_track(segments,"segmentCN_pre7")+ylim(0,4)
g4 <-get_TCN_track(segments,"segmentCN_pre8")+ylim(0,4)
grid.arrange(g1,g2,g3,g4)

g1 <-get_TCN_track(segments,"segmentCN_pre9")+ylim(0,4)
g2 <-get_TCN_track(segments,"segmentCN_pre10")+ylim(0,4)
g3 <-get_TCN_track(segments,"segmentCN_pre11")+ylim(0,4)
g4 <-get_TCN_track(segments,"segmentCN_pre12")+ylim(0,4)
grid.arrange(g1,g2,g3,g4)

g1 <-get_TCN_track(segments,"segmentCN_pre13")+ylim(0,4)
g2 <-get_TCN_track(segments,"segmentCN_pre14")+ylim(0,4)
g3 <-get_TCN_track(segments,"segmentCN_pre15")+ylim(0,4)
g4 <-get_TCN_track(segments,"segmentCN_pre16")+ylim(0,4)
grid.arrange(g1,g2,g3,g4)
