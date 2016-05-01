#############################
##### SNP array package #####
#############################
# Author: Lachlan McIntosh
rm(list=ls())
library(ggplot2)

#### unpaired test dataset: SAMPLE 13 ####
args=(commandArgs(TRUE))
sample <-  as.character(args[1])
#sample <- 9
PARENT_FOLDER = paste("/home/users/allstaff/lmcintosh/P2_LEON/ILO2.58-8359/","ascat_",as.character(sample),"/",sep="")
filename = paste(PARENT_FOLDER,"/raw_data.txt",sep="")
# load functions
source(paste(getwd(),"utils.R"))

library(data.table)
data <- read_raw_illumina_file(filename)
data <- rename_raw_columns(data)
data <- preprocess_raw_data(data)
library(PSCBS)
List[segments,fit,data] <- do_some_seg(data,0)
data <-  reformat(data,segments)

#### find het SNPs if the SNPs are not annotated ####
# not required for the illumina data
# find_all_het_SNPs(data)
#### find het SNPs if the SNPs are not annotated ####
data <-  renormalise(data,0)

#library(mgcv)
folder ="/home/users/allstaff/lmcintosh/ITH/SNP_stuff/data/"
GC_filenames <- lapply(1:50*50,function(x) paste(folder,"intervals_",as.character(x),".out",sep=""))
pGC <- lapply(GC_filenames, read_GC_file) # parallel was slower for loading files.

for (i in seq_along(pGC)) {
  setnames(pGC[[i]],c("chr","start","end","pAT","pGC","A","C","G","T","N","O","length"))
  pGC[[i]]$pos <- i*25+pGC[[i]]$start
  pGC[[i]]$loc <- as.numeric(pGC[[i]]$chr)+as.numeric(pGC[[i]]$pos)/max(pGC[[i]]$pos)
}
GC <- pGC[[1]]
GC$Position <- GC$start + GC$length/2
GC[,"GC50"] <- GC$pGC
for(i in 2:50){
  GC[,paste("GC",as.character(50*i),sep="")] <- pGC[[i]]$pGC
}
GC$name <- paste(as.character(GC[,c("chr")]),as.character(GC[,c("Position")]),sep=" ")
rm(pGC)
# need to find a faster way to do this.
data$name <- paste(as.character(data[,c("Chr")]),as.character(data[,c("Position")]),sep=" ")
data <- merge(data,GC,by=c("name","name"))
data$Position <- data$Position.x


library(mgcv)
temp_data <- data
data <- temp_data
segr <- 1
pr <- 0
old_segments <- list(segments)
while(segr < 5){ # could potentially do this until convergence.... but might just shrink all to zero.
  while(pr < 4*segr-1){
    old_segment_ids <- get_segment_ids(segments)
    data <-  renormalise(data,pr)
    #gam <- gam(R ~ 0+ reddyeAT_spec + greendyeAT_spec + reddyeGC_spec + greeendyeGC_spec + s(GC50)+s(GC150)+s(GC500)+s(GC2500),data=data)
    gam <- gam(R ~ 0+ reddyeAT + greendyeAT + reddyeGC + greeendyeGC + s(GC50)+s(GC150)+s(GC500)+s(GC2500),data=data)
    pr <- pr + 1
    name <- as.character(pr)
    print(name)
    summary(gam)
    data <- get_new_snp_estimates(gam=gam,data=data,
              new_snp_name = paste("CT_pre",name,sep=""),
              old_segment_name = paste("segmentCT_pre",as.character(pr-1),sep=""))
    segments <- get_new_seg_estimates(data=data,segments=segments,
                  new_snp_name = paste("CT_pre",name,sep=""),
                  new_seg_name = paste("segmentCT_pre",name,sep=""))
    data <- put_seg_estimates_into_data_array(data=data,segments=segments,
              segment_name = paste("segmentCT_pre",name,sep=""),
              new_snp_name = paste("segmentCT_pre",name,sep=""))
    }
    segr <- segr+1
    if(segr <5){
      old_segments[[segr]] <- segments
      pr <- pr + 1
      result <- do_some_seg(data,pr,paste("CT_pre",as.character(pr-1),sep=""))
      segments <- result$segs
      fit <- result$fit
      data <- result$data
    }
  }

           # > get_var(1,2,100,matrix(c(1,0,0,1),nrow=2))
           # Error: could not find function "get_var"
           # > get_var <- function(major,minor,n,sigma) sum(sigma)/(n*(major+minor)^2) - 2*(sigma[1,2]+sigma[2,2])/(n*minor*(major+minor))+sigma[2,2]/(n*minor^2)
           # >
           #   > get_var(1,2,100,matrix(c(1,0,0,1),nrow=2))

# sample <- 1
# PARENT_FOLDER = paste("/home/users/allstaff/lmcintosh/P2_LEON/ILO2.58-8359/","ascat_",as.character(sample),"/",sep="")PARENT_FOLDER
# load(paste(PARENT_FOLDER,"data.Rda",sep=""))
# load(paste(PARENT_FOLDER,"segments.Rda",sep=""))
for(pr in 1:16){
  name <- as.character(pr)
  if(pr%%4 == 0) {
    segments <- get_new_seg_estimates(data=data,segments=segments,
    new_snp_name = paste("segmentCT_pre",name,sep=""),
    new_seg_name = paste("segmentCT_pre",name,sep=""))
  }else{
    segments <- get_new_seg_estimates(data=data,segments=segments,
    new_snp_name = paste("CT_pre",name,sep=""),
    new_seg_name = paste("segmentCT_pre",name,sep=""))
  }
  print(pr)
}

get_TCN_track <- function(mat,seg_est_name){
  g <- ggplot(mat) +
    geom_segment(aes_string(x = "tcnStart", y = seg_est_name, xend = "tcnEnd", yend = seg_est_name),size=2)+
    facet_grid(.~chromosome,scales = "free_x", space = "free")
  return(g)
}

g1 <-get_TCN_track(old_segments[[3]],"segmentCT_pre0")+ylim(0,4)
g2 <-get_TCN_track(old_segments[[3]],"segmentCT_pre1")+ylim(0,4)
g3 <-get_TCN_track(old_segments[[3]],"segmentCT_pre2")+ylim(0,4)
g4 <-get_TCN_track(old_segments[[3]],"segmentCT_pre3")+ylim(0,4)
grid.arrange(g1,g2,g3,g4)

g1 <-get_TCN_track(old_segments[[3]],"segmentCT_pre4")+ylim(0,3)
g2 <-get_TCN_track(old_segments[[3]],"segmentCT_pre5")+ylim(0,3)
g3 <-get_TCN_track(old_segments[[3]],"segmentCT_pre6")+ylim(0,3)
g4 <-get_TCN_track(old_segments[[3]],"segmentCT_pre7")+ylim(0,3)
grid.arrange(g1,g2,g3,g4)

g1 <-get_TCN_track(old_segments[[4]],"segmentCT_pre8")+ylim(0,3)
g2 <-get_TCN_track(old_segments[[4]],"segmentCT_pre9")+ylim(0,3)
g3 <-get_TCN_track(old_segments[[4]],"segmentCT_pre10")+ylim(0,3)
g4 <-get_TCN_track(old_segments[[4]],"segmentCT_pre11")+ylim(0,3)
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
                                      new_snp_name = paste("segmentCT_pre",name,sep=""),
                                      new_seg_name = paste("segmentCT_pre",name,sep=""))
  }else{
    segments <- get_new_seg_estimates(data=data,segments=segments,
                                      new_snp_name = paste("CT_pre",name,sep=""),
                                      new_seg_name = paste("segmentCT_pre",name,sep=""))
  }
}

> for(pr in 1:16){
  +   name <- as.character(pr)
  +   if(pr%%4 == 0) {
    +     segments <- get_new_seg_estimates(data=data,segments=segments,
                                            +                                       new_snp_name = paste("segmentCT_pre",name,sep=""),
                                            +                                       new_seg_name = paste("segmentCT_pre",name,sep=""))
    +   }else{
      +   segments <- get_new_seg_estimates(data=data,segments=segments,
                                            +                                     new_snp_name = paste("CT_pre",name,sep=""),
                                            +                                     new_seg_name = paste("segmentCT_pre",name,sep=""))
      +   }
  +   print(pr)
  + }

g1 <-get_TCN_track(segments,"segmentCT_pre1")+ylim(0,4)
g2 <-get_TCN_track(segments,"segmentCT_pre2")+ylim(0,4)
g3 <-get_TCN_track(segments,"segmentCT_pre3")+ylim(0,4)
g4 <-get_TCN_track(segments,"segmentCT_pre4")+ylim(0,4)
grid.arrange(g1,g2,g3,g4)

g1 <-get_TCN_track(segments,"segmentCT_pre5")+ylim(0,4)
g2 <-get_TCN_track(segments,"segmentCT_pre6")+ylim(0,4)
g3 <-get_TCN_track(segments,"segmentCT_pre7")+ylim(0,4)
g4 <-get_TCN_track(segments,"segmentCT_pre8")+ylim(0,4)
grid.arrange(g1,g2,g3,g4)

g1 <-get_TCN_track(segments,"segmentCT_pre9")+ylim(0,4)
g2 <-get_TCN_track(segments,"segmentCT_pre10")+ylim(0,4)
g3 <-get_TCN_track(segments,"segmentCT_pre11")+ylim(0,4)
g4 <-get_TCN_track(segments,"segmentCT_pre12")+ylim(0,4)
grid.arrange(g1,g2,g3,g4)

g1 <-get_TCN_track(segments,"segmentCT_pre13")+ylim(0,4)
g2 <-get_TCN_track(segments,"segmentCT_pre14")+ylim(0,4)
g3 <-get_TCN_track(segments,"segmentCT_pre15")+ylim(0,4)
g4 <-get_TCN_track(segments,"segmentCT_pre16")+ylim(0,4)
grid.arrange(g1,g2,g3,g4)

renormalise
function(data,round){
  # should renormalise the homs and the hets seperately.
  if(round == 0){
    homs <- data$Allele1...Design == data$Allele1...Design
    med_homs <- median(data[which(homs),"reddye"]+data[which(homs),"greendye"],na.rm=TRUE)/2
    med_hets <- median(data[which(!homs),"reddye"]+data[which(!homs),"greendye"],na.rm=TRUE)/2
    data[which(homs),"reddye"] <- data[which(homs),"reddye"]/median(data[which(homs),"reddye"],na.rm=TRUE)*med_homs
    data[which(homs),"greendye"] <- data[which(homs),"greendye"]/median(data[which(homs),"greendye"],na.rm=TRUE)*med_homs
    data[which(!homs),"reddye"] <- data[which(!homs),"reddye"]/median(data[which(!homs),"reddye"],na.rm=TRUE)*med_hets
    data[which(!homs),"greendye"] <- data[which(!homs),"greendye"]/median(data[which(!homs),"greendye"],na.rm=TRUE)*med_hets
    data$CNsnp <- data$reddye+data$greendye

    data$reddyeAT <- data$reddye*(as.numeric(data$AT1))
    data$reddyeGC <- data$reddye*(1-as.numeric(data$AT1))
    data$greendyeAT <- data$greendye*(as.numeric(data$AT2))
    data$greeendyeGC <- data$greendye*(1-as.numeric(data$AT2))
  }
  data$R <- data$CNsnp/data[,paste("segmentCT_pre",as.character(round),sep='')]- median(data$CNsnp/data[,paste("segmentCT_pre",as.character(round),sep='')],na.rm=TRUE)
  data[,"reddyeAT_spec"] <- data$reddye*(as.numeric(data$AT1))/data[,paste("segmentCT_pre",as.character(round),sep='')]
  data[,"reddyeGC_spec"] <- data$reddye*(1-as.numeric(data$AT1))/data[,paste("segmentCT_pre",as.character(round),sep='')]
  data[,"greendyeAT_spec"] <- data$greendye*(as.numeric(data$AT2))/data[,paste("segmentCT_pre",as.character(round),sep='')]
  data[,"greeendyeGC_spec"] <- data$greendye*(1-as.numeric(data$AT2))/data[,paste("segmentCT_pre",as.character(round),sep='')]
  print("done")
  return(data)
}
