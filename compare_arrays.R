rm(list=ls())
require(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(manipulate)
library(cluster)
library(igraph)
library(knitr)

# the metadata
metadata3=read.table("/home/users/allstaff/lmcintosh/P2_LEON/clare-vincent/metadata/clonal_SNP_aggregated_textwrangeler.txt",sep="\t",header=TRUE)
metadata3$Server.Location <- gsub("O2_",'O2.',gsub(".txt","/",gsub("(_FinalReport)", "/ascat_", metadata3$Report.data.filename, perl=TRUE),perl=TRUE),perl=TRUE)
metadata3$Parent.Row <- match(metadata3$Source,metadata3$Sample.ID)
head(metadata3)

#files_without_parents <- is.na(metadata3$Parent.Row)
metadata3$is_root <- as.character(metadata3$Sample.ID) == as.character(metadata3$Source)

BASE="/home/users/allstaff/lmcintosh/run_pres/"
metadata3$Segmented.File.Exists <- file.exists(paste(BASE,metadata3$Server.Location,"raw_data.txt.segmented",sep=""))

# print the file names that need to be re run!
metadata3[which(!metadata3$Segmented.File.Exists),"Server.Location"]
metadata3[which(metadata3$Segmented.File.Exists),"Server.Location"]

datalist = list()
for(i in 1:nrow(metadata3)){
  if(metadata3[i,"Segmented.File.Exists"]){
    datalist[[i]] <-readRDS(file=paste(BASE,metadata3[i,"Server.Location"],"raw_data.txt.segmented" ,sep=""))
    datalist[[i]]$chr <- factor(datalist[[i]]$chr, levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y'))
    datalist[[i]]$TCN <- datalist[[i]]$tcnMean
    datalist[[i]]$TCN_new <- datalist[[i]]$segmentCN_pre_clustered13
    datalist[[i]]$start <- datalist[[i]]$tcnStart
    datalist[[i]]$end <- datalist[[i]]$tcnEnd

  }else{
    datalist[[i]] <- NULL
  }
}
names(datalist) <- metadata3$Sample.ID

setwd("/wehisan/home/allstaff/l/lmcintosh/SAP")
source(paste(getwd(),"/utils.R",sep=""))
centromeres <- load_centromeres()

i<- i+1
print_pair(i)




metadata3$diff25 <- rep(NA,nrow(metadata3))
metadata3$diff10 <- rep(NA,nrow(metadata3))

head(datalist[[1]])


for(i in 1:nrow(metadata3)){
  print(as.numeric(i)/nrow(metadata3))
  if(!metadata3[i,"is_root"] && metadata3[i,"Segmented.File.Exists"] && metadata3[metadata3[i,"Parent.Row"],"Segmented.File.Exists"]){
    if(metadata3[metadata3[i,"Parent.Row"],"Sample.ID"] %in% names(datalist)){
      metadata3[i,"diff10"] <- get_dual_TCN_track_seperated(datalist=datalist,sample1=as.character(metadata3[i,"Sample.ID"]),sample2=as.character(metadata3[i,"Source"]),thresh=0.1,segdiff=10^4,min_length=20,make_plots=FALSE)
      #function(datalist,sample1,sample2,thresh,segdiff,min_length,make_plots){
      metadata3[i,"diff25"] <- get_dual_TCN_track_seperated(datalist,as.character(metadata3[i,"Sample.ID"]),as.character(metadata3[i,"Source"]),thresh=0.1,segdiff=10^4,min_length=20,make_plots=FALSE)
    } else{
      print(paste("ERROR",as.character(metadata3[i,"Source"])))
    }
  }
}

#PUT PDX.GENERATION IN HERE
ggplot(data=metadata3[which(metadata3$Patient.ID =="M1310"),],aes(x=Patient.tumor.kind,y=diff10,position="jitter",label=Patient.ID,col=Patient.ID))+geom_text(position = position_dodge(.5))
ggplot(data=metadata3,aes(x=Patient.tumor.kind,y=diff25,position="jitter",label=Patient.ID,col=Patient.ID))+geom_text(position = position_dodge(.5))


metadata3$proportion_diff <- metadata3$num_segs_s1 <- metadata3$num_segs_s2 <- metadata3$change_in_num_segs <- metadata3$boths <- metadata3$starts <- metadata3$ends <- rep(NA,nrow(metadata3))
for(i in which(metadata3$Patient.ID == "M1310")){
  if(!metadata3[i,"is_root"] && metadata3[i,"Corrected.File.Exists"] && metadata3[metadata3[i,"Parent.Row"],"Corrected.File.Exists"]){
    if(metadata3[metadata3[i,"Parent.Row"],"Sample.ID"] %in% names(datalist)){
      print(i)
      metadata3[i,c("proportion_diff","num_segs_s1","num_segs_s2","change_in_num_segs","boths","starts","ends","ors")]  <- unlist(get_dual_CN_track_seperated(datalist,as.character(metadata3[i,"Sample.ID"]),as.character(metadata3[i,"Source"]),0.5, 7*10^7,min_length=100,make_plots=FALSE))
      #metadata3[i,"diff25"] <- get_dual_CN_track_seperated(datalist,as.character(metadata3[i,"Sample.ID"]),as.character(metadata3[i,"Source"]),0.25,make_plot=FALSE)
    } else{
      print(paste("ERROR",as.character(metadata3[i,"Source"])))
    }
  }
}

levels(metadata3$Patient.tumor.kind) <- c("primary","mets")
ggplot(data=metadata3[which(metadata3$Patient.ID == "M1310"),],aes(x=Patient.tumor.kind,y=boths/num_segs_s1,position="jitter",label=Patient.ID,col=Patient.ID))+labs(x="type",y="proportion of segments where both boundaries are consistent across generations")+geom_text()
ggplot(data=metadata3[which(metadata3$Patient.ID == "M1310"),],aes(x=Patient.tumor.kind,y=ors/num_segs_s1,position="jitter",label=Patient.ID,col=Patient.ID))+labs(x="type",y="proportion of segments where either of the boundaries are consistent across generations")+geom_text()


get_grid_plot <- function(dat2,MIN_LENGTH){
  g <- ggplot(dat2[which(dat2$length>MIN_LENGTH),],aes(major,minor))+geom_point(aes(col=chr,size=length))+geom_path(aes(col=chr))+coord_fixed()+xlab("CN 1")+ylab("CN 2")+ggtitle("Grid plot")
  return(g)
}
#
# head(datalist[[99]])

get_grid_plot(datalist[[1]],0)
get_grid_plot(datalist[[98]],0)
get_grid_plot(datalist[[99]],10)
get_grid_plot(datalist[[100]],10)
get_grid_plot(datalist[[99]],10)
get_grid_plot(datalist[[99]],10)
get_grid_plot(datalist[[99]],10)
get_grid_plot(datalist[[103]],10)
get_grid_plot(datalist[[105]],10)
get_grid_plot(datalist[[107]],10)
get_grid_plot(datalist[[110]],10)
get_grid_plot(datalist[[112]],10)

ggplot(data=metadata3[which(metadata3$Patient.ID == "M1310"),],aes(x=Patient.tumor.kind,y=proportion_diff,position="jitter",label=Patient.ID,col=Patient.ID))+geom_text()+labs(x="type",y="proportion of genome with inconsistent allele specific copy number between generations")


ggplot(data=metadata3[which(metadata3$Patient.ID == "M1310"),],aes(x=Patient.tumor.kind,y=change_in_num_segs,position="jitter",label=Patient.ID,col=Patient.ID))+geom_text(position = position_dodge(.5))

ggplot(data=metadata3[which(metadata3$Patient.ID == "M1310"),],aes(x=Patient.tumor.kind,y=boths/num_segs_s1,position="jitter",label=Patient.ID,col=Patient.ID))+labs(x="type",y="proportion of segment boundaries consistent across generations")+geom_text()

ggplot(data=metadata3[which(metadata3$Patient.ID == "M1310"),],aes(x=Patient.tumor.kind,y=starts/num_segs_s1,position="jitter",label=Patient.ID,col=Patient.ID))+labs(x="type",y="proportion of the begining of segment boundaries consistent across generations")+geom_text()
ggplot(data=metadata3[which(metadata3$Patient.ID == "M1310"),],aes(x=Patient.tumor.kind,y=ends/num_segs_s1,position="jitter",label=Patient.ID,col=Patient.ID))+labs(x="type",y="proportion of segment boundaries consistent across generations")+geom_text()
```
