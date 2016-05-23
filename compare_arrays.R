---
  title: "Compare_SNP_arrays"
output:
  html_document:
  pandoc_args:
  - +RTS
- -K1024m
- -RTS
pdf_document: default
---
  # Preprocessing
  ```{r, echo=TRUE,warnings=FALSE}
rm(list=ls())
require(data.table) # useful for comparing intervals
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
    datalist[[i]] <- unphase(data.frame(read.table(paste(BASE,metadata3[i,"Server.Location"],"raw_data.txt.segmented" ,sep=""),sep=",")))
    datalist[[i]]$chr <- factor(datalist[[i]]$chr, levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y'))
  }else{
    datalist[[i]] <- NULL
  }
}
names(datalist) <- metadata3$Sample.ID

file = "centromeres.txt"
BASE = '/home/users/allstaff/lmcintosh/ITH/SNP_stuff/data/'
centromeres <- read.table(paste(BASE,file,sep=''), sep="\t", header=FALSE)
colnames(centromeres) <- c("chr","x")
centromeres$x <- as.numeric(centromeres$x)*10^6 # is in megabase format.
#
similarities_in_boundaries <- function(s1,s2,which_values = "both",max_dist,min_length){
  ors <- starts <- ends <- boths <- 0
  for(i in 1:nrow(s1)){
    if(any(s2$chr == s1[i,"chr"])){
      closeststart = which(abs(s1[i,"start"]-s2[which(s2$chr == s1[i,"chr"]),"start"]) == min(abs(s1[i,"start"]-s2[which(s2$chr == s1[i,"chr"]),"start"])))
      closestend = which(abs(s1[i,"end"]-s2[which(s2$chr == s1[i,"chr"]),"end"]) == min(abs(s1[i,"end"]-s2[which(s2$chr == s1[i,"chr"]),"end"])))
      #closeststart
      #closestend

      if(abs(s1[i,"start"] - s1[closeststart,"start"]) < max_dist){
        starts <- starts + 1
        #print(c(i,closeststart))
      }

      if(abs(s1[i,"end"] - s1[closestend,"end"]) < max_dist){
        ends <- ends + 1
        #print(c(i,closestend))
      }

      if(abs(s1[i,"start"] - s1[closeststart,"start"]) < max_dist | abs(s1[i,"end"] - s1[closestend,"end"]) < max_dist){
        ors <- ors + 1
        #print(c(i,closeststart,closestend))
      }

      if(abs(s1[i,"start"] - s1[closeststart,"start"]) < max_dist & abs(s1[i,"end"] - s1[closestend,"end"]) < max_dist){
        boths <- boths + 1
        #print(c(i,closeststart,closestend))
      }
    }
  }
  if(which_values == "both")(return(boths))
  else if(which_values == "start")(return(starts))
  else if(which_values == "end")(return(ends))
  else if(which_values == "or")(return(ors))
}

get_dual_CN_track_seperated <- function(datalist,sample1,sample2,thresh,segdiff,min_length,make_plots){
  overall <- NULL
  for(chr in levels(datalist[[sample1]][,"chr"])){
    if(sum(datalist[[sample1]][,"chr"] == chr) > 0 && sum(datalist[[sample1]][,"chr"] == chr) > 0){
      # the intervals for sample 2
      intervals2 <- datalist[[sample1]][which(datalist[[sample1]][,"chr"] == chr),c("start","end","major","minor","chr")]
      names(intervals2) <- c("start","end","major2","minor2","chr")
      # the intervals for sample 1
      intervals1 <- datalist[[sample2]][which(datalist[[sample2]][,"chr"] == chr),c("start","end","major","minor","chr")]
      names(intervals1) <- c("start","end","major1","minor1","chr")

      setDT(intervals1)  ## convert to data.table without copy
      setDT(intervals2)

      setkey(intervals2, "start", "end")
      ans = foverlaps(intervals1, intervals2, type="any")
      ans = ans[, `:=`(start = pmax(start, i.start), end = pmin(end, i.end))]
      ans = ans[, `:=`(i.start=NULL, i.end=NULL)][start <= end]
      ans$majordiff <- abs(ans$major1-ans$major2)
      ans$minordiff <- abs(ans$minor1-ans$minor2)
      ans$length <- ans$end-ans$start
      if(chr == "1"){
        overall = ans
      }else{
        overall <- rbind(overall,ans)
      }
    }
  }
  overall <- data.frame(overall)

  # find the regions of the genome where the average copy number is different
  differentmajor <- overall[which(abs(overall$major1 - overall$major2) > thresh),]
  differentminor <- overall[which(abs(overall$minor1 - overall$minor2) > thresh),]
  differenteither <- overall[which(abs(overall$minor1 - overall$minor2) > thresh | abs(overall$major1 - overall$major2) > thresh),]

  # find the proporiton of the genome that is different
  proportion_diff <- sum(as.numeric(differenteither$length))/sum(as.numeric(overall$length))
  proportion_diff_major <- sum(as.numeric(differentmajor$length))/sum(as.numeric(overall$length))
  proportion_diff_minor <- sum(as.numeric(differentminor$length))/sum(as.numeric(overall$length))


  s1 <- datalist[[sample1]][,c("chr","start","end")]
  s2 <- datalist[[sample2]][,c("chr","start","end")]


  s1 <- s1[which(s1$end - s1$start > min_length),]
  s2 <- s2[which(s2$end - s2$start > min_length),]

  a <- abs_diff_in_approx_num_segments <- abs(nrow(s2) - nrow(s1))
  boths <- similarities_in_boundaries(s1,s2,"both",segdiff,min_length)
  starts <- similarities_in_boundaries(s1,s2,"start",segdiff,min_length)
  ends <- similarities_in_boundaries(s1,s2,"end",segdiff,min_length)
  ors <- similarities_in_boundaries(s1,s2,"or",segdiff,min_length)

  # calculate the maximum average copy number for plotting
  majormax <- 4
  #majormax <- max(c(overall$major1,overall$major2,overall$minor1,overall$minor2))+0.1
  #majormax <- majormax
  if(make_plots){
    g1 <- ggplot() +
      geom_rect(data=differentmajor, aes(xmin=start,xmax=end,ymin=0,ymax=4),alpha=0.2)+
      geom_segment(data=overall,aes_string(x = "start", y = "major1", xend = "end", yend = "major1"),size=2,alpha=0.5)+
      geom_segment(data=overall,aes_string(x = "start", y = "major2", xend = "end", yend = "major2"),size=1,alpha=0.5)+
      geom_vline(data=centromeres,aes(xintercept=x,col=chr))+
      facet_grid(.~chr,scales = "free_x", space = "free")+
      scale_x_continuous(breaks=seq(0,3*10^9,50*10^6))+ #make 50mb ticks...
      theme(axis.text.x = element_blank(),legend.position="none",legend.background = element_rect(fill = "white"))+
      #theme_bw()+
      xlab("50 MB ticks")+ylim(0,4)+
      #scale_y_continuous(breaks=number_ticks(20))+
      ylab("major CN")+ggtitle(paste("CN track, difference = ",proportion_diff_major,sep=""))

    g2 <- ggplot() +
      geom_rect(data=differentminor, aes(xmin=start,xmax=end,ymin=0,ymax=4),alpha=0.2)+
      geom_segment(data=overall,aes_string(x = "start", y = "minor1", xend = "end", yend = "minor1",position = position_jitter(w = 0.1, h = 0.1)),size=2,alpha=0.5)+
      geom_segment(data=overall,aes_string(x = "start", y = "minor2", xend = "end", yend = "minor2"),size=1,alpha=0.5)+
      geom_vline(data=centromeres,aes(xintercept=x,col=chr))+
      facet_grid(.~chr,scales = "free_x", space = "free")+
      #     scale_x_continuous(breaks=seq(0,3*10^9,50*10^6))+ #make 50mb ticks...
      theme(axis.text.x = element_blank(),legend.position="none",legend.background = element_rect(fill = "white"))+
      #theme_bw()+
      xlab("50 MB ticks")+ylim(0,4)+
      #scale_y_continuous(breaks=number_ticks(20))+
      ylab("minor CN")+ggtitle(paste("CN track, difference = ",proportion_diff_minor,sep=""))
    g <- grid.arrange(g1,g2,nrow=2,top=paste("Overall proportion different = ",proportion_diff,sep=""))
    print(g)

  }

  return(list(proportion_diff=proportion_diff, num_segs_s1 = nrow(s1), num_segs_s2 = nrow(s2), change_in_num_segs = nrow(s1)-nrow(s2), boths=boths, starts = starts, ends= ends, ors = ors))
}
```

```{r}
metadata3$diff25 <- rep(NA,nrow(metadata3))
metadata3$diff10 <- rep(NA,nrow(metadata3))


for(i in 1:nrow(metadata3)){
  print(as.numeric(i)/nrow(metadata3))
  if(!metadata3[i,"is_root"] && metadata3[i,"Corrected.File.Exists"] && metadata3[metadata3[i,"Parent.Row"],"Corrected.File.Exists"]){
    if(metadata3[metadata3[i,"Parent.Row"],"Sample.ID"] %in% names(datalist)){
      metadata3[i,"diff10"] <- get_dual_CN_track_seperated(datalist,as.character(metadata3[i,"Sample.ID"]),as.character(metadata3[i,"Source"]),thresh=0.1,segdiff=10^4,min_length=20,make_plots=FALSE)
      #function(datalist,sample1,sample2,thresh,segdiff,min_length,make_plots){
      metadata3[i,"diff25"] <- get_dual_CN_track_seperated(datalist,as.character(metadata3[i,"Sample.ID"]),as.character(metadata3[i,"Source"]),thresh=0.1,segdiff=10^4,min_length=20,make_plots=FALSE)
    } else{
      print(paste("ERROR",as.character(metadata3[i,"Source"])))
    }
  }
}

#levels(metadata3$Patient.tumor.kind) <- c("primary","secondary")

PUT PDX.GENERATION IN HERE
ggplot(data=metadata3[which(metadata3$Patient.ID =="M1310"),],aes(x=Patient.tumor.kind,y=diff10,position="jitter",label=Patient.ID,col=Patient.ID))+geom_text(position = position_dodge(.5))
ggplot(data=metadata3,aes(x=Patient.tumor.kind,y=diff25,position="jitter",label=Patient.ID,col=Patient.ID))+geom_text(position = position_dodge(.5))
```

# FOR THE GRANT PROPOSAL!
```{r}
#metadata3 <- metadata3[which(metadata3$Patient.ID == "M1310"),]
#metadata3$diff25 <- rep(NA,nrow(metadata3))
#metadata3$num_segments <- rep(NA,nrow(metadata3))
#metadata3$num_novel_segment_boundaries <- rep(NA,nrow(metadata3))
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

get_grid_plot(datalist[[97]],0)
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
