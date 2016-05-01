
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
[1] 1
[1] 2
[1] 3
[1] 4
[1] 5
[1] 6
[1] 7
[1] 8
[1] 9
[1] 10
[1] 11
[1] 12
[1] 13
[1] 14
[1] 15
[1] 16
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

#############################
  > ##### SNP array package #####
> #############################
  > # Author: Lachlan McIntosh
  > rm(list=ls())
>
  > library(ggplot2)
>
  > #### unpaired test dataset: SAMPLE 13 ####
> args=(commandArgs(TRUE))
> sample <-  as.character(args[1])
> sample <- 9
> print(sample)
[1] 9
> PARENT_FOLDER = paste("/home/users/allstaff/lmcintosh/P2_LEON/ILO2.58-8359/","ascat_",as.character(sample),"/",sep="")
> filename = paste(PARENT_FOLDER,"/raw_data.txt",sep="")
Error in names(frame)[names(frame) == "x"] <- name :
  names() applied to a non-vector
> # start
  > # load functions
  > source("/home/users/allstaff/lmcintosh/ITH/SNP_stuff/data/MAIN_accessory.R")
>
  >
  > library(data.table)
> data <- read_raw_illumina_file(filename)
Read 2379855 rows and 30 (of 30) columns from 0.299 GB file in 00:00:06
Warning message:
  In fread(filename, sep = "\t", sep2 = "auto", nrows = -1L, header = TRUE,  :
             Bumped column 6 to type character on data row 18, field contains 'Y'. Coercing previously read values in this column from logical, integer or numeric back to character which may not be lossless; e.g., if '00' and '000' occurred before they will now be just '0', and there may be inconsistencies with treatment of ',,' and ',NA,' too (if they occurred in this column before the bump). If this matters please rerun and set 'colClasses' to 'character' for this column. Please note that column type detection uses the first 5 rows, the middle 5 rows and the last 5 rows, so hopefully this message should be very rare. If reporting to datatable-help, please rerun and include the output from verbose=TRUE.
           > # rename columns to allow for the use of '$'
             > data <- rename_raw_columns(data)
           Warning messages:
             1: In `names<-.data.table`(`*tmp*`, value = value) :
             The colnames(x)<-value syntax copies the whole table. This is due to <- in R itself. Please change to setnames(x,old,new) which does not copy and is faster. See help('setnames'). You can safely ignore this warning if it is inconvenient to change right now. Setting options(warn=2) turns this warning into an error, so you can then use traceback() to find and change your colnames<- calls.
           2: In `names<-.data.table`(`*tmp*`, value = value) :
             The colnames(x)<-value syntax copies the whole table. This is due to <- in R itself. Please change to setnames(x,old,new) which does not copy and is faster. See help('setnames'). You can safely ignore this warning if it is inconvenient to change right now. Setting options(warn=2) turns this warning into an error, so you can then use traceback() to find and change your colnames<- calls.
           > data <- preprocess_raw_data(data)
           Warning message:
             In preprocess_raw_data(data) : NAs introduced by coercion
           >
             > # create a datastructure for PSCBS
             > library(PSCBS)
           > result <- do_some_seg(data,0)
           There were 18 warnings (use warnings() to see them)
           > segments <- result$segs
           > fit <- result$fit
           > data <- result$data
           >
             > data <-  reformat(data,segments)
           [1] "done"
           >
             > #### find het SNPs if the SNPs are not annotated ####
           > # not required for the illumina data
             > # find_all_het_SNPs(data)
             > #### find het SNPs if the SNPs are not annotated ####
           > data <-  renormalise(data,0)#, simplify = TRUE)
           [1] "done"
           >
             > #library(mgcv)
             > folder ="/home/users/allstaff/lmcintosh/ITH/SNP_stuff/data/"
           > GC_filenames <- lapply(1:50*50,function(x) paste(folder,"intervals_",as.character(x),".out",sep=""))
           >
             > pGC <- lapply(GC_filenames, read_GC_file) # parallel was slower for loading files.
           Read 2238691 rows and 12 (of 12) columns from 0.125 GB file in 00:00:03
           >
             > for (i in seq_along(pGC)) {
               +   setnames(pGC[[i]],c("chr","start","end","pAT","pGC","A","C","G","T","N","O","length"))
               +   pGC[[i]]$pos <- i*25+pGC[[i]]$start
               +   pGC[[i]]$loc <- as.numeric(pGC[[i]]$chr)+as.numeric(pGC[[i]]$pos)/max(pGC[[i]]$pos)
               + }
           >
             > GC <- pGC[[1]]
           > GC$Position <- GC$start + GC$length/2
           > GC[,"GC50"] <- GC$pGC
           > for(i in 2:50){
             +   GC[,paste("GC",as.character(50*i),sep="")] <- pGC[[i]]$pGC
             + }
           > head(GC)
           chr    start      end  pAT  pGC  A  C  G  T N O length      pos   loc Position GC50 GC100  GC150 GC200 GC250  GC300  GC350  GC400  GC450
           1   6 31239392 31239442 0.32 0.68  6 17 17 10 0 0     50 31239417 6.125 31239417 0.68  0.63 0.6667 0.705 0.712 0.6933 0.6943 0.7025 0.6756
           2   6 31323067 31323117 0.40 0.60  8 20 10 12 0 0     50 31323092 6.126 31323092 0.60  0.60 0.5800 0.590 0.592 0.5900 0.5943 0.5900 0.5933
           3   6 31964178 31964228 0.48 0.52 13 19  7 11 0 0     50 31964203 6.128 31964203 0.52  0.55 0.5933 0.595 0.588 0.5900 0.5943 0.5950 0.5978
           4   6 32548609 32548659 0.48 0.52 12 10 16 12 0 0     50 32548634 6.131 32548634 0.52  0.48 0.4733 0.490 0.472 0.4667 0.4629 0.4475 0.4467
           5   6 32632732 32632782 0.48 0.52 10 13 13 14 0 0     50 32632757 6.131 32632757 0.52  0.63 0.6267 0.595 0.632 0.6500 0.6543 0.6675 0.6756
           6   6 32634380 32634430 0.50 0.50 19 11 14  6 0 0     50 32634405 6.131 32634405 0.50  0.51 0.4933 0.515 0.524 0.5333 0.5314 0.5425 0.5267
           GC500  GC550  GC600  GC650  GC700  GC750  GC800  GC850  GC900  GC950 GC1000 GC1050 GC1100 GC1150 GC1200 GC1250 GC1300 GC1350 GC1400
           1 0.686 0.6964 0.7083 0.7092 0.7129 0.7107 0.7063 0.7082 0.7067 0.7063  0.703 0.6990 0.6936 0.6887 0.6875 0.6832 0.6815 0.6793 0.6757
           2 0.594 0.5964 0.5850 0.5815 0.5871 0.5827 0.5800 0.5824 0.5733 0.5737  0.574 0.5695 0.5673 0.5643 0.5617 0.5656 0.5623 0.5600 0.5600
           3 0.592 0.5909 0.5917 0.5985 0.6014 0.6040 0.5975 0.5941 0.5933 0.5968  0.601 0.5962 0.5982 0.5983 0.6050 0.6048 0.6069 0.6081 0.6079
           4 0.452 0.4455 0.4350 0.4462 0.4514 0.4467 0.4475 0.4576 0.4622 0.4663  0.470 0.4771 0.4745 0.4730 0.4725 0.4696 0.4685 0.4704 0.4729
           5 0.680 0.6745 0.6650 0.6677 0.6686 0.6653 0.6613 0.6565 0.6533 0.6389  0.633 0.6305 0.6327 0.6296 0.6300 0.6264 0.6254 0.6170 0.6050
           6 0.514 0.5000 0.4933 0.4846 0.4786 0.4707 0.4637 0.4612 0.4522 0.4537  0.458 0.4524 0.4527 0.4443 0.4442 0.4440 0.4446 0.4400 0.4386
           GC1450 GC1500 GC1550 GC1600 GC1650 GC1700 GC1750 GC1800 GC1850 GC1900 GC1950 GC2000 GC2050 GC2100 GC2150 GC2200 GC2250 GC2300 GC2350
           1 0.6724 0.6653 0.6606 0.6569 0.6564 0.6553 0.6526 0.6489 0.6411 0.6379 0.6344 0.6300 0.6298 0.6262 0.6242 0.6227 0.6169 0.6126 0.6089
           2 0.5614 0.5593 0.5600 0.5600 0.5612 0.5612 0.5651 0.5683 0.5686 0.5705 0.5708 0.5710 0.5707 0.5700 0.5712 0.5732 0.5720 0.5743 0.5783
           3 0.6069 0.6080 0.6065 0.6062 0.6067 0.6041 0.6069 0.6061 0.6070 0.6089 0.6108 0.6100 0.6098 0.6110 0.6093 0.6050 0.6044 0.6013 0.6004
           4 0.4724 0.4680 0.4684 0.4744 0.4770 0.4818 0.4857 0.4894 0.4865 0.4900 0.4892 0.4910 0.4927 0.4943 0.4935 0.4941 0.4933 0.4913 0.4894
           5 0.5959 0.5900 0.5832 0.5744 0.5703 0.5635 0.5577 0.5517 0.5470 0.5421 0.5369 0.5335 0.5322 0.5276 0.5256 0.5200 0.5178 0.5143 0.5119
           6 0.4372 0.4340 0.4323 0.4300 0.4255 0.4218 0.4160 0.4144 0.4146 0.4132 0.4118 0.4105 0.4127 0.4119 0.4149 0.4150 0.4151 0.4143 0.4140
           GC2400 GC2450 GC2500
           1 0.6071 0.6053 0.6032
           2 0.5792 0.5776 0.5784
           3 0.5992 0.5976 0.5960
           4 0.4858 0.4833 0.4804
           5 0.5112 0.5090 0.5056
           6 0.4100 0.4127 0.4108
           > GC$name <- paste(as.character(GC[,c("chr")]),as.character(GC[,c("Position")]),sep=" ")
           > rm(pGC)
           >
             > # need to find a faster way to do this.
             > data$name <- paste(as.character(data[,c("Chr")]),as.character(data[,c("Position")]),sep=" ")
           > data <- merge(data,GC,by=c("name","name"))
           > data$Position <- data$Position.x
           >
             > library(mgcv)
           > temp_data <- data
           > data <- temp_data
           > segr <- 1
           > pr <- 0
           > old_segments <- list(segments)
           > while(segr < 5){ # could potentially do this until convergence.... but might just shrink all to zero.
             +   while(pr < 4*segr-1){
               +     old_segment_ids <- get_segment_ids(segments)
               +     data <-  renormalise(data,pr)
               +     #gam <- gam(R ~ 0+ reddyeAT_spec + greendyeAT_spec + reddyeGC_spec + greeendyeGC_spec + s(GC50)+s(GC150)+s(GC500)+s(GC2500),data=data)
                 +     gam <- gam(R ~ 0+ reddyeAT + greendyeAT + reddyeGC + greeendyeGC + s(GC50)+s(GC150)+s(GC500)+s(GC2500),data=data)
                 +     pr <- pr + 1
                 +     name <- as.character(pr)
                 +     print(name)
                 +     summary(gam)
                 +     print(name)
                 +     data <- get_new_snp_estimates(gam=gam,data=data,
                                                     +                                   new_snp_name = paste("CT_pre",name,sep=""),
                                                     +                                   old_segment_name = paste("segmentCT_pre",as.character(pr-1),sep=""))
                 +     segments <- get_new_seg_estimates(data=data,segments=segments,
                                                         +                                       new_snp_name = paste("CT_pre",name,sep=""),
                                                         +                                       new_seg_name = paste("segmentCT_pre",name,sep=""))
                 +     data <- put_seg_estimates_into_data_array(data=data,segments=segments,
                                                                 +                                               segment_name = paste("segmentCT_pre",name,sep=""),
                                                                 +                                               new_snp_name = paste("segmentCT_pre",name,sep=""))
                 +   }
             +   segr <- segr+1
             +   if(segr <5){
               +     old_segments[[segr]] <- segments
               +     pr <- pr + 1
               +     result <- do_some_seg(data,pr,paste("CT_pre",as.character(pr-1),sep=""))
               +     segments <- result$segs
               +     fit <- result$fit
               +     data <- result$data
               +   }
             + }
           [1] "done"
           [1] "1"
           [1] "1"
           [1] "done"
           [1] "2"
           [1] "2"
           [1] "done"
           [1] "3"
           [1] "3"
           [1] "done"
           [1] "5"
           [1] "5"
           [1] "done"
           [1] "6"
           [1] "6"
           [1] "done"
           [1] "7"
           [1] "7"
           Error: NA/NaN/Inf in foreign function call (arg 7)
           In addition: Warning messages:
             1: In min(x) : no non-missing arguments to min; returning Inf
           2: In max(x) : no non-missing arguments to max; returning -Inf
           3: In min(x) : no non-missing arguments to min; returning Inf
           4: In max(x) : no non-missing arguments to max; returning -Inf
           5: In min(x) : no non-missing arguments to min; returning Inf
           6: In max(x) : no non-missing arguments to max; returning -Inf
           > get_var(1,2,100,matrix(c(1,0,0,1),nrow=2))
           Error: could not find function "get_var"
           > get_var <- function(major,minor,n,sigma) sum(sigma)/(n*(major+minor)^2) - 2*(sigma[1,2]+sigma[2,2])/(n*minor*(major+minor))+sigma[2,2]/(n*minor^2)
           >
             > get_var(1,2,100,matrix(c(1,0,0,1),nrow=2))
           [1] 0.001389
           > pr
           [1] 8
           > head(segments)
           chromosome tcnId dhId tcnStart  tcnEnd tcnNbrOfLoci   tcnMean tcnNbrOfSNPs tcnNbrOfHets dhNbrOfLoci dhStart   dhEnd dhMean    c1Mean
           1          1     1    1   534247  838496           58    0.7477           37           37          37  534247  838496 0.7606   0.08949
           2          1     2    1   838496 1538320          400   -0.3344          236          236         236  838496 1538320 0.7151  -0.04764
           3          1     3    1  1538320 2338822          583    0.5735          433          433         433 1538320 2338822 0.6272   0.10691
           4          1     4    1  2338822 2339356            3 -127.1763            3            3           3 2338822 2339356 0.4629 -34.15531
           5          1     5    1  2339356 5952768         4124    0.5652         2974         2974        2974 2339356 5952768 0.6485   0.09934
           6          1     6    1  5952768 6936324          830    0.8848          562          562         562 5952768 6936324 0.6928   0.13588
           c2Mean locationstart segmentCT_pre4 segmentCT_pre5 segmentCT_pre6 segmentCT_pre7
           1   0.6582         1.002         0.7477       -0.08164          1.749          4.151
           2  -0.2868         1.003        -0.3344        5.74270         90.733         64.211
           3   0.4666         1.006         0.5735        0.92724         13.186         12.983
           4 -93.0210         1.009      -127.1763      955.30826      64006.548      29167.848
           5   0.4659         1.009         0.5652        1.00560         16.069         15.583
           6   0.7489         1.024         0.8848       -0.12732          1.776          4.877
           > result <- do_some_seg(data,pr,paste("CT_pre",as.character(pr-1),sep=""))
           Error: NA/NaN/Inf in foreign function call (arg 7)
           > seg_name = paste("CT_pre",as.character(pr-1),sep="")
           > dataPSCBS <- data[,c("Chr","Position","B.Allele.Freq",seg_name)] # don't know why I am stuck here
           > setnames(dataPSCBS, c("chromosome","x","betaT","CT"))
           > # remove incomplete data
             > dataPSCBS <- dataPSCBS[complete.cases(dataPSCBS),]
           > #### end loading of unpaired test dataset: SAMPLE 13 ####
           > ##### END INITIAL DATA PROCESSING #####
           > ##### PSCBS MAIN BLOCK #####
           > # PSCBS stands for parent specific circular binary segmentation.
             > # The PSCBS manual: https://cAB2an.r-project.org/web/packages/PSCBS/vignettes/PairedPSCBS.pdf
             > # remove and install the PSCBS package to ensure the most up to date version:
             > # remove.packages(c("PSCBS"))
             > # install.packages(c("PSCBS"))
             > # PSCBS BEGIN the PSCBS method as per the specification at
             > # https://cran.r-project.org/web/packages/PSCBS/PSCBS.pdf
             > # Henrik Bengtsson's utilities library:
             > verbose <- R.utils::Arguments$getVerbose(-10*interactive(), timestamp=TRUE)
           > # is this required for PSCBS?
             > # (OPTIONAL)
             > # Drop single-locus "TCN" outliers
             > dataPSCBSwso <- dropSegmentationOutliers(dataPSCBS)
           Error: NA/NaN/Inf in foreign function call (arg 7)
           > # drop large gaps (greater than 1 mega base inbetween observations)
             > gaps <- findLargeGaps(dataPSCBSwso, minLength = 1e+06)
           Error in findLargeGaps(dataPSCBSwso, minLength = 1e+06) :
             object 'dataPSCBSwso' not found
           > knownSegments <- gapsToSegments(gaps)
           Error in UseMethod("gapsToSegments") :
             no applicable method for 'gapsToSegments' applied to an object of class "c('standardGeneric', 'genericFunction', 'function', 'OptionalFunction', 'PossibleMethod', 'Ufunction', 'expressionORfunction', 'functionORNULL', 'optionalMethod')"
           > # A segmentation (OPTIONS)
             > # USE IF NON PAIRED NORMAL:
             > fit <- segmentByNonPairedPSCBS(dataPSCBSwso, preserveScale = FALSE, verbose =FALSE)
           Error in segmentByNonPairedPSCBS(dataPSCBSwso, preserveScale = FALSE,  :
                                              object 'dataPSCBSwso' not found
                                            > # USE IF PAIRED NORMAL (can skip tumor boost normalisation by setting tbn = FALSE)
                                              > dataPSCBS
                                            chromosome         x  betaT         CT
                                            1               10 100000920 0.9981  5.965e+00
                                            2               10 100003304 0.5124  5.148e+00
                                            3               10 100004116 0.0000  5.913e+00
                                            4               10

                                            > dataPSCBS[!complete.cases(dataPSCBS),]
                                            [1] chromosome x          betaT      CT
                                            <0 rows> (or 0-length row.names)
                                            > nrow(dataPSCBS)
                                            [1] 2211371
                                            > apply(dataPSCBS,2,max)
                                            chromosome          x      betaT         CT
                                            22  249213967          1     430442
                                            > apply(dataPSCBS,2,miin)
                                            Error in match.fun(FUN) : object 'miin' not found
                                            > apply(dataPSCBS,2,min)
                                            chromosome          x      betaT         CT
                                            1       2220          0     -13659
                                            > apply(dataPSCBS,2,median)
                                            chromosome          x      betaT         CT
                                            8.000e+00  7.032e+07  5.699e-01  2.396e+01
                                            > ggplot(dataPSCBS,x=CT)+geom_histogram()
                                            `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
                                            Warning message:
                                              Computation failed in `stat_bin()`:
                                              attempt to apply non-function
                                            > ggplot(dataPSCBS,x=CT)+geom_histogram(aes(CT))
                                            `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
                                            > ggplot(dataPSCBS[which(dataPSCBS$CT<100),])+geom_histogram(aes(CT))
                                            `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
                                            > ggplot(dataPSCBS[which(dataPSCBS$CT<100& dataPSCBS$CT>0),])+geom_histogram(aes(CT))
                                            `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
                                            > head(datas)
                                            Error in head(datas) : object 'datas' not found
                                            > head(data)
                                            name Sample.Index Sample.ID Sample.Name SNP.Index    SNP.Name Chr Position.x GT.Score GC.Score Allele1...AB Allele2...AB
                                            1 10 100000920            9     #D300          NA   1069194 kgp22831654  10  100000920   0.8500   0.8270            B            B
                                            2 10 100003304            9     #D300          NA   2360019   kgp613501  10  100003304   0.8478   0.8238            A            B
                                            3 10 100004116            9     #D300          NA   1032032 kgp21561049  10  100004116   0.7898   0.7997            A            A
                                            4 10 100004360            9     #D300          NA   1147329 kgp10159511  10  100004360   0.7708   0.7685            A            B
                                            5 10 100004906            9     #D300          NA   2175633  kgp7849652  10  100004906   0.8832   0.9198            B            B
                                            6 10 100006683            9     #D300          NA   1681523 kgp21741920  10  100006683   0.8778   0.8640            B            B
                                            Allele1...Top Allele2...Top Allele1...Forward Allele2...Forward Allele1...Design Allele2...Design Theta      R X.Raw Y.Raw     X     Y
                                            1             G             G                 C                 C                G                G 0.949 0.5659   498  2323 0.087 1.076
                                            2             A             G                 A                 G                T                C 0.624 0.3232  1085   743 0.231 0.345
                                            3             A             A                 A                 A                A                A 0.043 0.5774  5378   185 1.070 0.073
                                            4             A             G                 A                 G                A                G 0.643 0.3451  1518   993 0.273 0.434
                                            5             C             C                 C                 C                G                G 0.962 0.2682   306  1719 0.046 0.769
                                            6             G             G                 G                 G                G                G 0.960 0.4540   338  1941 0.055 0.873
                                            B.Allele.Freq Log.R.Ratio SNP.Aux   SNP ILMN.Strand Top.Genomic.Sequence Customer.Strand     CT location segment resegmentCT_PSCBS0
                                            1        0.9981     -0.1463       0 [A/G]         TOP                   NA             BOT 1.6632     10.4     105               1.05
                                            2        0.5124     -0.5815       0 [T/C]         BOT                   NA             TOP 0.9611     10.4     105               1.05
                                            3        0.0000     -0.3044       0 [A/G]         TOP                   NA             TOP 1.3628     10.4     105               1.05
                                            4        0.5312     -0.5399       0 [A/G]         TOP                   NA             TOP 1.0128     10.4     105               1.05
                                            5        0.9933     -0.5999       0 [T/G]         BOT                   NA             TOP 0.9390     10.4     105               1.05
                                            6        0.9984     -0.2878       0 [A/G]         TOP                   NA             TOP 1.3916     10.4     105               1.05
                                            segmentCT_pre0      M        m ishet segmentCT_PSCBS reddye greendye AT1 AT2 segmentCT2 segmentCT segmentLog.R.Ratio LOG_RATIO  LOG_DIFF
                                            1           1.05 1.6601 0.003160 FALSE            1.05 2.6180 0.006188   0   0      1.051     1.051            -0.5947    0.2460  0.448441
                                            2           1.05 0.4925 0.468622  TRUE            1.05 0.7766 0.917578   1   0      1.051     1.051            -0.5947    0.9777  0.013241
                                            3           1.05 0.0000 1.362773 FALSE            1.05 0.0000 2.668355   1   1      1.051     1.051            -0.5947    0.5118  0.290341
                                            4           1.05 0.5380 0.474805  TRUE            1.05 0.8485 0.929685   1   0      1.051     1.051            -0.5947    0.9078  0.054841
                                            5           1.05 0.9328 0.006292 FALSE            1.05 1.4710 0.012319   0   0      1.051     1.051            -0.5947    1.0087 -0.005159
                                            6           1.05 1.3894 0.002227 FALSE            1.05 2.1911 0.004360   0   0      1.051     1.051            -0.5947    0.4839  0.306941
                                            TCN_RATIO TCN_DIFF reddyeAT reddyeGC greendyeAT greeendyeGC CNsnp reddyeAT_spec reddyeGC_spec greendyeAT_spec greeendyeGC_spec chr
                                            1    1.5826  0.61226    2.618   0.0000   0.006188       0.000 2.624        0.6830        0.0000        0.001614           0.0000  10
                                            2    0.9145 -0.08990    1.553  -0.7766   0.917578       0.000 1.694        0.4052       -0.2026        0.239377           0.0000  10
                                            3    1.2967  0.31179    0.000   0.0000   5.336711      -2.668 2.668        0.0000        0.0000        1.392237          -0.6961  10
                                            4    0.9637 -0.03817    1.697  -0.8485   0.929685       0.000 1.778        0.4427       -0.2213        0.242535           0.0000  10
                                            5    0.8935 -0.11193    1.471   0.0000   0.012319       0.000 1.483        0.3838        0.0000        0.003214           0.0000  10
                                            6    1.3241  0.34060    2.191   0.0000   0.004360       0.000 2.195        0.5716        0.0000        0.001137           0.0000  10
                                            start       end  pAT  pGC  A  C  G  T N O length   pos  loc Position.y GC50 GC100  GC150 GC200 GC250  GC300  GC350  GC400  GC450
                                            1 100000895 100000945 0.62 0.38 16 11  8 15 0 0     50 1e+08 10.4      1e+08 0.38  0.37 0.3333 0.305 0.292 0.3033 0.2943 0.3025 0.3133
                                            2 100003279 100003329 0.62 0.38 12  9 10 19 0 0     50 1e+08 10.4      1e+08 0.38  0.32 0.2667 0.285 0.312 0.3100 0.3229 0.3300 0.3356
                                            3 100004091 100004141 0.48 0.52 11 10 16 13 0 0     50 1e+08 10.4      1e+08 0.52  0.50 0.4933 0.480 0.456 0.4467 0.4429 0.4400 0.4422
                                            4 100004335 100004385 0.44 0.56 14  5 23  8 0 0     50 1e+08 10.4      1e+08 0.56  0.46 0.5067 0.485 0.464 0.4600 0.4800 0.4850 0.4867
                                            5 100004881 100004931 0.60 0.40  9 16  4 21 0 0     50 1e+08 10.4      1e+08 0.40  0.44 0.4667 0.450 0.448 0.4333 0.4371 0.4350 0.4289
                                            6 100006658 100006708 0.50 0.50 15 10 15 10 0 0     50 1e+08 10.4      1e+08 0.50  0.55 0.5067 0.520 0.504 0.4967 0.4829 0.4725 0.4711
                                            GC500  GC550  GC600  GC650  GC700  GC750  GC800  GC850  GC900  GC950 GC1000 GC1050 GC1100 GC1150 GC1200 GC1250 GC1300 GC1350 GC1400
                                            1 0.324 0.3255 0.3200 0.3246 0.3200 0.3133 0.3212 0.3200 0.3267 0.3295  0.329 0.3276 0.3218 0.3183 0.3258 0.3304 0.3354 0.3378 0.3336
                                            2 0.332 0.3345 0.3367 0.3369 0.3314 0.3293 0.3275 0.3306 0.3333 0.3347  0.337 0.3390 0.3382 0.3443 0.3500 0.3544 0.3585 0.3615 0.3629
                                            3 0.456 0.4545 0.4500 0.4477 0.4443 0.4373 0.4437 0.4471 0.4411 0.4411  0.434 0.4295 0.4273 0.4200 0.4167 0.4176 0.4100 0.4059 0.4036
                                            4 0.490 0.4800 0.4750 0.4708 0.4629 0.4600 0.4562 0.4459 0.4422 0.4421  0.448 0.4467 0.4436 0.4409 0.4400 0.4392 0.4392 0.4378 0.4343
                                            5 0.426 0.4182 0.4183 0.4231 0.4300 0.4320 0.4363 0.4353 0.4356 0.4368  0.438 0.4410 0.4418 0.4426 0.4417 0.4432 0.4446 0.4422 0.4429
                                            6 0.468 0.4582 0.4517 0.4477 0.4471 0.4413 0.4437 0.4459 0.4544 0.4547  0.453 0.4495 0.4518 0.4530 0.4483 0.4504 0.4469 0.4459 0.4486
                                            GC1450 GC1500 GC1550 GC1600 GC1650 GC1700 GC1750 GC1800 GC1850 GC1900 GC1950 GC2000 GC2050 GC2100 GC2150 GC2200 GC2250 GC2300 GC2350
                                            1 0.3303 0.3320 0.3310 0.3312 0.3309 0.3288 0.3274 0.3311 0.3400 0.3468 0.3533 0.3590 0.3585 0.3610 0.3647 0.3686 0.3724 0.3770 0.3796
                                            2 0.3641 0.3653 0.3684 0.3681 0.3709 0.3718 0.3714 0.3717 0.3730 0.3721 0.3713 0.3735 0.3741 0.3762 0.3753 0.3773 0.3782 0.3787 0.3783
                                            3 0.4055 0.4060 0.4045 0.4031 0.4036 0.4029 0.3989 0.3989 0.3989 0.3995 0.4000 0.4035 0.4039 0.4019 0.4028 0.4041 0.4027 0.4043 0.4034
                                            4 0.4303 0.4287 0.4265 0.4256 0.4218 0.4218 0.4211 0.4194 0.4162 0.4126 0.4149 0.4120 0.4117 0.4129 0.4121 0.4105 0.4080 0.4048 0.4051
                                            5 0.4434 0.4440 0.4458 0.4462 0.4448 0.4465 0.4474 0.4467 0.4459 0.4442 0.4405 0.4360 0.4371 0.4386 0.4381 0.4359 0.4333 0.4309 0.4311
                                            6 0.4503 0.4507 0.4503 0.4444 0.4424 0.4418 0.4383 0.4400 0.4400 0.4442 0.4426 0.4430 0.4439 0.4424 0.4419 0.4395 0.4396 0.4404 0.4413
                                            GC2400 GC2450 GC2500  Position CT_pre1 segmentCT_pre1 CT_pre2 segmentCT_pre2 CT_pre3 segmentCT_pre3 resegmentCT_PSCBS4 segmentCT_pre4
                                            1 0.3812 0.3824 0.3820 100000920  1.9815         0.9944  1.5081         0.4953  1.5649         0.6066             0.6031         0.6031
                                            2 0.3779 0.3833 0.3868 100003304  0.8974         0.9944  0.3944         0.4953  0.5000         0.6066             0.6031         0.6031
                                            3 0.4012 0.4012 0.4012 100004116  1.5048         0.9944  0.9123         0.4953  1.1459         0.6066             0.6031         0.6031
                                            4 0.4062 0.4073 0.4092 100004360  1.1293         0.9944  0.6715         0.4953  0.7234         0.6066             0.6031         0.6031
                                            5 0.4279 0.4265 0.4272 100004906  0.8198         0.9944  0.4048         0.4953  0.4054         0.6066             0.6031         0.6031
                                            6 0.4412 0.4420 0.4444 100006683  1.6493         0.9944  1.2565         0.4953  1.2337         0.6066             0.6031         0.6031
                                            CT_pre5 segmentCT_pre5 CT_pre6 segmentCT_pre6 CT_pre7 segmentCT_pre7
                                            1  0.89340        0.09103   5.143          3.833   5.965          5.571
                                            2 -0.22919        0.09103   3.094          3.833   5.148          5.571
                                            3  0.16921        0.09103   4.653          3.833   5.913          5.571
                                            4  0.04792        0.09103   2.899          3.833   5.281          5.571
                                            5 -0.52095        0.09103   2.397          3.833   4.902          5.571
                                            6  0.39075        0.09103   3.307          3.833   5.473          5.571
                                            > renormalise
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
                                            > get_new_snp_estimates
                                            function(gam,data,new_snp_name,old_segment_name){
                                              if(!is.null(gam$na.action)){
                                                data <- data[-gam$na.action,]
                                              }
                                              data[,new_snp_name] <- (gam$residuals+1)*data[,old_segment_name]
                                              return(data)
                                            }
                                            > PARENT_FOLDER
                                            [1] "/home/users/allstaff/lmcintosh/P2_LEON/ILO2.58-8359/ascat_9/"
                                            > sample <- 1
                                            > PARENT_FOLDER = paste("/home/users/allstaff/lmcintosh/P2_LEON/ILO2.58-8359/","ascat_",as.character(sample),"/",sep="")
                                            > load(paste(PARENT_FOLDER,"data.Rda",sep=""))
                                            > load(paste(PARENT_FOLDER,"segments.Rda",sep=""))
                                            > for(pr in 1:16){
                                              +   name <- as.character(pr)
                                              +   if(pr%%4 == 0) {
                                                +     segments <- get_new_seg_estimates(data=data,segments=segments,
                                                                                        +                                       new_snp_name = paste("segmentCT_pre",name,sep=""),
                                                                                        +                                       new_seg_name = paste("segmentCT_pre",name,sep=""))
                                                +   }else{
                                                  +     segments <- get_new_seg_estimates(data=data,segments=segments,
                                                                                          +                                       new_snp_name = paste("CT_pre",name,sep=""),
                                                                                          +                                       new_seg_name = paste("segmentCT_pre",name,sep=""))
                                                  +   }
                                              +   print(pr)
                                              + }
                                            [1] 1
                                            [1] 2
                                            [1] 3
                                            [1] 4
                                            [1] 5
                                            [1] 6
                                            [1] 7
                                            [1] 8
                                            [1] 9
                                            [1] 10
                                            [1] 11
                                            [1] 12
                                            [1] 13
                                            [1] 14
                                            [1] 15
                                            [1] 16
                                            There were 50 or more warnings (use warnings() to see the first 50)
                                            > # plot old segments and the new ones!
                                              > get_TCN_track <- function(mat,seg_est_name){
                                                +   g <- ggplot(mat) +
                                                  +     geom_segment(aes_string(x = "tcnStart", y = seg_est_name, xend = "tcnEnd", yend = seg_est_name),size=2)+
                                                  +     facet_grid(.~chromosome,scales = "free_x", space = "free")
                                                +   return(g)
                                                + }
                                            >
                                              > #g1 <- get_TCN_track(old_segments,"tcnMean")+ylim(0,4)
                                              > g1 <-get_TCN_track(segments,"segmentCT_pre1")+ylim(0,4)
                                            > g2 <-get_TCN_track(segments,"segmentCT_pre2")+ylim(0,4)
                                            > g3 <-get_TCN_track(segments,"segmentCT_pre3")+ylim(0,4)
                                            > g4 <-get_TCN_track(segments,"segmentCT_pre4")+ylim(0,4)
                                            > grid.arrange(g1,g2,g3,g4)
                                            Warning messages:
                                              1: Removed 19 rows containing missing values (geom_segment).
                                            2: Removed 24 rows containing missing values (geom_segment).
                                            3: Removed 19 rows containing missing values (geom_segment).
                                            4: Removed 16 rows containing missing values (geom_segment).
                                            > load(paste(PARENT_FOLDER,"old_segments.Rda",sep=""))
                                            > head(old_segments[[1]])
                                            chromosome tcnId dhId tcnStart  tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs tcnNbrOfHets dhNbrOfLoci dhStart   dhEnd dhMean c1Mean c2Mean
                                            32          1     1    1    82154  838755           66   1.955           35           35          35   82154  838755 0.6958 0.2973  1.657
                                            33          1     2    1   838755 1564268          439   2.663          140          140         140  838755 1564268 0.8212 0.2381  2.425
                                            34          1     3    1  1564268 2338724          588   2.320          288          288         288 1564268 2338724 0.5583 0.5124  1.808
                                            35          1     4    1  2338724 2339356            4   4.424            1            1           1 2338724 2339356 0.1314 1.9213  2.503
                                            36          1     5    1  2339356 3084750          563   2.316          263          263         263 2339356 3084750 0.5710 0.4968  1.819
                                            37          1     6    1  3084750 3300535          295   2.028          153          153         153 3084750 3300535 0.8410 0.1612  1.867
                                            locationstart segmentCT_pre0
                                            32         1.000          1.955
                                            33         1.003          2.663
                                            34         1.006          2.320
                                            35         1.009          4.424
                                            36         1.009          2.316
                                            37         1.012          2.028
                                            > head(old_segments[[2]])
                                            chromosome tcnId dhId tcnStart  tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs tcnNbrOfHets dhNbrOfLoci dhStart   dhEnd dhMean c1Mean c2Mean
                                            32          1     1    1    82154  838755           66   1.955           35           35          35   82154  838755 0.6958 0.2973  1.657
                                            33          1     2    1   838755 1564268          439   2.663          140          140         140  838755 1564268 0.8212 0.2381  2.425
                                            34          1     3    1  1564268 2338724          588   2.320          288          288         288 1564268 2338724 0.5583 0.5124  1.808
                                            35          1     4    1  2338724 2339356            4   4.424            1            1           1 2338724 2339356 0.1314 1.9213  2.503
                                            36          1     5    1  2339356 3084750          563   2.316          263          263         263 2339356 3084750 0.5710 0.4968  1.819
                                            37          1     6    1  3084750 3300535          295   2.028          153          153         153 3084750 3300535 0.8410 0.1612  1.867
                                            locationstart segmentCT_pre0 segmentCT_pre1 segmentCT_pre2 segmentCT_pre3
                                            32         1.000          1.955          1.946          1.871          1.789
                                            33         1.003          2.663          2.390          2.294          2.201
                                            34         1.006          2.320          2.172          2.083          1.998
                                            35         1.009          4.424          2.777          3.333          2.799
                                            36         1.009          2.316          2.145          2.057          1.973
                                            37         1.012          2.028          1.964          1.887          1.810
                                            > head(old_segments[[3]])
                                            chromosome tcnId dhId tcnStart  tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs tcnNbrOfHets dhNbrOfLoci dhStart   dhEnd dhMean  c1Mean c2Mean
                                            1          1     1    1   534247  819209           37   1.684           19           19          19  534247  819209 0.6760 0.27291  1.411
                                            2          1     2    1   819209 1120144          221   2.042           82           82          82  819209 1120144 0.7875 0.21698  1.825
                                            3          1     3    1  1120144 1563654          215   2.325           57           57          57 1120144 1563654 0.8326 0.19459  2.131
                                            4          1     4    1  1563654 1970092          232   1.906          122          122         122 1563654 1970092 0.5017 0.47492  1.431
                                            5          1     5    1  1970092 1997873           29   2.451            8            8           8 1970092 1997873 0.9754 0.03008  2.421
                                            6          1     6    1  1997873 2115914          107   1.796           37           37          37 1997873 2115914 0.9465 0.04800  1.748
                                            locationstart segmentCT_pre4 segmentCT_pre5 segmentCT_pre6 segmentCT_pre7
                                            1         1.002          1.684          1.624          1.536          1.452
                                            2         1.003          2.042          1.973          1.864          1.767
                                            3         1.004          2.325          2.195          2.089          1.974
                                            4         1.006          1.906          1.832          1.736          1.642
                                            5         1.008          2.451          2.216          2.142          2.011
                                            6         1.008          1.796          1.740          1.647          1.559
                                            > head(old_segments[[2]])
                                            chromosome tcnId dhId tcnStart  tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs tcnNbrOfHets dhNbrOfLoci dhStart   dhEnd dhMean c1Mean c2Mean
                                            32          1     1    1    82154  838755           66   1.955           35           35          35   82154  838755 0.6958 0.2973  1.657
                                            33          1     2    1   838755 1564268          439   2.663          140          140         140  838755 1564268 0.8212 0.2381  2.425
                                            34          1     3    1  1564268 2338724          588   2.320          288          288         288 1564268 2338724 0.5583 0.5124  1.808
                                            35          1     4    1  2338724 2339356            4   4.424            1            1           1 2338724 2339356 0.1314 1.9213  2.503
                                            36          1     5    1  2339356 3084750          563   2.316          263          263         263 2339356 3084750 0.5710 0.4968  1.819
                                            37          1     6    1  3084750 3300535          295   2.028          153          153         153 3084750 3300535 0.8410 0.1612  1.867
                                            locationstart segmentCT_pre0 segmentCT_pre1 segmentCT_pre2 segmentCT_pre3
                                            32         1.000          1.955          1.946          1.871          1.789
                                            33         1.003          2.663          2.390          2.294          2.201
                                            34         1.006          2.320          2.172          2.083          1.998
                                            35         1.009          4.424          2.777          3.333          2.799
                                            36         1.009          2.316          2.145          2.057          1.973
                                            37         1.012          2.028          1.964          1.887          1.810
                                            > head(old_segments[[1]])
                                            chromosome tcnId dhId tcnStart  tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs tcnNbrOfHets dhNbrOfLoci dhStart   dhEnd dhMean c1Mean c2Mean
                                            32          1     1    1    82154  838755           66   1.955           35           35          35   82154  838755 0.6958 0.2973  1.657
                                            33          1     2    1   838755 1564268          439   2.663          140          140         140  838755 1564268 0.8212 0.2381  2.425
                                            34          1     3    1  1564268 2338724          588   2.320          288          288         288 1564268 2338724 0.5583 0.5124  1.808
                                            35          1     4    1  2338724 2339356            4   4.424            1            1           1 2338724 2339356 0.1314 1.9213  2.503
                                            36          1     5    1  2339356 3084750          563   2.316          263          263         263 2339356 3084750 0.5710 0.4968  1.819
                                            37          1     6    1  3084750 3300535          295   2.028          153          153         153 3084750 3300535 0.8410 0.1612  1.867
                                            locationstart segmentCT_pre0
                                            32         1.000          1.955
                                            33         1.003          2.663
                                            34         1.006          2.320
                                            35         1.009          4.424
                                            36         1.009          2.316
                                            37         1.012          2.028
                                            > head(old_segments[[2]])
                                            chromosome tcnId dhId tcnStart  tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs tcnNbrOfHets dhNbrOfLoci dhStart   dhEnd dhMean c1Mean c2Mean
                                            32          1     1    1    82154  838755           66   1.955           35           35          35   82154  838755 0.6958 0.2973  1.657
                                            33          1     2    1   838755 1564268          439   2.663          140          140         140  838755 1564268 0.8212 0.2381  2.425
                                            34          1     3    1  1564268 2338724          588   2.320          288          288         288 1564268 2338724 0.5583 0.5124  1.808
                                            35          1     4    1  2338724 2339356            4   4.424            1            1           1 2338724 2339356 0.1314 1.9213  2.503
                                            36          1     5    1  2339356 3084750          563   2.316          263          263         263 2339356 3084750 0.5710 0.4968  1.819
                                            37          1     6    1  3084750 3300535          295   2.028          153          153         153 3084750 3300535 0.8410 0.1612  1.867
                                            locationstart segmentCT_pre0 segmentCT_pre1 segmentCT_pre2 segmentCT_pre3
                                            32         1.000          1.955          1.946          1.871          1.789
                                            33         1.003          2.663          2.390          2.294          2.201
                                            34         1.006          2.320          2.172          2.083          1.998
                                            35         1.009          4.424          2.777          3.333          2.799
                                            36         1.009          2.316          2.145          2.057          1.973
                                            37         1.012          2.028          1.964          1.887          1.810
                                            > head(old_segments[[3]])
                                            chromosome tcnId dhId tcnStart  tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs tcnNbrOfHets dhNbrOfLoci dhStart   dhEnd dhMean  c1Mean c2Mean
                                            1          1     1    1   534247  819209           37   1.684           19           19          19  534247  819209 0.6760 0.27291  1.411
                                            2          1     2    1   819209 1120144          221   2.042           82           82          82  819209 1120144 0.7875 0.21698  1.825
                                            3          1     3    1  1120144 1563654          215   2.325           57           57          57 1120144 1563654 0.8326 0.19459  2.131
                                            4          1     4    1  1563654 1970092          232   1.906          122          122         122 1563654 1970092 0.5017 0.47492  1.431
                                            5          1     5    1  1970092 1997873           29   2.451            8            8           8 1970092 1997873 0.9754 0.03008  2.421
                                            6          1     6    1  1997873 2115914          107   1.796           37           37          37 1997873 2115914 0.9465 0.04800  1.748
                                            locationstart segmentCT_pre4 segmentCT_pre5 segmentCT_pre6 segmentCT_pre7
                                            1         1.002          1.684          1.624          1.536          1.452
                                            2         1.003          2.042          1.973          1.864          1.767
                                            3         1.004          2.325          2.195          2.089          1.974
                                            4         1.006          1.906          1.832          1.736          1.642
                                            5         1.008          2.451          2.216          2.142          2.011
                                            6         1.008          1.796          1.740          1.647          1.559
                                            > head(segments)
                                            chromosome tcnId dhId tcnStart  tcnEnd tcnNbrOfLoci tcnMean tcnNbrOfSNPs tcnNbrOfHets dhNbrOfLoci dhStart   dhEnd dhMean  c1Mean c2Mean
                                            1          1     1    1   723918 1120144          256   1.397          101          101         101  723918 1120144 0.7665 0.16308 1.2336
                                            2          1     2    1  1120144 1563654          215   1.629           57           57          57 1120144 1563654 0.8326 0.13628 1.4922
                                            3          1     3    1  1563654 1970092          232   1.304          122          122         122 1563654 1970092 0.5017 0.32476 0.9787
                                            4          1     4    1  1970092 1997873           29   1.840            8            8           8 1970092 1997873 0.9754 0.02259 1.8176
                                            5          1     5    1  1997873 2115914          107   1.194           37           37          37 1997873 2115914 0.9465 0.03191 1.1619
                                            6          1     6    1  2115914 2171412           30   1.776           13           13          13 2115914 2171412 0.2700 0.64830 1.1278
                                            locationstart segmentCT_pre12 segmentCT_pre13 segmentCT_pre14 segmentCT_pre15 segmentCT_pre1 segmentCT_pre2 segmentCT_pre3
                                            1         1.003           1.397           1.329           1.170           1.123          2.179          2.085          1.995
                                            2         1.004           1.629           1.492           1.362           1.256          2.527          2.434          2.341
                                            3         1.006           1.304           1.280           1.103           1.077          2.079          1.993          1.906
                                            4         1.008           1.840           1.368           1.549           1.090          2.631          2.536          2.451
                                            5         1.008           1.194           1.246           1.027           1.044          1.963          1.880          1.796
                                            6         1.008           1.776           1.347           1.502           1.068          2.581          2.476          2.392
                                            segmentCT_pre4 segmentCT_pre5 segmentCT_pre6 segmentCT_pre7 segmentCT_pre8 segmentCT_pre9 segmentCT_pre10 segmentCT_pre11
                                            1          1.994          1.926          1.820          1.725          1.730          1.604           1.488           1.397
                                            2          2.325          2.195          2.089          1.974          1.927          1.866           1.755           1.645
                                            3          1.906          1.832          1.736          1.642          1.713          1.503           1.403           1.304
                                            4          2.451          2.216          2.142          2.011          1.713          2.041           1.938           1.840
                                            5          1.796          1.740          1.647          1.559          1.713          1.395           1.293           1.194
                                            6          2.392          2.163          2.093          1.964          1.713          1.977           1.873           1.776
                                            segmentCT_pre16
                                            1              NA
                                            2              NA
                                            3              NA
                                            4              NA
                                            5              NA
                                            6              NA
                                            > head(old_segments[[5]])
                                            Error in old_segments[[5]] : subscript out of bounds
                                            > #g1 <- get_TCN_track(old_segments,"tcnMean")+ylim(0,4)
                                              > g1 <-get_TCN_track(old_segments[[2]],"segmentCT_pre0")+ylim(0,4)
                                            > g2 <-get_TCN_track(old_segments[[2]],"segmentCT_pre1")+ylim(0,4)
                                            > g3 <-get_TCN_track(old_segments[[2]],"segmentCT_pre2")+ylim(0,4)
                                            > g4 <-get_TCN_track(old_segments[[2]],"segmentCT_pre3")+ylim(0,4)
                                            > grid.arrange(g1,g2,g3,g4)
                                            Warning messages:
                                              1: Removed 116 rows containing missing values (geom_segment).
                                            2: Removed 17 rows containing missing values (geom_segment).
                                            3: Removed 54 rows containing missing values (geom_segment).
                                            4: Removed 27 rows containing missing values (geom_segment).
                                            > renormalise
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

g1 <-get_TCN_track(old_segments[[3]],"segmentCT_pre4")+ylim(0,4)
g2 <-get_TCN_track(old_segments[[3]],"segmentCT_pre5")+ylim(0,4)
g3 <-get_TCN_track(old_segments[[3]],"segmentCT_pre6")+ylim(0,4)
g4 <-get_TCN_track(old_segments[[3]],"segmentCT_pre7")+ylim(0,4)
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

