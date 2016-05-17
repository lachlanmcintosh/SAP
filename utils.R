renormalise2 <- function(data,segments,pr){
  # start with the original data:
  data$CNsnp <- data$CT
  data$CNseg <- data[,paste("segmentCN_pre_clustered",as.character(pr-1),sep="")]

  data$reddye <- data$CNsnp*data$B.Allele.Freq
  data$greendye <- data$CNsnp*(1-data$B.Allele.Freq)
  data$homs <- data$Allele1...Design == data$Allele2...Design

  # gredinv <- gam(formula=CNseg~0+s(reddye),data=data[which(data$homs & data$reddye > data$greendye),])
  # ggreeninv <- gam(formula=CNseg~0+s(greendye),data=data[which(data$homs & data$reddye < data$greendye),])
  gredinv <- lm(formula=CNseg~0+reddye,data=data[which(data$homs & data$reddye > data$greendye),])
  ggreeninv <- lm(formula=CNseg~0+greendye,data=data[which(data$homs & data$reddye < data$greendye),])

  new_green = data.frame(data$greendye)
  colnames(new_green) <- "greendye"
  data$newgreen <- predict(ggreeninv,newdata=new_green)

  new_red = data.frame(data$reddye)
  colnames(new_red) <- "reddye"
  data$newred <- predict(gredinv,newdata=new_red)

  data$newbaf <- data$newgreen/(data$newred+data$newgreen)

  # ok now apply the tumor boost correction:
  med_baf <- median(data[which(!data$homs),"newbaf"])

  data$newbaf2 <- data$newbaf
  data[which(data$newbaf < med_baf),"newbaf2"] <- 0.5*data[which(data$newbaf < med_baf),"newbaf"]/med_baf
  data[which(data$newbaf > med_baf),"newbaf2"] <- 1-0.5*(1-data[which(data$newbaf > med_baf),"newbaf"])/(1-med_baf)

  data$CNsnp <- data$newgreen +data$newred
  data$reddye <- data$CNsnp*(1-data$newbaf2)
  data$greendye <- data$CNsnp*data$newbaf2
  data$baf <- data$newbaf2

  ghominv <- lm(CNseg~0+CNsnp,data=data[which(data$homs),])
  gnonhominv <- lm(CNseg~0+CNsnp,data=data[which(!data$homs),])
#
#   summary(ghominv)
#   summary(gnonhominv)

  data$CNsnp2 <- data$CNsnp

  new_hom = data.frame(data[which(data$homs),"CNsnp"])
  colnames(new_hom) <- "CNsnp"
  data[which(data$homs),"CNsnp2"] <- predict(ghominv,newdata=new_hom)
  new_nonhom = data.frame(data[which(!data$homs),"CNsnp"])
  colnames(new_nonhom) <- "CNsnp"
  data[which(!data$homs),"CNsnp2"] <- predict(gnonhominv,newdata=new_nonhom)

  data$CNsnp <- data$CNsnp2

  # data$CNsnp <- data$CNsnp*2/median(data$CNsnp,na.rm=TRUE)

  segments <- get_new_seg_estimates(data=data,segments=segments,
                                    new_snp_name = "CNsnp",
                                    new_seg_name = paste("segmentCN_pre_nogam",as.character(pr),sep=""))
  data <- put_seg_estimates_into_data_array(data=data,segments=segments,
                                            segment_name = paste("segmentCN_pre_nogam",as.character(pr),sep=""),
                                            new_snp_name = paste("segmentCN_pre_nogam",as.character(pr),sep=""))

  data$CNseg <- data[,paste("segmentCN_pre_nogam",as.character(pr),sep="")]
  return(list(data=data,segments=segments))
}

cluster_ACN2 <- function(dat,iterations,local_adj_thresh,global_adj_thresh,p_close_local,p_close_global,TCN,TCN_se,first=TRUE){
  # let the CN estimates converge locally if below the specified local threshold
  if(first==TRUE){
    dat[,"TCN"] <- dat[,TCN] #segmentCN_pre12
  }
  dat[,"length"] <- dat[,"tcnNbrOfLoci"]
  dat <- dat[complete.cases(dat$TCN) & complete.cases(dat$length) & complete.cases(dat$chromosome),]
  for (i in 2:nrow(dat)){
    if (abs(dat[i,"TCN"]-dat[i-1,"TCN"])<local_adj_thresh & dat[i-1,"chromosome"] == dat[i,"chromosome"]){
      p <- p_close_local
      mean_TCN <- (dat[i,"TCN"]*dat[i,"length"]+dat[i-1,"TCN"]*dat[i-1,"length"])/(dat[i,"length"]+dat[i-1,"length"])
      stopifnot(abs(dat[i-1,"TCN"] - mean_TCN) < local_adj_thresh & abs(dat[i,"TCN"] - mean_TCN) < local_adj_thresh)
      dat[i-1,"TCN"] <- mean_TCN*p + dat[i-1,"TCN"]*(1-p)
      dat[i,"TCN"] <- mean_TCN*p + dat[i,"TCN"]*(1-p)
    }
  }

  for(i in 3:nrow(dat)){
    if (abs(dat[i,"TCN"]-dat[i-2,"TCN"])<local_adj_thresh & dat[i-2,"chromosome"] == dat[i,"chromosome"]){
      mean_TCN <- (dat[i,"TCN"]*dat[i,"length"]+dat[i-2,"TCN"]*dat[i-2,"length"])/(dat[i,"length"]+dat[i-2,"length"])
      p <- p_close_local
      stopifnot(abs(dat[i-2,"TCN"] - mean_TCN) < local_adj_thresh & abs(dat[i,"TCN"] - mean_TCN) < local_adj_thresh)
      dat[i-2,"TCN"] <- mean_TCN*p + dat[i-2,"TCN"]*(1-p)
      dat[i,"TCN"] <- mean_TCN*p + dat[i,"TCN"]*(1-p)
    }
  }

  if(global_adj_thresh > 0){
    hc <- hclust(dist(dat$TCN),method="complete")                # apply hirarchical clustering
    dat$memb <- cutree(hc, h=global_adj_thresh)

    for(i in 1:max(dat$memb)){
      temp <- dat[which(dat$memb == i),]
      mean_TCN <- sum(temp$TCN*temp$length)/sum(temp$length)
      stopifnot(abs(dat[which(dat$memb == i),"TCN"] - mean_TCN) < global_adj_thresh)
      dat[which(dat$memb == i),"TCN"] <- p_close_global*mean_TCN + (1-p_close_global)*dat[which(dat$memb == i),"TCN"]
    }
  }

  #stopifnot(abs(dat[,TCN] - dat$TCN) < local_adj_thresh+global_adj_thresh)

  # dat[,TCN] <- dat$TCN
  if(iterations > 1){
    return(cluster_ACN2(dat,(iterations -1),local_adj_thresh,global_adj_thresh,p_close_local,p_close_global,TCN,TCN_se,first=FALSE))
  } else if(iterations == 1){
    return(dat)
  }
}

cluster_ACN <- function(dat,iterations,local_adj_thresh,global_adj_thresh,p_close_local,p_close_global,TCN,TCN_se){
  # let the CN estimates converge locally if below the specified local threshold
  dat[,"TCN"] <- dat[,TCN] #segmentCN_pre12
  #dat[,"TCN_se"] <- dat[,TCN_se] #segmentCN_pre12
  dat[,"length"] <- dat[,"tcnNbrOfLoci"]
  #dat <- dat[complete.cases(dat$TCN)&complete.cases(dat$TCN_se),]
  dat <- dat[complete.cases(dat$TCN) & complete.cases(dat$length) & complete.cases(dat$chromosome),]
  for (i in 2:nrow(dat)){
    if (abs(dat[i,"TCN"]-dat[i-1,"TCN"])<local_adj_thresh & dat[i-1,"chromosome"] == dat[i,"chromosome"]){
      # t = (dat[i,"TCN"]-dat[i-1,"TCN"])/(dat[i,"TCN_se"]^2/dat[i,"length"]+dat[i-1,"TCN_se"]^2/dat[i-1,"length"])^0.5
      # df = ((dat[i,"TCN_se"]^2/dat[i,"length"])^2+(dat[i-1,"TCN_se"]^2/dat[i-1,"length"])^2)^2/
      #         ((dat[i,"TCN_se"]^2/dat[i,"length"])^2*(dat[i,"length"]-1)+
      #          (dat[i-1,"TCN_se"]^2/dat[i-1,"length"])^2*(dat[i,"length"]-1))
      # p <- 2*pt(-abs(t),df=df)
      # print(p)
      # p <- max(p,2*p_close_local)/2
      # print(p)
      # rather than use the t-distribution here, assign probabilities based upon the ECDF of the differences between segments. Augment this ECDF with unobserved true changes.
      p <- p_close_local
      mean_TCN <- (dat[i,"TCN"]*dat[i,"length"]+dat[i-1,"TCN"]*dat[i-1,"length"])/(dat[i,"length"]+dat[i-1,"length"])
      stopifnot(abs(dat[i-1,"TCN"] - mean_TCN) < local_adj_thresh & abs(dat[i,"TCN"] - mean_TCN) < local_adj_thresh)
      dat[i-1,"TCN"] <- mean_TCN*p + dat[i-1,"TCN"]*(1-p)
      dat[i,"TCN"] <- mean_TCN*p + dat[i,"TCN"]*(1-p)
    }
  }

  for(i in 3:nrow(dat)){
    if (abs(dat[i,"TCN"]-dat[i-2,"TCN"])<local_adj_thresh & dat[i-2,"chromosome"] == dat[i,"chromosome"]){
      mean_TCN <- (dat[i,"TCN"]*dat[i,"length"]+dat[i-2,"TCN"]*dat[i-2,"length"])/(dat[i,"length"]+dat[i-2,"length"])
      # t = (dat[i,"TCN"]-dat[i-2,"TCN"])/(dat[i,"TCN_se"]^2/dat[i,"length"]+dat[i-2,"TCN_se"]^2/dat[i-2,"length"])^0.5
      # df = ((dat[i,"TCN_se"]^2/dat[i,"length"])^2+(dat[i-2,"TCN_se"]^2/dat[i-2,"length"])^2)^2/
      #   ((dat[i,"TCN_se"]^2/dat[i,"length"])^2*(dat[i,"length"]-1)+
      #      (dat[i-2,"TCN_se"]^2/dat[i-2,"length"])^2*(dat[i,"length"]-1))
      # p <- 2*pt(-abs(t),df=df)
      # p <- max(p,p_close_local)
      p <- p_close_local
      stopifnot(abs(dat[i-2,"TCN"] - mean_TCN) < local_adj_thresh & abs(dat[i,"TCN"] - mean_TCN) < local_adj_thresh)

      dat[i-2,"TCN"] <- mean_TCN*p + dat[i-2,"TCN"]*(1-p)
      dat[i,"TCN"] <- mean_TCN*p + dat[i,"TCN"]*(1-p)
    }
  }

  if(global_adj_thresh > 0){
    hc <- hclust(dist(dat$TCN),method="complete")                # apply hirarchical clustering
    dat$memb <- cutree(hc, h=global_adj_thresh)

    for(i in 1:max(dat$memb)){
      temp <- dat[which(dat$memb == i),]
      mean_TCN <- sum(temp$TCN*temp$length)/sum(temp$length)
      stopifnot(abs(dat[which(dat$memb == i),"TCN"] - mean_TCN) < global_adj_thresh)
      dat[which(dat$memb == i),"TCN"] <- p_close_global*mean_TCN + (1-p_close_global)*dat[which(dat$memb == i),"TCN"]
    }
  }

  stopifnot(abs(dat[,TCN] - dat$TCN) < local_adj_thresh+global_adj_thresh)

  dat[,TCN] <- dat$TCN
  if(iterations > 1){
    return(cluster_ACN(dat,iterations -1,local_adj_thresh,global_adj_thresh,p_close_local,p_close_global,TCN))
  } else if(iterations == 1){
    return(dat)
  }
}

merge_data_with_GC <- function(data,GC){
  # need to find a faster way to do this.
  data$name <- paste(as.character(data[,c("Chr")]),as.character(data[,c("Position")]),sep=" ")
  data <- merge(data,GC,by=c("name","name"))
  data$Position <- data$Position.x
  return(data)
}

get_GC_file <- function(){
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
  return(GC)
}

#ML estimation with a glm
get_var <- function(major,minor,n,s1,s2,s3) {
  var <- (s1+2*s2+s3)/(n*(major+minor)^2) - 2*(s2+s3)/(n*minor*(major+minor))+s3/(n*minor^2)
  return(var)
}

#get_var(1,2,100,matrix(c(1,0,0,1),nrow=2))+get_var(1,2,100,matrix(c(1,0,0,1),nrow=2))

ML_est_of_dye_CN_rel <- function(segments,max_ord){
  names <- c("major","minor","total","s1","s2","s3","n")
  df1 <- segments[1:(nrow(segments)-1),names]
  df2 <- segments[2:nrow(segments),names]
  df_low <- df1
  df_low[which(df1$major+df1$minor > df2$major+df2$minor),] <- df2[which(df1$major+df1$minor > df2$major+df2$minor),]
  df_high <- df1
  df_high[which(df1$major+df1$minor < df2$major+df2$minor),] <- df2[which(df1$major+df1$minor < df2$major+df2$minor),]

  CNp <- cbind(df_low,df_high)
  colnames(CNp) <- c("major1","minor1","total1","s11","s21","s31","n1","major2","minor2","total2","s12","s22","s32","n2")

  r1 <- log((CNp$major1/CNp$major2) * (CNp$total2/CNp$total1))

  r2 <- log((CNp$minor2/CNp$minor2) * (CNp$total2/CNp$total1))

  s1 <- get_var(CNp$minor1,CNp$major1,CNp$n1,CNp$s31,CNp$s21,CNp$s11)+
          get_var(CNp$minor2,CNp$major2,CNp$n2,CNp$s32,CNp$s22,CNp$s12)

  s2 <-  get_var(CNp$major1,CNp$minor1,CNp$n1,CNp$s11,CNp$s21,CNp$s31)+
          get_var(CNp$major2,CNp$minor2,CNp$n2,CNp$s12,CNp$s22,CNp$s32)

  X <- sapply((CNp$total2 - CNp$total1), function(x) all.moments(r,order.max=max_ord))
  alpha <- c(0.4,0.4,0.2)
  beta <- c(1,rep(0,max_ord-1))

  theta <- 1/(max(cbind(r1,r2)) - min(cbind(r1,r2)))

  for(i in 1:100){
    density <- cbind(dnorm(x=r1,mean=X%*%beta,s1)/theta, dnorm(x=r2,mean=X%*%beta,s2)/theta, 1/theta^2)
    p <- density%*%t(alpha)/sum(density%*%t(alpha))
    alpha <- apply(p,1,sum)/nrow(p)
    beta <- solve(((diag(p[1,])+diag(p[2,]))%*%t(X)%*%X)) %*% (diag(p[1,]) %*% t(X) %*% r1 + diag(p[2,]) %*% t(X) %*% r2)
  }
  return(beta = beta, alpha = alpha)
}



combine_het_and_hom_info <- function(segments){
  for(i in 1:nrow(segments)){
    # bt <- binom.test(segments$length-segments$NOHfreq*segments$length, segments$NOHfreq*segments$length, p = 0.95,alternative="less")
    # if(bt$p.value < 0.0001){
    #   #then probably not a homozygote segment?
    #   # what do you want to do with this information?
    # }
    if(!is.na(segments[i,"mux"] + segments[i,"mux"])){
      segments[i,"average_TCN"] <- (segments[i,"mux"] + segments[i,"muy"])*segments[i,"NOHfreq"] + segments[i,"NOHTCN"]*(1-segments[i,"NOHfreq"])
      segments[i,"M"] <- segments[i,"major"] <- segments[i,"mux"]*segments[i,"average_TCN"]/(segments[i,"mux"] + segments[i,"muy"])
      segments[i,"m"] <- segments[i,"minor"] <- segments[i,"muy"]*segments[i,"average_TCN"]/(segments[i,"mux"] + segments[i,"muy"])
    } else{
      # then declare it is a homozygote segment?
      # can always qc later and chuck these out anyway. what else can you do?
      segments[i,"M"] <- segments[i,"major"] <- segments[i,"NOHTCN"]
      segments[i,"m"] <- segments[i,"minor"] <- 0
    }
  }
  return(segments)
}


clean_segments <- function(mat){
  thresh <- median(c(abs(mat[1:(nrow(mat)-1),"major"]- mat[2:(nrow(mat)),"major"]),abs(mat[1:(nrow(mat)-1),"minor"]- mat[2:(nrow(mat)),"minor"])))

  mat <- cluster_PSACN(mat,iterations=20,min_length=25,local_adj_thresh=thresh/2,global_adj_thresh=0,false_to_true_thresh=0,p_straighten=0,p_close_local=0.25,p_close_global=0,allelic_balance_info="not_present",print_plots=FALSE,small_length=25)

  mat <- cluster_PSACN(mat,iterations=20,min_length=25,local_adj_thresh=thresh,global_adj_thresh=thresh/2,false_to_true_thresh=0,p_straighten=0,p_close_local=0.5,p_close_global=0,allelic_balance_info="not_present",print_plots=FALSE,small_length=25)

  return(mat)
}

load_centromeres <- function(){
  file = "centromeres.txt"
  BASE = '/home/users/allstaff/lmcintosh/ITH/SNP_stuff/data/'
  centromeres <- read.table(paste(BASE,file,sep=''), sep="\t", header=FALSE)
  colnames(centromeres) <- c("chr","x")
  centromeres$x <- as.numeric(centromeres$x)*10^6 # is in megabase format
  return(centromeres)
}

do_some_seg <- function(data,round,seg_name="CT"){
  dataPSCBS <- data[,c("Chr","Position","B.Allele.Freq",seg_name)] # don't know why I am stuck here

  setnames(dataPSCBS, c("chromosome","x","betaT","CT"))

  # remove incomplete data
  data <- data[complete.cases(dataPSCBS),]
  dataPSCBS <- dataPSCBS[complete.cases(dataPSCBS),]

  #### end loading of unpaired test dataset: SAMPLE 13 ####
  ##### END INITIAL DATA PROCESSING #####

  ##### PSCBS MAIN BLOCK #####
  # PSCBS stands for parent specific circular binary segmentation.
  # The PSCBS manual: https://cAB2an.r-project.org/web/packages/PSCBS/vignettes/PairedPSCBS.pdf
  # remove and install the PSCBS package to ensure the most up to date version:
  # remove.packages(c("PSCBS"))
  # install.packages(c("PSCBS"))

  # PSCBS BEGIN the PSCBS method as per the specification at
  # https://cran.r-project.org/web/packages/PSCBS/PSCBS.pdf

  # Henrik Bengtsson's utilities library:
  verbose <- R.utils::Arguments$getVerbose(-10*interactive(), timestamp=TRUE)
  # is this required for PSCBS?

  # (OPTIONAL)
  # Drop single-locus "TCN" outliers
  dataPSCBSwso <- dropSegmentationOutliers(dataPSCBS)

  # drop large gaps (greater than 1 mega base inbetween observations)
  gaps <- findLargeGaps(dataPSCBSwso, minLength = 1e+06)
  knownSegments <- gapsToSegments(gaps)

  # A segmentation (OPTIONS)
  # USE IF NON PAIRED NORMAL:
  fit <- segmentByNonPairedPSCBS(dataPSCBSwso, preserveScale = FALSE, verbose =FALSE)
  # USE IF PAIRED NORMAL (can skip tumor boost normalisation by setting tbn = FALSE)
  # fit <- segmentByPairedPSCBS(data, knownSegments = knownSegments, preserveScale = FALSE, seed = 48879, verbose = -10)

  # identifying segments:
  segments <-  getSegments(fit)
  segments <- segments[which(segments$chromosome != 0),]
  segments <- segments[which(complete.cases(segments$chromosome)),]

  #need to put the segment copy in here.
  max_pos = max(as.numeric(data$Position))+1
  segments$locationstart <- as.numeric(segments$chromosome)+as.numeric(segments$tcnStart)/max_pos
  data$location <- as.numeric(data$Chr)+as.numeric(data$Position)/max_pos
  data$segment <- cut(x=as.numeric(data$location),breaks=as.numeric(segments$locationstart),labels=FALSE)
  data[,paste("resegmentCN_PSCBS",as.character(round),sep="")] <- segments[data$segment,"tcnMean"]

  segments[,paste("segmentCN_pre",round,sep="")] <- sapply(1:nrow(segments),function(x) mean(data[which(data$segment == x),paste("resegmentCN_PSCBS",as.character(round),sep="")],na.rm=TRUE))
  data[,paste("segmentCN_pre",round,sep="")] <- segments[data$segment,paste("segmentCN_pre",round,sep="")]

  # make some plots:
  ##### END PSCBS MAIN BLOCK #####

  ##### SUPPLEMENTARY CALLS #####
  # For each segment, the test for ROH calculates the fraction of SNPs that are called heterozygous.
  # fit <- callROH(fit, verbose = -10) # [2016-02-23 20:08:17] Exception: Cannot call ROH from 'NonPairedPSCBS' data.

  # could do this with a binomial test on the sameness of the two alleles at any location. Could also do a test on BAF. Could combine both of these tests.

  # The AB caller tests whether a segment’s DH is zero or not, by comparing its DH level (or more
  # precisely, the 5% quantile of its bootstrapped DH mean level) to a threshold. This threshold will
  # be a function of the noise level, because the noiser the BAF signals (and hence the DH signals),
  # the greater the bias of the DH mean level for segments in AB will be
  # fit <- callAB(fit, verbose = -10)



  # pdeltaAB <- parLapply(cl,pfit,function(x) estimateDeltaAB(x,scale=1))
  # pfit <- parLapply(cl,merge_my_lists(pfit,pdeltaAB),function(x) callAB(x$x, delta=x$y))
  #
  # # tests if a segment is not in a LOH region
  # # fit <- callLOH(fit, verbose = -10)
  #
  # # could alt tune some parameters
  # pdeltaLOH <-  parLapply(cl,pfit,function(x) estimateDeltaLOH(x, midpoint=1/2))
  # pfit <- parLapply(cl,merge_my_lists(pfit,pdeltaLOH),function(x) callLOH(x$x, delta=x$y))
  #
  # # neutral total copy number
  # pkappa <- parLapply(cl,pfit,estimateKappa)
  # pdeltaCN <- parLapply(cl,merge_my_lists(pfit,pkappa),function(x) estimateDeltaCN(x$x, scale=1, kappa=x$y))
  # pfit <- parLapply(cl,merge_my_lists(pfit,pdeltaCN),function(x) callNTCN(x$x, delta=x$y, verbose=-10))
  ##### END SUPPLEMENTARY CALLS #####

  ##### SAVE PSCBS OUTPUT #####
  # save output!
  #pathname <- writeSegments(fit, name=paste(filename,"_PSCBS_segments",Sys.time(),sep=""), simplify=TRUE)
  # make a report!
  #report(fit, sampleName=paste(filename,"_PSCBS_report",Sys.time(),sep=""), studyName="PSCBS-Ex", verbose=FALSE)

  ##### END SAVE PSCCBS OUTPUT #####
  return(list(segments = segments,fit = fit,data=data))
}


merge_my_three_lists <- function(listx,listy,listz){
  finallist = list()
  for(i in 1:length(listx)){
    finallist[[i]] <- list(x=listx[[i]],y=listy[[i]],z=listz[[i]])
  }
  return(finallist)
}

get_segment_ids <- function(segments){
  return(segments[,c("chromosome","tcnStart","tcnEnd")])
}

renormalise <- function(data,segments,pr){
  # start with the original data:
  data$CNsnp <- data$CT
  data$CNseg <- data[,paste("segmentCN_pre",as.character(pr-1),sep="")]

  data$reddye <- data$CNsnp*data$B.Allele.Freq
  data$greendye <- data$CNsnp*(1-data$B.Allele.Freq)
  data$homs <- data$Allele1...Design == data$Allele2...Design

  # calculate relationships to segment value and their inverses:
  # gred <- lm(formula=reddye~0+CNseg,data=data[which(data$homs & data$reddye > data$greendye),])
  # ggreen <- lm(formula= greendye~0+CNseg,data=data[which(data$homs & data$reddye < data$greendye),])
  gredinv <- lm(formula=CNseg~0+reddye,data=data[which(data$homs & data$reddye > data$greendye),])
  ggreeninv <- lm(formula=CNseg~0+greendye,data=data[which(data$homs & data$reddye < data$greendye),])

  # gred <- lm(formula=reddye~0+CNseg + I(CNseg^3),data=data[which(data$homs & data$reddye > data$greendye),])
  # ggreen <- lm(formula= greendye~0+CNseg + I(CNseg^3),data=data[which(data$homs & data$reddye < data$greendye),])
  # gredinv <- lm(formula=CNseg~0+reddye+I(reddye^3),data=data[which(data$homs & data$reddye > data$greendye),])
  # ggreeninv <- lm(formula=CNseg~0+greendye+I(greendye^3),data=data[which(data$homs & data$reddye < data$greendye),])

  # plot(gred)
  # plot(ggreen)
  # plot(gredinv)
  # plot(ggreeninv)
  # summary(gred)
  # summary(ggreen)
  summary(gredinv)
  summary(ggreeninv)

  new_green = data.frame(data$greendye)
  colnames(new_green) <- "greendye"
  # new_green2 <- predict(ggreeninv,newdata=new_green)
  # new_green2 <- data.frame(new_green2)
  # colnames(new_green2) <- "CNseg"
  # data$newgreen <- (data$greendye+predict(gred,newdata=new_green2))/2
  data$newgreen <- predict(ggreeninv,newdata=new_green)

  # it looks like the gams are overfitting and just simply assigning whatever fucking value that they want.
  # use polynomials instead.
  # maybe only use odd powers to ensure the monotonicty idea.

  new_red = data.frame(data$reddye)
  colnames(new_red) <- "reddye"
  # new_red2 <- predict(gredinv,newdata=new_red)
  # new_red2 <- data.frame(new_red2)
  # colnames(new_red2) <- "CNseg"
  data$newred <- predict(gredinv,newdata=new_red)

  data$newbaf <- data$newgreen/(data$newred+data$newgreen)

  ggplot(data)+geom_segment(aes(x=reddye,y=greendye,xend=newred,yend=newgreen),alpha=0.01)+xlim(0,5)+ylim(0,5)

  # fit a single guassian on the middle values - this will tell you how to correct for baf at high levels.
  # allow for abberant data - use the uniform prior.


  #ggplot(data)+geom_segment(aes(x=reddye,y=greendye,xend=new_red2,yend=new_green2),alpha=0.01)+xlim(0,5)+ylim(0,5)

  # sd(data$B.Allele.Freq,na.rm=TRUE)
  # sd(data$newbaf,na.rm=TRUE)
  # ggplot()+geom_histogram(aes(data[which(data$homs),"newbaf"]),binwidth=0.01)
  # ggplot()+geom_histogram(aes(data[which(!data$homs),"newbaf"]),binwidth=0.01)
  # ggplot()+geom_histogram(aes(data[which(!data$homs),"B.Allele.Freq"]),binwidth=0.01)

  # median(data[which(data$homs & data$newred < data$newgreen),"newred"])
  # median(data[which(data$homs & data$newred > data$newgreen),"newgreen"])

  # ok now apply the tumor boost correction:
  med_baf <- median(data[which(!data$homs),"newbaf"])

  data$newbaf2 <- data$newbaf
  data[which(data$newbaf < med_baf),"newbaf2"] <- 0.5*data[which(data$newbaf < med_baf),"newbaf"]/med_baf
  data[which(data$newbaf > med_baf),"newbaf2"] <- 1-0.5*(1-data[which(data$newbaf > med_baf),"newbaf"])/(1-med_baf)

  data$CNsnp <- data$newgreen +data$newred
  data$reddye3 <- data$CNsnp*(1-data$newbaf2)
  data$greendye3 <- data$CNsnp*data$newbaf2
  data$baf <- data$newbaf2
  ggplot(data)+geom_segment(aes(x=reddye,y=greendye,xend=reddye3,yend=greendye3),alpha=0.01)+xlim(0,5)+ylim(0,5)

  data$reddye <- data$reddye3
  data$greendye <- data$greendye3

  # ggplot()+geom_histogram(aes(data[which(data$homs),"newbaf2"]),binwidth=0.01)
  # ggplot()+geom_histogram(aes(data[which(!data$homs),"newbaf2"]),binwidth=0.01)

  # should now do correction of homs to non homs in total copy number in the same way we just did red to non red!
  # ghom <- lm(CNsnp~0+CNseg+I(CNseg^3),data=data[which(data$homs),])
  # gnonhom <- lm(CNsnp~0+CNseg+I(CNseg^3),data=data[which(!data$homs),])
  ghominv <- lm(CNseg~0+CNsnp,data=data[which(data$homs),])
  gnonhominv <- lm(CNseg~0+CNsnp,data=data[which(!data$homs),])

  # ghom <- lm(CNsnp~0+CNseg+I(CNseg^3),data=data[which(data$homs),])
  # gnonhom <- lm(CNsnp~0+CNseg+I(CNseg^3),data=data[which(!data$homs),])
  # ghominv <- lm(CNseg~0+CNsnp+I(CNsnp^3),data=data[which(data$homs),])
  # gnonhominv <- lm(CNseg~0+CNsnp+I(CNsnp^3),data=data[which(!data$homs),])

  # summary(ghom)
  # summary(gnonhom)
  summary(ghominv)
  summary(gnonhominv)


  data$CNsnp2 <- data$CNsnp
  new_hom = data.frame(data[which(data$homs),"CNsnp"])
  colnames(new_hom) <- "CNsnp"
  data[which(data$homs),"CNsnp2"] <- predict(ghominv,newdata=new_hom)
  # new_hom2 <- data.frame(new_hom2)
  # colnames(new_hom2) <- "CNseg"
  # data[which(data$homs),"CNsnp2"] <- (data[which(data$homs),"CNsnp"] + predict(gnonhom,newdata=new_hom2))/2

  new_nonhom = data.frame(data[which(!data$homs),"CNsnp"])
  colnames(new_nonhom) <- "CNsnp"
  data[which(!data$homs),"CNsnp2"] <- predict(gnonhominv,newdata=new_nonhom)
  # new_nonhom2 <- predict(gnonhominv,newdata=new_nonhom)
  # new_nonhom2 <- data.frame(new_nonhom2)
  # colnames(new_nonhom2) <- "CNseg"
  # data[which(!data$homs),"CNsnp2"] <- (data[which(!data$homs),"CNsnp"] + predict(ghom,newdata=new_nonhom2))/2
  data$CNsnp <- data$CNsnp2

  data$CNsnp <- data$CNsnp*2/median(data$CNsnp,na.rm=TRUE)

  segments <- get_new_seg_estimates(data=data,segments=segments,
                                    new_snp_name = "CNsnp",
                                    new_seg_name = paste("segmentCN_pre_nogam",as.character(pr),sep=""))
  data <- put_seg_estimates_into_data_array(data=data,segments=segments,
                                            segment_name = paste("segmentCN_pre_nogam",as.character(pr),sep=""),
                                            new_snp_name = paste("segmentCN_pre_nogam",as.character(pr),sep=""))


  data$CNseg <- data[,paste("segmentCN_pre_nogam",as.character(pr),sep="")]
  return(list(data=data,segments=segments))
}

get_TCN_track_precision <- function(mat,seg_est_name){
  g <- ggplot(mat) +
    geom_segment(aes_string(x = "tcnStart", y = seg_est_name, xend = "tcnEnd", yend = seg_est_name),size=2)+
    facet_grid(.~chromosome,scales = "free_x", space = "free")
  return(g)
}

reformat <- function(data,segments){
  data = data[which(data$Chr != 0),]

  #Find which segment each thing belongs to
  max_pos = max(as.numeric(data$Position))+1
  data$location <- as.numeric(data$Chr)+as.numeric(data$Position)/max_pos
  segments$locationstart <- as.numeric(segments$chromosome)+as.numeric(segments$tcnStart)/max_pos

  data$segment <- cut(x=as.numeric(data$location),breaks=as.numeric(segments$locationstart),labels=FALSE)
  data$segmentCN_PSCBS <- segments[data$segment,"tcnMean"]

  return(data)
}

get_merged_data <- function(data,GC){
  #common_names <- intersect(rownames(GC),rownames(data))
  return(cbind(data,GC[data$name,]))
}

get_new_snp_estimates <- function(gam,data,new_snp_name,old_segment_name){
  if(!is.null(gam$na.action)){
    data <- data[-gam$na.action,]
  }
  data[,new_snp_name] <- (gam$residuals+1)*data[,old_segment_name]
  data[,new_snp_name] <- data[,new_snp_name]/median(data[,new_snp_name],na.rm=TRUE)*2
  return(data)
}

get_new_snp_estimates4 <- function(gam,data,new_snp_name,old_seg_name){
  if(!is.null(gam$na.action)){
    data <- data[-gam$na.action,]
  }
  data[,new_snp_name] <- gam$residuals-gam$fitted.values
  data[,new_snp_name] <- data[,new_snp_name]/median(data[,new_snp_name],na.rm=TRUE)*2
  return(data)
}

get_new_snp_estimates2 <- function(gam,data,new_snp_name,old_segment_name){
  if(!is.null(gam$na.action)){
    data <- data[-gam$na.action,]
  }
  data[,new_snp_name] <- (gam$residuals-gam$fitted.values)*data[,old_segment_name]
  data[,new_snp_name] <- data[,new_snp_name]/median(data[,new_snp_name],na.rm=TRUE)*2
  return(data)
}

get_new_seg_estimates <- function(data,segments,new_snp_name,new_seg_name){
  segments[,new_seg_name] <- sapply(1:nrow(segments),function(x) mean(data[which(data$segment == x),new_snp_name],na.rm=TRUE))
  segments[,paste(new_seg_name,"_stderr",sep="")] <- sapply(1:nrow(segments),function(x) sd(data[which(data$segment == x),new_snp_name],na.rm=TRUE))
  return(segments)
}

get_TCN_track <- function(mat,seg_est_name){
  g <- ggplot(mat) +
    geom_segment(aes_string(x = "tcnStart", y = seg_est_name, xend = "tcnEnd", yend = seg_est_name),size=2)+
    facet_grid(.~chromosome,scales = "free_x", space = "free")
  return(g)
}

put_seg_estimates_into_data_array <- function(data,segments,segment_name,new_snp_name){
  data[,new_snp_name] <- segments[data$segment,segment_name]
  return(data)
}

redo_baf <- function(mat,strech_baf=F){
  # BAF STRECH
  mat <- mat[which(complete.cases(mat)),]
  mat[,c("major","minor")] <- cbind(apply(mat[,c("major","minor")],1,min),apply(mat[,c("major","minor")],1,max))
  mat$baf_final <- mat$major/(mat$major+mat$minor)
  minbaf = quantile(mat$baf_final,probs=0.01,na.rm=T)
  maxbaf = 1-minbaf
  for(i in 1:nrow(mat)){
    if(is.na(mat[i,"baf_final"])){
      next
    }
    if(abs(mat[i,"baf_final"]-0.5)<0.025){
      mat[i,"baf_final"] = 0.5
    } else if(mat[i,"baf_final"] < 0.5 & strech_baf){
      mat[i,"baf_final"] = (mat[i,"baf_final"]-minbaf)*0.5/(0.5-minbaf)
    } else if(mat[i,"baf_final"] > 0.5 & strech_baf){
      mat[i,"baf_final"] = (mat[i,"baf_final"]-0.5)*0.5/(maxbaf-0.5)+0.5
    }
  }


  # BAF SNAP
  range = 0.05
  total <- length <- 0
  for(i in 1:nrow(mat)){
    if(is.na(mat[i,"baf_final"])){
      next
    }
    if((mat[i,"baf_final"]-1/3)>-range & (mat[i,"baf_final"]-1/3)<range){
      total <- total + (1/3 - mat[i,"baf_final"])*mat[i,"length"]
      length <- length + mat[i,"length"]
    }
    if((mat[i,"baf_final"]-2/3)> -range & (mat[i,"baf_final"]-2/3)<range){
      total <- total + (2/3-mat[i,"baf_final"])*mat[i,"length"]
      length <- length + mat[i,"length"]
    }
  }
  corr <- (total/length)
  range <- corr + 0.005*sign(corr)
  fun <- function (x) (1/3-corr)^x - 1/3
  uni <- (uniroot(fun, c(0, 2))$root)
  for(i in 1:nrow(mat)){
    total <-  sum(mat[i,c("major","minor")])
    baf <-  mat[i,c("baf_final")]
    length <- length + mat[i,"length"]
    mat[i,c("major")] <- min(max(sign(baf-0.5)*abs(baf-0.5)^uni+0.5,0),1)*total
    mat[i,c("minor")] <- min(max(0.5-sign(baf-0.5)*abs(baf-0.5)^uni,0),1)*total
  }

  mat$baf_final <- mat$major/(mat$major+mat$minor)
  mat$gmmT <- mat$major+mat$minor
  mat$major <- mat$baf_final*mat$gmmT
  mat$major <- mat$baf_final*mat$gmmT
  return(mat)
}


relabel_segments <- function(data,segments){
  max_pos = max(as.numeric(data$Position))+1
  data$location <- as.numeric(data$Chr)+as.numeric(data$Position)/max_pos
  segments$locationstart <- as.numeric(segments$chromosome)+as.numeric(segments$tcnStart)/max_pos
  data$segment2 <- cut(x=as.numeric(data$location),breaks=as.numeric(segments$locationstart),labels=FALSE)
  return(data)
}

normalise_data <- function(data){
  data$M <- data$B.Allele.Freq*data$CN_pre3
  data$m <- (1-data$B.Allele.Freq)*data$CN_pre3
  for(i in sort(unique(data$segment))){
    # run the em algorithm only on the het snps!
    seg_data = (data$segment2 == i)
    het_data = (data$Allele1...AB != data$Allele2...AB)

    het_x <- data[which(seg_data&het_data),"M"]
    het_y <- data[which(seg_data&het_data),"m"]

    if(length(het_x)<2 | length(het_y)<2){ next}

    data[which(seg_data&het_data),"m"] <- (data[which(seg_data&het_data),"m"]-mean(het_y))*(sd(het_x)/sd(het_y)) + mean(het_x)

    hom_x <- data[which(seg_data&!het_data),"M"]
    hom_y <- data[which(seg_data&!het_data),"m"]
    if(length(hom_x)<2 | length(hom_y)<2){ next}

    data[which(seg_data&!het_data),"m"] <- (data[which(seg_data&!het_data),"m"]-mean(hom_y))*(sd(hom_x)/sd(hom_y)) + mean(hom_x)
    print(i)
  }
  return(data)
}

mark_outliers <- function(data,segments){
  data$is_outlier = rep(FALSE,nrow(data))
  #### detect outliers and take a first pass at the datas####
  for(i in 1:nrow(segments)){
    #find the data just for this segment:
    similarity <- 1/(as.matrix(dist(as.matrix(cbind(data[which(data$segment2 == i),"M"]/sd(data[which(data$segment2 == i),"M"]),data[which(data$segment2 == i),"m"]/sd(data[which(data$segment2 == i),"m"])))))^2+0.01)
    # make the similarity of any two points of different catagories 0:
    similarity <- similarity*outer(data[which(data$segment2 == i),"ishet"],data[which(data$segment2 == i),"ishet"], "==")
    consistent <- data.frame(count = apply(similarity,1,function(x) sum(x)))
    outlier_thresh <- quantile(consistent[,"count"],0.01,na.rm = TRUE)
    data[which(data$segment2 == i),"is_outlier"] <- consistent < outlier_thresh
    print(i)
  }
  return(data)
}

infer_info <- function(data,segments){
  segments$mux <- segments$muy <- segments$sigma1 <- segments$sigma2 <- segments$sigma3 <- segments$sigma4 <- segments$alpha <- segments$ABmeanx <-  segments$ABmeany <- segments$sigmaAB1 <- segments$sigmaAB2 <- segments$sigmaAB3 <- segments$sigmaAB4 <- segments$pmardia <- segments$phztest<- segments$NOHfreq <- segments$NOHTCN <- rep(NA,nrow(segments))
  for(i in 1:nrow(segments)){
    # run the em algorithm only on the het snps!
    seg_data = (data$segment == i) & (data$Allele1...AB != data$Allele2...AB)
    seg_data = seg_data[complete.cases(seg_data)]
    NOH_data <- (data$segment == i) & (data$Allele1...AB == data$Allele2...AB)
    if(sum(seg_data,na.rm=TRUE)>2){
      #info = infer_PSACN(data[which(seg_data),c("M","m")],data$is_outlier,rem_iterations=1000,importance=50)
      info = try(infer_PSACN(data[which(seg_data),c("M","m")],rep(FALSE,sum(seg_data,na.rm=TRUE)),rem_iterations=500,importance=200),silent=TRUE)
      #return(c(new_mu,sigma,alpha,rem_iterations,nrow(x),"ABABABABABAB",ABmean,ABcov,mardiapvalue,hzpvalue))

      try(segments[i,c("mux","muy","sigma1","sigma2","sigma3","sigma4","alpha","ABmeanx","ABmeany","sigmaAB1","sigmaAB2","sigmaAB3","sigmaAB4","pmardia","phztest")] <- info[-c(8,9,10)],silent=TRUE)
      segments[i,"NOHfreq"] <-sum(seg_data)/(sum(seg_data)+sum(NOH_data))
      segments[i,"NOHTCN"] <- mean(apply(data[NOH_data,c("M","m")],1,sum))
      print("")
      print(i)
      print(sum(seg_data,na.rm=TRUE))
    } else{
      print("")
      print("insufficient data points for EM")
      print(sum(seg_data,na.rm=TRUE))
    }
  }
  return(segments)
}


infer_PSACN <- function(major_minor_pairs,outliers,rem_iterations,importance=10,restarts=0){
  major_minor_pairs <- as.matrix(major_minor_pairs[which(!outliers),])
  x <- major_minor_pairs[complete.cases(major_minor_pairs),]
  # should be a n by 2 matrix
  stopifnot(ncol(x) == 2)
  stopifnot(nrow(x) >= 2)
  # there should be no na's by this point!
  stopifnot(!is.na(x))

  nullsigma <- cov(x)

  # partition x into two sets, those greater than the line y=x and otherwise!
  those <- x[,2]-x[,1] > median(x[,2]-x[,1])
  # obtain initial estimates for our parameters
  old_mu <- c(mean(x[which(those),"M"])+mean(x[which(!those),"m"]),mean(x[which(!those),"M"])+mean(x[which(those),"m"]))
  new_mu <- c(old_mu[2],old_mu[1]) #
  meanTCN <- sum(new_mu)
  null_mean <- c(0.5,0.5)*meanTCN
  sigma <- matrix(apply(apply(x,1,function(xi) ((xi-old_mu)%*%t(xi-old_mu))),1,sum)/(nrow(x)-1),nrow=2)
  alpha <- 0.5

  stopifnot(eigen(sigma)$values>c(0,0))

  # set some parameters for this function
  restart = FALSE
  mag = 10
  while (rem_iterations > 0 & sum(abs(old_mu-new_mu))>10^(-30)) {
    if(restart){
      return(infer_PSACN(major_minor_pairs,rem_iterations,restarts+1))
    }
    # the likelihood of each observation under the current parameter estimates are:
    stopifnot(eigen(sigma)$values>c(0,0))
    comp <- c(0.5,0.5)
    comp <- try(cbind(alpha*dmvnorm(as.matrix(x), as.matrix(new_mu), as.matrix(sigma)), (1-alpha)*dmvnorm(as.matrix(cbind(x[,c(2,1)])),as.matrix(new_mu), as.matrix(sigma))),silent=TRUE)

    # the relative likelihood of each estimate being in each component is:
    p <- comp/apply(comp,1,sum)
    # then the updated estimate for alpha is:
    alpha=sum(p[,1])/nrow(x)
    sol <- alpha
    sol <- try(solve(polynomial(c(alpha,1,-6,4))),silent=TRUE)
    alpha <- try(sol[which(sol>0 & sol < 1)],silent=TRUE)
    # the updated estiamte for mu is:
    old_mu <- new_mu
    new_mu <- apply((p[,2]*cbind(x[,2],x[,1]) + p[,1]*x),2,sum)/nrow(x)

    new_mu <- try(solve(nrow(x)*nullsigma%*%solve(sigma)+diag(c(1,1)))%*%(nrow(x)*nullsigma%*%solve(sigma)%*%new_mu+null_mean),silent=TRUE)

    N1 <- apply(x,1,function(xi) ((xi-new_mu)%*%t(xi-new_mu)))
    N2 <- apply(x,1,function(xi) ((c(xi[2],xi[1])-new_mu)%*%t(c(xi[2],xi[1])-new_mu)))
    sigma <- matrix(apply(rep(p[,1],each=4)*N1 +rep(p[,2],each=4)*N2,1,sum)/(nrow(x)-1),nrow=2)
    sigma <- (importance*nullsigma+nrow(x)*sigma)/(importance+nrow(x))

    stopifnot(eigen(sigma)$values>c(0,0))
    rem_iterations <- rem_iterations - 1
  }

  #print(mardiaTest)
  #print(hzTest)
  mardiapvalue = -1
  mardiapvalue = try(mardiaTest(x)$p.value.small,silent=TRUE)
  hzpvalue = -1
  hzpvalue = try(hzTest(x)$p.value,silent=TRUE)

  ABmean = apply(x,2,mean)
  ABcov = cov(x)

  #return(c("mu",new_mu,"sigma",sigma,"alpha",alpha,"it",rem_iterations,"N",nrow(x),"ABABABABABAB","mu",ABmean,"sigma",ABcov,"p1",mardiapvalue,"p2",hzpvalue))
  return(c(new_mu,sigma,alpha,rem_iterations,nrow(x),"ABABABABABAB",ABmean,ABcov,mardiapvalue,hzpvalue))
}

# Want to write a funciton to maximise alpha over its domain:
maximise_alpha <- function(p,n){
  max_likelihood = 0
  max_val = -1
  for(k in 0:n){
    if(k>1 & k<n-1){
      new_like = sum(sapply(1:n,function(x){lgamma(n-1) + (- lgamma(k) - lgamma(n-1-k)) * p[x,1] +( - lgamma(k-1) - lgamma(n-k)) * p[x,2]}))
    }else if(k>0 & k<n-1){
      new_like = sum(sapply(1:n,function(x){lgamma(n-1) + (-lgamma(k) - lgamma(n-1-k)) * p[x,1]}))
    } else if(k==0) {
      new_like = 0
    }
    if(new_like > max_likelihood){
      max_likelihood = new_like
      max_val = k
    }
  }
  return(max_val/n)
}

library(polynom)

rescale <- function(dat,rescale_to_one=TRUE){
  main_pts <- get_main_points3(dat,"major","minor",0.05,0.001)
  if(rescale_to_one){
    main_pts <- main_pts[which(abs(main_pts$major-1)<0.3&abs(main_pts$major-1)<0.3),]
  }
  most_freq <- main_pts[which(main_pts$freq == max(main_pts$freq)),c("major","minor")]
  scaling <- as.numeric(most_freq[which(apply(most_freq-1,1,abs)==min(apply(most_freq-1,1,abs)))]) # select the scaling of the two that is already closes to 1.
  dat$major <- dat$major/scaling
  dat$minor <- dat$minor/scaling
  return(dat)
}

EMA <- function (x, n=10){stats::filter(x,rep(1/n,n), sides=2)}
SMA <- function (x, n=10){runMean( x, n )}


# the closest function finds the row of the closest point in a data frame to a point provided given column headings.
closest <- function(x,major_name,minor_name,main_pts){
  final=which(abs(main_pts[,c(major_name)]-x[1])+abs(main_pts[,c(minor_name)]-x[2])==min(abs(main_pts[,c(major_name)]-x[1])+abs(main_pts[,c(minor_name)]-x[2])))
  return(final)
}

# the indicator function. I suspect that the TRUE FALSE type in R is of this version anyway - i.e. you don't need this function.
indicator<-function(condition) ifelse(condition,1,0)

# the number_ticks functions helps to display xticks when producing CN tracks:
number_ticks <- function(n) {function(limits) pretty(limits, n)}

# a function to calculate the distance between points on a grid plot. This is not the euclidean distance since we are calculating the distance between adjacent genomic states which we expect to move in horizontal and verticle paths. Using this distance function over the euclidean distance seemed to give better results in the algorithms that use it.
distance <- function(X1,X2){
  X1 <- unname(X1)
  X2 <- unname(X2)
  x1 =X1[1]
  y1=X1[2]
  x2 = X2[1]
  y2 = X2[2]
  return ((abs(x1-x2)+abs(y1-y2)))
}

# # get_main_points3 finds a crude approximation to the locations of the 'main points' on the grid plot. The crudeness of the approximation is set by two parameters, merge_pts_thresh which is the upper bound on precision of CN estimates, and freq_thresh which is a lower bound to how many SNPs can constitute one of these main points. Please seek examples of this function being used to have an idea of typical parameters.
# get_main_points3 <- function(dat,major_name,minor_name,merge_pts_thresh,freq_thresh){
#   dat$major <- dat[,major_name]
#   dat$minor <- dat[,minor_name]
#
#   main_pts <- dat[,c(major_name,minor_name,"length","length2","chr","start","end")]
#   main_pts$cluster <- main_pts$segment+1
#   main_pts["freq"] <- main_pts$length/sum(main_pts$length)
#   main_pts$TCN <- main_pts[,major_name]+main_pts[,minor_name]
#
#   # calculate all pairwise distances efficiently:
#   d <- as.matrix(dist(main_pts[,c(major_name,minor_name)]))
#   to_delete = order(d)
#
#   for (k in 1:length(to_delete)){
#     # calculate the array location given the 1d location in that array:
#     A = (to_delete[k]-1)%/%nrow(main_pts)+1
#     B = (to_delete[k]-1)%%nrow(main_pts)+1
#     i=min(A,B)
#     j=max(A,B)
#     if ((main_pts[i,"cluster"] != "0") & (main_pts[j,"cluster"] != "0") & (i!=j) & (!is.na(main_pts[i,major_name]) & !is.na(main_pts[i,minor_name]) & !is.na(main_pts[j,major_name]) & !is.na(main_pts[j,minor_name]))){
#       if (d[i,j] > merge_pts_thresh){
#         break
#       }
#       # then merge these two data points as we have guessed that they are the same allelic states in truth!
#       main_pts[j,major_name] = (main_pts[i,major_name]*main_pts[i,"length"]+main_pts[j,major_name]*main_pts[j,"length"])/(main_pts[i,"length"]+main_pts[j,"length"])
#       main_pts[j,minor_name] = (main_pts[i,minor_name]*main_pts[i,"length"]+main_pts[j,minor_name]*main_pts[j,"length"])/(main_pts[i,"length"]+main_pts[j,"length"])
#       main_pts[j,"length"]<- main_pts[i,"length"] + main_pts[j,"length"]
#       main_pts[j,"freq"]<- main_pts[i,"freq"] + main_pts[j,"freq"]
#       main_pts[i,"cluster"] <- "0"
#     }
#   }
#   main_pts <- subset(main_pts, cluster!= 0)
#   main_pts <- main_pts[which(main_pts$freq>freq_thresh),]
#   main_pts <- main_pts[order(main_pts$TCN),]
#
#   # Take the location of the main points to be the median - weighted by the number of SNPs in each segment. This prevents outliers from unnesecarily affecting our ability to remove bias:
#   for (i in 1:nrow(main_pts)){
#     close = dat[which(abs(dat$major-main_pts[i,"major"])+abs(dat$minor-main_pts[i,"minor"])<merge_pts_thresh),]
#     close_rep <- close[rep(1:nrow(close), times = max(10,ceiling(close$length/10))),]
#     main_pts[i,"major_med"] <- median(close_rep$major)
#     main_pts[i,"minor_med"] <- median(close_rep$minor)
#   }
#
#   main_pts$major <- main_pts$major_med
#   main_pts$minor <- main_pts$minor_med
#   keep_going = TRUE
#   while(keep_going){
#     break_out = FALSE
#     for (i in 1:(nrow(main_pts)-1)){
#       for (j in (i+1):(nrow(main_pts))){
#         if (abs(main_pts[i,"major"]-main_pts[j,"major"]) < merge_pts_thresh/2 & abs(main_pts[i,"minor"]-main_pts[j,"minor"]) < merge_pts_thresh/2){
#           if (main_pts[i,"freq"] < main_pts[j,"freq"]){
#             main_pts <- main_pts[-i,]
#             break_out = TRUE
#             break
#           } else{
#             main_pts <- main_pts[-j,]
#             break_out = TRUE
#             break
#           }
#         }
#       }
#       if (break_out){
#         break
#       }
#     }
#     if (!break_out){
#       keep_going = FALSE
#     }
#   }
#   return(main_pts)
# }


# get_main_points3 finds a crude approximation to the locations of the 'main points' on the grid plot. The crudeness of the approximation is set by two parameters, merge_pts_thresh which is the upper bound on precision of CN estimates, and freq_thresh which is a lower bound to how many SNPs can constitute one of these main points. Please seek examples of this function being used to have an idea of typical parameters.
get_main_points3 <- function(dat,major_name,minor_name,merge_pts_thresh,freq_thresh){
  # dat<- mat
  dat$major <- dat[,major_name]
  dat$minor <- dat[,minor_name]

  d <- daisy(dat[,c("major","minor")],weights=dat$length,metric="euclidean")
  #d <- dist(dat[,c("major","minor")])#,weights=dat$length,metric="euclidean")
  hc <- hclust(d,method="complete")                # apply hirarchical clustering
  #plot(hc)
  #rect.hclust(hc, h=merge_pts_thresh, border="red")
  memb <- cutree(hc, h=merge_pts_thresh)
  dat$memb <- memb
  #dat <- dat[rep(1:nrow(dat), times = max(10,ceiling(dat$length/50))),]
  #main_pts <- dat[,c("segment",major_name,minor_name,"length","length2","chr","start","end")]
  for(i in levels(factor(dat$memb))){
    major <- rep(weightedMedian(dat[which(dat$memb == i),c("major")],dat[which(dat$memb == i),c("length")]),nrow=length(dat[which(dat$memb == i),c("major")]))
    dat[which(dat$memb == i),c("major")] <- major
    minor <- rep(weightedMedian(dat[which(dat$memb == i),c("minor")],dat[which(dat$memb == i),c("length")]),nrow=length(dat[which(dat$memb == i),c("minor")]))
    dat[which(dat$memb == i),c("minor")] <- minor
    dat[which(dat$memb == i),c("TCN")] <- minor+major
    len <- rep(sum(dat[which(dat$memb == i),c("length")]),nrow=length(dat[which(dat$memb == i),c("length")]))
    dat[which(dat$memb == i),c("length")] <- len

  }
  main_pts <- dat[,c(major_name,minor_name,"major","minor","length","length2","chr","start","end","memb","TCN")]
  main_pts <- main_pts[order(main_pts$memb),]
  for( i in 2:nrow(main_pts)){
    if(main_pts[i,"memb"] == main_pts[i-1,"memb"]){
      main_pts[i-1,"memb"] = -1
    }
  }
  main_pts <- main_pts[which(main_pts$memb != -1),]
  main_pts$freq <- main_pts$length/sum(main_pts$length)
  main_pts <- main_pts[which(main_pts$freq>freq_thresh),]
  #
  #   main_pts$cluster <- main_pts$segment+1
  #   main_pts["freq"] <- main_pts$length/sum(main_pts$length)
  #   main_pts$TCN <- main_pts[,major_name]+main_pts[,minor_name]
  #
  #   # calculate all pairwise distances efficiently:
  #   d <- as.matrix(dist(main_pts[,c(major_name,minor_name)]))
  #   to_delete = order(d)
  #
  #   for (k in 1:length(to_delete)){
  #     # calculate the array location given the 1d location in that array:
  #     A = (to_delete[k]-1)%/%nrow(main_pts)+1
  #     B = (to_delete[k]-1)%%nrow(main_pts)+1
  #     i=min(A,B)
  #     j=max(A,B)
  #     if ((main_pts[i,"cluster"] != "0") & (main_pts[j,"cluster"] != "0") & (i!=j) & (!is.na(main_pts[i,major_name]) & !is.na(main_pts[i,minor_name]) & !is.na(main_pts[j,major_name]) & !is.na(main_pts[j,minor_name]))){
  #       if (d[i,j] > merge_pts_thresh){
  #         break
  #       }
  #       # then merge these two data points as we have guessed that they are the same allelic states in truth!
  #       main_pts[j,major_name] = (main_pts[i,major_name]*main_pts[i,"length"]+main_pts[j,major_name]*main_pts[j,"length"])/(main_pts[i,"length"]+main_pts[j,"length"])
  #       main_pts[j,minor_name] = (main_pts[i,minor_name]*main_pts[i,"length"]+main_pts[j,minor_name]*main_pts[j,"length"])/(main_pts[i,"length"]+main_pts[j,"length"])
  #       main_pts[j,"length"]<- main_pts[i,"length"] + main_pts[j,"length"]
  #       main_pts[j,"freq"]<- main_pts[i,"freq"] + main_pts[j,"freq"]
  #       main_pts[i,"cluster"] <- "0"
  #     }
  #   }
  #   main_pts <- subset(main_pts, cluster!= 0)
  #   main_pts <- main_pts[which(main_pts$freq>freq_thresh),]
  #   main_pts <- main_pts[order(main_pts$TCN),]
  #
  #   # Take the location of the main points to be the median - weighted by the number of SNPs in each segment. This prevents outliers from unnesecarily affecting our ability to remove bias:
  #   for (i in 1:nrow(main_pts)){
  #     close = dat[which(abs(dat[,major_name]-main_pts[i,major_name])+abs(dat[,minor_name]-main_pts[i,minor_name])<merge_pts_thresh),]
  #     close_rep <- close[rep(1:nrow(close), times = max(10,ceiling(close$length/10))),]
  #     main_pts[i,"major_med"] <- median(close_rep[,major_name])
  #     main_pts[i,"minor_med"] <- median(close_rep[,minor_name])
  #   }
  keep_going = TRUE
  while(keep_going){
    break_out = FALSE
    for (i in 1:(nrow(main_pts)-1)){
      for (j in (i+1):(nrow(main_pts))){
        if (abs(main_pts[i,"major"]-main_pts[j,"major"]) < merge_pts_thresh & abs(main_pts[i,"minor"]-main_pts[j,"minor"]) < merge_pts_thresh){
          if (main_pts[i,"freq"] < main_pts[j,"freq"]){
            main_pts <- main_pts[-i,]
            break_out = TRUE
            break
          } else{
            main_pts <- main_pts[-j,]
            break_out = TRUE
            break
          }
        }
      }
      if (break_out){
        break
      }
    }
    if (!break_out){
      keep_going = FALSE
    }
  }
  return(main_pts)
}

# rAB is a function used infer the bias in total copy number - see methods for the explanation and derivation.
rAB <- function(major,minor,major_next,minor_next){
  total <- major+minor
  total_next <- major_next+minor_next
  if (abs(major-major_next)>abs(minor-minor_next)){
    return((minor/minor_next)*(total_next/total))
  } else {
    return((major/major_next)*(total_next/total))
  }
}

#rAB2 is the same as rAB except that the ratio is always calculated moving from low total copy number to high total copy number.
# the idea behind this is that the bias is monotonic. rAB is deprecated.
rAB2 <- function(major,minor,major_next,minor_next){
  if(major+minor > major_next+minor_next){
    tempmajor = major
    tempminor = minor
    major = major_next
    minor = minor_next
    major_next = tempmajor
    minor_next = tempminor
  }
  total <- major+minor
  total_next <- major_next+minor_next
  if (abs(major-major_next)>abs(minor-minor_next)){
    return((minor/minor_next)*(total_next/total))
  } else {
    return((major/major_next)*(total_next/total))
  }
}

# get_edges find evidence for paths in the genome given a set of main_points and the adjacency information of the segments in main_paths. The path_thresh is an upper bound on precision when searching for paths between the main points.
# as far as I can see there is no use for get_edges anymore.
get_edges <- function(main_paths,main_points,path_thresh){
  relationships <- matrix(0,nrow(main_points),nrow(main_points))
  relationship_count <- matrix(0,nrow(main_points),nrow(main_points))
  for (i in 1:nrow(main_paths)){
    for (j in 1:nrow(main_points)){
      for (k in 1:nrow(main_points)){
        # if the main path goes from main point j to main point k add the weight!
        dist1 <- distance(data.matrix(main_paths[i,c("major1","minor1")]),data.matrix(main_points[j,c("major","minor")]))
        dist2 <- distance(data.matrix(main_paths[i,c("major2","minor2")]),data.matrix(main_points[k,c("major","minor")]))
        # to prevent information travelling across centromeres uncommen this line - however it did not provide much use:
        # cent = c[which(main_paths[nrow(main_paths),"Chromosome1"]==c$chr),"x"]
        # if (!((main_paths[i,"end1"]<cent & main_paths[i,"start2"]<cent) | (main_paths[i,"end2"]>cent & main_paths[i,"start1"]>cent))){ next }
        if (dist1 < path_thresh & dist2 <path_thresh & main_paths[i,"Chromosome1"] == main_paths[i,"Chromosome2"]){
          relationships[j,k] <- relationships[j,k] + main_paths[i,"weights"]
          relationship_count[j,k] <- relationship_count[j,k] +1
        }
      }
    }
  }
  # put all the info in the upper triangle only:
  relationships <- relationships+t(relationships)
  relationship_count <- relationship_count + t(relationship_count)
  relationships[which(lower.tri(relationships,diag=TRUE))] <- 0
  relationship_count[which(lower.tri(relationships,diag=TRUE))] <- 0
  relationships <- relationships+t(relationships)
  relationship_count <- relationship_count + t(relationship_count)
  return(list(relationships=relationships,relationship_count=relationship_count))
}

# cluster_PSACN is a nessecary function to reduce the noise in the dataset. The aim of this script in general is to make the data more accurate, however increasing acuracy is confounded by a a lack of precision is is more often a typical problem. Segments which constitute few SNPs and are hence the least precise can add so much confusion to the grid plot that our bias removal algorithm can get confused. In practise the bias removal algorithm is substantially better on some data when applied to the data that has been smoothed by cluster_PSACN.
# Furthermore, we can use this function not just to increase precision prior to bias removal, we can apply it following bias removal to greatly improve the appearence of our CN plots. Hence you will see this function being used at multiple stages through the script.
# How does it work? We have many pieces of information on how to smooth the data. The idea of the following function is to nudge the data towards what we expect to be the true state given the information at our disposal by sharing infromation across the genome allowing the CN track to converge on a simplified state. We know, because of the way we preprocessed the data, when a segment is in allelic balance and nudge it toward allelic balance. We know, or at least we can guess well, when there is a CN change in only one of the two chromosomes in truth and we can straighten out the respective line on the grid plot. Segments that are adjacent to each other, and are within a really small distance of each other on the grid plot as determined by our distance function are probably the same state in truth and we nudge them toward each other. We also expect and observe the same allelic states to repeat across the genome and hence we can share the local infromation - allelic balance, straightness, and closeness at a global level by imparting the small adjustment required in one location of the grid plot to another that shares a similar situation.
# The use of the following funciton should however be used with caution. Setting the parameters to be too loose on the data will allow the data to collapse and converge on an allelic state that is much more simplified than what exists in truth.
# separate is a modern addition which allows the convergence onto both copy number states at the same time.
cluster_PSACN <- function(dat,iterations,min_length,local_adj_thresh,global_adj_thresh,false_to_true_thresh,p_straighten,p_close_local,p_close_global,print_plots,small_length,allelic_balance_info="present",seperate=FALSE){
  print(paste("Smoothing to take place on the data for ",as.character(iterations)," iterations with the following parame
              ters:",sep=""))
  print(paste("    Minimum number of SNPs: ",as.character(min_length),sep=""))
  print(paste("    CN change ratio of true change to hypothetical false change: ",as.character(false_to_true_thresh),sep
              =""))
  print(paste("    Local noise threshold: ",as.character(local_adj_thresh),sep=""))
  print(paste("    Global noise threshold: ",as.character(global_adj_thresh),sep=""))
  print(paste("    Rate of convergence to straight CN changes: ",as.character(p_straighten),sep=""))
  print(paste("    Rate of convergence to a locally smooth state: ",as.character(p_close_local),sep=""))
  print(paste("    Rate of convergence to a globally smooth state: ",as.character(p_close_global),sep=""))

  # take only the subset of data that constitute more than min_length SNPs. In practise there is no need to take such a subset of data and a good min_length is 0.
  dat2 <- dat[which(dat$length>min_length),]
  dat2 <- unphase(dat2)
  # label each segment
  dat2$cluster <- 1:nrow(dat2)
  # iteratively increase precision:
  for (jiter in 1:iterations){
    #print(paste("Smoothing iteration ",as.character(jiter),sep=""))
    if (allelic_balance_info == "present"){
      # force points that we know to be in allelic balance to actually be so
      for (i in 1:nrow(dat2)){
        if (dat2[i,"allelic_balance"]){
          dat2[i,c("major")] <- dat2[i,c("minor")] <- mean(dat2[i,c("major")],dat2[i,c("minor")])
        }
      }
    }
  }

  # let the CN estimates converge locally if below the specified local threshold
  for (i in 2:nrow(dat2)){
    if (abs(dat2[i,"major"]-dat2[i-1,"major"])<local_adj_thresh/2 & abs(dat2[i,"minor"]-dat2[i-1,"minor"])<local_adj_thresh/2 & dat2[i-1,"chr"] == dat2[i,"chr"]){
      dat2[i,"major"] <- mean(dat2[i,"major"],dat2[i-1,"major"])
      dat2[i,"minor"] <- mean(dat2[i,"minor"],dat2[i-1,"minor"])
      dat2[i,"length"] <- dat2[i,"length"]+dat2[i-1,"length"]
      dat2[i,"length2"] <- dat2[i,"length2"]+dat2[i-1,"length2"]
      dat2[i,"start"] <- min(dat2[i,"start"],dat2[i-1,"start"])
      dat2[i,"end"] <- max(dat2[i,"end"],dat2[i-1,"end"])
      dat2[i-1,"cluster"]<- 0
    }
  }
  dat2 <- subset(dat2, cluster!= 0)

  # phase the SNPs:
  dat2 <- phase(dat2)

  # make some plots to observe the progress if desired:
  if (print_plots){
    print(paste("    The grid plot prior to iteration ",as.character(jiter)," is:",sep=""))
    print(get_grid_plot(dat2,MIN_LENGTH))
    print(paste("    The CN track prior to iteration ",as.character(jiter)," is:",sep=""))
    print(get_CN_track(dat2))
  }

  #force adjacent segments that have a true to false CN change less than the specified false_to_true_thresh to have a CN change on only one chromosome:
  for (i in 2:nrow(dat2)){
    # add a little bit to the jump if it is percieved to be zero to prevent division by zero:
    major_jump = max(abs(dat2[i-1,"major"]-dat2[i,"major"]),0.000001)
    minor_jump = max(abs(dat2[i-1,"minor"]-dat2[i,"minor"]),0.000001)
    # the ratio of the jumps in the major and minor axis is given by r1:
    r1 <- abs(major_jump/minor_jump)

    # if the adjacent segments are on differen't chromosomes, move on:
    if(dat2[i,"chr"] != dat2[i-1,"chr"]){ next }

    if (r1 > false_to_true_thresh){
      # then the change in the minor axis is false
      minor <- (dat2[i-1,"minor"]*dat2[i-1,"length"]+dat2[i,"minor"]*dat2[i,"length"])/(dat2[i,"length"]+dat2[i-1,"length"])
      dat2[i-1,"minor"] <- dat2[i-1,"minor"]*(1-p_straighten)+minor*p_straighten
      dat2[i,"minor"] <- dat2[i,"minor"]*(1-p_straighten)+minor*p_straighten
    }

    if (1/r1 > false_to_true_thresh) {
      # then the change in the major axis is false
      major <- (dat2[i-1,"major"]*dat2[i-1,"length"]+dat2[i,"major"]*dat2[i,"length"])/(dat2[i,"length"]+dat2[i-1,"length"])
      dat2[i-1,"major"] <- dat2[i-1,"major"]*(1-p_straighten)+major*p_straighten
      dat2[i,"major"] <- dat2[i,"major"]*(1-p_straighten)+major*p_straighten
    }
  }

  #find major/minor estimates again:
  dat2 <- unphase(dat2)
  dat2$TCN <- dat2$major + dat2$minor
  a3 <- as.matrix(dist(dat2$TCN))
  d3 <- which(a3<local_adj_thresh)-1
  for (d in d3){
    #print(c(i,j))
    i <- d%%nrow(a3)+1
    j <- d%/%nrow(a3)+1
    diff = 1#exp(-abs(dat2[j,"minor"]-dat2[i,"minor"])/adj_lim)
    TCN <- (dat2[j,"TCN"]*dat2[j,"length"]+dat2[i,"TCN"]*dat2[i,"length"])/(dat2[i,"length"]+dat2[j,"length"])
    dat2[j,"TCN"] <- dat2[j,"TCN"]*(1-p_close_local)+p_close_local*TCN
    dat2[i,"TCN"] <- dat2[i,"TCN"]*(1-p_close_local)+p_close_local*TCN
    dat2[i,c("major","minor")] <- dat2[i,c("major","minor")]*TCN/sum(dat2[i,c("major","minor")])
  }

  # calculate the global distances between points such that we can allow convergence globally after the application of the local rules allowing the sharing of information of what is learnt in these rules.
  a1 <- as.matrix(dist(dat2$major))
  a2 <- as.matrix(dist(dat2$minor))
  a3 <- as.matrix(dist(dat2[,c("major","minor")]))
  d1 <- which(a1<global_adj_thresh)-1
  d2 <- which(a2<global_adj_thresh)-1
  d3 <- which(a3<global_adj_thresh)-1
  M <- dat2$major-dat2$minor

  # let the adjacent CN estimates converge locally if below the specified local threshold
  for (i in 2:nrow(dat2)){
    local_adj_thresh2 <- local_adj_thresh/min(1,dat2[i,"length"]/small_length,dat2[i-1,"length"]/small_length)
    if (abs(dat2[i,"major"]-dat2[i-1,"major"])<local_adj_thresh2 & abs(dat2[i,"minor"]-dat2[i-1,"minor"])<local_adj_thresh2 & dat2[i-1,"chr"] == dat2[i,"chr"]){
      minor <- (dat2[i-1,"minor"]*dat2[i-1,"length"]+dat2[i,"minor"]*dat2[i,"length"])/(dat2[i,"length"]+dat2[i-1,"length"])
      dat2[i,"minor"] <- minor*p_close_local+(1-p_close_local)*dat2[i,"minor"]
      dat2[i-1,"minor"] <- minor*p_close_local+(1-p_close_local)*dat2[i-1,"minor"]

      major <- (dat2[i-1,"major"]*dat2[i-1,"length"]+dat2[i,"major"]*dat2[i,"length"])/(dat2[i,"length"]+dat2[i-1,"length"])
      dat2[i,"major"] <- major*p_close_local+(1-p_close_local)*dat2[i,"major"]
      dat2[i-1,"major"] <- major*p_close_local+(1-p_close_local)*dat2[i-1,"major"]

    }
  }

  # let the off by one semi adjacent CN estimates converge locally if below the specified local threshold. I added this step in since it was common to observe a CN state on adjacent sides of a segment to be very similar. This makes sense since everything starts at a CN state of 1,1 in truth and then some window of the genome undergoes a CN change. Unless this CN change hits a region of LOH or the end of the chromosome it is liekly that we observe the same CN state on either side of the window that changed!
  for(i in 3:nrow(dat2)){
    local_adj_thresh2 <- local_adj_thresh/min(1,dat2[i,"length"]/small_length,dat2[i-1,"length"]/small_length)
    if (abs(dat2[i,"major"]-dat2[i-2,"major"])<local_adj_thresh2 & abs(dat2[i,"minor"]-dat2[i-2,"minor"])<local_adj_thresh2 & dat2[i-2,"chr"] == dat2[i,"chr"]){
      minor <- (dat2[i-2,"minor"]*dat2[i-2,"length"]+dat2[i,"minor"]*dat2[i,"length"])/(dat2[i,"length"]+dat2[i-2,"length"])
      dat2[i-2,"minor"] <- dat2[i-2,"minor"]*(1-p_close_local)+p_close_local*minor
      dat2[i,"minor"] <- dat2[i,"minor"]*(1-p_close_local)+p_close_local*minor

      major <- (dat2[i-2,"major"]*dat2[i-2,"length"]+dat2[i,"major"]*dat2[i,"length"])/(dat2[i,"length"]+dat2[i-2,"length"])
      dat2[i-2,"major"] <- dat2[i-2,"major"]*(1-p_close_local)+p_close_local*major
      dat2[i,"major"] <- dat2[i,"major"]*(1-p_close_local)+p_close_local*major
    }
  }

  # let the CN estimates converge globally if below the estimates prior to this round of nudging were below the specified global threshold on each axis:
  if(seperate){
    for (d in d1){
      i <- d%%nrow(a1)+1
      j <- d%/%nrow(a1)+1
      major <- (dat2[j,"major"]*dat2[j,"length"]+dat2[i,"major"]*dat2[i,"length"])/(dat2[i,"length"]+dat2[j,"length"])
      dat2[j,"major"] <- dat2[j,"major"]*(1-p_close_global)+p_close_global*major
      dat2[i,"major"] <- dat2[i,"major"]*(1-p_close_global)+p_close_global*major
    }
    for (d in d2){
      i <- d%%nrow(a2)+1
      j <- d%/%nrow(a2)+1
      diff = 1#exp(-abs(dat2[j,"minor"]-dat2[i,"minor"])/adj_lim)
      minor <- (dat2[j,"minor"]*dat2[j,"length"]+dat2[i,"minor"]*dat2[i,"length"])/(dat2[i,"length"]+dat2[j,"length"])
      dat2[j,"minor"] <- dat2[j,"minor"]*(1-p_close_global)+p_close_global*minor
      dat2[i,"minor"] <- dat2[i,"minor"]*(1-p_close_global)+p_close_global*minor
    }
  } else{
    for (d in d3){
      i <- d%%nrow(a3)+1
      j <- d%/%nrow(a3)+1
      diff = 1#exp(-abs(dat2[j,"minor"]-dat2[i,"minor"])/adj_lim)
      Mm <- (dat2[j,c("major","minor")]*dat2[j,"length"]+dat2[i,c("major","minor")]*dat2[i,"length"])/(dat2[i,"length"]+dat2[j,"length"])
      dat2[i,c("major","minor")] <- dat2[j,c("major","minor")]*(1-p_close_global)+p_close_global*Mm
      dat2[i,c("major","minor")] <- dat2[j,c("major","minor")]*(1-p_close_global)+p_close_global*Mm
    }

    # let the CN estimates across chromosomes converge if the estimates prior to this rround of nudging were below the specified local threshold:
    for (i in 1:length(M)){
      if (abs(M[i])<local_adj_thresh){
        avg <- (dat2[i,c("minor")]+dat2[i,c("major")])/2*p_close_local
        dat2[i,c("major","minor")]<-cbind(avg,avg)+dat2[i,c("major","minor")]*(1-p_close_local)
      }
    }
  }
  dat2$t <- dat2$major + dat2$minor
  return(dat2)
}



  # the phase function, has a best guess at phasing the SNPs based upon the simple idea of making the lines in the grid plot as straight as possible.
phase <- function(mat){
  for (i in 2:nrow(mat)){
    major_jump = max(abs(mat[i-1,"major"]-mat[i,"major"]),0.0001)
    minor_jump = max(abs(mat[i-1,"minor"]-mat[i,"minor"]),0.0001)
    major_jump_alt = max(abs(mat[i-1,"major"]-mat[i,"minor"]),0.0001)
    minor_jump_alt = max(abs(mat[i-1,"minor"]-mat[i,"major"]),0.0001)
    if(mat[i,"chr"] != mat[i-1,"chr"]){ next }
    if (min(major_jump_alt,minor_jump_alt) < min(major_jump,minor_jump) & mat[i-1,"chr"] == mat[i,"chr"]){
      mat[i,c("major","minor")]<-mat[i,c("minor","major")]
    }
  }
  return(mat)
}


  # get_CN_track plots a copy number track of the data along with vertical lines for the centromeres.
get_CN_track <- function(mat){
  g <- ggplot(mat) +
    geom_segment(aes_string(x = "start", y = "minor", xend = "end", yend = "minor",col="chr"),size=2)+
    geom_segment(aes_string(x = "start", y = "major", xend = "end", yend = "major"),size=1)+
    geom_vline(data=centromeres,aes(xintercept=x,col=chr))+
    facet_grid(.~chr,scales = "free_x", space = "free")+
    scale_x_continuous(breaks=seq(0,3*10^9,50*10^6))+ #make 50mb ticks...
    theme(axis.text.x = element_blank(),legend.position="none",legend.background = element_rect(fill = "white"))+
    #theme_bw()+
    xlab("50 MB ticks")+
    scale_y_continuous(breaks=number_ticks(20))+
    ylab("major/minor CN")+ggtitle("CN track")
  return(g)
}

# a second, perhaps easier to view
get_CN_track_grid <- function(qat){
  qat <- droplevels(qat)
  g <- ggplot(qat) +
    geom_segment(aes_string(x = "start", y = "minor", xend = "end", yend = "minor",col="chr"),size=2)+
    geom_segment(aes_string(x = "start", y = "major", xend = "end", yend = "major"),size=1)+
    geom_vline(data=centromeres,aes(xintercept=x,col=chr))+
    facet_wrap(~chr,scales = "free_x")+
    scale_x_continuous(breaks=seq(0,3*10^9,50*10^6))+ #make 50mb ticks...
    scale_y_continuous(breaks=number_ticks(5))+theme_bw()+
    #theme(axis.text.x = element_blank(),legend.position="none",legend.background = element_rect(fill = "white"))+
    xlab("50 MB ticks")+
    ylab("major/minor CN")+ggtitle("CN track")
  return(g)
}

get_TCN_track <- function(mat){
  mat$total <- mat$major + mat$minor
  g <- ggplot(mat) +
    geom_segment(aes_string(x = "start", y = "total", xend = "end", yend = "total",col="chr"),size=2)+
    geom_vline(data=centromeres,aes(xintercept=x,col=chr))+
    facet_grid(.~chr,scales = "free_x", space = "free")+
    scale_x_continuous(breaks=seq(0,3*10^9,50*10^6))+ #make 50mb ticks...
    theme(axis.text.x = element_blank(),legend.position="none",legend.background = element_rect(fill = "white"))+    xlab("50 MB ticks")+
    scale_y_continuous(breaks=number_ticks(20))+
    ylab("TCN")+ggtitle("TCN track")
  return(g)
}

# get_grid_plot plots the data as a grid plot...
# the grid plot appearance is dramatically thrown off by the presence of small subclones.
get_grid_plot <- function(dat2,MIN_LENGTH){
  g <- ggplot(dat2[which(dat2$length>MIN_LENGTH),],aes(major,minor))+geom_point(aes(col=chr,size=length))+geom_path(aes(col=chr))+coord_fixed()+xlab("CN 1")+ylab("CN 2")+ggtitle("Grid plot")
  return(g)
}

get_grid_plot_no_lines <- function(dat2,MIN_LENGTH){
  g <- ggplot(dat2[which(dat2$length>MIN_LENGTH),],aes(major,minor))+geom_point(aes(col=chr,size=length))+coord_fixed()+xlab("CN 1")+ylab("CN 2")
  return(g)
}


# get_next_node is the function used iteratively to grow our tree of true NSEW paths through the genome. The aim of this script is to resotre a grid structure to the genome. To do so this script iteratively grows a tree on top of the grid plot for which there is evidence that movements only this tree are in truth verticle or horizontal. Growing this tree is easily confounded by a lack of precision in the average cn calls and a lack of phasing information.
  # This function takes a node that is currently on our tree, with the information of all the main points and their relationships (by relationship I mean the prescence of an observation which indicates that two such points on a grid plot differ in allelic state in only one of the two chromosomes) and tries to find another node to add to the tree which differs in the allelic state on only one of the two chromosomes in truth. You will observe many heuristics to prevent bad edges from beign added (a bad edge would be a node for which the alleic state of itselfg and its parent on the tree differed in both chromosonal arms.) The main idea behind these heuristics were the observations that the bias funciton was convex and that the average cn calls were good in the region of 1,1 and 2,0 (the range of allelic state that the SNP arrays were designed for). To further improve the chance of corrctly guesisng an underlying lattice on the grid plot that was actually a grid in truth the algorithm only allows for one movement in each NSEW direction from any node. The idea behind this is that you can be pretty confident in general about adding the first adjacent node to another since you add the node you are most confident about, but by the time you are adding a second, third or fourth you are probably adding in nonsense information.
get_next_node <- function(node,rel,pts,visited,radius_thresh){
  lower_bound <- min(pts[visited,"TCN"])
  upper_bound <- max(pts[visited,"TCN"])
  if (sum(rel[node,]) == 0){
    return(c(node,node,0))
  }

  # first check to see if there relationships funciton found any evidence of further nodes emanating from this point:
  potential <- as.numeric(rel[node,]>0)
  if (sum(potential) == 0){
    return(c(node,node,0))
  }
  # now go through all the potential nodes that we haven't added to our tree and add the next best node in a greedy fashion:
  for (i in 1:nrow(great_pts)){
    # Supposing that the data we observe was a good representation of the underlying data we could suppose that the straightness of a line apriori to bias removal was a good indication of a line that was in truth straight. This helped some datasets a little bit and ruinedthe chance at bias removal on others since data that was badly phased was often observed to be fairly straight. Adding in a badly phased node is dangerous because by definition you have an inverted set of allelic states for the node you are potentially adding in. Also, on some data sets were aided by encouraging the next node to be as far away as possible from the current node in any one given direction since CN changes are composed with themselves and searchign for the greatest jumps at any one stage of the algorithm should allow you to get a broader range of allelic states from which to infer bias. This actually cuased some problems in some datasets and was also hence removed from guiding our choice. The choice of the node to add in next was simply found to be a node which consisted of a large number of SNPs for which there was a lot of evidence of allelic paths between that node and one already on our tree.
    # gradient <- atan(abs((great_pts[node,"minor"]-great_pts[i,"minor"])/(abs(great_pts[node,"major"]-great_pts[i,"major"])+0.0000001)))
    # closeness <- 1 #abs(abs(gradient/pi*2)-0.5)*2
    potential[i] <- potential[i]*great_pts[i,"freq"]*great_pts[node,"freq"]*relationships[min(node,i),max(node,i)]#*closeness*relationship_count[min(node,i),max(node,i)]*(distance(as.matrix(great_pts[node,c("major","minor")]),as.matrix(great_pts[i,c("major","minor")]))-2*radius_thresh)#
    total_i <- great_pts[i,"major"]+great_pts[i,"minor"]
    total_node <- great_pts[node,"major"]+great_pts[node,"minor"]
    if (lower_bound < total_i & total_i < upper_bound){
      # then we have picked a point within the range of points currently explored, go find something different!
      potential[i] = 0
    }
    # the next four if statements ensure we are adding in a node that under a node currenlty on our tree that is in a direction we have not already explored.
    if (great_pts[node,"N"] & abs(great_pts[i,"major"]-great_pts[node,"major"]) < abs(great_pts[i,"minor"]-great_pts[node,"minor"]) & (great_pts[i,"minor"]-great_pts[node,"minor"]) > 0 ){
      potential[i] = 0
    }
    if (great_pts[node,"S"] & abs(great_pts[i,"major"]-great_pts[node,"major"]) < abs(great_pts[i,"minor"]-great_pts[node,"minor"]) & (great_pts[i,"minor"]-great_pts[node,"minor"]) < 0 ){
      potential[i] = 0
    }
    if (great_pts[node,"E"] & abs(great_pts[i,"major"]-great_pts[node,"major"]) > abs(great_pts[i,"minor"]-great_pts[node,"minor"]) & (great_pts[i,"major"]-great_pts[node,"major"])>0){
      potential[i] = 0
    }
    if (great_pts[node,"W"] & abs(great_pts[i,"major"]-great_pts[node,"major"]) > abs(great_pts[i,"minor"]-great_pts[node,"minor"]) & (great_pts[i,"major"]-great_pts[node,"major"])<0){
      potential[i] = 0
    }

    if (i%in%visited){
      # don't visit the same node twice!
      potential[i] = 0
    }
    # don't allow for paths along the either the x or y axis (when one of the two allelic states is zero or close to zero), since it appears that values close ot the axis are snapped there by the illumina preprocessing and that really messes up estimation of bias...
    if (great_pts[node,"minor"] < 0.1 & great_pts[i,"minor"] < 0.1){
      potential[i] = 0
    }
    # as we observe that the bias is convex, ensure that we have this assumption prior to adding a node to our tree:
    if(sign(great_pts[i,"major"]-great_pts[node,"major"]) == sign(great_pts[i,"minor"]-great_pts[node,"minor"]) ){
      potential[i] = 0
    }
  }

  # if there were no potential nodes to add in, finish here.
  if(sum(potential)<=0){
    return(c(node,node,0))
  }
  next_step <- which(potential==max(potential))
  weight <- potential[next_step]
  return(c(node,next_step,weight))
}

jump <- function(great_pts,node,i){
  return(distance(great_pts[node,c("major","minor")],great_pts[i,c("major","minor")]))
}

# unphase merely calculates the major and minor cn values for any segment. By unphasing we can better share information in both the bias removal algorithm and the cluster_PSACN function between symmetrical states (two states such as (x,y) and (y,x)) on the grid plot.
unphase <- function(mat2){
  M <- apply(mat2[,c("major","minor")],1,max)
  m <- apply(mat2[,c("major","minor")],1,min)
  mat2[,c("major")] <- M
  mat2[,c("minor")] <- m
  return(mat2)
}

# some funcitons for traversing the tree in the best path... havn't thouroughly tested yet - should come back to later. The basic principle is to use dynamic programming to recursively build a data structure that tells you how best to traverse the tree such that you move in a path from low TCN to high TCN via the most nodes.
get_path_pts <- function(graph_adj){
  graph_adj2 <- build(graph_adj)
  all_largest = which(graph_adj2 == max(graph_adj2))%/%nrow(graph_adj)+1
  path = c()
  # hence the node we want to start at is largest_j.
  while(TRUE){
    next_node <- all_largest[which(great_pts[all_largest,"freq"]==max(great_pts[all_largest,"freq"]))]
    path = cbind(next_node,path)
    # but how do you get up to next_node?
    if(sum(graph_adj2[,next_node]) == 0){
      break
    }
    all_largest <- which(graph_adj2[,next_node]==max(graph_adj2[,next_node]))
  }
  return(unname(path))
}

# build a data structure to find the best way to tranverse the tree - a dynamic programming approach to tree tranversal.
build <- function(graph_adj){
  for( i in 1:nrow(graph_adj)){
    for (j in 1:nrow(graph_adj)){
      if(graph_adj[i,j]>0){
        for (k in 1:nrow(graph_adj)){
          if (graph_adj[j,k]>0){
            graph_adj[j,k] = graph_adj[j,k]+graph_adj[i,j]
          }
        }
      }
    }
  }
  return(graph_adj)
}



remove_bias <- function(matta,mattn,widthsstart,widthsend,multi){
  mattn <- matta #<- mat
  # widthsstart = 3
  # widthsend = 30
  # multi = 3
  mattab <- matta[which(matta$length >1),]
  pts <- cbind(mattab[1:(nrow(mattab)-1),c("major","minor","chr","length")],mattab[2:nrow(mattab),c("major","minor","chr","length")])
  colnames(pts) <- c("x1","y1","chr1","length1","x2","y2","chr2","length2")
  rownames(pts) <- 1:nrow(pts)
  pts$compareable <- pts$chr1 == pts$chr2

  # pts$grad1 <- (pts$y1 - pts$y2)/(pts$x1-pts$x2)
  # pts$grad2 <- (pts$y1 - pts$x2)/(pts$x1-pts$y2)
  # pts$grad3 <- (pts$x1 - pts$x2)/(pts$y1-pts$y2)
  # pts$grad4 <- (pts$x1 - pts$y2)/(pts$y1-pts$x2)

  #pts <- pts[which(complete.cases())]
  # of all the gradients, take only the ones that have a magnitude of less than one - since there are pairs of gradients that are each others reciprocals
  # grads <- unlist(apply(pts[,c("grad1","grad2","grad3","grad4")],1,function(x) sort(x[which(abs(x)<=1)])[c(1,2)]))
  # pts$g1 <- grads[seq(1,length(grads),2)]
  # pts$g2 <- grads[seq(2,length(grads),2)]
  # convert these gradients to an angle:
  # pts$angle1 <- atan((pts$g1))
  # pts$angle2 <- atan((pts$g2))

  pts$tcnlow <- apply(cbind(pts$x1+pts$y1,pts$x2+pts$y2),1,min)
  pts$tcnhigh <- apply(cbind(pts$x1+pts$y1,pts$x2+pts$y2),1,max)
  pts$center <- apply(cbind(pts$x1+pts$y1,pts$x2+pts$y2),1,mean)
  # these 'rad' values show what the correction would be given the rAB
  pts$rad1 <- apply(pts[,c("x1","y1","x2","y2")],1,function(x) rAB2(x[1],x[2],x[3],x[4]))
  pts$rad2 <- apply(pts[,c("x1","y1","y2","x2")],1,function(x) rAB2(x[1],x[2],x[3],x[4]))
  # pts$radeff1 <- 1/(pts$tcnhigh - pts$tcnlow) * log(1+pts$rad1 - pts$tcnhigh + pts$tcnlow)
  # pts$radeff2 <- 1/(pts$tcnhigh - pts$tcnlow) * log(1+pts$rad2 - pts$tcnhigh + pts$tcnlow)

  pts$logradeff1 <- 1/(pts$tcnhigh - pts$tcnlow) * log(pts$rad1)# - pts$tcnhigh + pts$tcnlow)
  pts$logradeff2 <- 1/(pts$tcnhigh - pts$tcnlow) * log(pts$rad2)# - pts$tcnhigh + pts$tcnlow)
  pts$radeff1 <- exp(pts$logradeff1)
  pts$radeff2 <- exp(pts$logradeff1)

  # the changes are only interesting if they are of the same chromosome:
  pts <- pts[which(pts$chr1 == pts$chr2),]
  #arbritrarily rearrange the values such that
  for(i in 1:nrow(pts)){
    if((!is.nan(pts[i,"grad1"])  & pts[i,"g1"] == pts[i,"grad1"]) | (!is.nan(pts[i,"grad3"])  & pts[i,"g1"] == pts[i,"grad3"])){
      next
      #pts[i,"r1"] = pts[i,"rad1"]
      #pts[i,"r2"] = pts[i,"rad2"]
      #print("A")
    }
    else{
      temp <- pts[i,"rad2"]
      pts[i,"rad2"] <- pts[i,"rad1"]
      pts[i,"rad1"] <- temp
      temp <- pts[i,"radeff2"]
      pts[i,"radeff2"] <- pts[i,"radeff1"]
      pts[i,"radeff1"] <- temp
    }
  }



  pts$size = 1/(1/pts$length1 +1/pts$length2)

  pts$sizelog = -log(1/pts$length1 +1/pts$length2)#^2
  pts$angle_per1 <- abs(pi/4-pts$angle1)
  pts$angle_per2 <- abs(pi/4-pts$angle2)
  #pts$size7 <- order(pts$size)
  g1 <- ggplot() + geom_segment(data=pts,aes(x=tcnlow,xend=tcnhigh,y=(angle1),yend=(angle1)),alpha=0.1,size=5) + geom_segment(data=pts,aes(x=tcnlow,xend=tcnhigh,y=(angle2),yend=(angle2)),alpha=0.1,size=5)+ theme(legend.position="none")
  g1b <- ggplot() + geom_segment(data=pts,aes(x=tcnlow,xend=tcnhigh,y=(angle_per1),yend=(angle_per1)),alpha=0.1,size=5) + geom_segment(data=pts,aes(x=tcnlow,xend=tcnhigh,y=(angle_per2),yend=(angle_per2)),alpha=0.1,size=5)+ theme(legend.position="none")
  g2 <- ggplot() + geom_segment(data=pts,aes(x=tcnlow,xend=tcnhigh,y=(g1),yend=(g1)),alpha=0.1,size=5) + geom_segment(data=pts,aes(x=tcnlow,xend=tcnhigh,y=g2,yend=g2),alpha=0.1,size=5)+ theme(legend.position="none")
  g3 <- ggplot() + geom_segment(data=pts,aes(x=tcnlow,xend=tcnhigh,y=(rad1),yend=(rad1)),alpha=0.1,size=5) + geom_segment(data=pts,aes(x=tcnlow,xend=tcnhigh,y=(rad2),yend=(rad2)),alpha=0.1,size=5)+ylim(0,5)+ theme(legend.position="none")

  g3b <- ggplot() + geom_segment(data=pts,aes(x=tcnlow,xend=tcnhigh,y=(radeff1),yend=(radeff1)),alpha=0.1,size=5) + geom_segment(data=pts,aes(x=tcnlow,xend=tcnhigh,y=(radeff2),yend=(radeff2)),alpha=0.1,size=5)+ylim(0,2)+theme(legend.position="none")

  g4 <- get_grid_plot(matta,0)+ theme(legend.position="none")
  g5 <- get_CN_track(matta)+ theme(legend.position="none",axis.text.x=element_blank(),axis.text.y=element_blank())
  g6 <- grid.arrange(g4,g5,ncol=2)
  grid.arrange(g1,g1b,g2,g3,g3b,g6,nrow=6)


  ggplot(mat) +
    # geom_segment(aes_string(x = "start", y = "minor", xend = "end", yend = "minor",col="chr"),size=2)+
    # geom_segment(aes_string(x = "start", y = "major", xend = "end", yend = "major"),size=1)+
    geom_segment(aes_string(x = "start", y = "gmmT", xend = "end", yend = "gmmT",col="chr"),size=0.5)+
    geom_segment(aes_string(x = "start", y = "tcnMean", xend = "end", yend = "tcnMean"),size=3)+
    geom_vline(data=centromeres,aes(xintercept=x,col=chr))+
    facet_grid(.~chr,scales = "free_x", space = "free")+
    scale_x_continuous(breaks=seq(0,3*10^9,50*10^6))+ #make 50mb ticks...
    theme(axis.text.x = element_blank(),legend.position="none",legend.background = element_rect(fill = "white"))+
    #theme_bw()+
    xlab("50 MB ticks")+
    scale_y_continuous(breaks=number_ticks(20))+
    ylab("major/minor CN")+ggtitle("CN track")


  pts$jump = pts$tcnhigh - pts$tcnlow
  pts2 <- pts[,c("rad1","radeff1","angle1","size","sizelog","jump","center")]
  pts3 <- pts[,c("rad2","radeff2","angle2","size","sizelog","jump","center")]
  colnames(pts3) <- colnames(pts2)
  pts2 <- rbind(pts2,pts3)

  # should try different rotations, if you rotate
  interval_widths = seq(30,100,10)
  interval_means = seq(0,2,0.01)
  interval_alphas = seq(-pi/4,pi/4,0.01)

  counts <- matrix(rep(0,length(interval_means)*length(interval_alphas)),ncol=length(interval_alphas))
  rownames(counts) <- as.character(interval_means)
  colnames(counts) <- as.character(interval_alphas)

  for(i in 1:length(interval_means)){
    for(j in 1:length(interval_alphas)){
      for(width in 1:length(interval_widths)){
        counts[i,j] = counts[i,j] + sum(pts2[,c("sizelog")] * exp(-((width/pts2[,"jump"]^2)*(1+tan(interval_alphas[j])^2)*abs(pts2[,"radeff1"]-pts2[,"center"]*tan(interval_alphas[j])-interval_means[i]))),na.rm=TRUE)
      }
    }
    print(i/length(interval_means))
  }
  count_temp <- counts
  counts <- count_temp

  double_down <- function(mat){
    return( (mat[1:(nrow(mat)-2),1:(ncol(mat)-2)]+
               mat[2:(nrow(mat)-1),1:(ncol(mat)-2)]+
               mat[3:(nrow(mat)),1:(ncol(mat)-2)]+
               mat[1:(nrow(mat)-2),2:(ncol(mat)-1)]+
               mat[1:(nrow(mat)-2),3:(ncol(mat))]+
               mat[2:(nrow(mat)-1),2:(ncol(mat)-1)]+
               mat[2:(nrow(mat)-1),3:(ncol(mat))]+
               mat[3:(nrow(mat)),2:(ncol(mat)-1)]+
               mat[3:(nrow(mat)),3:(ncol(mat))]))}

  n_down <- function(mat,count){
    if(count <=0){
      return(mat)
    }
    return(n_down(double_down((mat/max(mat))),count-1))
  }
  times=40
  counts <- n_down(counts,times)
  plot_ly(z = counts, x=interval_alphas,y=interval_means,type = "heatmap")

  best = which(counts == max(counts), arr.ind = TRUE)
  best_rad <- interval_means[best[1]+times]
  best_alpha <- interval_alphas[best[2]+times]

  best_rad
  best_alpha

  ggplot() + geom_segment(data=pts,aes(x=tcnlow,xend=tcnhigh,y=(radeff1),yend=(radeff1),col=factor(angle1>0.5 | angle1< -0.5)),alpha=0.1,size=5) + geom_segment(data=pts,aes(x=tcnlow,xend=tcnhigh,y=(radeff2),yend=(radeff2),col=factor(angle1>0.5 | angle1< -0.5)),alpha=0.1,size=5)+ylim(0,5)+xlim(0,5)+theme(legend.position="none")+geom_abline(intercept=best_rad,slope=tan(best_alpha))

  pts$information1 <- abs(pts$radeff1-pts$center*tan(best_alpha)-best_rad)#*exp(-pts$jump^2))
  pts$information2 <- abs(pts$radeff2-pts$center*tan(best_alpha)-best_rad)#*exp(-pts$jump^2))
  #plot(ecdf(log(c(pts$information1,pts$information2))))
  #values <- sort(log(c(pts$information1,pts$information2)))
  #sapply(values,function(x) sum())
  #median(count(as.integer(x))$freq)
  #hs <- sort(E(x))
  #gradient <- (hs[1:(length(hs)-1)]- hs[2:(length(hs))])/(x[1:(length(x)-1)]- x[2:(length(x))])
  #gradient <- gradient[which(!is.nan(gradient))]
  #gradient <- SMA(gradient,10)

  #ggplot()+geom_point(aes(x=1:length(gradient),y=gradient))
  ggplot() + geom_segment(data=pts,aes(x=tcnlow,xend=tcnhigh,y=(radeff1),yend=(radeff1),col=information1),alpha=0.1,size=5) + geom_segment(data=pts,aes(x=tcnlow,xend=tcnhigh,y=(radeff2),yend=(radeff2),col=information2),alpha=0.1,size=5)+ylim(0,5)+xlim(0,5)+theme(legend.position="none")+geom_abline(intercept=best_rad,slope=tan(best_alpha))+ scale_colour_gradient(low="red")

  pts$inf <- log(apply(pts[,c("information1","information2")],1,max))

  pts <- pts[order(-pts$inf),]
  ggplot(pts)+geom_segment(aes(x=x1,y=y1,xend=x2,yend=y2,col=inf,size=sizelog))+coord_fixed()+xlab("CN 1")+ylab("CN 2")+ggtitle("Grid plot")+ scale_colour_gradient(low="red")



  # ok so for each segment there were two angles that could possibly be true, and in truth only at most one of these is true.
  # so remove the furthest angle from the data and continue
  # this shouldnt change anything because it can't make a second largest cluster larger than the first.
  pts2 <- pts[,!names(pts) %in% c("radeff2","information2","angle2")]
  pts3 <- pts[,!names(pts) %in% c("radeff1","information1","angle1")]
  colnames(pts3) <- colnames(pts2)

  # ok continue with this idea - i think we are on to something!

  removefrompts2 <- removefrompts3 <- rep(FALSE,nrow(pts2))
  for( i in 1:nrow(pts2)){
    if(!complete.cases(pts3[i,]) | !complete.cases(pts2[i,])){
      if(!complete.cases(pts3[i,])){
        removefrompts3[i] = TRUE
      }
      if(!complete.cases(pts2[i,])){
        removefrompts2[i] = TRUE
      }
      next
    }
    if(pts2[i,"radeff1"] == pts3[i,"radeff1"]){
      removefrompts3[i] = TRUE
    } else if(pts2[i,"information1"] < pts3[i,"information1"]){
      removefrompts3[i] = TRUE
    } else {
      removefrompts2[i] = TRUE
    }
    # if(is.na(pts2[i, "angle1"]) | pts2[i,"angle1"] < -0.5 | pts2[i,"angle1"] > 0.5){
    #   removefrompts2[i] = TRUE
    # }
    # if(is.na(pts3[i, "angle1"]) | pts3[i,"angle1"] < -0.5 | pts3[i,"angle1"] > 0.5){
    #   removefrompts2[i] = TRUE
    # }
  }
  pts <- rbind(pts2[which(removefrompts2),],pts3[which(removefrompts3),])
  #pts <- pts2

  #ggplot()+geom_point(aes(x=1:nrow(counts),y=counts[,which(sums$score == max(sums$score))]))
  pts$T1 <- pts$x1+pts$y1
  pts$T2 <- pts$x2+pts$y2
  pts$TCNend = apply(pts[,c("T1","T2")],1,max)
  pts$TCN = apply(pts[,c("T1","T2")],1,min)
  #pts$reff = ((1+(pts$rad1-1))^(1/(pts$TCNend-pts$TCN))-1)
  #pts <-pts[-which(pts$information1<1),]
  #stopifnot(pts$reff > 0)
  pts <- pts[complete.cases(pts),]
  ys = sort(c(pts$tcnlow,pts$tcnhigh))

  ggplot()+geom_point(aes(x=1:length(ys),y=ys))+ylab("TCN")+xlab("independent segments")
  intervals <- data.frame("TCN" = unique(round(ys,2)))
  intervals$r <- intervals$TCN*0
  intervals$weights <- intervals$TCN*0
  intervals$TCNend <- intervals$TCN*0

  for(i in 1:(nrow(intervals)-1)){
    #rvalues = c()
    reffvalues = c()
    weights = c()
    number = 0
    for(j in 1:nrow(pts)){
      if((pts[j,"TCN"] <= intervals[i+1,"TCN"] &
          pts[j,"TCNend"] >= intervals[i,"TCN"] &
          pts[j,"TCNend"]-pts[j,"TCN"] > 0.2 &
          abs(pts[j,"rad1"])<5)){
        #rvalues[number+1] = ((pts[j,"radeff1"]+1-pts[j,"tcnhigh"] +pts[j,"tcnlow"])^(intervals[i+1,"TCN"]-intervals[i,"TCN"])-1)
        #1/(pts$tcnhigh - pts$tcnlow) * log(1+pts$rad1 - pts$tcnhigh + pts$tcnlow)
        reffvalues[number+1] = pts[j,"radeff1"]
        weights[number+1] = sum(pts[j,"sizelog"] *exp(-(60:200/pts2[,"jump"]^2*pts[j,"information1"])^2))
        number = number + 1
      }
    }
    #intervals[i,"r"] = sum(rvalues*weights)/sum(weights)
    intervals[i,"reff"] = sum(reffvalues*weights)/sum(weights)
    intervals[i,"weights"] = sum(weights)
    intervals[i,"TCNend"] = intervals[i+1,"TCN"]
    intervals[i,"probes"] = sum(mattab[which(mattab$t < intervals[i+1,"TCN"] & mattab$t > intervals[i,"TCN"]),"length"])
  }
  #intervals$r
  intervals <- intervals[complete.cases(intervals),]
  intervals <- intervals[which(intervals$TCNend >0),]
  intervals$reffbeta <- ((1+intervals$r)^(1/(intervals$TCNend-intervals$TCN))-1)
  ggplot(data=intervals) + geom_segment(aes(x=TCN,xend=TCNend, y=reff,yend=reff,size=probes))

  #reffbeta is clearly better.

  tempintervals <- intervals
  #intervals <- tempintervals

  for(i in 1:(nrow(intervals))){
    # if(sum(abs(intervals[i,"angle"] - intervals[,"angle"]) <0.075)<nrow(intervals)/3){
    #   intervals[i,"weights"] = 0
    # }
    #
    # if(sum(abs(intervals[i,"r"] - intervals[,"r"]) <1)<nrow(intervals)/3){
    #   intervals[i,"weights"] = 0
    # }
    if(sum(abs(intervals[i,"reff"] - intervals[,"reff"]) <0.5)<nrow(intervals)/3){
      intervals[i,"weights"] = 0
    }
  }

  newintervals <- intervals[which(intervals$weights>0),]

  # want to break up newintervals into smaller pieces and then apply some smoothing!
  max_interval_width = 0.01
  count = nrow(newintervals)
  initsize =  nrow(newintervals)
  for(i in 1:(nrow(newintervals)-1)){
    diff = newintervals[i+1,"TCN"]- newintervals[i,"TCN"]

    for(j in 1:as.integer(diff/max_interval_width)){
      newintervals[count+1,]= newintervals[i,]
      newintervals[count+1,"TCN"] = newintervals[i,"TCN"]+j*max_interval_width
      newintervals[count+1,"TCNend"] = min(newintervals[i+1,"TCN"],newintervals[i,"TCN"]+(j+1)*max_interval_width)
      count = count + 1
    }
    newintervals[i,"TCNend"] = newintervals[i,"TCN"]+max_interval_width
  }
  row.names(newintervals) <- 1:nrow(newintervals)

  # break up the last one on its own:
  diff = newintervals[initsize,"TCNend"]- newintervals[initsize,"TCN"]
  for(j in 1:as.integer(diff/max_interval_width)){
    newintervals[count+1,]= newintervals[initsize,]
    newintervals[count+1,"TCN"] = newintervals[initsize,"TCN"]+j*max_interval_width
    newintervals[count+1,"TCNend"] = min(newintervals[initsize,"TCNend"],newintervals[initsize,"TCN"]+(j+1)*max_interval_width)
    count = count + 1
  }
  newintervals[initsize,"TCNend"] = newintervals[initsize,"TCN"]+max_interval_width
  # finished breaking up the intervals!
  # now sort them:
  newintervals <- newintervals[order(newintervals$TCN),]

  goodintervals = newintervals
  goodintervals$reffEMA <- EMA(goodintervals$reff,n=as.integer(0.05*nrow(goodintervals)))
  goodintervals[which(is.na(goodintervals$reffEMA)),"reffEMA"] <- goodintervals[which(is.na(goodintervals$reffEMA)),"reff"]
  goodintervals$reffSMA <- SMA(goodintervals$reff,n=as.integer(0.05*nrow(goodintervals)))
  goodintervals[which(is.na(goodintervals$reffSMA)),"reffSMA"] <- goodintervals[which(is.na(goodintervals$reffSMA)),"reff"]

  g1 <- ggplot(data=goodintervals[which(goodintervals$weights>0),]) + geom_segment(aes(x=TCN,xend=TCNend, y=reff,yend=reff,size=weights,col=log(probes)))+ scale_colour_gradient(low="red")
  g3 <- ggplot(data=goodintervals[which(goodintervals$weights>0),]) + geom_segment(aes(x=TCN,xend=TCNend, y=reffEMA,yend=reffEMA,size=weights,col=log(probes)))+ scale_colour_gradient(low="red")
  g5 <- ggplot(data=goodintervals[which(goodintervals$weights>0),]) + geom_segment(aes(x=TCN,xend=TCNend, y=reffSMA,yend=reffSMA,size=weights,col=log(probes)))+ scale_colour_gradient(low="red")
  grid.arrange(g1,g3,g5)

  #goodintervals <- goodintervals[which(goodintervals$weights > 1),]
  # extend beyond the start and to the end.
  low_mean <- apply(goodintervals[1:(nrow(goodintervals)/4),],2,mean)
  low_mean["TCN"] <- 0
  low_mean["TCNend"] <- goodintervals[1,"TCN"]

  high_mean <-  apply(goodintervals[(nrow(goodintervals)/4):(nrow(goodintervals)),],2,mean)
  high_mean["TCN"] <- goodintervals[nrow(goodintervals),"TCN"]#,error= function(x){message(x);stop(filename)}
  high_mean["TCNend"] <- 10

  # want to break up newintervals into smaller pieces and then apply some smoothing!
  max_interval_width = 0.01
  diff = high_mean["TCNend"]- high_mean["TCN"]
  count = 1
  for(j in 1:as.integer(diff/max_interval_width)){
    if(j==1){
      high_mean = rbind(high_mean,high_mean)
    } else{
      high_mean = rbind(high_mean,high_mean[1,])
    }
    high_mean[nrow(high_mean),"TCN"] = high_mean[1,"TCN"]+j*max_interval_width
    high_mean[nrow(high_mean),"TCNend"] = min(high_mean[1,"TCNend"],high_mean[1,"TCN"]+(j+1)*max_interval_width)
  }
  high_mean[1,"TCNend"] = high_mean[1,"TCN"]+max_interval_width

  goodintervals$fake = FALSE
  low_mean = as.data.frame(t(low_mean))
  low_mean["fake"] = TRUE
  high_mean = as.data.frame(high_mean)
  high_mean$fake = TRUE
  goodintervals <- rbind(data.frame(low_mean),data.frame(goodintervals),data.frame(high_mean))
  row.names(goodintervals) <- 1:nrow(goodintervals)

  #g2 <- ggplot(data=goodintervals[which(goodintervals$weights>0),]) + geom_segment(aes(x=TCN,xend=TCNend, y=reffEMA,yend=reffEMA,size=weights,col=weights>1))+ylim(0,6)
  #grid.arrange(g2,nrow=5)

  # g1 <- ggplot(data=goodintervals[which(!goodintervals$fake),]) + geom_segment(aes(x=TCN,xend=TCNend, y=tan(angle),yend=tan(angle),size=weights,col=weights>1))
  # g2 <- ggplot(data=goodintervals[which(!goodintervals$fake),]) + geom_segment(aes(x=TCN,xend=TCNend, y=r,yend=r,size=weights,col=weights>1))
  # g3 <- ggplot(data=goodintervals[which(!goodintervals$fake),]) + geom_segment(aes(x=TCN,xend=TCNend, y=reff,yend=reff,size=weights,col=weights>1))
  # grid.arrange(g1,g2,g3)
  #((pts[j,"radeff1"]+1-pts[j,"tcnhigh"] +pts[j,"tcnlow"])^(intervals[i+1,"TCN"]-intervals[i,"TCN"])-1)


  # ((pts[j,"radeff1"]+1-pts[j,"tcnhigh"] +pts[j,"tcnlow"])^(intervals[i+1,"TCN"]-intervals[i,"TCN"])-1)
  #
  # pts$radeff1 <- 1/(pts$tcnhigh - pts$tcnlow) * log(1+pts$rad1 - pts$tcnhigh + pts$tcnlow)
  #
  # goodintervals$r <- (goodintervals$reff+1-goodintervals$TCNend + goodintervals$TCN)^(goodintervals$TCNend - goodintervals$TCN) -1
  # goodintervals$rbeta <- (goodintervals$reffbeta+1-goodintervals$TCNend + goodintervals$TCN)^(goodintervals$TCNend - goodintervals$TCN) -1
  # goodintervals$rEMA <- (goodintervals$reffEMA+1-goodintervals$TCNend + goodintervals$TCN)^(goodintervals$TCNend - goodintervals$TCN) -1
  # goodintervals$rbetaEMA <- (goodintervals$reffbetaEMA+1-goodintervals$TCNend + goodintervals$TCN)^(goodintervals$TCNend - goodintervals$TCN) -1
  # goodintervals$rSMA <- (goodintervals$reffSMA+1-goodintervals$TCNend + goodintervals$TCN)^(goodintervals$TCNend - goodintervals$TCN) -1
  # goodintervals$rbetaSMA <- (goodintervals$reffbetaSMA+1-goodintervals$TCNend + goodintervals$TCN)^(goodintervals$TCNend - goodintervals$TCN) -1

  # this part is what is wrong the compounding shouldnt be over these intervals...... yes it should.

  goodintervals$r <- exp(goodintervals$reff*(goodintervals$TCNend -goodintervals$TCN))
  goodintervals$rEMA <- exp(goodintervals$reffEMA*(goodintervals$TCNend -goodintervals$TCN))
  goodintervals$rSMA <- exp(goodintervals$reffSMA*(goodintervals$TCNend -goodintervals$TCN))

  #### UP TO HERE ####
  goodintervals$CNr = goodintervals$CNrEMA  = goodintervals$CNrSMA  = rep(1,nrow(goodintervals))
  for(i in 2:nrow(goodintervals)){
    goodintervals[i,"CNr"] = goodintervals[i-1,"CNr"]*(goodintervals[i,"r"])
    goodintervals[i,"CNrEMA"] = goodintervals[i-1,"CNrEMA"]*(goodintervals[i,"rEMA"])
    goodintervals[i,"CNrSMA"] = goodintervals[i-1,"CNrSMA"]*(goodintervals[i,"rSMA"])
  }
  g1<-ggplot(data=goodintervals) + geom_segment(aes(x=TCN,xend=TCNend, y=CNr,yend=CNr,size=weights,col=log(probes)))+xlim(0,4)+ylim(0,10)+ scale_colour_gradient(low="red")
  g2<-ggplot(data=goodintervals) + geom_segment(aes(x=TCN,xend=TCNend, y=CNrEMA,yend=CNrEMA,size=weights,col=log(probes)))+xlim(0,4)+ylim(0,10)+ scale_colour_gradient(low="red")
  g3<-ggplot(data=goodintervals) + geom_segment(aes(x=TCN,xend=TCNend, y=CNrSMA,yend=CNrSMA,size=weights,col=log(probes)))+xlim(0,4)+ylim(0,10)+ scale_colour_gradient(low="red")
  grid.arrange(g1,g2,g3)
  # work out which of these relationships to use.

  mattn$t <- mattn$major + mattn$minor
  mattn$newtr  <-mattn$newtrEMA  <-mattn$newtrSMA <- mattn$t*0
  for(i in 1:nrow(mattn)){
    CT <- newTr <- newTrEMA <- newTrSMA <- mattn[i,"t"]
    for(j in 1:(nrow(goodintervals)-1)){
      if( goodintervals[j,"TCN"] < CT & goodintervals[j,"TCNend"]>=CT){
        if(j == 1){
          #linear interpolation to zero.
          newTr <- newTrEMA <- newTrSMA <- goodintervals[1,"TCNend"]/CT
          break
        }
        alpha = (CT-goodintervals[j,"TCN"])/(goodintervals[j,"TCNend"]-goodintervals[j,"TCN"])
        newTr <- goodintervals[j,"CNr"]*alpha+goodintervals[j+1,"CNr"]*(1-alpha)
        newTrEMA <- goodintervals[j,"CNrEMA"]*alpha+goodintervals[j+1,"CNrEMA"]*(1-alpha)
        newTrSMA <- goodintervals[j,"CNrSMA"]*alpha+goodintervals[j+1,"CNrSMA"]*(1-alpha)
        break
      }
    }
    mattn[i,"newtr"] <-newTr
    mattn[i,"newtrEMA"] <-newTrEMA
    mattn[i,"newtrSMA"] <-newTrSMA
  }

  matte <- mattn
  matte$major <- matte$major*matte$newtrEMA/matte$t
  matte$minor <- matte$minor*matte$newtrEMA/matte$t
  g3 <- get_grid_plot(phase(rescale(matte,rescale_to_one = FALSE)),200)+ theme(legend.position="none")+xlim(0,3)+ylim(0,3)+ggtitle("rEMA")

  matte <- mattn
  matte$major <- matte$major*matte$newtrSMA/matte$t
  matte$minor <- matte$minor*matte$newtrSMA/matte$t
  g5 <- get_grid_plot(phase(rescale(matte,rescale_to_one = FALSE)),200)+ theme(legend.position="none")+xlim(0,3)+ylim(0,3)+ggtitle("rSMA")

  matte <- mattn
  matte$major <- matte$major*matte$newtr/matte$t
  matte$minor <- matte$minor*matte$newtr/matte$t
  #ggplot(matte) + geom_point(aes(x=jitter((matte$major+matte$minor),factor=50),y=jitter(round(major+minor),factor=1.05)))+ geom_abline(intercept = 0, slope = 1)+xlim(0,5)+ylim(0,5)
  #m1 <- lm(matte$major+matte$minor ~ round(matte$major+matte$minor))
  g1 <- get_grid_plot(phase(rescale(matte,rescale_to_one = FALSE)),200)+ theme(legend.position="none")+xlim(0,3)+ylim(0,3)+ggtitle("r")
  get_CN_track(phase(rescale(matte,rescale_to_one = FALSE)))

  grid.arrange(g1,g3,g5,ncol=3)
  # now we have fitted the gradients - can we fit to have the most points on integers?
  # how would you even do that? you would start by expecting that you have done a good job already with the gradients.
  # actually that should just be added to the clustering algorithm.

  return(list("mat" = matte,"pts" = pts,"counts" = counts,"intervals" = intervals,"pts2" = pts2))
}




















# remove_bias_grid <- function(){
#   # begin our algorithm:
#   final = FALSE
#   prev_range = 0
#   for (k in 0:11){
#     print(paste("Iteration ",as.character(k+1)))
#     if (final == TRUE){
#       k = k-2
#     }
#     radius_thresh <- 0.015*(10-k) # this describes the precision that we expect - the scale that we are observeign the grid plot from. We go through a range of values -  from having blurry eyes to sucessively sharper eyes when looking at the grid plot.
#     dat<- mat2[which(mat2$length>0),]
#     main_pts <- get_main_points3(dat,"major","minor",radius_thresh/3,0.005)
#     main_pts_ones <- main_pts[which(abs(main_pts$major - main_pts$minor) < 0.5),]
#     ones <- main_pts_ones[which(main_pts_ones$freq == max(main_pts_ones$freq)),]
#     main_pts <- main_pts[order(main_pts$TCN),]
#
#     # now create a data structure to analyse the main paths:
#     main_segs <-dat[which(dat$length>10),c("chr","segment","length","major","minor","start","end")]
#     rownames(main_segs)<- NULL
#     main_paths <- cbind(main_segs[1:(nrow(main_segs)-1),],main_segs[2:(nrow(main_segs)),])
#     colnames(main_paths)<- c("Chromosome1","segment1","length1","major1","minor1","start1","end1","Chromosome2","segment2","length2","major2","minor2","start2","end2")
#     # weight the paths such that paths between well established points are given better values...
#     main_paths["weights"] <- apply(main_paths[,c("length1","length2")],1,function(x) return(1/(1/x[1]^0.5+1/x[2]^0.5)))
#
#     # take a confident subset of the main points:
#     great_pts <- main_pts[which(main_pts$length>100),]
#     great_pts$name <- 1:nrow(great_pts)
#     great_pts$N <- rep(FALSE,nrow(great_pts))
#     great_pts$E <- rep(FALSE,nrow(great_pts))
#     great_pts$S <- rep(FALSE,nrow(great_pts))
#     great_pts$W <- rep(FALSE,nrow(great_pts))
#     great_pts$from <- rep(0,nrow(great_pts))
#     great_pts$to <- rep(0,nrow(great_pts))
#     great_pts$TCN_true = great_pts$TCN*0
#
#     pair <- get_edges(main_paths,great_pts,radius_thresh)
#     relationships <- pair$relationships
#     relationship_count <- pair$relationship_count
#     for(i in 1:nrow(relationships)){
#       for(j in 1:nrow(relationships)){
#         if(abs(great_pts[i,"major"] - great_pts[j,"major"]) + abs(great_pts[i,"major"] - great_pts[j,"major"]) <0.2){
#           relationships[i,j] = relationships[j,i] = 0
#           relationship_count[i,j] = relationship_count[j,i] = 0
#         }
#       }
#     }
#
#     # relationships describe how much evidence there is for these two approximate points to be either in a N,S,E,W direction by viewing how much evidence there is for paths between these points in our genome...
#
#     relationships <- relationships+t(relationships)
#     relationship_count <- relationship_count + t(relationship_count)
#     relationships[which(lower.tri(relationships,diag=TRUE))] <- 0
#     relationship_count[which(lower.tri(relationships,diag=TRUE))] <- 0
#     relationships <- relationships+t(relationships)
#     relationship_count <- relationship_count + t(relationship_count)
#
#     # now we have found the main paths and the main points we automatically remove the inferred bias between these points.
#     great_pts_B <- great_pts[which(abs(great_pts$major - great_pts$minor)<0.5),]
#     start_B = which(great_pts_B$freq == max(great_pts_B$freq))
#     start = which(great_pts_B[start_B,"major"] == great_pts$major)
#     curr_node = start
#     great_pts[start,"TCN_true"]<- 2
#     visited=c(start)
#     pairs = c()
#     repeat{
#       next_node <- curr_node
#       max_path <- 0
#       for (node in visited){
#         #print(node)
#         triple = get_next_node(node,relationships,great_pts,visited,radius_thresh)
#         pot_curr_node <- triple[1]
#         pot_next_node <- triple[2]
#         weight <- triple[3]
#         #print(triple)
#         if (weight>max_path){
#           next_node <- pot_next_node
#           curr_node <- pot_curr_node
#           max_path = weight
#         }
#       }
#       if(max_path==0){
#         break
#       }
#       next_node <- as.numeric(next_node)
#       #print(curr_node)
#       #print(next_node)
#       r <- rAB(great_pts[curr_node,"major"],great_pts[curr_node,"minor"],great_pts[next_node,"major"],great_pts[next_node,"minor"])
#       if (abs(great_pts[curr_node,"major"]-great_pts[next_node,"major"]) < abs(great_pts[curr_node,"minor"]-great_pts[next_node,"minor"])){
#         if ((great_pts[curr_node,"minor"]-great_pts[next_node,"minor"]) > 0){
#           great_pts[curr_node,"S"] = TRUE
#           great_pts[next_node,"N"] = TRUE
#         }
#         else {
#           great_pts[curr_node,"N"] = TRUE
#           great_pts[next_node,"S"] = TRUE
#         }
#       }
#       if (abs(great_pts[curr_node,"major"]-great_pts[next_node,"major"]) > abs(great_pts[curr_node,"minor"]-great_pts[next_node,"minor"])){
#         if ((great_pts[curr_node,"major"]-great_pts[next_node,"major"]) > 0){
#           great_pts[curr_node,"W"]  = TRUE
#           great_pts[next_node,"E"]  = TRUE
#         }
#         else {
#           great_pts[curr_node,"E"]  = TRUE
#           great_pts[next_node,"W"]  = TRUE
#         }
#       }
#       great_pts[next_node,"from"] <- curr_node
#       great_pts[curr_node,"to"] <- next_node
#       curr_obs_TCN <- great_pts[curr_node,"TCN"]
#       next_obs_TCN <- great_pts[next_node,"TCN"]
#       curr_true_TCN <- great_pts[curr_node,"TCN_true"]
#       great_pts[next_node,"TCN_true"]<- r*curr_true_TCN
#       pairs <- rbind(pairs,cbind(curr_node,next_node))
#       (curr_node <- next_node)
#       visited <- cbind(visited,as.numeric(next_node))
#     }
#
#     graph_adj = matrix(rep(0,nrow(great_pts)^2),ncol=nrow(great_pts))
#
#     for (i in 1:nrow(great_pts)){
#       if (great_pts[i,"from"]!=0){
#         graph_adj[min(great_pts[i,"from"],i),max(great_pts[i,"from"],i)] = 1
#       }
#       if (great_pts[i,"to"]!=0){
#         graph_adj[min(great_pts[i,"from"],i),max(great_pts[i,"from"],i)] = 1
#       }
#     }
#
#     g <- graph.adjacency(graph_adj, mode="directed")
#     # begin our algorithm:
#     final = FALSE
#     prev_range = 0
#     for (k in 1:11){
#       print(paste("Iteration ",as.character(k+1)))
#       if (final == TRUE){
#         k = k-2
#       }
#       radius_thresh <- 0.02*(10-k) # this describes the precision that we expect - the scale that we are observeign the grid plot from. We go through a range of values -  from having blurry eyes to sucessively sharper eyes when looking at the grid plot.
#       dat<- mat2[which(mat2$length>0),]
#       main_pts <- get_main_points3(dat,"major","minor",radius_thresh,0.001)
#       main_pts_ones <- main_pts[which(abs(main_pts$major - main_pts$minor) < 0.5),]
#       ones <- main_pts_ones[which(main_pts_ones$freq == max(main_pts_ones$freq)),]
#       main_pts <- main_pts[order(main_pts$TCN),]
#
#       # now create a data structure to analyse the main paths:
#       main_segs <-dat[,c("chr","segment","length","major","minor","start","end")]
#       rownames(main_segs)<- NULL
#       main_paths <- cbind(main_segs[1:(nrow(main_segs)-1),],main_segs[2:(nrow(main_segs)),])
#       colnames(main_paths)<- c("Chromosome1","segment1","length1","major1","minor1","start1","end1","Chromosome2","segment2","length2","major2","minor2","start2","end2")
#       # weight the paths such that paths between well established points are given better values...
#       main_paths["weights"] <- apply(main_paths[,c("length1","length2")],1,function(x) return(1/(1/x[1]^0.5+1/x[2]^0.5)))
#
#       # take a confident subset of the main points:
#       great_pts <- main_pts[which(main_pts$length>100),]
#       great_pts$name <- 1:nrow(great_pts)
#       great_pts$N <- rep(FALSE,nrow(great_pts))
#       great_pts$E <- rep(FALSE,nrow(great_pts))
#       great_pts$S <- rep(FALSE,nrow(great_pts))
#       great_pts$W <- rep(FALSE,nrow(great_pts))
#       great_pts$from <- rep(0,nrow(great_pts))
#       great_pts$to <- rep(0,nrow(great_pts))
#       great_pts$TCN_true = great_pts$TCN*0
#
#       pair <- get_edges(main_paths,great_pts,radius_thresh)
#       relationships <- pair$relationships
#       relationship_count <- pair$relationship_count
#
#       # relationships describe how much evidence there is for these two approximate points to be either in a N,S,E,W direction by viewing how much evidence there is for paths between these points in our genome...
#
#       relationships <- relationships+t(relationships)
#       relationship_count <- relationship_count + t(relationship_count)
#       relationships[which(lower.tri(relationships,diag=TRUE))] <- 0
#       relationship_count[which(lower.tri(relationships,diag=TRUE))] <- 0
#       relationships <- relationships+t(relationships)
#       relationship_count <- relationship_count + t(relationship_count)
#
#       # now we have found the main paths and the main points we automatically remove the inferred bias between these points.
#       great_pts_B <- great_pts[which(abs(great_pts$major - great_pts$minor)<0.5),]
#       start_B = which(great_pts_B$freq == max(great_pts_B$freq))
#       start = which(great_pts_B[start_B,"major"] == great_pts$major)
#       curr_node = start
#       great_pts[start,"TCN_true"]<- 2
#       visited=c(start)
#       pairs = c()
#       repeat{
#         next_node <- curr_node
#         max_path <- 0
#         for (node in visited){
#           #print(node)
#           triple = get_next_node(node,relationships,great_pts,visited,radius_thresh)
#           pot_curr_node <- triple[1]
#           pot_next_node <- triple[2]
#           weight <- triple[3]
#           #print(triple)
#           if (weight>max_path){
#             next_node <- pot_next_node
#             curr_node <- pot_curr_node
#             max_path = weight
#           }
#         }
#         if(max_path==0){
#           break
#         }
#         next_node <- as.numeric(next_node)
#         #print(curr_node)
#         #print(next_node)
#         r <- rAB(great_pts[curr_node,"major"],great_pts[curr_node,"minor"],great_pts[next_node,"major"],great_pts[next_node,"minor"])
#         if (abs(great_pts[curr_node,"major"]-great_pts[next_node,"major"]) < abs(great_pts[curr_node,"minor"]-great_pts[next_node,"minor"])){
#           if ((great_pts[curr_node,"minor"]-great_pts[next_node,"minor"]) > 0){
#             great_pts[curr_node,"S"] = TRUE
#             great_pts[next_node,"N"] = TRUE
#           }
#           else {
#             great_pts[curr_node,"N"] = TRUE
#             great_pts[next_node,"S"] = TRUE
#           }
#         }
#         if (abs(great_pts[curr_node,"major"]-great_pts[next_node,"major"]) > abs(great_pts[curr_node,"minor"]-great_pts[next_node,"minor"])){
#           if ((great_pts[curr_node,"major"]-great_pts[next_node,"major"]) > 0){
#             great_pts[curr_node,"W"]  = TRUE
#             great_pts[next_node,"E"]  = TRUE
#           }
#           else {
#             great_pts[curr_node,"E"]  = TRUE
#             great_pts[next_node,"W"]  = TRUE
#           }
#         }
#         great_pts[next_node,"from"] <- curr_node
#         great_pts[curr_node,"to"] <- next_node
#         curr_obs_TCN <- great_pts[curr_node,"TCN"]
#         next_obs_TCN <- great_pts[next_node,"TCN"]
#         curr_true_TCN <- great_pts[curr_node,"TCN_true"]
#         great_pts[next_node,"TCN_true"]<- r*curr_true_TCN
#         pairs <- rbind(pairs,cbind(curr_node,next_node))
#         (curr_node <- next_node)
#         visited <- cbind(visited,as.numeric(next_node))
#       }
#       graph_adj = matrix(rep(0,nrow(great_pts)^2),ncol=nrow(great_pts))
#
#       for (i in 1:nrow(great_pts)){
#         if (great_pts[i,"from"]!=0){
#           graph_adj[min(great_pts[i,"from"],i),max(great_pts[i,"from"],i)] = 1
#         }
#         if (great_pts[i,"to"]!=0){
#           graph_adj[min(great_pts[i,"from"],i),max(great_pts[i,"from"],i)] = 1
#         }
#       }
#
#       g <- graph.adjacency(graph_adj, mode="directed")
#
#
#       #  dat3 <- remove_bias(mat3,mat3,3,60,3)
#       # dat <- rescale(phase(dat3$mat),rescale_to_one=FALSE)
#       #get_grid_plot(dat,100)
#       #get_grid_plot(mat,100)
#       #get_CN_track(dat)
#       # iteratively apply remove_bias, stopping at the round which has the smallest summed square distance to the nearest integer.
#       # METHOD 3, my own dynamic programming approach:
#       print("The inferred tree for which edges  represent horizontal and verticle steps in truth:")
#       plot(g)
#       path <- get_path_pts(graph_adj)
#       path_pts <- great_pts[path,]
#       table <- NULL
#       table <- great_pts[path,c("TCN","TCN_true")]
#
#       g <- ggplot(dat[which(dat$length > MIN_LENGTH),],aes(major,minor))+geom_path(aes(col=chr))+geom_point(data=great_pts,aes(major,minor,label=name))+geom_path(data=path_pts,aes(major,minor),size=3,col="#FFFFFF")+geom_text(data=great_pts,aes(major,minor,label=name),size=6)+coord_fixed()
#       print("The path taken through the main points found is:")
#       print(path)
#       print("The main points on our grid plot (the path from low TCN to high TCN  composed of (what is in truth) entirely verticle and horizontal steps is in white):")
#       print(g)
#       print("The weighted adjacency matrix of paths between the main points in our grid plot is:")
#       print(signif(relationships,digits=3))
#       print("The adjacency matrix of counts of independent observations of paths between the main points in our grid plot is:")
#       print(signif(relationship_count,digits=3))
#       print("Information on the main points:")
#       print(great_pts[,c("name","major","minor","TCN","TCN_true","freq")])
#       range <- max(table[,"TCN_true"])- min(table[which(table$TCN_true > 0),"TCN_true"])
#       if (final == TRUE){
#         break
#       }
#       if (range > 0.9*prev_range){
#         prev_range = range
#       } else {
#         final = TRUE
#       }
#       if (nrow(table)<2){
#         next
#       }
#       for (i in 2:nrow(table)){
#         if (table[i-1,"TCN_true"]>table[i,"TCN_true"]){
#           final = TRUE
#         }
#       }
#     }
#   }
# }


mardiaTest <- function (data, cov = TRUE, qqplot = FALSE){
  if (dim(data)[2] < 2 || is.null(dim(data))) {stop("number of variables must be equal or greater than 2")}
  if (dim(data)[1] < 2 || is.null(dim(data))) {stop("number of observations must be equal or greater than 2")}

  if (!is.data.frame(data) && !is.matrix(data)) stop('Input must be one of classes \"data frame\" or \"matrix\"')

  dataframe=as.data.frame(data)
  dname <- deparse(substitute(data))
  data <- as.matrix(data)
  n <- dim(data)[1]
  p <- dim(data)[2]
  data.org <- data
  data <- scale(data, scale = FALSE)
  if (cov) {
    S <- ((n - 1)/n) * cov(data)
  }
  else {
    S <- cov(data)
  }
  D <- data %*% solve(S) %*% t(data)
  g1p <- sum(D^3)/n^2
  g2p <- sum(diag((D^2)))/n
  df <- p * (p + 1) * (p + 2)/6
  k <- (p + 1) * (n + 1) * (n + 3)/(n * ((n + 1) * (p + 1) -
                                           6))
  small.skew <- n * k * g1p/6
  skew <- n * g1p/6
  kurt <- (g2p - p * (p + 2)) * sqrt(n/(8 * p * (p + 2)))
  p.skew <- pchisq(skew, df, lower.tail = FALSE)
  p.small <- pchisq(small.skew, df, lower.tail = FALSE)
  p.kurt <- 2 * (1 - pnorm(abs(kurt)))
  #     if (qqplot) {
  #       d <- diag(D)
  #       r <- rank(d)
  #       chi2q <- qchisq((r - 0.5)/n, p)
  #       plot(d, chi2q, pch = 19, main = "Chi-Square Q-Q Plot",
  #            xlab = "Squared Mahalanobis Distance", ylab = "Chi-Square Quantile")
  #       abline(0, 1, lwd = 2, col = "black")
  #     }
  result <- list("mardia", g1p = g1p, chi.skew = skew, p.value.skew = p.skew,
                 chi.small.skew = small.skew, p.value.small = p.small, g2p = g2p,
                 z.kurtosis = kurt, p.value.kurt = p.kurt, dname = dname, dataframe = dataframe)
  return(result)
}


hzTest <-function (data, cov = TRUE, qqplot = FALSE){
  if (dim(data)[2] < 2 || is.null(dim(data))) {stop("number of variables must be equal or greater than 2")}
  if (!is.data.frame(data) && !is.matrix(data)) stop('Input must be one of classes \"data frame\" or \"matrix\"')

  dataframe=as.data.frame(data)
  dname <- deparse(substitute(data))
  data <- as.matrix(data)
  n <- dim(data)[1]
  p <- dim(data)[2]
  data.org = data

  if (cov){
    S <- ((n-1)/n)*cov(data)
  }
  else    {
    S <- cov(data)
  }

  dif <- scale(data, scale = FALSE)


  Dj <- diag(dif%*%solve(S)%*%t(dif))  #squared-Mahalanobis' distances

  Y <- data%*%solve(S)%*%t(data)


  Djk <- - 2*t(Y) + matrix(diag(t(Y)))%*%matrix(c(rep(1,n)),1,n) + matrix(c(rep(1,n)),n,1)%*%diag(t(Y))

  b <- 1/(sqrt(2))*((2*p + 1)/4)^(1/(p + 4))*(n^(1/(p + 4))) #smoothing
  {                                                                 #parameter
    if (qr(S)$rank == p){
      HZ = n * (1/(n^2) * sum(sum(exp( - (b^2)/2 * Djk))) - 2 *
                  ((1 + (b^2))^( - p/2)) * (1/n) * (sum(exp( - ((b^2)/(2 *
                                                                         (1 + (b^2)))) * Dj))) + ((1 + (2 * (b^2)))^( - p/2)))
    }
    else {
      HZ = n*4
    }

  }
  wb <- (1 + b^2)*(1 + 3*b^2)

  a <- 1 + 2*b^2

  mu <- 1 - a^(- p/2)*(1 + p*b^2/a + (p*(p + 2)*(b^4))/(2*a^2)) #HZ mean

  si2 <- 2*(1 + 4*b^2)^(- p/2) + 2*a^( - p)*(1 + (2*p*b^4)/a^2 + (3*p*
                                                                    (p + 2)*b^8)/(4*a^4)) - 4*wb^( - p/2)*(1 + (3*p*b^4)/(2*wb) + (p*
                                                                                                                                     (p + 2)*b^8)/(2*wb^2)) #HZ variance

  pmu <- log(sqrt(mu^4/(si2 + mu^2))) #lognormal HZ mean
  psi <- sqrt(log((si2 + mu^2)/mu^2)) #lognormal HZ variance

  P <- 1 - plnorm(HZ,pmu,psi) #P-value associated to the HZ statistic



  if (qqplot){
    d <- Dj
    r <- rank(d)
    chi2q <- qchisq((r-0.5)/n,p)
    plot(d, chi2q, pch = 19, main = "Chi-Square Q-Q Plot",
         xlab = "Squared Mahalanobis Distance",ylab="Chi-Square Quantile")
    abline(0, 1,lwd = 2, col = "black")
  }


  result <- list("hz", HZ = HZ, p.value = P, dname = dname, dataframe = dataframe)

  return(result)
}

mvnormalmixEM = function (x, lambda = NULL, mu = NULL, sigma = NULL, k = 2, arbmean = TRUE, arbvar = TRUE,
                          epsilon = 1e-08, maxit = 10000, verb = FALSE)
{
  if(arbmean == FALSE && arbvar == FALSE){
    stop(paste("Must change constraints on mu and/or sigma!","\n"))
  }
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  tmp <- mvnormalmix.init(x = x, lambda = lambda, mu = mu,
                          sigma = sigma, k = k, arbmean=arbmean, arbvar = arbvar)
  lambda <- tmp$lambda
  mu<-tmp$mu
  sigma <- tmp$sigma
  k = tmp$k
  diff <- 1
  iter <- 0
  if (arbmean==FALSE){
    comp <- lapply(1:k, function(i) lambda[i] * dmvnorm(x,
                                                        mu, sigma[[i]]))
  } else{
    if (arbvar==FALSE) {
      comp <- lapply(1:k, function(i) lambda[i] * dmvnorm(x,
                                                          mu[[i]], sigma))
    }
    else comp <- lapply(1:k, function(i) lambda[i] * dmvnorm(x,
                                                             mu[[i]], sigma[[i]]))
  }
  comp <- sapply(comp, cbind)
  compsum <- apply(comp, 1, sum)
  obsloglik <- sum(log(compsum))
  ll <- obsloglik
  restarts <- 0
  while (diff > epsilon & iter < maxit) {
    if (arbvar) {
      z = matrix(nrow = n, ncol = k)
      for (i in 1:n) {
        for (j in 1:k) {
          z.denom = c()
          for (m in 1:k) {
            z.denom = c(z.denom, lambda[m]/lambda[j] *
                          (det(sigma[[j]])/det(sigma[[m]]))^(0.5) *
                          exp(-0.5 * ((x[i, ] - mu[[m]]) %*% solve(sigma[[m]]) %*% t(t(x[i, ] - mu[[m]])) -
                                      (x[i, ] - mu[[j]]) %*% solve(sigma[[j]]) %*% t(t(x[i, ] - mu[[j]])) )))
          }
          z[i, j] = 1/sum(z.denom)
        }
      }
      z = z/apply(z,1,sum)
      #	  z[,k]=1-apply(as.matrix(z[,(1:(k-1))]),1,sum)
      sing <- sum(is.nan(z))
      lambda.new <- apply(z, 2, mean)
      if (sum(lambda.new < 1e-08)>0 || is.na(sum(lambda.new))) {
        sing <- 1
      }
      else {
        if(arbmean==FALSE) {
          mu.new <- lapply(1:k, function(j) sapply(1:p,
                                                   function(i) apply(z * x[, i], 2, sum))[j, ])
          mu.new <- apply(sapply(mu.new,as.vector),1,sum)/n
          mu.new <- lapply(1:k, function(j) mu.new)
        } else{
          mu.new <- lapply(1:k, function(j) sapply(1:p,
                                                   function(i) apply(z * x[, i], 2, sum))[j, ]/sum(z[,
                                                                                                     j]))
        }
        sigma.new <- lapply(1:k, function(j) matrix(apply(sapply(1:n,
                                                                 function(i) z[i, j] * (x[i, ] - mu.new[[j]]) %*%
                                                                   t(x[i, ] - mu.new[[j]])), 1, sum), p, p)/sum(z[,
                                                                                                                  j]))
        lambda <- lambda.new
        mu <- mu.new
        sigma <- sigma.new
        comp <- lapply(1:k, function(i) lambda[i] * dmvnorm(x,
                                                            mu[[i]], sigma[[i]]))
        comp <- sapply(comp, cbind)
        compsum <- apply(comp, 1, sum)
        newobsloglik <- sum(log(compsum))
      }
    }
    else {
      z = matrix(nrow = n, ncol = k)
      sigma.inv = solve(sigma)
      for (i in 1:n) {
        for (j in 1:k) {
          z.denom = c()
          for (m in 1:k) {
            z.denom = c(z.denom, lambda[m]/lambda[j] *
                          (det(sigma.inv)/det(sigma.inv))^(0.5) *
                          exp(-0.5 * ((x[i, ] - mu[[m]]) %*% sigma.inv %*%
                                        t(t(x[i, ] - mu[[m]])) - (x[i, ] - mu[[j]]) %*%
                                        sigma.inv %*% t(t(x[i, ] - mu[[j]])))))
          }
          z[i, j] = 1/sum(z.denom)
        }
      }
      #	  z[,k]=1-apply(as.matrix(z[,(1:(k-1))]),1,sum)
      z = z/apply(z,1,sum)

      sing <- sum(is.nan(z))
      lambda.new <- apply(z, 2, mean)
      if (sum(lambda.new < 1e-08)>0 || is.na(sum(lambda.new))) {
        sing <- 1
      }
      else {
        if(arbmean==FALSE) {
          mu.new <- lapply(1:k, function(j) sapply(1:p,
                                                   function(i) apply(z * x[, i], 2, sum))[j, ])
          mu.new <- apply(sapply(mu.new,as.vector),1,sum)/n
          mu.new <- lapply(1:k, function(j) mu.new)
        } else{
          mu.new <- lapply(1:k, function(j) sapply(1:p,
                                                   function(i) apply(z * x[, i], 2, sum))[j, ]/sum(z[,
                                                                                                     j]))
        }
        temp.sig <- lapply(1:k, function(j) matrix(apply(sapply(1:n,
                                                                function(i) z[i, j] * (x[i, ] - mu.new[[j]]) %*%
                                                                  t(x[i, ] - mu.new[[j]])), 1, sum), p, p))
        sigma.new <- matrix(apply(sapply(temp.sig, as.vector),
                                  1, sum), p, p)/n
        lambda <- lambda.new
        mu <- mu.new
        sigma <- sigma.new
        comp <- lapply(1:k, function(i) lambda[i] * dmvnorm(x,
                                                            mu[[i]], sigma))
        comp <- sapply(comp, cbind)
        compsum <- apply(comp, 1, sum)
        newobsloglik <- sum(log(compsum))
      }
    }
    if (sing > 0 || is.na(newobsloglik) || abs(newobsloglik) == Inf){# || sum(z) != n) {
      cat("Need new starting values due to singularity...",
          "\n")
      restarts <- restarts + 1
      if(restarts>15) stop("Too many tries!")
      tmp <- mvnormalmix.init(x = x, k = k, arbmean=arbmean, arbvar = arbvar)
      lambda <- tmp$lambda
      mu <- tmp$mu
      sigma <- tmp$sigma
      k = tmp$k
      diff <- 1
      iter <- 0
      if (arbvar) {
        comp <- lapply(1:k, function(i) lambda[i] * dmvnorm(x,
                                                            mu[[i]], sigma[[i]]))
      }
      else comp <- lapply(1:k, function(i) lambda[i] *
                            dmvnorm(x, mu[[i]], sigma))
      comp <- sapply(comp, cbind)
      compsum <- apply(comp, 1, sum)
      obsloglik <- sum(log(compsum))
      ll <- obsloglik
    }
    else {
      diff <- newobsloglik - obsloglik
      obsloglik <- newobsloglik
      ll <- c(ll, obsloglik)
      iter <- iter + 1
      if (verb) {
        cat("iteration=", iter, "diff=", diff, "log-likelihood",
            obsloglik, "\n")
      }
    }
  }
  if(arbmean==FALSE) {
    mu = mu[[1]]
  }
  if (iter == maxit) {
    cat("WARNING! NOT CONVERGENT!", "\n")
  }
  colnames(z) <- c(paste("comp", ".", 1:k, sep = ""))
  cat("number of iterations=", iter, "\n")
  a=list(x=x, lambda = lambda, mu = mu, sigma = sigma,
         loglik = obsloglik, posterior = z, all.loglik=ll, restarts=restarts, ft="mvnormalmixEM")
  #class(a) = "mixEM"
  return(a)
}

# have to be careful about using other peoples packages - what if they are not supported in later versions of R?
#packageurl <- "https://cran.r-project.org/src/contrib/MVN_4.0.tar.gz" #http://cran.r-project.org/src/contrib/Archive/XXXX/XXXX_A.B.C.tar.gz"
#install.packages(packageurl, contriburl=NULL, type="source")


find_het_SNPs <- function(data){
  all_bafs <- data.frame(baf=data$B.Allele.Freq,chr=data$Chr)
  lower_bafs <- all_bafs[which(all_bafs[,"baf"]<0.2),]
  upper_bafs <- all_bafs[which(all_bafs[,"baf"]>0.8),]
  lower_thresh <- quantile(aggregate(.~chr,data=lower_bafs,function(x) quantile(x,0.99))[,"baf"],0.5)
  upper_thresh <- quantile(aggregate(.~chr,data=upper_bafs,function(x) quantile(x,0.01))[,"baf"],0.5)
  print(c(lower_thresh,upper_thresh))
  ggplot(mixdata[which(mixdata$B.Allele.Freq<0.1),])+geom_histogram(aes(B.Allele.Freq),binwidth=0.001)
  ggplot(data.frame(mixdata))+geom_point(aes(M,m,col=factor(B.Allele.Freq<lower_thresh | B.Allele.Freq > upper_thresh)))+geom_density2d(aes(M,m))
}

get_mixtures<- function(data, lambda = NULL, mu = NULL, sigma = NULL, k = 2){
  if(nrow(data) <=1){
    return(NA)
  }
  values <- NULL
  attempt <- 0
  while( is.null(values) && attempt <= 100 ) {
    attempt <- attempt + 1
    try(
      values <- mvnormalmixEM(data, lambda = NULL, mu = NULL, sigma = NULL, k = 2,arbmean = TRUE, arbvar = TRUE, epsilon = 1e-04,maxit = 10000, verb = FALSE)
    )
  }
  return(values)
}

read_raw_illumina_file <- function(filename){
  return(as.data.table(fread(filename, sep="\t", sep2="auto", nrows=-1L, header=TRUE, na.strings="NA",
                             stringsAsFactors=FALSE, verbose=getOption("datatable.verbose"), autostart=1L,
                             skip=10, select=NULL, drop=NULL, colClasses=NULL,
                             integer64=getOption("datatable.integer64"),         # default: "integer64"
                             #dec=if (sep!=".") "." else ",", col.names,
                             #check.names=FALSE, encoding="unknown", strip.white=TRUE,
                             showProgress=getOption("datatable.showProgress"),   # default: TRUE
                             data.table=TRUE))) # default: TRUE))
}

read_GC_score_file <- function(filename,max_pos){
  dt <- as.data.table(fread(filename, sep="\t", sep2="auto", nrows=-1L, header=TRUE, na.strings="NA",
                            stringsAsFactors=FALSE, verbose=getOption("datatable.verbose"), autostart=1L,
                            skip=0, select=NULL, drop=NULL, colClasses=NULL,
                            integer64=getOption("datatable.integer64"),         # default: "integer64"
                            #dec=if (sep!=".") "." else ",", col.names,
                            #check.names=FALSE, encoding="unknown", strip.white=TRUE,
                            showProgress=TRUE, #getOption("datatable.showProgress"),   # default: TRUE
                            data.table=TRUE))
  setnames(dt,colnames(dt),c("chr","base","GC"))
  dt$chr <- gsub("chr","",dt$chr)
  dt$chr <- gsub("X","23",dt$chr)
  #data$location <- as.numeric(data$Chr)+as.numeric(data$Position)/max_pos
  dt$location <- as.numeric(dt$chr)+as.numeric(dt$base)/max_pos
  #data$segment <- findInterval(data$location, segments$locationstart)
  return(as.data.frame(dt))
}

read_illumina_annotation_file <- function(filename,max_pos){
  dt <- as.data.table(fread(filename, sep=",", sep2="auto", nrows=-1L, header=TRUE, na.strings="NA",
                            stringsAsFactors=FALSE, verbose=getOption("datatable.verbose"), autostart=1L,
                            skip=7, select=NULL, drop=NULL, colClasses=NULL,
                            integer64=getOption("datatable.integer64"),         # default: "integer64"
                            #dec=if (sep!=".") "." else ",", col.names,
                            #check.names=FALSE, encoding="unknown", strip.white=TRUE,
                            showProgress=TRUE, #getOption("datatable.showProgress"),   # default: TRUE
                            data.table=TRUE))
  return(as.data.frame(dt))
}

read_GC_file <- function(filename,max_pos){
  dt <- as.data.table(fread(filename, sep="auto", sep2="auto", nrows=-1L, header=FALSE, na.strings="NA",
                            stringsAsFactors=FALSE, verbose=getOption("datatable.verbose"), autostart=1L,
                            skip=1, select=NULL, drop=NULL, colClasses=NULL,
                            integer64=getOption("datatable.integer64"),         # default: "integer64"
                            #dec=if (sep!=".") "." else ",", col.names,
                            #check.names=FALSE, encoding="unknown", strip.white=TRUE,
                            showProgress=TRUE, #getOption("datatable.showProgress"),   # default: TRUE
                            data.table=TRUE))
  return(as.data.frame(dt))
}

reinfer_segments <- function(data,segments){
  #data = pdata[[1]]
  #fit = pfit[[1]]
  #segments = psegments[[1]]

  data = data[which(data$Chr != 0),]
  segments = segments[which(segments$chromosome != 0),]

  data$M <- data$CT * data$B.Allele.Freq
  data$m <- data$CT * (1-data$B.Allele.Freq)

  #### find het SNPs if the SNPs are annotated by their expected value ####
  data$ishet <- data$Allele1...AB != data$Allele2...AB
  #Find which segment each fucking things belongs to
  data$segment <- cbind(data$location, findInterval(data$location, segments$locationstart))[,2]

  data$is_outlier = rep(FALSE,nrow(data))
  segments$t1 <- segments$t2 <- segments$t3 <- segments$t4 <- segments$t5 <- segments$t6 <- segments$t7 <- rep(TRUE,nrow(segments))
  #### detect outliers and take a first pass at the datas####
  for(i in 1:nrow(segments)){
    print(i)
    #find the data just for this segment:
    those = which(data$segment == i & complete.cases(data[,c("M","m")]))

    tempdata <- data[those,]
    if(nrow(tempdata) == 0){
      segments[i,"t1"] <- segments[i,"t2"]<- segments[i,"t3"]<- segments[i,"t4"]<- segments[i,"t5"]<- segments[i,"t6"] <- segments[i,"t7"] <- FALSE
      next
    }
    tempdata$M <- tempdata$M/sd(tempdata$M)
    tempdata$m <- tempdata$m/sd(tempdata$m)
    similarity <- 1/(as.matrix(dist(as.matrix(tempdata[,c("M","m")])))^2+0.01)
    # make the similarity of any two points of different catagories 0:
    similarity <- similarity*outer(tempdata$ishet,tempdata$ishet, "==")
    consistent <- data.frame(count = apply(similarity,1,function(x) sum(x)))
    outlier_thresh <- quantile(consistent[,"count"],0.01,na.rm = TRUE)
    tempdata$is_outlier <- data[those,"is_outlier"] <- consistent < outlier_thresh
    #### finish throwing away outliers ####

    #ggplot(tempdata)+geom_point(aes(M,m))

    # infer some coordinates:
    if(nrow(tempdata[which(tempdata$ishet & !tempdata$is_outlier),]) <= 2){
      segments[i,"t1"] <- segments[i,"t2"]<- segments[i,"t3"]<- segments[i,"t4"]<- segments[i,"t5"]<- segments[i,"t6"] <- segments[i,"t7"] <- FALSE
      next
    }
    values <- get_mixtures(tempdata[which(tempdata$ishet & !tempdata$is_outlier),c("M","m")], lambda = NULL, mu = NULL, sigma = NULL, k = 2)

    # lambda is the value of the mixture weights.
    # mu ater the locations of the mixtures
    # sigma is the covariance matrix of each mixture.

    mtest <- mardiaTest(tempdata[which(tempdata$ishet),c("M","m")], qqplot = TRUE)
    hztest <- hzTest(tempdata[which(tempdata$ishet),c("M","m")], cov = TRUE, qqplot = FALSE)

    # 1) the tests that we have is: are the het total copy numbers the same as the ones inferred on the x or y axis?
    # 2) is the theoretical symmetry broken?
    # 3) mtest$p.value.skew > alpha
    # 4) mtest$p.value.kurt > alpha
    # 5) hztest$p.value > alpha
    # 6) are the mixture weights significantly different from 0.5?

    alpha = 0.001
    segments[i,"t1"] <- min(sum(values$mu[[1]]),sum(values$mu[[2]]))/sum(tempdata[which(!tempdata$ishet),c("M","m")])*nrow(tempdata[which(!tempdata$ishet),]) #this should test only one of the two components.
    segments[i,"t2"] <- !((values$mu[[1]][1] > values$mu[[2]][1] & values$mu[[1]][2] < values$mu[[2]][2]) |
                            (values$mu[[1]][1] < values$mu[[2]][1] & values$mu[[1]][2] > values$mu[[2]][2]))
    segments[i,"t3"] <- mtest$p.value.skew
    segments[i,"t4"] <- mtest$p.value.kurt
    segments[i,"t5"] <- hztest$p.value
    segments[i,"t6"] <- abs(min(values$lambda[1]/values$lambda[2],values$lambda[2]/values$lambda[1])-1)
    segments[i,"lambda"] <- values$lambda[1]
    segments[i,"mux1"] <- values$mu[[1]][1]
    segments[i,"muy1"] <- values$mu[[1]][2]
    segments[i,"mux2"] <- values$mu[[2]][1]
    segments[i,"muy2"] <- values$mu[[2]][2]
    segments[i,"s11"] <- values$sigma[[1]][1,1]
    segments[i,"s12"] <- values$sigma[[1]][1,2]
    segments[i,"s21"] <- values$sigma[[1]][2,1]
    segments[i,"s22"] <- values$sigma[[1]][2,2]
    segments[i,"ss11"] <- values$sigma[[2]][1,1]
    segments[i,"ss12"] <- values$sigma[[2]][1,2]
    segments[i,"ss21"] <- values$sigma[[2]][2,1]
    segments[i,"ss22"] <- values$sigma[[2]][2,2]


    segments[i,"mux"] <- mean(tempdata[which(tempdata$ishet),c("M")])
    segments[i,"muy"] <- mean(tempdata[which(tempdata$ishet),c("m")])
    segments[i,"c11"] <- cov(tempdata[which(tempdata$ishet),c("M","m")])[1,1]
    segments[i,"c12"] <- cov(tempdata[which(tempdata$ishet),c("M","m")])[1,2]
    segments[i,"c21"] <- cov(tempdata[which(tempdata$ishet),c("M","m")])[2,1]
    segments[i,"c22"] <- cov(tempdata[which(tempdata$ishet),c("M","m")])[2,2]
    #segments[i,"mux"] <- mean(tempdata[which(tempdata$ishet),c("M")])
  }
  # t1 == true => mixture of two guassians. if we don't have a mixture t1 will fail - it shouldn't be close to zero, if we do have a mixture we hope t1 to suceed
  # t2 == true => no mixture of two guassians. likewise if we don't have amixture we would expect t2 to fail - it should only suceed if the data is hetrogenous, if we do have a mixture we hope t2 to
  # t3,t4,t5 == true => no mixture. t3,t4,t5 will be large and will not fail only if we dont have a mixture
  # t6 == true => no mixture. is more likely to suceed on larger segments, and will only suceed if the data is a mixture
  # t7 being true implies that there was enough data to infer a mixture from - i.e. two points.
  return(segments)
}

preprocess_raw_data <- function(x){
  colnames(x) <- gsub(' ', '.', colnames(x))
  colnames(x) <- gsub('-', '.', colnames(x))

  # force the type that we want:
  x$Chr <- as.integer(as.character(x$Chr)) # I don't think that we want this, come back to it later.
  x$Position <- as.integer(as.character(x$Position))

  # order by genomic location
  x <- x[order(x$Chr,x$Position),]

  # name the rows
  row.names(x) <- x$SNP.Name

  # initial dye-copy number transformation (as used by ascat)
  x[,"CT"] <- 2*2^(x$Log.R.Ratio/0.55)
  return(as.data.frame(x))
}

merge_my_lists <- function(listx,listy){
  listz = list()
  for(i in 1:length(listx)){
    listz[[i]] <- list(x=listx[[i]],y=listy[[i]])
  }
  return(listz)
}