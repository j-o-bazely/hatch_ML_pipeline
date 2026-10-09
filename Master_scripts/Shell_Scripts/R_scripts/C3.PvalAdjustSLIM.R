###          C3. PQLseq P Value adjustment with SLIM                ###

## Date: 06.10.2026

## Author: James B, adapted from Charley Y

## Description:
##  This code runs multiple testing correction on PQLseq-generated P values using SLIM  ##



## Pull iteration, available cores and chromosome ##
args <- commandArgs(trailingOnly = TRUE)
iter <- as.numeric(args[1])
numCore <- as.numeric(args[2])

## Set paths
workdir <- Sys.getenv(
  "CH2_WORKDIR",
  unset = normalizePath(getwd(), winslash = "/", mustWork = FALSE)
)
tmpdir <- file.path(workdir, "tmp")



####    Merge all chromosomes   ####

## Create chromosome vector to loop through ##
chroms <- strsplit(
  Sys.getenv("CH2_CHROMOSOMES", unset = paste(c("chr0", paste0("chr", 1:28)), collapse = ",")),
  ","
)[[1]]

## Create emopty lists to store results ##
List_fit <- list()
List_dfs_fit <- list() 

## Load in data for each chromosome ##
for (chr in chroms) {
  print(chr)
  List_fit[[ chr ]] <- readRDS(file.path(tmpdir, paste0("hatch_trainingSamples_PQLseqOutput_Fit_iter", iter, "_SLK063_ragtag_", chr, ".RDS")))
  List_dfs_fit[[ chr ]] <- List_fit[[ chr ]]
}

##  Bind chromosomes into a single dataframe ##
df_AllChrms <- do.call(rbind, List_dfs_fit)

##  Check the % of NaNs for h2  ##
length(df_AllChrms$h2[grep("NaN", df_AllChrms$h2)])/nrow(df_AllChrms)*100

##  Save object  ##
saveRDS(
  df_AllChrms,
  file.path(tmpdir, paste0("hatch_trainingSamples_PQLseqFit_iter", iter, "_wholeGenome.RDS"))
)


####    Run SLIM correction   ####  


# SLIMcodes.R downloaded from: https://www.aspendb.org/downloads 

#### Define SLIM functions #### 

# Nothing changed from default in SLIMcodes.R (methylKit also uses default settings)
SLIMfunc<-function(rawp,STA=.1,Divi=10,Pz=0.05,B=100,Bplot=FALSE)
{
  ####################
  m <- length(rawp) 
  
  ########################
  alpha_mtx=NULL;#
  pi0s_est_COM=NULL;
  SzCluster_mtx=NULL;
  P_pi1_mtx=NULL;
  pi1_act_mtx=NULL;
  Dist_group=NULL;
  Num_mtx=NULL;
  Gamma=NULL;
  PI0=NULL;
  
  #############
  ##observed points
  lambda_ga=seq(0,1,0.001);
  gamma_ga=sapply(lambda_ga,f1,rawp=rawp);
  Gamma=c(Gamma,gamma_ga);
  alpha_mtx=c(alpha_mtx,gamma_ga[which(lambda_ga==0.05)]);
  
  
  ###esimation
  pi0_mtx=NULL;
  x.axis=NULL;
  y.axis=NULL;
  itv=(1-STA)/Divi;
  for (i in 1:Divi)##10 for uniform data
  {
    cutoff=STA+(i/Divi)*(1-STA);##10 for uniform data
    lambda=seq(cutoff-itv,cutoff,itv/10);
    gamma_mtx=sapply(lambda,f1,rawp=rawp);
    LModel=lm(gamma_mtx~lambda);
    pi0_mtx=c(pi0_mtx,coefficients(LModel)[2]);
  }
  
  
  ##################################
  ########searching
  N_COM=NULL;
  N_rawp=NULL;
  maxFDR_mtx=NULL;
  quapoint_mtx=NULL; 
  if (B<=1) B=100;
  quapoint_mtx=seq(0.01,0.99,1/B);
  for (k in 1:length(quapoint_mtx))
  {
    qua_point=quapoint_mtx[k];
    
    pi0_combLR=min(quantile(pi0_mtx,qua_point),1);#mean(pi0_mtx);#median();# qua_point=0.78 for desreasing distribution;
    ##0.4 for uniform or normal distribution;
    pi0_est=pi0_combLR;
    
    ###########Calculate independent index of raw p vlaues
    PI0=rbind(PI0,pi0_mtx);
    
    pi0s_est_COM=c(pi0s_est_COM,pi0_est);
    ##Condition1
    P_pi1=sort(rawp)[max(length(rawp)*(1-pi0_est),1)];##
    P_pi1_mtx=c(P_pi1_mtx,P_pi1);
    
    pi0=pi0_est;
    if (is.null(Pz)) Pz=0.05;
    maxFDR=Pz*pi0/(1-(1-Pz)*pi0);
    maxFDR_mtx=c(maxFDR_mtx,maxFDR);
    
    qvalues_combLR=QValuesfun(rawp,pi0);
    qvalue_cf=maxFDR;
    selected=which(qvalues_combLR<qvalue_cf);
    Sel_qvalues_combLR=selected;
    
    pi1_act_mtx=c(pi1_act_mtx,length(Sel_qvalues_combLR)/length(rawp));
    N_COM=c(N_COM,list(Sel_qvalues_combLR));
    Num_mtx=c(Num_mtx,length(Sel_qvalues_combLR));
  }
  length(N_COM)
  length(quapoint_mtx)
  
  ####doing judging
  ##by max FDR
  pi1s_est_COM=1-pi0s_est_COM;
  Diff=sum(rawp<=Pz)/length(rawp)-pi1_act_mtx;
  
  ###
  loc=which.min(abs(Diff));
  Diff.loc=Diff[loc];
  selQuantile=quapoint_mtx[loc];
  pi0_Est=min(1,pi0s_est_COM[loc]);
  maxFDR.Pz=Pz*pi0_Est/(1-(1-Pz)*pi0_Est);
  
  if(Bplot)
  {
    windows();
    par(mfrow=c(1,2));
    hist(rawp,main="Histogram of p-value");
    gamma_ga=sapply(lambda_ga,f1,rawp=rawp);
    plot(lambda_ga,gamma_ga,type="l",main="Relationship of p- and q value",xlab=expression(lambda),ylab=expression(gamma),cex.lab=1.45,cex.axis=1.42)
    #par(xaxp=c(0,1,10));
    #axis(1);
    ##qvalues
    qValues=QValuesfun(rawp,pi0=pi0_Est);
    gammaq_ga=sapply(lambda_ga,f1,rawp=qValues);
    lines(lambda_ga,gammaq_ga,col="blue",lwd=2,lty="dashed")
    abline(v=Pz,col="black",lwd=2,lty="dotdash")
    abline(v=maxFDR.Pz,col="blue",lwd=2,lty="dotdash")
    text(0.75,0.6,labels=paste("L=",round(abs(Diff.loc),4),sep=""));
    leg=list(bquote("CPD of p-value"),bquote("CPD of q-value"),bquote("Pmax"==.(Pz)),bquote("FDRmax"==.(round(maxFDR.Pz,2))));
    legend("bottomright",legend=as.expression(leg),lwd=2,lty=c("solid","dashed","dotdash","dotdash"),col=c("black","blue","black","blue"));
  }
  
  return(list(pi0_Est=pi0_Est,selQuantile=selQuantile));
}
f1 <- function(cutoff,rawp){sum(rawp<cutoff)/length(rawp)*100};
QValuesfun<-function(rawp,pi0)
{
  order_rawp=sort(rawp);
  qvalues=pi0*length(order_rawp)*order_rawp/c(1:length(order_rawp));
  temp=cummin(qvalues[seq(length(qvalues),1,-1)])
  qvalues=temp[seq(length(temp),1,-1)];
  qvalues=qvalues[order(order(rawp))]
}

#### Run SLIM on PQLseq output ####

## Load PQLseq whole genome object generated previously ##
df_AllChrms <- readRDS(file.path(tmpdir, paste0("hatch_trainingSamples_PQLseqFit_iter", iter, "_wholeGenome.RDS")))

## Pull p values ##
pvals <- df_AllChrms$pvalue

# Run SLIM on these p values ##
qvals = {QValuesfun(pvals,
                 SLIMfunc(pvals,STA=.1,Divi=10,Pz=0.05,B=100,Bplot=FALSE)$pi0_Est)
}

## Add calculated Q values to PQLseq output ##
slim_output <- df_AllChrms
slim_output$qvalue <- qvals

####    Checks    ####
## Check number of sites which failed to converge ##
algorithmFailure <- nrow(slim_output[slim_output$qvalue == 0 & slim_output$converged == "FALSE",])  # -> All sites that are q=0 also failed to converge -> represent algorithm failure and should be removed

message(paste0(algorithmFailure, ' sites failed to converge and will be removed'))

## Filter out failed convergences
slim_output_filtered <- slim_output[slim_output$converged == "TRUE", ]

## Save filtered object
saveRDS(
  slim_output_filtered,
  file.path(tmpdir, paste0("hatch_trainingSamples_PQLseqFit_iter", iter, "_wholeGenome_corrected.RDS"))
)


