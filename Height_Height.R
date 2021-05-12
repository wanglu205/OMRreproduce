#########Height-Height#######################
## .. indicates the work directory
#####load packages
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(AnnotationHub)

########################################################################
##
##    Step 1: Data qc and intersect
##
########################################################################
height_male <- fread("GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.txt")
height_female <- fread("GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_WOMEN_N.txt")
height_male$Freq.Hapmap.Ceu <- as.numeric(height_male$Freq.Hapmap.Ceu)
height_female$Freq.Hapmap.Ceu <- as.numeric(height_female$Freq.Hapmap.Ceu)

height_male$min_allele_frequency <- ifelse(height_male$Freq.Hapmap.Ceu > 0.5,1-height_male$Freq.Hapmap.Ceu,height_male$Freq.Hapmap.Ceu)
height_female$min_allele_frequency <- ifelse(height_female$Freq.Hapmap.Ceu > 0.5,1-height_female$Freq.Hapmap.Ceu,height_female$Freq.Hapmap.Ceu)
height_male_maf <- height_male[which(height_male$min_allele_frequency > 0.05),]
height_female_maf <- height_female[which(height_female$min_allele_frequency > 0.05),]
intersect_RSID<-intersect(height_male_maf$MarkerName,height_female_maf$MarkerName)

bim_all<-c()
for(chr_num in 1:22){
  chr_bim_name<-paste0("../eur_chr_",chr_num,".bim")
  bim_file<-fread(chr_bim_name)
  bim_all<-rbind(bim_all,bim_file)
}
bim_all<-data.frame(bim_all)
colnames(bim_all)<-c("chr","rsid_1000","sex","bp","A1_1000","A2_1000")

merged_file<-merge(height_male_maf, bim_all, by.x = "MarkerName", by.y = "rsid_1000")
merged_file_all <- merge(height_female_maf, merged_file, by = "MarkerName")
merged_file_MHC<-merged_file_all[!(merged_file_all$chr==6&merged_file_all$bp>2*10^7&merged_file_all$bp<3*10^7),]

for(chr_num in 1:22){
  merged_file_MHC_subset<-merged_file_MHC$MarkerName[merged_file_MHC$chr==chr_num]
  rsid_name<-paste0("../chr_",chr_num,".txt")
  write.table(merged_file_MHC_subset,rsid_name,col.names = F,row.names = F,quote = F)
  plink_cmd<-paste0("plink --bfile ../eur_chr_",chr_num," --extract ", rsid_name," --make-bed --out ../eur_chr_",chr_num)
  system(plink_cmd) 
}
plink_cmd<-paste0("plink --bfile eur_chr_1 --merge-list ../merge_list.txt --make-bed --out eur_chr_all")
system(plink_cmd)
ldsc_cmd<-paste0("python ldsc.py --bfile ../eur_chr_all --l2 --ld-wind-kb 10000 --out ../eur_chr_all")
system(ldsc_cmd)
########################################################################
##
##    Step 2: Match with ref_panel on chr bp and Allele Frequency
##
########################################################################
library(data.table)
bim_all<-c()
for(chr_num in 1:22){
  chr_bim_name<-paste0("../eur_chr_",chr_num,".bim")
  bim_file<-fread(chr_bim_name)
  bim_all<-rbind(bim_all,bim_file)
}
bim_all<-data.frame(bim_all)
colnames(bim_all)<-c("chr","rsid_1000","sex","bp","A1_1000","A2_1000")

#--height_male
height_male <- fread("GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_MEN_N.txt")
combined_file<- unique(merge(bim_all,height_male,by.x=c("rsid_1000"),by.y=c("MarkerName")))
combined_file<- combined_file[combined_file$Freq.Hapmap.Ceu!=".",]
combined_file<- combined_file[!duplicated(combined_file$rsid_1000),]
combined_file$A1<-toupper(combined_file$A1)
combined_file$BETA[combined_file$A1!=combined_file$A1_1000]<- (-1)*combined_file$BETA[combined_file$A1!=combined_file$A1_1000]
combined_file$Z<-combined_file$BETA/combined_file$SE.2gc
output_file<-data.frame("chr" = combined_file$chr,"bp" = combined_file$bp,"SNP" = combined_file$rsid_1000,"A1" = combined_file$A1_1000,"A2" = combined_file$A2_1000,"Z" = combined_file$Z,"P" = combined_file$P.2gc,"N" = combined_file$N,"beta"=combined_file$BETA,"se"=combined_file$SE.2gc)
output_file_name<-paste0("../height.male.processed.sumstats.txt")
fwrite(output_file,output_file_name,sep=" ")

#--height_female
height_female <- fread("GIANT_Randall2013PlosGenet_stage1_publicrelease_HapMapCeuFreq_HEIGHT_WOMEN_N.txt")
combined_file<- unique(merge(bim_all,height_female,by.x=c("rsid_1000"),by.y=c("MarkerName")))
combined_file<- combined_file[combined_file$Freq.Hapmap.Ceu!=".",]
combined_file<- combined_file[!duplicated(combined_file$rsid_1000),]
combined_file$A1<-toupper(combined_file$A1)
combined_file$BETA[combined_file$A1!=combined_file$A1_1000]<- (-1)*combined_file$BETA[combined_file$A1!=combined_file$A1_1000]
combined_file$Z<-combined_file$BETA/combined_file$SE.2gc
output_file<-data.frame("chr" = combined_file$chr,"bp" = combined_file$bp,"SNP" = combined_file$rsid_1000,"A1" = combined_file$A1_1000,"A2" = combined_file$A2_1000,"Z" = combined_file$Z,"P" = combined_file$P.2gc,"N" = combined_file$N,"beta"=combined_file$BETA,"se"=combined_file$SE.2gc)
output_file_name<-paste0("/net/mulan/home/yuef/realdata/Height/height.female.processed.sumstats.txt")
fwrite(output_file,output_file_name,sep=" ")
########################################################################
##
##    Step 3: MR analysis
##
########################################################################
library(data.table)
library(OMR)
library(MendelianRandomization)
library(MRMix)
library(BWMR)
library(ggplot2)
library(readr)
library(dplyr)
library(cause)
Z1<-fread("../height.male.processed.sumstats.txt")
Z2 <- fread("../height.female.processed.sumstats.txt")
##OMR
LDscorein<-read.table("../eur_chr_all.l2.ldscore.gz",header = T)
LDscore<-LDscorein[LDscorein$SNP %in% Z1$SNP,]
LDscore<-LDscore[match(Z1$SNP,LDscore$SNP),]
Z2<-Z2[match(Z1$SNP,Z2$SNP),]
l.j <- LDscore$L2
Z <- cbind(Z2$Z,Z1$Z)
n2 <- round(mean(Z1$N))
n1 <- round(mean(Z2$N))
num.per <- 503
num.snp <- nrow(Z)
numCore = 25
OMR_out <- omr(n1,n2,num.per,l.j,Z,coreNum)

##SNP clumping 
clump <- data.frame("SNP"=Z1$SNP,"P"=Z1$P)
clump_name <- paste0("../clump_data.txt")
fwrite(clump,file=clump_name, quote=F,sep=" ")
out_name <- paste0("../clump_data")
cmd <- paste0("plink --bfile ../eur_chr_all --clump ",clump_name," --clump-p1 5e-08 --clump-r2 0.1 --clump-kb 10000 --out ",out_name)
system(cmd)
system(paste0("rm -f ",clump_name))
ind_sig_snp <- fread(paste0(out_name,".clumped"))#change name
dim(ind_sig_snp)
numIV = nrow(ind_sig_snp)

##IVW
f <- match(ind_sig_snp$SNP,Z1$SNP)
Zx <- Z1$beta[f]
Zy <-  Z2$beta[f]
mr.obj = mr_input(Zx, Z1$se[f],Zy , Z2$se[f],snps = ind_sig_snp$SNP)
IVW = mr_ivw(mr.obj)

##Egger
egger = mr_egger(mr.obj)

##MRMix
est = MRMix(Zx, Zy ,Z1$se[f], Z2$se[f])

##BWMR
bwmr_est <- BWMR(Zx, Zy ,Z1$se[f], Z2$se[f])

##CAUSE
cause_data <- data.frame("snp"=Z1$SNP,"beta_hat_1"=Z1$beta,"seb1"=Z1$se,"beta_hat_2"=Z2$beta,"seb2"=Z2$se,"A1"=Z1$A1,"A2"=Z1$A2)
pruned_SNP <- c()
for (chr in 1:22){
  name1 <- paste0("chr",chr,"_AF0.05_0.1.RDS")
  name2 <- paste0("chr",chr,"_AF0.05_snpdata.RDS")
  ld <- readRDS(paste0("LD/",name1))
  snp_info <- readRDS(paste0("LD/",name2))
  variants <- cause_data %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))
  variants$snp <- as.character(variants$snp)
  pruned <- ld_prune(variants = variants, 
                     ld = ld, total_ld_variants = snp_info$SNP, 
                     pval_cols = c("pval1"), 
                     pval_thresh = c(1e-3))
  pruned_SNP <- c(pruned_SNP,pruned)
  print(chr)
}
X <- new_cause_data(cause_data)
set.seed(100)
varlist <- with(X, sample(snp, size=300000, replace=FALSE))
params <- est_cause_params(X, varlist)
num_cause <- length(pruned_SNP)
res <- cause(X=X, variants = pruned_SNP, param_ests = params)
cause_p <- summary(res)$p
ci_size=0.95
fit <- res[["causal"]]
ix <- which(fit$params == "gamma")
qs <- with(fit$marge_post[[ix]], step_quantile(c(0.5, (1-ci_size)/2, 1-((1-ci_size)/2)),
                                               begin, end, post))
cause_est <- qs[1]
cause_est_lower <- qs[2]
cause_est_upper <- qs[3]
cause_num = num_cause
   
#---------------save result-------------------------#
MRres = list()
MRres$OMR_est = OMR_out$alpha
MRres$OMR_sd = OMR_out$se
MRres$OMR_P = OMR_out$pvalue
MRres$IVW_est = IVW$Estimate
MRres$IVW_sd = IVW$StdError
MRres$IVW_p = IVW$Pvalue
MRres$egger_est = egger$Estimate
MRres$egger_sd = egger$StdError.Est
MRres$egg_p = egger$Pvalue.Est
MRres$mix_est = est$theta
MRres$mix_sd = est$SE_theta
MRres$mix_p = data.frame(est)$pvalue_theta
MRres$bwmr_est = bwmr_est$beta
MRres$bwmr_sd = bwmr_est$se_beta
MRres$bwmr_p = bwmr_est$P_value
MRres$numIV = numIV
MRres$cause_est = cause_est
MRres$cause_est_lower = cause_est_lower
MRres$cause_est_upper = cause_est_upper
MRres$cause_p = cause_p
MRres$cause_num = num_cause

