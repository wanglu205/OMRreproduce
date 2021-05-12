#########CAD-CAD#############
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
#---------maf 0.05 intersect-------------------------#
c4d <- fread("../C4D_CAD_DISCOVERY_METAANALYSIS_UPDATE.TXT")
cad <- fread("../CARDIoGRAM_GWAS_RESULTS.txt")
c4d$min_allele_frequency <- ifelse(c4d$REFERENCE_ALLELE_FREQ > 0.5,1-c4d$REFERENCE_ALLELE_FREQ,c4d$REFERENCE_ALLELE_FREQ)
cad$min_allele_frequency <- ifelse(cad$ref_allele_frequency > 0.5,1-cad$ref_allele_frequency,cad$ref_allele_frequency)
c4d_maf <- c4d[which(c4d$min_allele_frequency > 0.05),]
cad_maf <- cad[which(cad$min_allele_frequency > 0.05),]
intersect_RSID<-intersect(c4d_maf$SNP,cad_maf$SNP)

#--------- change hg18 to hg19
c4d_subset<-c4d[c4d$SNP%in%intersect_RSID,]
df <- data.frame(cbind(matrix(unlist(strsplit(as.character(c4d_subset$CHR_POSB36),"_")),byrow=T,ncol=2)[,1], matrix(unlist(strsplit(as.character(c4d_subset$CHR_POSB36),"_")),byrow=T,ncol=2)[,2],
matrix(unlist(strsplit(as.character(c4d_subset$CHR_POSB36),"_")),byrow=T,ncol=2)[,2],c4d_subset$LOG_ODDS,c4d_subset$SNP))
colnames(df) <- c('chr', 'start', 'end',"beta","rsid")
gr <- makeGRangesFromDataFrame(df, TRUE)
ahub <- AnnotationHub()
ahub.chain <- subset(ahub, rdataclass == "ChainFile" & species == "Homo sapiens")
chain <- ahub.chain[ahub.chain$title == "hg18ToHg19.over.chain.gz"]
chain <- chain[[1]]

pos_hg18_to_hg19 <- liftOver(gr, chain)
pos_hg18_to_hg19 <- as.data.frame(unlist(pos_hg18_to_hg19))
pos_hg18_to_hg19$chr<-matrix(unlist(strsplit(as.character(pos_hg18_to_hg19$seqnames),"chr")),byrow=T,ncol=2)[,2]
pos_hg18_to_hg19$bp<-pos_hg18_to_hg19$start

library(data.table)
bim_all<-c()
for(chr_num in 1:22){
  chr_bim_name<-paste0("../eur_chr_",chr_num,".bim")
  bim_file<-fread(chr_bim_name)
  bim_all<-rbind(bim_all,bim_file)
}
bim_all<-data.frame(bim_all)
colnames(bim_all)<-c("chr","rsid_1000","sex","bp","A1_1000","A2_1000")

merged_file<-merge(pos_hg18_to_hg19, bim_all, by.x = c("chr", "bp"), by.y = c("chr", "bp"))
cad_subset<-cad[cad$SNP%in%intersect_RSID,]
df <- data.frame(cbind(matrix(unlist(strsplit(as.character(cad_subset$chr_pos_),":")),byrow=T,ncol=2)[,1], matrix(unlist(strsplit(as.character(cad_subset$chr_pos_),":")),byrow=T,ncol=2)[,2],
matrix(unlist(strsplit(as.character(cad_subset$chr_pos_),":")),byrow=T,ncol=2)[,2],cad_subset$log_odds,cad_subset$SNP))
colnames(df) <- c('chr', 'start', 'end',"beta","rsid")

gr <- makeGRangesFromDataFrame(df, TRUE)
pos_hg18_to_hg19 <- liftOver(gr, chain)
pos_hg18_to_hg19 <- as.data.frame(unlist(pos_hg18_to_hg19))
pos_hg18_to_hg19$chr<-matrix(unlist(strsplit(as.character(pos_hg18_to_hg19$seqnames),"chr")),byrow=T,ncol=2)[,2]
pos_hg18_to_hg19$bp<-pos_hg18_to_hg19$start

merged_file_all <- merge(pos_hg18_to_hg19, merged_file, by.x = c("chr", "bp"), by.y = c("chr", "bp"))
merged_file_MHC<-merged_file_all[!(merged_file_all$chr==6&merged_file_all$bp>28477797&merged_file_all$bp<33448354),]

for(chr_num in 1:22){
  merged_file_MHC_subset<-merged_file_MHC$rsid_1000[merged_file_MHC$chr==chr_num]
  rsid_name<-paste0("../chr_",chr_num,".txt")
  write.table(merged_file_MHC_subset,rsid_name,col.names = F,row.names = F,quote = F)
  plink_cmd<-paste0("../plink --bfile ../eur_chr_",chr_num," --extract ", rsid_name," --make-bed --out ../eur_chr_",chr_num)
  system(plink_cmd) 
}

merge_list<-cbind(paste0(rep("eur_chr_",21),seq(2,22,1),rep(".bed",22)),paste0(rep("eur_chr_",21),seq(2,22,1),rep(".bim",22)),paste0(rep("eur_chr_",21),seq(2,22,1),rep(".fam",22)))
write.table(merge_list,"../merge_list.txt",col.names = F,row.names = F,quote=F)
plink_cmd<-paste0("plink --bfile ../eur_chr_1 --merge-list ../merge_list.txt --make-bed --out ../eur_chr_all")
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

#--c4d
c4d <- fread("../C4D_CAD_DISCOVERY_METAANALYSIS_UPDATE.TXT")
c4d$start <-matrix(unlist(strsplit(as.character(c4d$CHR_POSB36),"_")),byrow=T,ncol=2)[,2]
c4d$end <-matrix(unlist(strsplit(as.character(c4d$CHR_POSB36),"_")),byrow=T,ncol=2)[,2]
c4d$CHR<-matrix(unlist(strsplit(as.character(c4d$CHR_POSB36),"_")),byrow=T,ncol=2)[,1]
gr <- makeGRangesFromDataFrame(c4d, TRUE)
ahub <- AnnotationHub()
ahub.chain <- subset(ahub, rdataclass == "ChainFile" & species == "Homo sapiens")
chain <- ahub.chain[ahub.chain$title == "hg18ToHg19.over.chain.gz"]
chain <- chain[[1]]

pos_hg18_to_hg19 <- liftOver(gr, chain)
pos_hg18_to_hg19 <- as.data.frame(unlist(pos_hg18_to_hg19))
pos_hg18_to_hg19$chr<-matrix(unlist(strsplit(as.character(pos_hg18_to_hg19$seqnames),"chr")),byrow=T,ncol=2)[,2]
pos_hg18_to_hg19$bp<-pos_hg18_to_hg19$start

combined_file<-merge(bim_all,pos_hg18_to_hg19,by.x=c("chr","bp" ),by.y=c("chr","bp"))
combined_file$A1<-toupper(combined_file$OTHER_ALLELE)

combined_file$LOG_ODDS[combined_file$A1!=combined_file$A1_1000]<- (-1)*combined_file$LOG_ODDS[combined_file$A1!=combined_file$A1_1000]
combined_file$Z<-combined_file$LOG_ODDS/combined_file$LOG_ODDS_SE

output_file<-data.frame("chr" = combined_file$chr,"bp" = combined_file$bp,"SNP" = combined_file$rsid_1000,"A1" = combined_file$A1_1000,"A2" = combined_file$A2_1000,"Z" = combined_file$Z,"P" = combined_file$PVALUE,"N" = combined_file$N_CASE + combined_file$N_CONTROL,"beta"=combined_file$LOG_ODDS,"se"=combined_file$LOG_ODDS_SE)
output_file_name<-paste0("../C4D.processed.sumstats.txt")
fwrite(output_file,output_file_name,sep=" ")

#--cad
cad <- fread("../CARDIoGRAM_GWAS_RESULTS.txt")
cad_subset<-cad[cad$SNP%in%c4d$SNP,]
cad_subset$start <-matrix(unlist(strsplit(as.character(cad_subset$chr_pos),":")),byrow=T,ncol=2)[,2]
cad_subset$end <-matrix(unlist(strsplit(as.character(cad_subset$chr_pos),":")),byrow=T,ncol=2)[,2]
cad_subset$CHR<-matrix(unlist(strsplit(as.character(cad_subset$chr_pos),":")),byrow=T,ncol=2)[,1]

gr <- makeGRangesFromDataFrame(cad_subset, TRUE)
pos_hg18_to_hg19 <- liftOver(gr, chain)
pos_hg18_to_hg19 <- as.data.frame(unlist(pos_hg18_to_hg19))
pos_hg18_to_hg19$chr<-matrix(unlist(strsplit(as.character(pos_hg18_to_hg19$seqnames),"chr")),byrow=T,ncol=2)[,2]
pos_hg18_to_hg19$bp<-pos_hg18_to_hg19$start

combined_file<-merge(bim_all,pos_hg18_to_hg19,by.x=c("chr","bp" ),by.y=c("chr","bp"))
combined_file$A1<-toupper(combined_file$other_allele)

combined_file$log_odds[combined_file$A1!=combined_file$A1_1000]<- (-1)*combined_file$log_odds[combined_file$A1!=combined_file$A1_1000]
combined_file$Z<-combined_file$log_odds/combined_file$log_odds_se  

output_file<-data.frame("chr" = combined_file$chr,"bp" = combined_file$bp,"SNP" = combined_file$rsid_1000,"A1" = combined_file$A1_1000,"A2" = combined_file$A2_1000,"Z" = combined_file$Z,"P" = combined_file$pvalue,"N" = combined_file$N_case + combined_file$N_control,"beta"=combined_file$log_odds,"se"=combined_file$log_odds_se)
output_file_name<-paste0("../CAD.processed.sumstats.txt")
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
Z1<-fread("../C4D.processed.sumstats.txt")
Z2 <- fread("../CAD.processed.sumstats.txt")
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
ind_sig_snp <- fread(paste0(out_name,".clumped"))
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

