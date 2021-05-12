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
cad <- fread("../CARDIoGRAM_GWAS_RESULTS.txt")
cad$min_allele_frequency <- ifelse(cad$ref_allele_frequency > 0.5,1-cad$ref_allele_frequency,cad$ref_allele_frequency)
cad_maf <- cad[which(cad$min_allele_frequency > 0.05),]
df <- data.frame(cbind(matrix(unlist(strsplit(as.character(cad_maf$chr_pos_),":")),byrow=T,ncol=2)[,1], matrix(unlist(strsplit(as.character(cad_maf$chr_pos_),":")),byrow=T,ncol=2)[,2],
matrix(unlist(strsplit(as.character(cad_maf$chr_pos_),":")),byrow=T,ncol=2)[,2],cad_maf$log_odds,cad_maf$SNP))
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

bim_all<-fread("../eur_chr_all.bim")
colnames(bim_all)<-c("chr","rsid_1000","sex","bp","A1_1000","A2_1000")

merged_file<-merge(pos_hg18_to_hg19, bim_all, by.x = c("chr", "bp"), by.y = c("chr", "bp"))#611132 SNP

asthma <- fread("../TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv")
colnames(asthma)[3] <- "bp"
asthma$chr <- as.character(asthma$chr)
merged_file_all <- merge(asthma, merged_file, by.x = c("chr", "bp"), by.y = c("chr", "bp"))

merged_file_MHC<-merged_file_all[!(merged_file_all$chr==6&merged_file_all$bp>2*10^7&merged_file_all$bp<3*10^7),]#610651 SNP

for(chr_num in 1:22){  
  merged_file_MHC_subset<-merged_file_MHC$rsid_1000[merged_file_MHC$chr==chr_num]
  rsid_name<-paste0("../chr_",chr_num,".txt")
  write.table(merged_file_MHC_subset,rsid_name,col.names = F,row.names = F,quote = F)
  plink_cmd<-paste0("../plink --bfile ../eur_chr_",chr_num," --extract ", rsid_name," --make-bed --out ../eur_chr_",chr_num)
  system(plink_cmd) 
}

merge_list<-cbind(paste0(rep("eur_chr_",21),seq(2,22,1),rep(".bed",22)),paste0(rep("eur_chr_",21),seq(2,22,1),rep(".bim",22)),paste0(rep("eur_chr_",21),seq(2,22,1),rep(".fam",22)))
write.table(merge_list,"../merge_list.txt",col.names = F,row.names = F,quote=F)
plink_cmd<-paste0("../plink --bfile ../eur_chr_1 --merge-list ../merge_list.txt --make-bed --out ../eur_chr_all")
system(plink_cmd)
ldsc_cmd<-paste0("python ../ldsc.py --bfile ../eur_chr_all --l2 --ld-wind-kb 10000 --out ../eur_chr_all")
system(ldsc_cmd)
########################################################################
##
##    Step 2: Match with ref_panel on chr bp and Allele Frequency
##
########################################################################
bim_all<-c()
for(chr_num in 1:22){
  chr_bim_name<-paste0("../eur_chr_",chr_num,".bim")
  bim_file<-fread(chr_bim_name)
  bim_all<-rbind(bim_all,bim_file)
}
bim_all<-data.frame(bim_all)
colnames(bim_all)<-c("chr","rsid_1000","sex","bp","A1_1000","A2_1000")

#--cad
cad <- fread("../CARDIoGRAM_GWAS_RESULTS.txt")
cad_subset<-cad
cad_subset$start <-matrix(unlist(strsplit(as.character(cad_subset$chr_pos),":")),byrow=T,ncol=2)[,2]
cad_subset$end <-matrix(unlist(strsplit(as.character(cad_subset$chr_pos),":")),byrow=T,ncol=2)[,2]
cad_subset$CHR<-matrix(unlist(strsplit(as.character(cad_subset$chr_pos),":")),byrow=T,ncol=2)[,1]

gr <- makeGRangesFromDataFrame(cad_subset, TRUE)
pos_hg18_to_hg19 <- liftOver(gr, chain)
pos_hg18_to_hg19 <- as.data.frame(unlist(pos_hg18_to_hg19))
pos_hg18_to_hg19$chr<-matrix(unlist(strsplit(as.character(pos_hg18_to_hg19$seqnames),"chr")),byrow=T,ncol=2)[,2]
pos_hg18_to_hg19$bp<-pos_hg18_to_hg19$start

combined_file<-merge(bim_all,pos_hg18_to_hg19,by.x=c("chr","bp" ),by.y=c("chr","bp"))
combined_file$A1<-toupper(combined_file$reference_allele)

combined_file$log_odds[combined_file$A1!=combined_file$A1_1000]<- (-1)*combined_file$log_odds[combined_file$A1!=combined_file$A1_1000]
combined_file$Z<-combined_file$log_odds/combined_file$log_odds_se
  
output_file<-data.frame("chr" = combined_file$chr,"bp" = combined_file$bp,"SNP" = combined_file$rsid_1000,"A1" = combined_file$A1_1000,"A2" = combined_file$A2_1000,"Z" = combined_file$Z,"P" = combined_file$pvalue,"N" = combined_file$N_case + combined_file$N_control,"beta"=combined_file$log_odds,"se"=combined_file$log_odds_se)
output_file_name<-paste0("../CAD_test.processed.sumstats.txt")
fwrite(output_file,output_file_name,sep=" ")

#--ASTHMA
asthma <- fread("../TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv")
colnames(asthma)[3] <- "bp"
asthma$chr <- as.character(asthma$chr)

combined_file<-merge(bim_all,asthma,by.x=c("chr","bp" ),by.y=c("chr","bp"))
combined_file$A1<-toupper(combined_file$alternate_allele)

combined_file$European_ancestry_beta_rand[combined_file$A1!=combined_file$A1_1000]<- (-1)*combined_file$European_ancestry_beta_rand[combined_file$A1!=combined_file$A1_1000]
combined_file$Z<-combined_file$European_ancestry_beta_rand/combined_file$European_ancestry_se_rand
  
output_file<-data.frame("chr" = combined_file$chr,"bp" = combined_file$bp,"SNP" = combined_file$rsid_1000,"A1" = combined_file$A1_1000,"A2" = combined_file$A2_1000,"Z" = combined_file$Z,"P" = combined_file$European_ancestry_pval_rand
,"N" = 19954+107715,"beta"=combined_file$European_ancestry_beta_rand,"se"=combined_file$European_ancestry_se_rand)
output_file_name<-paste0("../ASTHMA.processed.sumstats.txt")
fwrite(output_file,output_file_name,sep=" ")

#---TRAIT X
trait_X = list.files()
for(i in 1:length(trait_X)){
exp <- fread(paste0("../",trait_X[i]))
combined_file<-merge(bim_all,exp,by.x=c("chr","bp" ),by.y=c("chr","bp"))
print(dim(combined_file))
output_file_name<-paste0("../",trait_X[i])
fwrite(combined_file,output_file_name,sep=" ")
print(i)
}
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
args<-as.numeric(commandArgs(TRUE))
index_X<-args[1]
index_Y<-args[2]

trait_X = list.files("/net/fantasia/home/borang/data2/real_data/processed_sumstat_less_SNP_4")
asthma = "../ASTHMA.processed.sumstats.txt"
cad = "../CAD_test.processed.sumstats.txt"
trait_Y = c(asthma,cad)
x_name = unlist(strsplit(trait_X[index_X],".", fixed= T ))[1]
   
summary_name_1<-trait_X[index_X]
summary_name_2<-trait_Y[index_Y]

Z1<-fread(summary_name_1)
Z2<-fread(summary_name_2)

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
two_study=T
numCore = 25
OMR_out <- omr(n1,n2,num.per,l.j,Z,coreNum)
##SNP clumping 
bimfile <- read.table("../eur_chr_all.bim",stringsAsFactors = F)
clump <- data.frame("SNP"=bimfile$V2,"P"=px)
clump_name <- paste0(trait_X[index_X],".clumped.txt")
fwrite(clump,file=clump_name, quote=F,sep=" ")
out_name <- paste0(trait_X[index_X],".clumped")
cmd <- paste0("../plink --bfile ../eur_chr_all --allow-no-sex --clump ",clump_name," --clump-p1 5e-08 --clump-r2 0.1 --clump-kb 10000 --out ",out_name)
system(cmd)
system(paste0("rm -f ",clump_name))

##IVW
ind_sig_snp <- fread(out_name)
numIV = nrow(ind_sig_snp)
f <- match(ind_sig_snp$SNP,bimfile$V2)
Zx <- Z1$Z[f]/sqrt(Z1$N[f])
Zy <-  Z2$beta[f]
mr.obj = mr_input(Zx, 1/sqrt(Z1$N[f]),Zy , Z2$se[f],snps = ind_sig_snp$SNP)
IVW = mr_ivw(mr.obj)
  
##Egger
egger = mr_egger(mr.obj)

##MRMix
est = MRMix(Zx, Zy ,1/sqrt(Z1$N[f]), Z2$se[f])

##BWMR
bwmr_est <- BWMR(Zx, Zy ,1/sqrt(Z1$N[f]), Z2$se[f])

##CAUSE
X <- data.frame("snp"=Z1$rsid_1000,"beta_hat_1"=Z1$Z/sqrt(Z1$N),"seb1"=1/sqrt(Z1$N),"beta_hat_2"=Z2$beta,"seb2"=Z2$se,"A1"=Z1$A1_1000,"A2"=Z1$A2_1000,stringsAsFactors=F)
X <- new_cause_data(X)
set.seed(100)
varlist <- with(X, sample(snp, size=300000, replace=FALSE))
params <- est_cause_params(X, varlist)
pruned_SNP <- c()
for (chr in 1:22){
name1 <- paste0("chr",chr,"_AF0.05_0.1.RDS")
name2 <- paste0("chr",chr,"_AF0.05_snpdata.RDS")
ld <- readRDS(paste0("LD/",name1))
snp_info <- readRDS(paste0("LD/",name2))
variants <- X %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))
variants$snp <- as.character(variants$snp)
pruned <- ld_prune(variants = variants, 
                            ld = ld, total_ld_variants = snp_info$SNP, 
                            pval_cols = c("pval1"), 
                            pval_thresh = c(1e-3))
pruned_SNP <- c(pruned_SNP,pruned)
print(chr)
}
num_cause <- length(pruned_SNP)
res <- cause(X=X, variants = pruned_SNP, param_ests = params)
cause_z <- with(res$elpd, z[model1=="sharing" & model2=="causal"])
cause_p <- pnorm(cause_z)
ci_size=0.95
fit <- res[["causal"]]
ix <- which(fit$params == "gamma")
qs <- with(fit$marge_post[[ix]], step_quantile(c(0.5, (1-ci_size)/2, 1-((1-ci_size)/2)),
                                                       begin, end, post))
cause_est <- qs[1]
cause_est_lower <- qs[2]
cause_est_upper <- qs[3]

#---------------save result-------------------------#
MRres = list()
MRres$OMR_est = OMR_out$alpha
MRres$OMR_sd = OMR_out$se
MRres$OMR_P = OMR_out$pvalue
MRres$IVW_est = IVW$Estimate
MRres$IVW_sd = IVW$StdError
MRres$IVW_p = IVW$Pvalue
res$egger_est = egger$Estimate
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
  
  
 

