########LOAD PACKAGES################
args<- as.numeric(commandArgs(TRUE))
library(mvtnorm)
library(MASS)
library(snpStats)
library(Rcpp)
library(RcppArmadillo)
library(inline)
library(data.table)
library(OMR)
library(MendelianRandomization)
library(MRMix)
library(BWMR)
library(ggplot2)
library(readr)
library(dplyr)
library(cause) 
sourceCpp("./funs.cpp")

##set simulation parameters
sigma <- args[1]
gamma <- args[2]
rho_idx<- args[3]
alpha_idx<- args[4]
j <- args[5]
sigma_beta_index <- c(0.25,0.5,0.125)
sigma_gamma_index <- c(0.25,0.5,0)
alpha_true_index <- c(0,0.25,0.15,0.35)
rho_index <- c(0)
num.per<-1000
num.snp <-345373
sigma_beta_true<-sigma_beta_index[sigma]/num.snp ###The scale is at heritability/p level
sigma_gamma_true<-sigma_gamma_index[gamma]/num.snp
rho<-rho_index[rho_idx]
alpha_true<-alpha_true_index[alpha_idx]
n1<-30000
n2<-30000
omega_matrix<-matrix(c(1,rho,rho,1),ncol=2,nrow=2)

##import genotype matirx
sample <- read.plink("Kaiser_chr1_1.bed","Kaiser_chr1_1.bim", "Kaiser_chr1_1.fam")
gen1 <- as(sample$genotypes, "numeric")
gen1=-gen1+2
sample <- read.plink("Kaiser_chr1_2.bed", "Kaiser_chr1_2.bim", "Kaiser_chr1_2.fam")
gen2 <- as(sample$genotypes, "numeric")
gen2=-gen2+2

##simulate summary statitics
load("/net/mulan/home/walu/Kaiser/data/Kaiser_chr_1_total.bim.Rdata")
pos <- bimfile[,4]/1e6
B <- as.matrix(cbind(rnorm(num.snp,0,sqrt(sigma_beta_true)),rnorm(num.snp,0,sqrt(sigma_gamma_true))))
C <- matrix(1,length(pos),2)
beta_exp1 <- calExp2(pos,B,C,1,0.1,12,3)
beta_exp2 <- calExp2(pos,B,C,1,0.1,12,2)
eu <- rmvnorm(num.per,mean = rep(0,2),omega_matrix)
ex <- lmRcpp(gen1, as.matrix(eu[,1]))
ey <- lmRcpp(gen2, as.matrix(eu[,2]))
betax <- beta_exp1[,1]+ex*sqrt(num.per)/sqrt(n2)
sex <- 1/sqrt(n2)
betay <- alpha_true*beta_exp2[,1] + beta_exp2[,2]+ey*sqrt(num.per)/sqrt(n2)
sey <- 1/sqrt(n1)
px <- 2*(1-pnorm(sqrt(n2)*abs(betax)))
rm(B,C,beta_exp1,beta_exp2,eu,ex,ey,betax,betay,sex,sey,px,name)  

##perform MR analysis
##OMR
zx <- betax/sex
zy <- betay/sey
Z <- cbind(zy,zx)
load("/net/mulan/home/walu/Kaiser/data/Kaiser_chr_1_total.l2.ldscore.Rdata")
two_study=T;coreNum = 1
OMR_out <- omr(n1,n2,num.per,l.j,Z,coreNum)

##SNP clumping 
clump <- data.frame("SNP"=bimfile$V2,"P"=px)
clump_name <- paste0("z_two_sigma",sigma_beta_index[sigma],"_gamma",sigma_gamma_index[gamma],"_rho",rho,"_alpha",alpha_true,"_",s,".txt")
fwrite(clump,file=clump_name, quote=F,sep=" ")
out_name <- paste0("z_two_sigma",sigma_beta_index[sigma],"_gamma",sigma_gamma_index[gamma],"_rho",rho,"_alpha",alpha_true,"_",s)
cmd <- paste0("./plink --bfile Kaiser_chr_1_total --allow-no-sex --clump ",clump_name," --clump-p1 5e-08 --clump-r2 0.1 --clump-kb 10000 --out ",out_name)
system(cmd)
system(paste0("rm -f ",clump_name))

##IVW
ind_sig_snp <- fread(paste0(out_name,".clumped"))#change name
numIV = nrow(ind_sig_snp)
f <- match(ind_sig_snp$SNP,bimfile$V2)
bx <- betax[f]
bxse <- rep(sex, length(bx))
by <-  betay[f]
byse <- rep(sey, length(by))
mr.obj = mr_input(bx, bxse,by , byse,snps = ind_sig_snp$SNP)
IVW = mr_ivw(mr.obj)

##Egger
egger = mr_egger(mr.obj)

##MRMix
est = MRMix(bx, by, bxse, byse)

##BWMR
bwmr_est <- BWMR(bx, by, bxse, byse)

##CAUSE
X <- data.frame("snp"=bimfile$V2,"beta_hat_1"=betax,"seb1"=sex,"beta_hat_2"=betay,"seb2"=sey,"A1"=bimfile$V5,"A2"=bimfile$V6,stringsAsFactors=F)
X <- new_cause_data(X)
set.seed(100)
varlist <- with(X, sample(snp, size=300000, replace=FALSE))
params <- est_cause_params(X, varlist)
ld <- readRDS("LD/chr1_AF0.05_0.1.RDS")
snp_info <- readRDS("LD/chr1_AF0.05_snpdata.RDS")
variants <- X %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))
variants$snp <- as.character(variants$snp)
pruned <- ld_prune(variants = variants, 
                            ld = ld, total_ld_variants = snp_info$SNP, 
                            pval_cols = c("pval1"), 
                            pval_thresh = c(1e-3))
num_cause <- length(pruned)
res <- cause(X=X, variants = pruned, param_ests = params)
cause_p <- summary(res)$p
ci_size=0.95
fit <- res[["causal"]]
ix <- which(fit$params == "gamma")
qs <- with(fit$marge_post[[ix]], step_quantile(c(0.5, (1-ci_size)/2, 1-((1-ci_size)/2)),
                                                       begin, end, post))
cause_est <- qs[1]
cause_est_lower <- qs[2]
cause_est_upper <- qs[3]

##MR results
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

 
