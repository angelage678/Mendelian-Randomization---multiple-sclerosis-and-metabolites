#packages
library(TwoSampleMR)
library(dplyr)
library(MendelianRandomization)

#met-a
data_met_a=ao
data_met_a$new_name_1=sapply(strsplit(data_met_a$id, split= "-"),"[",1)
data_met_a$new_name_2=sapply(strsplit(data_met_a$id, split= "-"),"[",2)
data_met_a=data_met_a[data_met_a$new_name_1=="met",]
data_met_a=data_met_a[data_met_a$new_name_2=="a",]

#met-c
data_met_c=ao
data_met_c$new_name_1=sapply(strsplit(data_met_c$id, split= "-"),"[",1)
data_met_c$new_name_2=sapply(strsplit(data_met_c$id, split= "-"),"[",2)
data_met_c=data_met_c[data_met_c$new_name_1=="met",]
data_met_c=data_met_c[data_met_c$new_name_2=="c",]

#met-d
data_met_d=ao
data_met_d$new_name_1=sapply(strsplit(data_met_d$id, split= "-"),"[",1)
data_met_d$new_name_2=sapply(strsplit(data_met_d$id, split= "-"),"[",2)
data_met_d=data_met_d[data_met_d$new_name_1=="met",]
data_met_d=data_met_d[data_met_d$new_name_2=="d",]

#ao <- available_outcomes()

Pathway_out="/Users/yinhuige/Desktop/"

final_res=data.frame("id.exposure","exposure","id.outcome","outcome", "b_IVW_MRE", "se_IVW_MRE", "pval_IVW_MRE", "b_IVW_FE", "se_IVW_FE", "pval_IVW_FE", "b_Egger", "se_Egger", "pval_Egger", "Egger_intercept", "pval_intercept", "Het_IVW_pval", "Het_Egger_pval", "b_W_Med", "se_W_Med", "pval_W_Med", "b_W_Mod", "se_W_Mod", "pval_W_Mod", "nsnps")

write.table(final_res, file= paste(Pathway_out,"MendelianRandomizationResults.txt", sep=""), col.names = FALSE, append = TRUE, row.names = F, quote = FALSE, na = "-",sep='\t')

all_data_met <- c(data_met_a$id, data_met_c$id, data_met_d$id)
options(scipen = 5)
#running MendelianRandomization package
for (e in all_data_met) {
  for (n in c("ieu-b-18")) {
    try(exp_dat <- extract_instruments(outcomes=e),silent = T)
    #alternate p-value:
    #try(exp_dat <- extract_instruments(outcomes=e),silent = T, p1=1e-6)
    if (exists("exp_dat")==TRUE){
      if (length(exp_dat$SNP)>0) {
        exp_dat$exposure=sapply(strsplit(exp_dat$exposure,fixed = TRUE, split= " ||"),"[",1)
        outcome_dat <-extract_outcome_data(
          snps = exp_dat$SNP,
          outcomes=n)
        if (length(outcome_dat$SNP)>0) {
          outcome_dat$outcome=sapply(strsplit(outcome_dat$outcome,fixed = TRUE, split= " ||"),"[",1)
          #My_MR(exp_dat,outcome_dat)
          rm(dat)
          try(dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = outcome_dat), silent=TRUE)
          if (exists("exp_dat")==TRUE){
            dat=dat[dat$mr_keep==TRUE,]
            dat=dat[is.na(dat$beta.outcome)==FALSE ,]
            if (length(which(dat$mr_keep=='TRUE'))>2) {
              MR_IVW_MRE = MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure, by = dat$beta.outcome, byse = dat$se.outcome), model = "random")
              MR_IVW_FE = MendelianRandomization::mr_ivw(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure, by = dat$beta.outcome, byse = dat$se.outcome), model = "fixed")
              MR_Egger = MendelianRandomization::mr_egger(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                                   by = dat$beta.outcome, byse = dat$se.outcome))
              MR_W_Med = MendelianRandomization::mr_median(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                                    by = dat$beta.outcome, byse = dat$se.outcome), weighting = "weighted")
              MR_W_Mod = mr_mbe(mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                         by = dat$beta.outcome, byse = dat$se.outcome), weighting = "weighted")
              final_res=data.frame(dat$id.exposure[1],dat$exposure[1],dat$id.outcome[1],dat$outcome[1],
                                   MR_IVW_MRE@Estimate,MR_IVW_MRE@StdError,MR_IVW_MRE@Pvalue,
                                   MR_IVW_FE@Estimate,MR_IVW_FE@StdError,MR_IVW_FE@Pvalue,
                                   MR_Egger@Estimate,MR_Egger@StdError.Est,MR_Egger@Pvalue.Est,
                                   MR_Egger@Intercept,MR_Egger@Pvalue.Int,
                                   MR_IVW_FE@Heter.Stat[2],MR_Egger@Heter.Stat[2],
                                   MR_W_Med@Estimate,MR_W_Med@StdError,MR_W_Med@Pvalue,
                                   MR_W_Mod@Estimate,MR_W_Mod@StdError,MR_W_Mod@Pvalue,
                                   MR_IVW_FE@SNPs)
              write.table(final_res, file= paste(Pathway_out,"MendelianRandomizationResults.txt", sep=""), col.names = FALSE, append = TRUE,
                          row.names = F, quote = FALSE, na = "-",sep='\t')
            }
          }
        }
      }
    }
  }
}