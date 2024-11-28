#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# WORKFLOW FOR SERUM JOINT MODELLING

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(JMbayes2))
suppressPackageStartupMessages(library(nlme))

set.seed(2222)

########################################################################### PREP

args <- commandArgs(trailingOnly = TRUE)
serum <- read.csv(args[1])
demog <- read.csv(args[2])
protein <- args[3]

long <- merge(serum,demog,by.x="SampleID",by.y="SERUM_OLINK_MANIFEST",all.x=T) %>%
  filter(GROUP=="AMYOTROPHIC LATERAL SCLEROSIS") %>%
  group_by(IMCM_ID) %>%
  arrange(IMCM_ID,AGE_AT_SAMPLING) %>%
  mutate(Time=AGE_AT_DEAD_OR_CENSORED-first(AGE_AT_SAMPLING)) %>%
  mutate(death=case_when(DEAD_OR_CENSORED=="DEAD"~1,DEAD_OR_CENSORED=="CENSORED"~0,TRUE~NA)) %>%
  mutate(DELTA=AGE_AT_SAMPLING-first(AGE_AT_SAMPLING)) %>%
  select(c(1,3:5418,IMCM_ID,AGE_AT_SAMPLING,Time,death,DELTA))

surv <- long %>%
  group_by(IMCM_ID) %>%
  arrange(IMCM_ID,AGE_AT_SAMPLING) %>%
  filter(row_number() == 1)

################################################################ EXAMPLE PROTEIN

res <- tryCatch(
  { 
    survFit <- survival::coxph(Surv(Time, death) ~ 1, data = surv, x = TRUE)
    lme_protein <- nlme::lme(as.formula(paste(protein, "~ DELTA")), 
                             random = ~ DELTA | IMCM_ID, 
                             data = long,
                             control = lmeControl(opt = 'optim'))
    joint_protein <- JMbayes2::jm(survFit, lme_protein, time_var = "DELTA",
                                  id_var="IMCM_ID", n_iter = 10000L, n_burnin = 1000L)
    joint_protein_results <- data.frame(row.names = NULL,
                                        Protein = protein,
                                        marginal_DIC = joint_protein$fit_stats$marginal$DIC,
                                        marginal_WAIC = joint_protein$fit_stats$marginal$WAIC,
                                        marginal_LPML = joint_protein$fit_stats$marginal$LPML,
                                        Alphas_Mean = joint_protein$statistics$Mean$alphas,
                                        Alphas_StDev = joint_protein$statistics$SD$alphas,
                                        Alphas_P = joint_protein$statistics$P$alphas,
                                        Rhat = joint_protein$statistics$Rhat$alphas[,"Point est."])
    write.csv(joint_protein_results,paste0(protein,"_joint_results.csv"),row.names=FALSE)
    saveRDS(joint_protein,paste0(protein,"_joint_results.RDS"),compress = "gzip")
  },
  error=function(err) {
    empty_obj <- list()
    file.create(paste0(protein,"_joint_results.csv"))
    saveRDS(empty_obj,paste0(protein,"_joint_results.RDS"),compress = "gzip")
  })
