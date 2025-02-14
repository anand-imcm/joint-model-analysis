#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-# WORKFLOW FOR SERUM JOINT MODELLING
print(paste0("START script ", format(Sys.time(), "%H:%M:%S"), " on ", Sys.Date()))
required_packages <- c("dplyr", "JMbayes2", "survival", "nlme")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(missing_packages)) {
  install.packages(missing_packages, dependencies = TRUE)
}
suppressPackageStartupMessages(lapply(required_packages, require, character.only = TRUE))
set.seed(2222)
########################################################################### PREP
print(paste0("START loading inputs ", format(Sys.time(), "%H:%M:%S"), " on ", Sys.Date()))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  if (interactive()) {
    message("Using default test file paths.")
    serum <- read.csv("test/serum_test.csv")
    demog <- read.csv("test/demographics_test.csv")
    protein <- "A"
  } else {
    stop("Insufficient arguments supplied. Please provide: serum(csv file), demographics(csv file), and protein name")
  }
} else {
  message("Using commandline arguments.")
  serum <- read.csv(args[1])
  demog <- read.csv(args[2])
  protein <- args[3]
}
print(paste0("END loading inputs ", format(Sys.time(), "%H:%M:%S"), " on ", Sys.Date()))
print(paste0("START merging data ", format(Sys.time(), "%H:%M:%S"), " on ", Sys.Date()))
long <- merge(serum,demog,by.x="SampleID",by.y="SERUM_OLINK_MANIFEST",all.x=T) %>% 
  filter(GROUP=="Group_2") %>% 
  group_by(IMCM_ID) %>%
  arrange(IMCM_ID,AGE_AT_SAMPLING) %>%
  mutate(Time=AGE_AT_A_OR_B-first(AGE_AT_SAMPLING)) %>%
  mutate(death=case_when(A_OR_B=="A"~1,A_OR_B=="B"~0,TRUE~NA)) %>%
  mutate(DELTA=AGE_AT_SAMPLING-first(AGE_AT_SAMPLING)) %>%
  select(c(1,3:5418,IMCM_ID,AGE_AT_SAMPLING,Time,death,DELTA))
print(paste0("END merging data ", format(Sys.time(), "%H:%M:%S"), " on ", Sys.Date()))
print(paste0("START filter ", format(Sys.time(), "%H:%M:%S"), " on ", Sys.Date()))
surv <- long %>%
  group_by(IMCM_ID) %>%
  arrange(IMCM_ID,AGE_AT_SAMPLING) %>%
  filter(row_number() == 1)
print(paste0("END filter ", format(Sys.time(), "%H:%M:%S"), " on ", Sys.Date()))
################################################################ EXAMPLE PROTEIN
print(paste0("START survival::coxph ", format(Sys.time(), "%H:%M:%S"), " on ", Sys.Date()))
survFit <- survival::coxph(Surv(Time, death) ~ 1, data = surv, x = TRUE)
print(paste0("END survival::coxph ", format(Sys.time(), "%H:%M:%S"), " on ", Sys.Date()))
print(paste0("START nlme::lme ", format(Sys.time(), "%H:%M:%S"), " on ", Sys.Date()))
lme_protein <- nlme::lme(as.formula(paste(protein, "~ DELTA")), 
                         random = ~ DELTA | IMCM_ID, 
                         data = long,
                         control = lmeControl(opt = 'optim'))
print(paste0("END nlme::lme ", format(Sys.time(), "%H:%M:%S"), " on ", Sys.Date()))
print(paste0("START JMbayes2::jm ", format(Sys.time(), "%H:%M:%S"), " on ", Sys.Date()))
joint_protein <- JMbayes2::jm(survFit, lme_protein, time_var = "DELTA",
                           id_var="IMCM_ID", n_iter = 10000L, n_burnin = 1000L)
print(paste0("END JMbayes2::jm ", format(Sys.time(), "%H:%M:%S"), " on ", Sys.Date()))
print(paste0("START joint_protein_results ", format(Sys.time(), "%H:%M:%S"), " on ", Sys.Date()))
joint_protein_results <- data.frame(row.names = NULL,
                                    Protein = protein,
                                    marginal_DIC = joint_protein$fit_stats$marginal$DIC,
                                    marginal_WAIC = joint_protein$fit_stats$marginal$WAIC,
                                    marginal_LPML = joint_protein$fit_stats$marginal$LPML,
                                    Alphas_Mean = joint_protein$statistics$Mean$alphas,
                                    Alphas_StDev = joint_protein$statistics$SD$alphas,
                                    Alphas_P = joint_protein$statistics$P$alphas,
                                    Rhat = joint_protein$statistics$Rhat$alphas[,"Point est."])
print(paste0("END joint_protein_results ", format(Sys.time(), "%H:%M:%S"), " on ", Sys.Date()))
print(paste0("START write.csv ", format(Sys.time(), "%H:%M:%S"), " on ", Sys.Date()))
write.csv(joint_protein_results,paste0(protein,"_joint_results.csv"),row.names=FALSE)
print(paste0("END write.csv ", format(Sys.time(), "%H:%M:%S"), " on ", Sys.Date()))
print(paste0("START saveRDS ", format(Sys.time(), "%H:%M:%S"), " on ", Sys.Date()))
saveRDS(joint_protein,paste0(protein,"_joint_results.RDS"),compress = "gzip")
print(paste0("END saveRDS ", format(Sys.time(), "%H:%M:%S"), " on ", Sys.Date()))
print(paste0("END script ", format(Sys.time(), "%H:%M:%S"), " on ", Sys.Date()))
