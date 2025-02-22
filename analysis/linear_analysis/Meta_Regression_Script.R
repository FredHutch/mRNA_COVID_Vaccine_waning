################################################################################
######################Meta Regression Script####################################
################################################################################

output_folder = ifelse(
  .Platform$OS.type == "windows",
  paste0("data/Meta_Analysis_Results/Results_", Sys.Date()),
  paste0("/fh/fast/schiffer_j/COVIDmodeling/Meta_Analysis_Results/Results_", Sys.Date())
)

library(R.utils)
library(purrr)
args = commandArgs(trailingOnly=TRUE, asValues = TRUE, 
                   defaults = c(MULTI_LEVEL_FLAG = TRUE,
                                AGE_ADJUSTMENT = FALSE,
                                BY_STUDY_TYPE = FALSE))

args$input_data = "X:\\fast\\schiffer_j\\COVIDmodeling\\data\\Waning_Table_02_21_25.rds"
args$AGE_ADJUSTMENT = TRUE
args$BY_STUDY_TYPE = FALSE

if(args$BY_STUDY_TYPE){
output_folder = ifelse(
  .Platform$OS.type == "windows",
  paste0("data/Meta_Analysis_Results/Results_BY_TYPE_", Sys.Date()),
  paste0("/fh/fast/schiffer_j/COVIDmodeling/Meta_Analysis_Results/Results_BY_TYPE_", Sys.Date())
)
}


source("analysis/linear_analysis/Meta_Regression_Functions.R")
if(args$AGE_ADJUSTMENT){
  source("analysis/linear_analysis/Meta_Regression_Functions_Age.R")
  output_folder = paste0(output_folder, "_age_adjust")  
}

library(data.table)

MULTI_LEVEL_FLAG = args$MULTI_LEVEL_FLAG
BY_STUDY_TYPE = args$BY_STUDY_TYPE
LOG_TRUNCATION = log(1 - 99/100)


#Waning_Table_Reduced = readRDS("data/Waning_Table_01_18_25.rds")

Waning_Table_Reduced = data.table(readRDS(args$input_data))
if(!BY_STUDY_TYPE){
  Waning_Table_Reduced[,STUDY_TYPE_SIMPLE:="All"]
}

Waning_Table_Reduced = Waning_Table_Reduced[!is.na(STUDY_TYPE_SIMPLE)]



dir.create(output_folder)
saveRDS(Waning_Table_Reduced, paste0(output_folder, "/Meta_Regression_Data.rds"))
#source("Parameter_Estimation.R")
Waning_Table_Reduced[,IMM := factor(IMM, levels = c("V", "I", "H", "B", "HB", "BvV", "BB", "BvB", "HvI", "HBvI", "Hvf", "HI", "HBvH")),]
Waning_Table_Reduced[,STRAIN := factor(STRAIN, levels = c("Pre-Omicron", "Omicron")),]
Waning_Table_Reduced[,PSTRAIN := factor(PSTRAIN, levels = c("Omicron", "Pre-Omicron")),]

Waning_Table_Reduced[
  TYPE!="SD" & IMM%in%c("V", "B", "I", "H", "HB") & !is.na(WEIGHT),
  .(studies = 
      paste0(length(unique(STUDY)), " (",.N, ")")),
  .(IMM, 
    TYPE, 
    PSTRAIN, 
    STRAIN)
]->STUDY_SUMMARY
STUDY_SUMMARY[is.na(PSTRAIN), PSTRAIN:=""]
STUDY_SUMMARY[,TYPE:=factor(TYPE, levels = c("S", "D", "M")),]
STUDY_SUMMARY[,IMM:=factor(IMM, levels = c("V", "B", "I", "H", "HB")),]
STUDY_SUMMARY[,PSTRAIN:=factor(PSTRAIN, levels = c("","Pre-Omicron", "Omicron")),]
STUDY_SUMMARY[,STRAIN:=factor(STRAIN, levels = c("Pre-Omicron", "Omicron")),]

STUDY_SUMMARY[TYPE %in% c("S", "D", "M") & 
                IMM %in% c("V", "B", "H", "HB"),
              strain_model:=T,]
STUDY_SUMMARY[TYPE %in% c("S", "D", "M") & 
                IMM %in% c("V", "B"),
              booster_model:=T,]
STUDY_SUMMARY[TYPE %in% c("S", "D") & 
                IMM %in% c("B", "HB"),
              hybrid_boosted_model:=T,]
STUDY_SUMMARY[TYPE %in% c("S", "D") & 
                IMM %in% c("I", "V", "H"),
              infected_model:=T,]
STUDY_SUMMARY[TYPE %in% c("S", "D") & 
                IMM %in% c("I", "V", "H"),
              pstrain_model:=T,]

saveRDS(STUDY_SUMMARY, paste0(output_folder, "/STUDY_SUMMARY.rds"))

# Model 1 ####
RFIG_STRAIN_DATA = 
  Waning_Table_Reduced[TYPE %in% c("S", "D", "M") &
                         IMM %in% c("V", "B", "H", "HB") &
                         !is.na(LOWER) & LOWER!=UPPER,
                       fit_rma_strain_wrap(.SD, multilevel = MULTI_LEVEL_FLAG),
                       .(TYPE, STUDY_TYPE_SIMPLE)]

Waning_Table_Reduced[TYPE %in% c("S", "D", "M") &
                       IMM %in% c("V", "B", "H", "HB") & !is.na(STUDY_TYPE_SIMPLE) &
                       !is.na(LOWER) & LOWER!=UPPER,.(M = length(unique(STUDY)),.N),
                     STUDY_TYPE_SIMPLE]->Model_1


saveRDS(RFIG_STRAIN_DATA, paste0(output_folder, "/RFIG_STRAIN_DATA.rds"))

## Log Version ####

RFIG_STRAIN_DATA_LOG = 
  Waning_Table_Reduced[TYPE %in% c("S", "D", "M") &
                         IMM %in% c("V", "B", "H", "HB") &
                         !is.na(LOWER) & LOWER!=UPPER,
                       .(VALUE = ifelse(VALUE<100, LHR, LOG_TRUNCATION),
                         TYPE,
                         STUDY_TYPE_SIMPLE,
                         STRAIN,
                         AGE_GROUP,
                         STUDY,
                         effect_id,
                         MIN_TIME,
                         LOWER = ifelse(UPPER<100, LHR_LO, LOG_TRUNCATION), 
                         UPPER = LHR_HI,
                         STDERROR = ifelse(UPPER<100,
                                           (LHR_HI - LHR_LO)/2,
                                           ifelse(VALUE<100,
                                           LHR_HI - LHR,
                                           LHR_HI - LOG_TRUNCATION))/qnorm(.975))
  ][,
    fit_rma_strain_wrap(.SD, multilevel = MULTI_LEVEL_FLAG),
    .(TYPE, STUDY_TYPE_SIMPLE)]



Waning_Table_Reduced[TYPE %in% c("S", "D", "M") &
                       IMM %in% c("V", "B", "H", "HB") &
                       !is.na(LOWER) & LOWER!=UPPER & UPPER<100,
                     .(M = length(unique(STUDY)),.N),STUDY_TYPE_SIMPLE
                     ]->Model_1_log


saveRDS(RFIG_STRAIN_DATA_LOG, paste0(output_folder, "/RFIG_STRAIN_DATA_LOG.rds"))


# Model 2 ####
RFIG_BOOSTER_DATA = 
  Waning_Table_Reduced[TYPE %in% c("S", "M", "D") &
                         IMM %in% c("V", "B", "H", "HB") &
                         !is.na(LOWER) & LOWER!=UPPER,
                       fit_rma_booster_wrap(.SD, multilevel = MULTI_LEVEL_FLAG),
                       .(TYPE, STRAIN, STUDY_TYPE_SIMPLE)]

saveRDS(RFIG_BOOSTER_DATA, paste0(output_folder, "//RFIG_BOOSTER_DATA.rds"))

## Log Version ####
RFIG_BOOSTER_DATA_LOG = 
  Waning_Table_Reduced[TYPE %in% c("S", "D", "M") &
                         IMM %in% c("V", "B", "H", "HB") &
                         !is.na(LOWER) & LOWER!=UPPER,
                       .(VALUE = ifelse(VALUE<100, LHR, LOG_TRUNCATION),
                         TYPE,
                         STUDY_TYPE_SIMPLE,
                         STRAIN,
                         AGE_GROUP,
                         Boosted,
                         STUDY,
                         effect_id,
                         MIN_TIME,
                         LOWER = ifelse(UPPER<100, LHR_LO, LOG_TRUNCATION), 
                         UPPER = LHR_HI,
                         STDERROR = ifelse(UPPER<100,
                                           (LHR_HI - LHR_LO)/2,
                                           ifelse(VALUE<100,
                                                  LHR_HI - LHR,
                                                  LHR_HI - LOG_TRUNCATION))/qnorm(.975))
  ][,
    fit_rma_booster_wrap(.SD, multilevel = MULTI_LEVEL_FLAG),
    .(TYPE, STRAIN, STUDY_TYPE_SIMPLE)]

saveRDS(RFIG_BOOSTER_DATA_LOG, paste0(output_folder, "//RFIG_BOOSTER_DATA_LOG.rds"))


Waning_Table_Reduced[TYPE %in% c("S", "D", "M") &
                       IMM %in% c("V", "B") &
                       !is.na(LOWER) & LOWER!=UPPER,.(M = length(unique(STUDY)),.N),
                     STUDY_TYPE_SIMPLE]->Model_2

# Model 2 alternative version ####

RFIG_BOOSTER_DATA_mod = 
  Waning_Table_Reduced[TYPE %in% c("S", "M", "D") &
                         IMM %in% c("V", "B") &
                         !is.na(LOWER) & LOWER!=UPPER,
                       fit_rma_booster_wrap(.SD, multilevel = MULTI_LEVEL_FLAG),
                       .(TYPE, STRAIN, STUDY_TYPE_SIMPLE)]


saveRDS(RFIG_BOOSTER_DATA_mod, paste0(output_folder, "/RFIG_BOOSTER_DATA_mod.rds"))

## Log Version ####
RFIG_BOOSTER_DATA_LOG_mod = 
  Waning_Table_Reduced[TYPE %in% c("S", "M", "D") &
                         IMM %in% c("V", "B") &
                         !is.na(LOWER) & LOWER!=UPPER,
                       .(VALUE = ifelse(VALUE<100, LHR, LOG_TRUNCATION),
                         TYPE,
                         STUDY_TYPE_SIMPLE,
                         STRAIN,
                         AGE_GROUP,
                         Boosted,
                         STUDY,
                         effect_id,
                         MIN_TIME,
                         LOWER = ifelse(UPPER<100, LHR_LO, LOG_TRUNCATION), 
                         UPPER = LHR_HI,
                         STDERROR = ifelse(UPPER<100,
                                           (LHR_HI - LHR_LO)/2,
                                           ifelse(VALUE<100,
                                                  LHR_HI - LHR,
                                                  LHR_HI - LOG_TRUNCATION))/qnorm(.975))
  ][,
    fit_rma_booster_wrap(.SD, multilevel = MULTI_LEVEL_FLAG),
    .(TYPE, STRAIN, STUDY_TYPE_SIMPLE)]


saveRDS(RFIG_BOOSTER_DATA_LOG_mod, paste0(output_folder, "/RFIG_BOOSTER_DATA_LOG_mod.rds"))


Waning_Table_Reduced[TYPE %in% c("S", "D", "M") &
                       IMM %in% c("V", "B") &
                       !is.na(LOWER) & LWEIGHT!= 0 & is.finite(LWEIGHT) & VALUE<100,
                     .(M = length(unique(STUDY)),.N),STUDY_TYPE_SIMPLE]->Model_2_log

# Model 3 ####
RFIG_INF_DATA = 
  Waning_Table_Reduced[TYPE %in% c("S", "D") & 
                         !(TYPE != "S" & STRAIN == "Pre-Omicron") &
                         IMM %in% c("I", "V", "H") &
                         !is.na(LOWER) & LOWER!=UPPER,
                       fit_rma_inf_wrap(.SD, multilevel = MULTI_LEVEL_FLAG),
                       .(TYPE, STRAIN, STUDY_TYPE_SIMPLE)]



Waning_Table_Reduced[TYPE %in% c("D") &
                       IMM %in% c("I", "V") & STRAIN == "Pre-Omicron" &
                       !is.na(LOWER) & LOWER!=UPPER & STUDY_TYPE_SIMPLE != "Cohort",
                     fit_rma_inf_alt_wrap(.SD, multilevel = MULTI_LEVEL_FLAG),
                     .(TYPE, STRAIN, STUDY_TYPE_SIMPLE)]->.out

.out[,TYPE:="D",]
.out[,STRAIN:="Pre-Omicron",]

RFIG_INF_DATA = rbind(RFIG_INF_DATA, .out)

saveRDS(RFIG_INF_DATA, paste0(output_folder, "/RFIG_INF_DATA.rds"))

Waning_Table_Reduced[TYPE %in% c("S", "D") &
                       IMM %in% c("I", "V", "H") &
                       !is.na(LOWER) & LOWER!=UPPER,.(M = length(unique(STUDY)),.N), 
                     STUDY_TYPE_SIMPLE]->Model_3

## Log Version ####
RFIG_INF_DATA_LOG = 
  Waning_Table_Reduced[TYPE %in% c("S", "D") & 
                         !(TYPE != "S" & STRAIN == "Pre-Omicron") &
                         IMM %in% c("I", "V", "H") &
                         !is.na(LOWER) & LOWER!=UPPER,
                       .(VALUE = ifelse(VALUE<100, LHR, LOG_TRUNCATION),
                         TYPE,
                         STUDY_TYPE_SIMPLE,
                         STRAIN,
                         AGE_GROUP,
                         IMM,
                         STUDY,
                         effect_id,
                         MIN_TIME,
                         LOWER = ifelse(UPPER<100, LHR_LO, LOG_TRUNCATION), 
                         UPPER = LHR_HI,
                         STDERROR = ifelse(UPPER<100,
                                           (LHR_HI - LHR_LO)/2,
                                           ifelse(VALUE<100,
                                                  LHR_HI - LHR,
                                                  LHR_HI - LOG_TRUNCATION))/qnorm(.975))
  ][,
    fit_rma_inf_wrap(.SD, multilevel = MULTI_LEVEL_FLAG),
    .(TYPE, STRAIN, STUDY_TYPE_SIMPLE)]



Waning_Table_Reduced[TYPE %in% c("D") &
                       IMM %in% c("I", "V") & STRAIN == "Pre-Omicron" &
                       !is.na(LOWER) & LOWER!=UPPER & STUDY_TYPE_SIMPLE != "Cohort",
                     .(VALUE = ifelse(VALUE<100, LHR, LOG_TRUNCATION),
                       TYPE, 
                       STUDY_TYPE_SIMPLE,
                       STRAIN,
                       AGE_GROUP,
                       IMM,
                       STUDY,
                       effect_id,
                       MIN_TIME,
                       LOWER = ifelse(UPPER<100, LHR_LO, LOG_TRUNCATION), 
                       UPPER = LHR_HI,
                       STDERROR = ifelse(UPPER<100,
                                         (LHR_HI - LHR_LO)/2,
                                         ifelse(VALUE<100,
                                                LHR_HI - LHR,
                                                LHR_HI - LOG_TRUNCATION))/qnorm(.975))
][,
                     fit_rma_inf_alt_wrap(.SD, multilevel = MULTI_LEVEL_FLAG),
                     .(TYPE, STRAIN, STUDY_TYPE_SIMPLE)]->.out

.out[,TYPE:="D",]
.out[,STRAIN:="Pre-Omicron",]

RFIG_INF_DATA_LOG = rbind(RFIG_INF_DATA_LOG, .out)

saveRDS(RFIG_INF_DATA_LOG, paste0(output_folder, "/RFIG_INF_DATA_LOG.rds"))

Waning_Table_Reduced[TYPE %in% c("S", "D") &
                       IMM %in% c("I", "V", "H") &
                       !is.na(LOWER) & LWEIGHT!= 0 & is.finite(LWEIGHT) & VALUE<100,
                     .(M = length(unique(STUDY)),.N), STUDY_TYPE_SIMPLE]->Model_3_log


# Model 4 ####
Waning_Table_Reduced[TYPE %in% c("S", "D") & !(TYPE != "S" & STRAIN == "Pre-Omicron") &
                       IMM %in% c("I", "V", "H") &
                       !is.na(LOWER) & LOWER!=UPPER & STUDY_TYPE_SIMPLE!="Cohort"]->.tmp

.tmp[!is.na(PSTRAIN) & PSTRAIN == "Pre-Omicron", IMMP:=paste0(IMM,PSTRAIN),]
.tmp[PSTRAIN %in% c("Omicron", "Mix"), IMMP:=paste0(IMM,"Post-Omicron"),]
.tmp[is.na(PSTRAIN) & IMM == "V", IMMP:="V",]

RFIG_INF_DATA_pstrain = .tmp[TYPE %in% c("S", "D") & STRAIN == "Omicron" &
                               IMM %in% c("I", "V", "H") &
                               !is.na(LOWER) & LOWER!=UPPER
                             ,fit_rma_pstrain_wrap(.SD, multilevel = MULTI_LEVEL_FLAG),
                             .(TYPE, STRAIN, STUDY_TYPE_SIMPLE)]


Waning_Table_Reduced[TYPE %in% c("S", "D") & STRAIN == "Omicron" &
                       IMM %in% c("I", "V", "H") &
                       !is.na(LOWER) & LOWER!=UPPER,.(M = length(unique(STUDY)),.N), 
                     STUDY_TYPE_SIMPLE]->Model_4


RFIG_INF_DATA_pstrain = merge(RFIG_INF_DATA_pstrain,
                              data.table(
                                IMMP = c("V", "V", "HPre-Omicron", "HPost-Omicron","IPre-Omicron", "IPost-Omicron"),
                                PSTRAIN = rep(c("Pre-Omicron", "Post-Omicron"), 3),
                                IMM = rep(c("V", "H", "I"), each = 2)
                              ), all = T, allow.cartesian = T)


saveRDS(RFIG_INF_DATA_pstrain, paste0(output_folder, "/RFIG_INF_DATA_pstrain.rds"))

## Log Version ####

RFIG_INF_DATA_pstrain_LOG = .tmp[TYPE %in% c("S", "D") & STRAIN == "Omicron" &
                               IMM %in% c("I", "V", "H") &
                                 !is.na(LOWER) & LOWER!=UPPER,
                               .(VALUE = ifelse(VALUE<100, LHR, LOG_TRUNCATION),
                                 TYPE,
                                 STUDY_TYPE_SIMPLE,
                                 STRAIN,
                                 AGE_GROUP,
                                 IMMP,
                                 STUDY,
                                 effect_id,
                                 MIN_TIME,
                                 LOWER = ifelse(UPPER<100, LHR_LO, LOG_TRUNCATION), 
                                 UPPER = LHR_HI,
                                 STDERROR = ifelse(UPPER<100,
                                                   (LHR_HI - LHR_LO)/2,
                                                   ifelse(VALUE<100,
                                                          LHR_HI - LHR,
                                                          LHR_HI - LOG_TRUNCATION))/qnorm(.975))
][,fit_rma_pstrain_wrap(.SD, multilevel = MULTI_LEVEL_FLAG),
                             .(TYPE, STRAIN, STUDY_TYPE_SIMPLE)]


Waning_Table_Reduced[TYPE %in% c("S", "D") & STRAIN == "Omicron" &
                       IMM %in% c("I", "V", "H") &
                       !is.na(LOWER) & LWEIGHT!= 0 & is.finite(LWEIGHT) & VALUE<100,
                     .(M = length(unique(STUDY)),.N),STUDY_TYPE_SIMPLE]->Model_4_log


RFIG_INF_DATA_pstrain_LOG = merge(RFIG_INF_DATA_pstrain_LOG,
                              data.table(
                                IMMP = c("V", "V", "HPre-Omicron", "HPost-Omicron","IPre-Omicron", "IPost-Omicron"),
                                PSTRAIN = rep(c("Pre-Omicron", "Post-Omicron"), 3),
                                IMM = rep(c("V", "H", "I"), each = 2)
                              ), all = T, allow.cartesian = T)


saveRDS(RFIG_INF_DATA_pstrain_LOG, paste0(output_folder, "/RFIG_INF_DATA_pstrain_LOG.rds"))

# Model 5 ####

RFIG_HYBRID_BOOSTED_DATA = 
  Waning_Table_Reduced[TYPE %in% c("S", "D") &
                         IMM %in% c("B", "HB") & STRAIN == "Omicron" &
                         !is.na(LOWER) & LOWER!=UPPER & STUDY_TYPE_SIMPLE != "Cohort",
                       fit_rma_hboosted_wrap(.SD, multilevel = MULTI_LEVEL_FLAG),
                       .(TYPE, STRAIN, STUDY_TYPE_SIMPLE)]



Waning_Table_Reduced[TYPE %in% c("S", "D") & STRAIN == "Omicron" &
                       IMM %in% c("B", "HB") & STUDY_TYPE_SIMPLE != "Cohort" &
                       !is.na(LOWER) & LOWER!=UPPER,.(M = length(unique(STUDY)),.N),
                     STUDY_TYPE_SIMPLE]->Model_5

saveRDS(RFIG_HYBRID_BOOSTED_DATA, paste0(output_folder, "/RFIG_HYBRID_BOOSTED_DATA.rds"))

## Log Version ####

RFIG_HYBRID_BOOSTED_DATA_LOG = 
  Waning_Table_Reduced[TYPE %in% c("S", "D") &
                         IMM %in% c("B", "HB") & STRAIN == "Omicron" &
                         !is.na(LOWER) & LOWER!=UPPER & STUDY_TYPE_SIMPLE != "Cohort",
                       .(VALUE = ifelse(VALUE<100, LHR, LOG_TRUNCATION),
                         TYPE,
                         STUDY_TYPE_SIMPLE,
                         STRAIN,
                         AGE_GROUP,
                         Infected,
                         STUDY,
                         effect_id,
                         MIN_TIME,
                         LOWER = ifelse(UPPER<100, LHR_LO, LOG_TRUNCATION), 
                         UPPER = LHR_HI,
                         STDERROR = ifelse(UPPER<100,
                                           (LHR_HI - LHR_LO)/2,
                                           ifelse(VALUE<100,
                                                  LHR_HI - LHR,
                                                  LHR_HI - LOG_TRUNCATION))/qnorm(.975))
  ][,fit_rma_hboosted_wrap(.SD, multilevel = MULTI_LEVEL_FLAG),
                       .(TYPE, STRAIN, STUDY_TYPE_SIMPLE)]


Waning_Table_Reduced[TYPE %in% c("S", "D") & STRAIN == "Omicron" &
                       IMM %in% c("B", "HB") &
                       !is.na(LOWER) & LOWER!=UPPER& STUDY_TYPE_SIMPLE != "Cohort",.(M = length(unique(STUDY)),.N),
                     STUDY_TYPE_SIMPLE]->Model_5_log

saveRDS(RFIG_HYBRID_BOOSTED_DATA_LOG, paste0(output_folder, "/RFIG_HYBRID_BOOSTED_DATA_LOG.rds"))

saveRDS(rbindlist(list(Model_1, Model_2, Model_3, Model_4, Model_5), idcol = "Model"),
        paste0(output_folder, "/MODEL_SUMMARY.rds"))


saveRDS(rbindlist(list(Model_1_log, Model_2_log, Model_3_log, Model_4_log, Model_5_log), idcol = "Model"),
        paste0(output_folder, "/MODEL_SUMMARY_LOG.rds"))
