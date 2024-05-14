################################################################################
######################Meta Regression Script####################################
################################################################################
MULTI_LEVEL_FLAG = T
library(data.table)
library(metafor)
Waning_Table_Reduced = readRDS("Waning_Table_1_30.rds")
output_folder = paste0("Meta_Analysis_Results/Results_", Sys.Date())

dir.create(output_folder)
saveRDS(Waning_Table_Reduced, file = paste0(output_folder, "/Meta_Regression_Data.rds"))
#source("Parameter_Estimation.R")
Waning_Table_Reduced[,IMM := factor(IMM, levels = c("V", "I", "H", "B", "HB", "BvV", 
                                                    "BB", "BvB", "HvI", "HBvI", "Hvf", "HI", "HBvH")),]
Waning_Table_Reduced[,STRAIN := factor(STRAIN, levels = c("Pre-Omicron", "Omicron")),]
Waning_Table_Reduced[,PSTRAIN := factor(PSTRAIN, levels = c("Omicron", "Pre-Omicron", "None")),]


NONLINEAR_CURVEFITS = 
  Waning_Table_Reduced[TYPE %in% c("S", "M", "D", "SD") &
                         IMM %in% c("V", "B", "I", "H", "HB") &
                         PSTRAIN != "Omicron" & !is.na(LOWER) & LOWER!=UPPER,
                       cbind(
                         MODEL = "Main Fit",
                         fit_rma_nonlinear_wrap(
                           .SD, multilevel = MULTI_LEVEL_FLAG)
                         ),
                       .(TYPE, STRAIN)]


NONLINEAR_CURVEFITS_REDUCED = 
  Waning_Table_Reduced[TYPE %in% c("S", "M", "D", "SD") &
                         IMM %in% c("V", "B", "HI", "HB") &
                         PSTRAIN != "Omicron" & !is.na(LOWER) & LOWER!=UPPER,
                       cbind(
                         MODEL = "Combine Hybrid and Natural",
                         fit_rma_nonlinear_wrap(
                           .SD, multilevel = MULTI_LEVEL_FLAG)
                         ),
                       .(TYPE, STRAIN)]


NONLINEAR_CURVEFITS_MATCHED = 
  Waning_Table_Reduced[TYPE %in% c("S", "M", "D", "SD") &
                         IMM %in% c("V", "B", "H", "I", "HB") &
                         ((STRAIN == "Omicron" & PSTRAIN== "Omicron") |
                            (STRAIN == "Pre-Omicron" & PSTRAIN %in% c("Pre-Omicron", "None"))) &
                         !is.na(LOWER) & LOWER!=UPPER,
                       cbind(
                         MODEL = "Use Omicron v Omicron data",
                         fit_rma_nonlinear_wrap(.SD, multilevel = MULTI_LEVEL_FLAG)
                         ),
                       .(TYPE)]

NONLINEAR_CURVEFITS_MATCHED_REDUCED = 
  Waning_Table_Reduced[TYPE %in% c("S", "M", "D", "SD") &
                         IMM %in% c("V", "B", "HI", "HB") &
                         ((STRAIN == "Omicron" & PSTRAIN== "Omicron") |
                            (STRAIN == "Pre-Omicron" & PSTRAIN %in% c("Pre-Omicron", "None"))) &
                         !is.na(LOWER) & LOWER!=UPPER,
                       cbind(
                         MODEL = "Use Omicron v Omicron data and combine Hybrid and Natural",
                         fit_rma_nonlinear_wrap(.SD, multilevel = MULTI_LEVEL_FLAG)
                         ),
                       .(TYPE)]
Non_Linear_Fit = rbind(NONLINEAR_CURVEFITS, 
      NONLINEAR_CURVEFITS_MATCHED, 
      NONLINEAR_CURVEFITS_MATCHED_REDUCED, 
      NONLINEAR_CURVEFITS_REDUCED, 
      fill = T)



saveRDS(Non_Linear_Fit, file = paste0(output_folder, "/Non_Linear_Fit.rds"))

setnames(Non_Linear_Fit, old = c("se", "STRAIN", "IMM"), 
         new = c("STDERROR", "STRAINmodel", "IMMmodel"), skip_absent = T)
Non_Linear_Fit[is.na(STRAINmodel),STRAINmodel:="Matched",]


# Efficacy measures not compared to Naive ####
Waning_Table_Reduced[TYPE %in% c("S", "M", "D", "SD") &
                       !(IMM %in% c("V", "B", "HI", "HB", "I", "H")) &
                       !is.na(LOWER) & LOWER!=UPPER
                     ]->OUT_OF_SAMPLE_OBS



OUT_OF_SAMPLE_OBS[IMM == "BvB",IMM_BASE:="B",]
OUT_OF_SAMPLE_OBS[IMM %in% c("HBvV", "BvV"), IMM_BASE:="V",]
OUT_OF_SAMPLE_OBS[IMM %in% c("HvI", "HBvI"), IMM_BASE:="I",]
OUT_OF_SAMPLE_OBS[IMM %in% c("BvB", "BvV"), IMM_MAIN:="B",]
OUT_OF_SAMPLE_OBS[IMM %in% c("HBvV", "HBvI"), IMM_MAIN:="HB",]
OUT_OF_SAMPLE_OBS[IMM %in% c("HvI", "Hvf"), IMM_MAIN:="H",]


matching.table = data.table(
  STRAIN = rep(c("Pre-Omicron", "Omicron"), each = 4),
  PSTRAIN = c("None", "None", "Pre-Omicron", "Pre-Omicron",
              "Pre-Omicron", "None", "Omicron", "Omicron"),
  STRAINmodel = c("Pre-Omicron", "Matched", "Pre-Omicron", "Matched",
                  "Omicron", "Omicron", "Pre-Omicron", "Matched")
)


matching.table2 = data.table(
  IMMmodel = c("V", "I", "H", "HI", "HI", "HB", "B"),
  IMM = c("V", "I", "H", "H", "I", "HB", "B")
)


Non_Linear_Fit[
  matching.table, 
  on = .(STRAINmodel), allow.cartesian = T
][
  matching.table2,
  on = .(IMMmodel), allow.cartesian = T
]->Non_Linear_Fit_expand



merge(merge(merge(OUT_OF_SAMPLE_OBS, Non_Linear_Fit_expand[WANING == 1], 
      by.x = c("TYPE", "PSTRAIN", "STRAIN", "IMM_BASE"), 
      by.y = c("TYPE", "PSTRAIN", "STRAIN", "IMM"), 
      allow.cartesian = T, suffixes = c(".data", ".base")),
      Non_Linear_Fit_expand[WANING == 1],
      by.x = c("TYPE", "PSTRAIN", "STRAIN", "IMM_MAIN", "MODEL", "method", "STRAINmodel", "OUTPUT"), 
      by.y = c("TYPE", "PSTRAIN", "STRAIN", "IMM", "MODEL", "method", "STRAINmodel", "OUTPUT"),
      suffixes = c("", ".waned")
      ),
      Non_Linear_Fit_expand[WANING == 0],
      by.x = c("TYPE", "PSTRAIN", "STRAIN", "IMM_MAIN", "MODEL", "method", "STRAINmodel", "OUTPUT"), 
      by.y = c("TYPE", "PSTRAIN", "STRAIN", "IMM", "MODEL", "method", "STRAINmodel", "OUTPUT"),
      suffixes = c("", ".initial")
)->OUT_OF_SAMPLE_OBS




OUT_OF_SAMPLE_OBS[,VALUE.predict.model:=VALUE.initial * exp(-decay_rate * MIN_TIME) + 
                    VALUE * (1 - exp(-decay_rate * MIN_TIME)),]


OUT_OF_SAMPLE_OBS[,STDERROR.predict.model:=sqrt(STDERROR.initial^2 * exp(-2 *decay_rate * MIN_TIME)+ 
                                            STDERROR^2 * (1 - exp(-decay_rate * MIN_TIME))^2),]
OUT_OF_SAMPLE_OBS[,VALUE.predict.data:=100 - (1 - VALUE.data/100) * (100 - VALUE.base),]
OUT_OF_SAMPLE_OBS[,STDERROR.predict.data:=sqrt((STDERROR.data/(100 - VALUE.data))^2  + 
                                            (STDERROR.base/(100 - VALUE.base))^2) * (100 - VALUE.predict.data),]


OUT_OF_SAMPLE_OBS[,t.value := (VALUE.predict.model - VALUE.predict.data)/sqrt(STDERROR.predict.data^2 + STDERROR.predict.model^2),]

OUT_OF_SAMPLE_OBS[,log_lik := dt(t.value, Inf, log = T),]
OUT_OF_SAMPLE_OBS[,sum(log_lik),.(MODEL, method)]


# Efficacy measures not compared to Naive ####

Waning_Table_Reduced[TYPE %in% c("S", "M", "D", "SD") &
                       (IMM %in% c("V", "B", "HI", "HB", "I", "H")) &
                       !is.na(LOWER) & LOWER!=UPPER&
                       (STRAIN == "Omicron" & PSTRAIN== "Omicron")
]->OUT_OF_SAMPLE_OBS2

merge(merge(
  OUT_OF_SAMPLE_OBS2,
Non_Linear_Fit[
  matching.table2, 
  on = .(IMMmodel), allow.cartesian = T
][WANING == 1 & STRAINmodel == "Pre-Omicron"],
by.x = c("TYPE","IMM"), 
by.y = c("TYPE","IMM"),
suffixes = c(".data", ".waned"), allow.cartesian = T
),
Non_Linear_Fit[
  matching.table2, 
  on = .(IMMmodel), allow.cartesian = T
][WANING == 0 & STRAINmodel == "Pre-Omicron"],
by.x = c("TYPE","IMM", "MODEL", "method"), 
by.y = c("TYPE","IMM", "MODEL", "method"),
suffixes = c(".waned", ".initial")
)->OUT_OF_SAMPLE_OBS2_expand


OUT_OF_SAMPLE_OBS2_expand[,VALUE.predict.model:=VALUE * exp(-decay_rate * MIN_TIME) + 
                    VALUE * (1 - exp(-decay_rate * MIN_TIME)),]


OUT_OF_SAMPLE_OBS2_expand[,STDERROR.predict.model:=sqrt(STDERROR^2 * exp(-2 *decay_rate * MIN_TIME)+ 
                                                  STDERROR^2 * (1 - exp(-decay_rate * MIN_TIME))^2),]

OUT_OF_SAMPLE_OBS2_expand[,t.value := (VALUE.predict.model - VALUE.data)/sqrt(STDERROR.data^2 + STDERROR.predict.model^2),]

OUT_OF_SAMPLE_OBS2_expand[,log_lik := dt(t.value, Inf, log = T),]

merge(
merge(
  OUT_OF_SAMPLE_OBS[,.(nLL1 = -sum(log_lik)),.(MODEL, method, TYPE, STRAINmodel)],
  OUT_OF_SAMPLE_OBS2_expand[,.(nLL2 = -sum(log_lik), STRAINmodel = "Pre-Omicron"),.(MODEL, method, TYPE)],
  all = TRUE
),
  Non_Linear_Fit[OUTPUT == "BIC",.(BIC = VALUE),.(MODEL, method, TYPE, STRAINmodel)]
)->LL_summary