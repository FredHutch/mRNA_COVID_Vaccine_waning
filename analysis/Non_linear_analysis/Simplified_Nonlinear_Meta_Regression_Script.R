################################################################################
######################Simplified Meta Regression Script#########################
################################################################################

#' The purpose of this script is to run a non-linear meta-regression making
#' several simplifying assumptions informed both by biology, previous linear
#' and non-linear meta-regressions.
#' We take advantage of previous code by creating new variables IMM_FIT and
#' STRAIN_SIMPLE

source("analysis/linear_analysis/Meta_Regression_Functions.R")
MULTI_LEVEL_FLAG = T
decay_rate = 0.05 #Twenty weeks
library(data.table)
library(metafor)
Waning_Table_Reduced = readRDS("data/Waning_Table_5_9_24.rds")
output_folder = paste0("data/Meta_Analysis_Results/Results_", Sys.Date())

dir.create(output_folder)
saveRDS(Waning_Table_Reduced, file = paste0(output_folder, "/Meta_Regression_Data.rds"))
#source("Parameter_Estimation.R")
Waning_Table_Reduced[,IMM := factor(IMM, levels = c("V", "I", "H", "B", "HB", "BvV", 
                                                    "BB", "BvB", "HvI", "HBvI", "Hvf", "HI", "HBvH")),]
Waning_Table_Reduced[,STRAIN := factor(STRAIN, levels = c("Pre-Omicron", "Omicron")),]
Waning_Table_Reduced[,PSTRAIN := factor(PSTRAIN, levels = c("Omicron", "Pre-Omicron", "None")),]
Waning_Table_Reduced[((STRAIN == "Omicron" & PSTRAIN== "Omicron") |
                        (STRAIN == "Pre-Omicron" & PSTRAIN %in% c("Pre-Omicron", "None"))),
                     STRAIN_SIMPLE:="Matched",]
Waning_Table_Reduced[STRAIN == "Omicron" & PSTRAIN %in% c("Pre-Omicron", "None"),
                     STRAIN_SIMPLE:="Omicron",]



# Severe Disease/Omicron Variants ####

# Lump together V, B and I/H/HB
Waning_Table_Reduced[STRAIN_SIMPLE=="Omicron" & TYPE == "M" & IMM == "V",
                     IMM_FIT:="V",]
Waning_Table_Reduced[STRAIN_SIMPLE=="Omicron" & TYPE == "M" & IMM == "B",
                     IMM_FIT:="B",]
Waning_Table_Reduced[STRAIN_SIMPLE=="Omicron" & TYPE == "M" & 
                       IMM %in% c("I", "H", "HB"),
                     IMM_FIT:="I_H_HB",]

RMA_M_omicron = rma.mv(VALUE, 
       V = STDERROR^2, 
       mods = ~WANING:(IMM_FIT == "B") + WANING + (IMM_FIT == "V"),
       struct = "GEN",
       data = Waning_Table_Reduced[STRAIN_SIMPLE == "Omicron" & 
                                     TYPE == "M" & 
                                     IMM %in% c("H", "I", "V", "B", "HB") & 
                                     !is.na(LOWER) & LOWER!=UPPER,
                                   .(VALUE, STDERROR, 
                                     IMM_FIT = factor(IMM_FIT, levels = c("V", "B", "I_H_HB")), 
                                     STUDY, effect_id, 
                   WANING = 1 - exp(-decay_rate * MIN_TIME))], 
       random = list(~WANING|STUDY, 
                     ~WANING|effect_id)
)

template_M_omicron = data.table(
  expand.grid(
    IMM_FIT = c("V", "B", "I_H_HB"),
    WANING = c(0, 1)
  )
)
parameters_M_omicron = extract_parameters(RMA_M_omicron, template_M_omicron, 
                                          formula(~WANING:(IMM_FIT == "B") + WANING + (IMM_FIT == "V")))

parameters_M_omicron[,TYPE_FIT:="M",]
parameters_M_omicron[,STRAIN_FIT:="Omicron",]
# Severe Disease/Pre-Omicron Variants ####

# Lump together V/H/I and B/HB
Waning_Table_Reduced[STRAIN_SIMPLE=="Matched" & TYPE == "M" & IMM %in% c("V", "H", "I"),
                     IMM_FIT:=c("V_H_I"),]
Waning_Table_Reduced[STRAIN_SIMPLE=="Matched" & TYPE == "M" & 
                       IMM %in% c("B", "HB"),
                     IMM_FIT:="B_HB",]

RMA_M_matched = rma.mv(VALUE, 
       V = STDERROR^2,
       mods = ~WANING,
       struct = "GEN",
       data = Waning_Table_Reduced[STRAIN_SIMPLE == "Matched" & TYPE == "M" & 
                                     IMM %in% c("H", "I", "V", "B", "HB")& 
                                     !is.na(LOWER) & LOWER!=UPPER,
                                   .(VALUE, STDERROR, IMM_FIT, STUDY, effect_id,
                                     WANING = 1 - exp(-decay_rate * MIN_TIME))], 
       random = list(~WANING|STUDY, 
                     ~WANING|effect_id)
)

template_M_matched = data.table(
  expand.grid(
    IMM_FIT = c("B_HB", "V_H_I"),
    WANING = c(0, 1)
  )
)
parameters_M_matched = extract_parameters(RMA_M_matched, template_M_matched, 
                                          formula(~WANING)
)

parameters_M_matched[,TYPE_FIT:="M",]
parameters_M_matched[,STRAIN_FIT:="Matched",]

# Infection/Omicron Variants ####


Waning_Table_Reduced[STRAIN_SIMPLE=="Omicron" & TYPE == "SD",
                     IMM_FIT:=IMM,]
Waning_Table_Reduced[STRAIN_SIMPLE=="Omicron" & TYPE == "SD",
                     IMM_FIT2:=IMM_FIT,]
Waning_Table_Reduced[STRAIN_SIMPLE=="Omicron" & TYPE == "SD" & 
                       IMM %in% c("H", "HB"),
                     IMM_FIT2:="H_HB",]

RMA_SD_omicron = rma.mv(VALUE, 
       V = STDERROR^2, 
       mods = ~WANING*(IMM_FIT2) + I(1 - WANING):IMM_FIT + 0,
       struct = "GEN",
       data = Waning_Table_Reduced[STRAIN_SIMPLE == "Omicron" & TYPE == "SD" & 
                                     IMM %in% c("H", "I", "V", "B", "HB") &
                                     !is.na(LOWER) & LOWER!=UPPER,
                                   .(VALUE, STDERROR, IMM_FIT, IMM_FIT2, STUDY, effect_id,
                                     WANING = 1 - exp(-decay_rate * MIN_TIME))], 
       random = list(~WANING|STUDY, 
                     ~WANING|effect_id)
)

template_SD_omicron = data.table(
  expand.grid(
    IMM_FIT = c("B", "HB", "V", "H", "I"),
    WANING = c(0, 1)
  )
)

template_SD_omicron[,IMM_FIT2:=IMM_FIT,]
template_SD_omicron[IMM_FIT %in% c("H", "HB"),IMM_FIT2:="H_HB",]


prediction_grid = 
  model.matrix( ~WANING*(IMM_FIT2) + I(1 - WANING):IMM_FIT + 0, template_SD_omicron)[,row.names(RMA_SD_omicron$b)]

model_prediction = as.data.table(predict.rma(RMA_SD_omicron, newmods = prediction_grid))

parameters_SD_omicron = cbind(template_SD_omicron, model_prediction[,.(VALUE = pred, se, 
                                    LOWER = ci.lb, 
                                    UPPER = ci.ub)])


parameters_SD_omicron[,TYPE_FIT:="SD",]
parameters_SD_omicron[,STRAIN_FIT:="Omicron",]

# Infection/Pre-Omicron Variants ####

Waning_Table_Reduced[STRAIN_SIMPLE=="Matched" & TYPE == "SD" & IMM == "V",
                     IMM_FIT:="V",]
Waning_Table_Reduced[STRAIN_SIMPLE=="Matched" & TYPE == "SD" & IMM == "B",
                     IMM_FIT:="B",]
Waning_Table_Reduced[STRAIN_SIMPLE=="Matched" & TYPE == "SD" & IMM == "I",
                     IMM_FIT:="I",]
Waning_Table_Reduced[STRAIN_SIMPLE=="Matched" & TYPE == "SD" & 
                       IMM %in% c("H", "HB"),
                     IMM_FIT:="H_HB",]

RMA_SD_matched = rma.mv(VALUE, 
       V = STDERROR^2, 
       mods = ~WANING:(IMM_FIT == "V") + WANING,
       struct = "GEN",
       data = Waning_Table_Reduced[STRAIN_SIMPLE == "Matched" & TYPE == "SD" & 
                                     IMM %in% c("H", "I", "V", "B", "HB") & 
                                     !is.na(LOWER) & LOWER!=UPPER,
                                   .(VALUE, STDERROR, IMM_FIT, STUDY, effect_id,
                                     WANING = 1 - exp(-decay_rate * MIN_TIME))], 
       random = list(~WANING|STUDY, 
                     ~WANING|effect_id)
)


template_SD_matched = data.table(
  expand.grid(
    IMM_FIT = c("B", "V", "I", "H_HB"),
    WANING = c(0, 1)
  )
)

parameters_SD_matched = extract_parameters(RMA_SD_matched, template_SD_matched, 
                                           formula(~WANING:(IMM_FIT == "V") + WANING)
)

parameters_SD_matched[,TYPE_FIT:="SD",]
parameters_SD_matched[,STRAIN_FIT:="Matched",]

Non_linear_parameters_simplified = rbind(parameters_M_matched, parameters_M_omicron, parameters_SD_matched, parameters_SD_omicron, fill = TRUE)
Non_linear_parameters_simplified[,CONSTRAINED := TRUE,]
Non_linear_parameters_simplified[,IMM_FIT2:=NULL,]
saveRDS(Non_linear_parameters_simplified, file = paste0(output_folder, "/Non_Linear_Fit_simplified.rds"))