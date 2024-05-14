################################################################################
######################Meta Regression Script####################################
################################################################################
MULTI_LEVEL_FLAG = T

Waning_Table_Reduced = readRDS("Waning_Table_5_24_24.rds")
output_folder = paste0("Meta_Analysis_Results/Results_", Sys.Date())

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

RFIG_STRAIN_DATA = 
  Waning_Table_Reduced[TYPE %in% c("S", "D", "M") &
                         IMM %in% c("V", "B", "H", "HB") &
                         !is.na(LOWER) & LOWER!=UPPER,
                       fit_rma_strain_wrap(.SD, multilevel = MULTI_LEVEL_FLAG),
                       .(TYPE)]

saveRDS(RFIG_STRAIN_DATA, paste0(output_folder, "/RFIG_STRAIN_DATA.rds"))


RFIG_BOOSTER_DATA = 
  Waning_Table_Reduced[TYPE %in% c("S", "M", "D") &
                         IMM %in% c("V", "B", "H", "HB") &
                         !is.na(LOWER) & LOWER!=UPPER,
                       fit_rma_booster_wrap(.SD, multilevel = MULTI_LEVEL_FLAG),
                       .(TYPE, STRAIN)]

saveRDS(RFIG_BOOSTER_DATA, paste0(output_folder, "//RFIG_BOOSTER_DATA.rds"))


RFIG_BOOSTER_DATA_mod = 
  Waning_Table_Reduced[TYPE %in% c("S", "M", "D") &
                         IMM %in% c("V", "B") &
                         !is.na(LOWER) & LOWER!=UPPER,
                       fit_rma_booster_wrap(.SD, multilevel = MULTI_LEVEL_FLAG),
                                             .(TYPE, STRAIN)]

saveRDS(RFIG_BOOSTER_DATA_mod, paste0(output_folder, "/RFIG_BOOSTER_DATA_mod.rds"))




RFIG_INF_DATA = 
  Waning_Table_Reduced[TYPE %in% c("S", "D") & 
                         !(TYPE != "S" & STRAIN == "Pre-Omicron") &
                         IMM %in% c("I", "V", "H") &
                         !is.na(LOWER) & LOWER!=UPPER,
                       fit_rma_inf_wrap(.SD, multilevel = MULTI_LEVEL_FLAG),
                       .(TYPE, STRAIN)]



Waning_Table_Reduced[TYPE %in% c("D") &
                       IMM %in% c("I", "V") & STRAIN == "Pre-Omicron" &
                       !is.na(LOWER) & LOWER!=UPPER,
                     fit_rma_inf_alt_wrap(.SD, multilevel = MULTI_LEVEL_FLAG),
                     .(TYPE, STRAIN)]->.out

.out[,TYPE:="D",]
.out[,STRAIN:="Pre-Omicron",]

RFIG_INF_DATA = rbind(RFIG_INF_DATA, .out)

saveRDS(RFIG_INF_DATA, paste0(output_folder, "/RFIG_INF_DATA.rds"))


Waning_Table_Reduced[TYPE %in% c("S", "D") & !(TYPE != "S" & STRAIN == "Pre-Omicron") &
                       IMM %in% c("I", "V", "H") &
                       !is.na(LOWER) & LOWER!=UPPER]->.tmp

.tmp[!is.na(PSTRAIN) & PSTRAIN == "Pre-Omicron", IMMP:=paste0(IMM,PSTRAIN),]
.tmp[PSTRAIN %in% c("Omicron", "Mix"), IMMP:=paste0(IMM,"Post-Omicron"),]
.tmp[is.na(PSTRAIN) & IMM == "V", IMMP:="V",]

RFIG_INF_DATA_pstrain = .tmp[TYPE %in% c("S", "D") & STRAIN == "Omicron" &
                               IMM %in% c("I", "V", "H") &
                               !is.na(LOWER) & LOWER!=UPPER
                             ,fit_rma_pstrain_wrap(.SD, multilevel = MULTI_LEVEL_FLAG),
                             .(TYPE, STRAIN)]




RFIG_INF_DATA_pstrain = merge(RFIG_INF_DATA_pstrain,
                              data.table(
                                IMMP = c("V", "V", "HPre-Omicron", "HPost-Omicron","IPre-Omicron", "IPost-Omicron"),
                                PSTRAIN = rep(c("Pre-Omicron", "Post-Omicron"), 3),
                                IMM = rep(c("V", "H", "I"), each = 2)
                              ), all = T, allow.cartesian = T)


saveRDS(RFIG_INF_DATA_pstrain, paste0(output_folder, "/RFIG_INF_DATA_pstrain.rds"))


RFIG_HYBRID_BOOSTED_DATA = 
  Waning_Table_Reduced[TYPE %in% c("S", "D") &
                         IMM %in% c("B", "HB") & STRAIN == "Omicron" &
                         !is.na(LOWER) & LOWER!=UPPER,
                       fit_rma_hboosted_wrap(.SD, multilevel = MULTI_LEVEL_FLAG),
                       .(TYPE, STRAIN)]


saveRDS(RFIG_HYBRID_BOOSTED_DATA, paste0(output_folder, "/RFIG_HYBRID_BOOSTED_DATA.rds"))
