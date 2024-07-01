output_folder = "data/Meta_Analysis_Results/Results_2024-06-24/"


library(data.table)
Meta_Regression_Data = readRDS(paste0(output_folder, "Meta_Regression_Data.rds"))
Non_Linear_Fit <- readRDS(paste0(output_folder, "Non_Linear_Fit.rds"))

setnames(Non_Linear_Fit, old = c("se", "STRAIN", "IMM"), 
         new = c("STDERROR", "STRAIN_FIT", "IMM_FIT"), skip_absent = T)
Non_Linear_Fit[is.na(STRAIN_FIT),STRAIN_FIT:="Matched",]

decay_rate = 1/20
# Filter Parameters Based on Constraints ####

## Generate Samples ####
N_samples = 1000

sample.parameter = 
  Non_Linear_Fit[OUTPUT == "parameters",.(iter = seq(N_samples), 
                                          SAMPLE = rnorm(N_samples, mean = VALUE, sd = STDERROR)),
                 .(IMM_FIT, WANING, MODEL, method, STRAIN_FIT, TYPE)]


## Filter out values >100 ####

#Only check SD (lump ES and ED into single estimate) and EM
filter_SD_only = sample.parameter[TYPE %in% c("SD", "M"),
                 mean(SAMPLE>100) == 0,
                 .(MODEL, method, iter, STRAIN_FIT, TYPE)][
                   (V1),
                   .(MODEL, method, iter, STRAIN_FIT, TYPE),]

#Check ES, ED and EM
filter_S_AND_D = sample.parameter[TYPE %in% c("S", "D", "M"),
                                  mean(SAMPLE>100) == 0,
                                  .(MODEL, method, iter, STRAIN_FIT, TYPE)][
                                    (V1),
                                    .(MODEL, method, iter, STRAIN_FIT, TYPE),]



## Filter based on additional constraints ####
# Assume hybrid better than primary AND
# Assume boosted hybrid better than hybrid

check_IMM_full = function(X){
  X[,
    (sum(SAMPLE * (IMM_FIT == "V"))>sum(SAMPLE * (IMM_FIT == "H"))) +
      (sum(SAMPLE * (IMM_FIT == "V"))>sum(SAMPLE * (IMM_FIT == "B"))) +
      (sum(SAMPLE * (IMM_FIT == "H"))>sum(SAMPLE * (IMM_FIT == "HB"))) * 
      sum(IMM_FIT == "HB")+
      (sum(SAMPLE * (IMM_FIT == "B"))>sum(SAMPLE * (IMM_FIT == "HB"))) * 
      sum(IMM_FIT == "HB"),
    .(WANING)][,sum(V1)==0,]
}





check_IMM_red = function(X){
  X[,
    (sum(SAMPLE * (IMM_FIT == "V"))>sum(SAMPLE * (IMM_FIT == "HI"))) +
      (sum(SAMPLE * (IMM_FIT == "V"))>sum(SAMPLE * (IMM_FIT == "B"))) +
      (sum(SAMPLE * (IMM_FIT == "HI"))>sum(SAMPLE * (IMM_FIT == "HB"))) * 
      sum(IMM_FIT == "HB")+
    (sum(SAMPLE * (IMM_FIT == "B"))>sum(SAMPLE * (IMM_FIT == "HB"))) * 
      sum(IMM_FIT == "HB"),,
    .(WANING)][,sum(V1)==0,]
}

### Fit S and D together, H and I together ####

sample.parameter[
  #Filter out negative values or those above 100
  filter_SD_only[
  MODEL %in% c("Combine Hybrid and Natural",
               "Use Omicron v Omicron data and combine Hybrid and Natural")
], 
                 on = .(MODEL, method, iter, TYPE, STRAIN_FIT)][
                   #Check that V>HI and HI>HB
                   TYPE %in% c("SD"),
                 check_IMM_red(.SD),
                 .(MODEL,method,iter, STRAIN_FIT, TYPE)
                 ][(V1)][
                   #Compute number satisfying all conditions by model and outcome
                   ,.N,.(MODEL, method, STRAIN_FIT, TYPE)][
                     #Merge with BIC data
                   Non_Linear_Fit[OUTPUT == "BIC" &
                                    MODEL %in% c("Combine Hybrid and Natural",
                                                 "Use Omicron v Omicron data and combine Hybrid and Natural") & 
                                    TYPE %in% c("SD"),
                                  .(BIC = VALUE),
                                  .(MODEL, method, STRAIN_FIT, TYPE)], 
                   on = .(MODEL, method, STRAIN_FIT, TYPE)
                   ] ->filter_SD_only_reduced



### Fit S and D together, H and I separately ####

sample.parameter[
  #Filter out negative values or those above 100
  filter_SD_only[
    MODEL %in% c("Main Fit",
                 "Use Omicron v Omicron data")
  ], 
  on = .(MODEL, method, iter, TYPE, STRAIN_FIT)][
    #Check that V>HI and HI>HB
    TYPE %in% c("SD"),
    check_IMM_full(.SD),
    .(MODEL,method,iter, STRAIN_FIT, TYPE)
  ][(V1)][
    #Compute number satisfying all conditions by model and outcome
    ,.N,.(MODEL, method, STRAIN_FIT, TYPE)][
      #Merge with BIC data
      Non_Linear_Fit[OUTPUT == "BIC" & 
                       TYPE %in% c("SD") &
                     MODEL %in% c("Main Fit",
                                  "Use Omicron v Omicron data"),
                     .(BIC = VALUE),
                     .(MODEL, method, STRAIN_FIT, TYPE)], 
      on = .(MODEL, method, STRAIN_FIT, TYPE)
    ] ->filter_SD_only_full

### Fit S and D separately, H and I together ####

sample.parameter[
  #Filter out negative values or those above 100
  filter_S_AND_D[
    MODEL %in% c("Combine Hybrid and Natural",
                 "Use Omicron v Omicron data and combine Hybrid and Natural")
  ], 
  on = .(MODEL, method, iter, TYPE, STRAIN_FIT)][
    #Check that V>HI and HI>HB
    TYPE %in% c("S", "D", "M"),
    check_IMM_red(.SD),
    .(MODEL,method,iter, STRAIN_FIT, TYPE)
  ][(V1)][
    #Compute number satisfying all conditions by model and outcome
    ,.N,.(MODEL, method, STRAIN_FIT, TYPE)][
      #Merge with BIC data
      Non_Linear_Fit[OUTPUT == "BIC" & 
                       TYPE %in% c("S", "D", "M") & 
                       MODEL %in% c("Combine Hybrid and Natural",
                                    "Use Omicron v Omicron data and combine Hybrid and Natural"),
                     .(BIC = VALUE),
                     .(MODEL, method, STRAIN_FIT, TYPE)], 
      on = .(MODEL, method, STRAIN_FIT, TYPE)
    ] ->filter_S_AND_D_reduced

### Fit S and D separately, H and I separately ####

sample.parameter[
  #Filter out negative values or those above 100
  filter_S_AND_D[
    MODEL %in% c("Main Fit",
                 "Use Omicron v Omicron data")
  ], 
  on = .(MODEL, method, iter, TYPE, STRAIN_FIT)][
    #Check that V>HI and HI>HB
    TYPE %in% c("S", "D", "M"),
    check_IMM_full(.SD),
    .(MODEL,method,iter, STRAIN_FIT, TYPE)
  ][(V1)][
    #Compute number satisfying all conditions by model and outcome
    ,.N,.(MODEL, method, STRAIN_FIT, TYPE)][
      #Merge with BIC data
      Non_Linear_Fit[OUTPUT == "BIC" & 
                       TYPE %in% c("S", "D", "M") &
                       MODEL %in% c("Main Fit",
                                    "Use Omicron v Omicron data"),
                     .(BIC = VALUE),
                     .(MODEL, method, STRAIN_FIT, TYPE)], 
      on = .(MODEL, method, STRAIN_FIT, TYPE)
    ] ->filter_S_AND_D_full

## Combine filters and calculate proportion of samples kept ####
filter_all = rbind(
  filter_SD_only_reduced,
  filter_SD_only_full,
  filter_S_AND_D_reduced,
  filter_S_AND_D_full
)

filter_all[,nLL3:= - log(N/N_samples),]

# Out of sample comparisons ####

#Mapping between fitted quantities and literature values
#Note that different models map differently. Some models
# model "Pre-Omicron" immunity specifically, others model
# "matched" which incorporates both Pre-Omicron and Omicron-induced
# protection against Omicron.

matching.table = data.table(
  #Strain column refers to the "true" value from the meta-regression data
  STRAIN = rep(c("Pre-Omicron", "Omicron"), each = 4),
  
  #PSTRAIN = prior infecting strain, if any
  PSTRAIN = c("None", "None", "Pre-Omicron", "Pre-Omicron",
              "Pre-Omicron", "None", "Omicron", "Omicron"),
  
  #How this is modeled
  STRAIN_FIT = c("Pre-Omicron", "Matched", "Pre-Omicron", "Matched",
                 "Omicron", "Omicron", "Pre-Omicron", "Matched")
)


matching.table2 = data.table(
  IMM_FIT = c("V", "I", "H", "HI", "HI", "HB", "B"),
  IMM = c("V", "I", "H", "H", "I", "HB", "B")
)

## Efficacy measures not compared to Naive ####
Meta_Regression_Data[TYPE %in% c("S", "M", "D", "SD") &
                       !(IMM %in% c("V", "B", "HI", "HB", "I", "H")) &
                       !is.na(LOWER) & LOWER!=UPPER
]->OUT_OF_SAMPLE_OBS


## Categorize by what immunity type the measure is compared to ####
OUT_OF_SAMPLE_OBS[IMM == "BvB",IMM_BASE:="B",]
OUT_OF_SAMPLE_OBS[IMM %in% c("HBvV", "BvV"), IMM_BASE:="V",]
OUT_OF_SAMPLE_OBS[IMM %in% c("HvI", "HBvI"), IMM_BASE:="I",]
OUT_OF_SAMPLE_OBS[IMM %in% c("BvB", "BvV"), IMM_MAIN:="B",]
OUT_OF_SAMPLE_OBS[IMM %in% c("HBvV", "HBvI"), IMM_MAIN:="HB",]
OUT_OF_SAMPLE_OBS[IMM %in% c("HvI", "Hvf"), IMM_MAIN:="H",]


Non_Linear_Fit[
  matching.table, 
  on = .(STRAIN_FIT), allow.cartesian = T
][
  matching.table2,
  on = .(IMM_FIT), allow.cartesian = T
]->Non_Linear_Fit_expand


### Merge Data Sets ####
#Note that this assumes that the base immunity is always "waned"
merge(merge(merge(OUT_OF_SAMPLE_OBS, Non_Linear_Fit_expand[WANING == 1], 
                  by.x = c("TYPE", "PSTRAIN", "STRAIN", "IMM_BASE"), 
                  by.y = c("TYPE", "PSTRAIN", "STRAIN", "IMM"), 
                  allow.cartesian = T, suffixes = c(".data", ".base")),
            Non_Linear_Fit_expand[WANING == 1],
            by.x = c("TYPE", "PSTRAIN", "STRAIN", "IMM_MAIN", "MODEL", "method", "STRAIN_FIT", "OUTPUT"), 
            by.y = c("TYPE", "PSTRAIN", "STRAIN", "IMM", "MODEL", "method", "STRAIN_FIT", "OUTPUT"),
            suffixes = c("", ".waned")
),
Non_Linear_Fit_expand[WANING == 0],
by.x = c("TYPE", "PSTRAIN", "STRAIN", "IMM_MAIN", "MODEL", "method", "STRAIN_FIT", "OUTPUT"), 
by.y = c("TYPE", "PSTRAIN", "STRAIN", "IMM", "MODEL", "method", "STRAIN_FIT", "OUTPUT"),
suffixes = c("", ".initial")
)->OUT_OF_SAMPLE_OBS


### Compare prediction with data ####

# Model prediction (incorporates waning)
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


## Efficacy measures of Omicron vs Omicron ####

#Subset of meta regression data
Meta_Regression_Data[TYPE %in% c("S", "M", "D", "SD") &
                       (IMM %in% c("V", "B", "HI", "HB", "I", "H")) &
                       !is.na(LOWER) & LOWER!=UPPER&
                       (STRAIN == "Omicron" & PSTRAIN== "Omicron")
]->OUT_OF_SAMPLE_OBS2

### Merge data ####
#Note that we only merge models that fit "Pre-Omicron" specifically and
#ignored these Omicron vs Omicron observations.
merge(merge(
  OUT_OF_SAMPLE_OBS2,
  Non_Linear_Fit[
    matching.table2, 
    on = .(IMM_FIT), allow.cartesian = T
  ][WANING == 1 & STRAIN_FIT == "Pre-Omicron"],
  by.x = c("TYPE","IMM"), 
  by.y = c("TYPE","IMM"),
  suffixes = c(".data", ".waned"), allow.cartesian = T
),
Non_Linear_Fit[
  matching.table2, 
  on = .(IMM_FIT), allow.cartesian = T
][WANING == 0 & STRAIN_FIT == "Pre-Omicron"],
by.x = c("TYPE","IMM", "MODEL", "method"), 
by.y = c("TYPE","IMM", "MODEL", "method"),
suffixes = c(".waned", ".initial")
)->OUT_OF_SAMPLE_OBS2_expand

### Compare Prediction and Data ####
OUT_OF_SAMPLE_OBS2_expand[,VALUE.predict.model:=VALUE * exp(-decay_rate * MIN_TIME) + 
                            VALUE * (1 - exp(-decay_rate * MIN_TIME)),]


OUT_OF_SAMPLE_OBS2_expand[,STDERROR.predict.model:=sqrt(STDERROR^2 * exp(-2 *decay_rate * MIN_TIME)+ 
                                                          STDERROR^2 * (1 - exp(-decay_rate * MIN_TIME))^2),]

OUT_OF_SAMPLE_OBS2_expand[,t.value := (VALUE.predict.model - VALUE.data)/sqrt(STDERROR.data^2 + STDERROR.predict.model^2),]

OUT_OF_SAMPLE_OBS2_expand[,log_lik := dt(t.value, Inf, log = T),]

# Merge Results of Out of Sample and Filtering Analysis ####
merge(filter_all, merge(
  merge(
    OUT_OF_SAMPLE_OBS[,.(nLL1 = -sum(log_lik)),.(MODEL, method, TYPE, STRAIN_FIT)],
    OUT_OF_SAMPLE_OBS2_expand[,.(nLL2 = -sum(log_lik), STRAIN_FIT = "Pre-Omicron"),.(MODEL, method, TYPE)],
    all = TRUE
  ),
  Non_Linear_Fit[OUTPUT == "BIC",.(BIC = VALUE),.(MODEL, method, TYPE, STRAIN_FIT)]
))->LL_summary

# Find best fits models for each outcome ####

## Score models with BIC ####
LL_summary[is.na(nLL2),nLL2 := 0,]
LL_summary[is.na(nLL3),nLL3 := Inf,]
LL_summary[,TOTAL_BIC_constrained := BIC + 2 * nLL1 + 2 * nLL2 + 2 * nLL3,]
LL_summary[,TOTAL_BIC_unconstrained := BIC + 2 * nLL1 + 2 * nLL2,]

## Find best unconstrained fits (no filtering) ####
LL_summary[
  LL_summary[,.(TOTAL_BIC_unconstrained = min(TOTAL_BIC_unconstrained, na.rm = T),
                method = method[which.min(TOTAL_BIC_unconstrained)], 
                MODEL = MODEL[which.min(TOTAL_BIC_unconstrained)]),
             .(STRAIN_FIT, TYPE)],
  on = .(MODEL, method, STRAIN_FIT, TYPE, TOTAL_BIC_unconstrained)]->Best_Fits_unconstrained

## Find best constrained fits (with filtering) ####
LL_summary[
  LL_summary[,.(TOTAL_BIC_constrained = min(TOTAL_BIC_constrained, na.rm = T),
                method = method[which.min(TOTAL_BIC_constrained)], 
                MODEL = MODEL[which.min(TOTAL_BIC_constrained)]),
             .(STRAIN_FIT, TYPE)],
  on = .(MODEL, method, STRAIN_FIT, TYPE, TOTAL_BIC_constrained)]->Best_Fits_constrained

## Find best fits from regression alone (no out of sample) ####
LL_summary[
  LL_summary[,.(BIC = min(BIC, na.rm = T),
                method = method[which.min(BIC)], 
                MODEL = MODEL[which.min(BIC)]),
             .(STRAIN_FIT, TYPE)],
  on = .(MODEL, method, STRAIN_FIT, TYPE, BIC)]->Best_Fits_regression

Best_Fits_constrained[,TOTAL_BIC:=TOTAL_BIC_constrained,]
Best_Fits_unconstrained[,TOTAL_BIC:=TOTAL_BIC_unconstrained,]
Best_Fits_regression[,TOTAL_BIC:=BIC,]

Best_Fits_unconstrained[,constrained:=FALSE,]
Best_Fits_constrained[,constrained:=TRUE,]
Best_Fits_regression[,constrained:=FALSE,]
Best_Fits_unconstrained[,out_of_sample:=TRUE,]
Best_Fits_constrained[,out_of_sample:=TRUE,]
Best_Fits_regression[,out_of_sample:=FALSE,]

Best_Fits = rbind(Best_Fits_constrained, Best_Fits_unconstrained, Best_Fits_regression, fill = TRUE)

## Compare modeling ES and ED alone vs together ####
Best_Fits[TYPE!="M",
          USE_SD:=sum(TOTAL_BIC * (TYPE == "SD"))<sum(TOTAL_BIC * (TYPE!="SD")),
          .(STRAIN_FIT, constrained)]

## Evaluate choice of using vs not using Omicron data in fitting ####
Best_Fits[STRAIN_FIT!="Omicron",
          USE_matched:=sum(TOTAL_BIC * (STRAIN_FIT == "Matched"))<sum(TOTAL_BIC * (STRAIN_FIT!="Matched")),
          .(TYPE, constrained)]

Best_Fits[
  (TYPE == "SD" & USE_SD == TRUE)|(TYPE%in% c("S", "D") & USE_SD == FALSE)|
    (is.na(USE_SD))
  ][
    (STRAIN_FIT == "Matched" & USE_matched== TRUE)|(STRAIN_FIT == c("Pre-Omicron") & USE_matched == FALSE)|
      (is.na(USE_matched))
  ]->Best_Fits_reduced

Non_Linear_Fit[Best_Fits, on = .(MODEL, method, TYPE, STRAIN_FIT)
               ][OUTPUT == "parameters",.(VALUE, STDERROR, LOWER, UPPER),
                 .(IMM_FIT, WANING_FIT = WANING, TYPE_FIT = TYPE, STRAIN_FIT, 
                   CONSTRAINED = constrained, USE_NON_STANDARD = out_of_sample)]->Non_Linear_Fit_reduced

saveRDS(Non_Linear_Fit_reduced, file = paste0(output_folder, "/Non_Linear_Fit_reduced.rds"))
saveRDS(LL_summary, file = paste0(output_folder, "/Fitting_Summary.rds"))

check_S_D_M = function(X){
  X[,(SAMPLE[3]>SAMPLE[1]) + #S>D
      (SAMPLE[1]>SAMPLE[2]), #D>M
    .(IMM, WANING)][,sum(V1)==0,STRAIN_FIT]
}


check_SD_M = function(X){
  X[,(SAMPLE[4]>SAMPLE[2]), #SD>M
    .(IMM, STRAIN_FIT, WANING)][,sum(V1)==0,STRAIN_FIT]
}

check_IMM_full_individual = function(X){
  X[,
    .(V_greater_than_H = (sum(VALUE * (IMM_FIT == "V"))>sum(VALUE * (IMM_FIT == "H"))),
      V_greater_than_B = (sum(VALUE * (IMM_FIT == "V"))>sum(VALUE * (IMM_FIT == "B"))),
      H_greater_than_HB = (sum(VALUE * (IMM_FIT == "H"))>sum(VALUE * (IMM_FIT == "HB"))) * sum(IMM_FIT == "HB"),
      B_greater_than_HB = (sum(VALUE * (IMM_FIT == "B"))>sum(VALUE * (IMM_FIT == "HB"))) * sum(IMM_FIT == "HB")),
  ]
}

Non_Linear_Fit_reduced[,check_IMM_full_individual(.SD),.(WANING_FIT, TYPE_FIT, STRAIN_FIT, CONSTRAINED, USE_NON_STANDARD)]->Problems