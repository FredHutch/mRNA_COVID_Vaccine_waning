################################################################################
######################Meta Regression Functions#################################
################################################################################

#Purpose: Defining the meta regression functions used in waning meta-analysis

# Model One: Efficacy ~ TIME * STRAIN ####

# Wrapping function loops over all methods
fit_rma_strain_wrap = function(X, multilevel = F){
  if(!multilevel){
    fit_rma_booster(X, "single_random", multilevel)
  }
  
  .DT = data.table(method = c("intercept only: study and estimate",
                              "slope-intercept correlation: study and estimate",
                              "slope-intercept-modulator correlation: study and estimate",
                              "intercept only: study only",
                              "slope-intercept correlation: study only",
                              "slope-intercept-modulator correlation: study only",
                              "intercept only: estimate only",
                              "slope-intercept correlation: estimate only",
                              "slope-intercept-modulator correlation: estimate only"))
  .DT[,fit_rma_strain(X,method, multilevel = multilevel),method]
  
}


fit_rma_strain = function(X, Method, multilevel = F){
  #Old method, not multilevel
  out.RMA = rma(VALUE, 
                sei=STDERROR, 
                mods = ~MIN_TIME * STRAIN, 
                data = X
  )
  
  if(multilevel){
    if(Method == "intercept only: study and estimate"){
      #Single Random effect
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * STRAIN, 
                       data = X, random = list(~1|STUDY, 
                                               ~1|effect_id)
      )
    }
    
    if(Method == "slope-intercept correlation: study and estimate"){
      #Auto-correlated random effects following exp(-(time_1 - time_2)/rho)
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * STRAIN, struct = "GEN",
                       data = X, random = list(~MIN_TIME|STUDY,
                                               ~MIN_TIME|effect_id)
      )
    }
    
    if(Method == "slope-intercept-modulator correlation: study and estimate"){
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * STRAIN, struct = "GEN",
                       data = X, random = list(~STRAIN + MIN_TIME|STUDY,
                                               ~STRAIN + MIN_TIME|effect_id)
      )
    }
    
    if(Method == "intercept only: study only"){
      #Single Random effect
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * STRAIN, 
                       data = X, random = ~1|STUDY
      )
    }
    
    if(Method == "slope-intercept correlation: study only"){
      #Auto-correlated random effects following exp(-(time_1 - time_2)/rho)
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * STRAIN, struct = "GEN",
                       data = X, random = ~MIN_TIME|STUDY
      )
    }
    
    if(Method == "slope-intercept-modulator correlation: study only"){
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * STRAIN, struct = "GEN",
                       data = X, random = ~STRAIN + MIN_TIME|STUDY
      )
    }
    
    if(Method == "intercept only: estimate only"){
      #Single Random effect
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * STRAIN, 
                       data = X, random = ~1|effect_id
      )
    }
    
    if(Method == "slope-intercept correlation: estimate only"){
      #Auto-correlated random effects following exp(-(time_1 - time_2)/rho)
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * STRAIN, struct = "GEN",
                       data = X, random = ~MIN_TIME|effect_id
      )
    }
    
    if(Method == "slope-intercept-modulator correlation: estimate only"){
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * STRAIN, struct = "GEN",
                       data = X, random = ~STRAIN + MIN_TIME|effect_id
      )
    }
  }
  
  #Model predictions
  .out = as.data.table(predict.rma(out.RMA))
  
  #Model weights
  .out$WEIGHT = weights(out.RMA)
  
  #Combine model predictions, weights, and input data
  .out = cbind(.out[,.(pred, se, ci.lb, ci.ub, WEIGHT),], X[,.(VALUE, LOWER, UPPER, STRAIN, MIN_TIME),])
  
  .out$OUTPUT = "Predicted"
  
  
  #Computing of global means
  template_a = data.table(
    expand.grid(MIN_TIME = c(2, 30),
                STRAIN = c("Pre-Omicron", "Omicron")
    )
  )
  prediction.grid_a = model.matrix(~MIN_TIME * STRAIN, template_a)[,-1]
  prediction.grid_b = cbind(
    MIN_TIME = c(-28, -28),
    `STRAINOmicron` = c(0,0),
    `MIN_TIME:STRAINOmicron` = c(-28, 0)
  )
  
  .out2a = as.data.table(predict.rma(out.RMA, newmods = prediction.grid_a))
  .out2b = as.data.table(predict.rma(out.RMA, newmods = prediction.grid_b, intercept = FALSE))
  
  .out2 = cbind(rbind(.out2a[,.(pred, se, ci.lb, ci.ub),], 
                      .out2b[,.(pred, se, ci.lb, ci.ub),]), 
                rbind(template_a, 
                      data.frame(
                        STRAIN = c("Omicron", 
                                   "Pre-Omicron")), 
                      fill = T)
  )
  .out2$OUTPUT = "Standards"
  
  effect.grid = cbind(
    MIN_TIME = c(0, 0, 0),
    `STRAINOmicron` = c(1, 1, 0),
    `MIN_TIME:STRAINOmicron` = c(2, 30, 1)
  )
  
  .out3 = as.data.table(predict.rma(out.RMA, newmods = effect.grid, intercept = FALSE))
  
  .out3 = cbind(.out3[,.(pred, se, ci.lb, ci.ub),], MIN_TIME = c(2, 30, NA), STRAIN = "p-value")
  .out3$OUTPUT = "Contrasts"
  
  .out4 = data.table(
    VALUE = out.RMA$fit.stats$REML,
    OUTPUT = row.names(out.RMA$fit.stats)
  )
  
  .out = rbind(.out, .out2, .out3, .out4, fill = T)
  
}

# Model Two: Efficacy ~ TIME * Dose ####

fit_rma_booster_wrap = function(X, multilevel = F){
  if(!multilevel){
    fit_rma_booster(X, "single_random", multilevel)
  }
  
  .DT = data.table(method = c("intercept only: study and estimate",
                              "slope-intercept correlation: study and estimate",
                              "slope-intercept-modulator correlation: study and estimate",
                              "intercept only: study only",
                              "slope-intercept correlation: study only",
                              "slope-intercept-modulator correlation: study only",
                              "intercept only: estimate only",
                              "slope-intercept correlation: estimate only",
                              "slope-intercept-modulator correlation: estimate only"))
  .DT[,fit_rma_booster(X,method, multilevel = multilevel),method]
  
}


fit_rma_booster = function(X, Method, multilevel = F){
  out.RMA = rma(VALUE, 
                sei=STDERROR, 
                mods = ~MIN_TIME * Boosted, 
                data = X
  )
  
  if(multilevel){
    if(Method == "intercept only: study and estimate"){
      #Single Random effect
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * Boosted, 
                       data = X, random = list(~1|STUDY, 
                                               ~1|effect_id)
      )
    }
    
    if(Method == "slope-intercept correlation: study and estimate"){
      #Auto-correlated random effects following exp(-(time_1 - time_2)/rho)
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * Boosted, struct = "GEN",
                       data = X, random = list(~MIN_TIME|STUDY,
                                               ~MIN_TIME|effect_id)
      )
    }
    
    if(Method == "slope-intercept-modulator correlation: study and estimate"){
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * Boosted, struct = "GEN",
                       data = X, random = list(~Boosted + MIN_TIME|STUDY,
                                               ~Boosted + MIN_TIME|effect_id)
      )
    }
    
    if(Method == "intercept only: study only"){
      #Single Random effect
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * Boosted, 
                       data = X, random = ~1|STUDY
      )
    }
    
    if(Method == "slope-intercept correlation: study only"){
      #Auto-correlated random effects following exp(-(time_1 - time_2)/rho)
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * Boosted, struct = "GEN",
                       data = X, random = ~MIN_TIME|STUDY
      )
    }
    
    if(Method == "slope-intercept-modulator correlation: study only"){
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * Boosted, struct = "GEN",
                       data = X, random = ~Boosted + MIN_TIME|STUDY
      )
    }
    
    if(Method == "intercept only: study only"){
      #Single Random effect
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * Boosted, 
                       data = X, random = ~1|effect_id
      )
    }
    
    if(Method == "slope-intercept correlation: study only"){
      #Auto-correlated random effects following exp(-(time_1 - time_2)/rho)
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * Boosted, struct = "GEN",
                       data = X, random = ~MIN_TIME|effect_id
      )
    }
    
    if(Method == "slope-intercept-modulator correlation: study only"){
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * Boosted, struct = "GEN",
                       data = X, random = ~Boosted + MIN_TIME|effect_id
      )
    }
  }
  .out = as.data.table(predict.rma(out.RMA))
  
  .out$WEIGHT = weights(out.RMA)
  .out = cbind(.out[,.(pred, se, ci.lb, ci.ub, WEIGHT)], X[,.(VALUE, LOWER, UPPER, Boosted, MIN_TIME),])
  
  .out$OUTPUT = "Predicted"
  template_a = data.table(
    expand.grid(MIN_TIME = c(2, 30),
                #STRAIN = c("Omicron", "Pre-Omicron"),
                Boosted = c(T, F)
    )
  )
  
  prediction.grid_a = model.matrix(~MIN_TIME * Boosted, template_a)[,-1]
  prediction.grid_b = cbind(
    MIN_TIME = c(-28, -28),
    `BoostedTRUE` = c(0,0),
    `MIN_TIME:BoostedTRUE` = c(0, -28)
  )
  
  .out2a = as.data.table(predict.rma(out.RMA, newmods = prediction.grid_a))
  .out2b = as.data.table(predict.rma(out.RMA, newmods = prediction.grid_b, intercept = FALSE))
  
  .out2 = cbind(rbind(.out2a[,.(pred, se, ci.lb, ci.ub),], .out2b[,.(pred, se, ci.lb, ci.ub),]), 
                rbind(template_a, 
                      data.frame(
                        Boosted = c(F, 
                                    T)), 
                      fill = T)
  )
  
  .out2$OUTPUT = "Standards"
  
  effect.grid = cbind(
    MIN_TIME = c(0, 0, 0),
    `BoostedTRUE` = c(1,1,0),
    `MIN_TIME:BoostedTRUE` = c(2, 30, 1)
  )
  
  .out3 = as.data.table(predict.rma(out.RMA, newmods = effect.grid, intercept = FALSE))
  
  .out3 = cbind(.out3[,.(pred, se, ci.lb, ci.ub),], MIN_TIME = c(2, 30, NA), Boosted = "p-value")
  .out3$OUTPUT = "Contrasts"
  
  
  .out4 = data.table(
    VALUE = out.RMA$fit.stats$REML,
    OUTPUT = row.names(out.RMA$fit.stats)
  )
  rbind(.out, .out2, .out3, .out4, fill = T)
  
}


# Model Three: Efficacy ~ TIME * Infection ####

fit_rma_inf_wrap = function(X, multilevel = F){
  if(!multilevel){
    fit_rma_inf(X, "single_random", multilevel)
  }
  
  .DT = data.table(method = c("intercept only: study and estimate",
                              "slope-intercept correlation: study and estimate",
                              "slope-intercept-modulator correlation: study and estimate",
                              "intercept only: study only",
                              "slope-intercept correlation: study only",
                              "slope-intercept-modulator correlation: study only",
                              "intercept only: estimate only",
                              "slope-intercept correlation: estimate only",
                              "slope-intercept-modulator correlation: estimate only"))
  .DT[,fit_rma_inf_error_catch(X,method, multilevel = multilevel),method]
  
}

fit_rma_inf_error_catch = function(X, Method, multilevel = F){
  result = tryCatch({
    fit_rma_inf(X, Method, multilevel = multilevel)
  },
  error = function(e){
    print(e)
    return(data.table(
      pred = numeric(0),
      se = numeric(0),
      ci.lb = numeric(0),
      ci.ub = numeric(0),
      WEIGHT = numeric(0),
      VALUE = numeric(0),
      LOWER = numeric(0),
      UPPER = numeric(0),
      IMM = integer(0),
      MIN_TIME = numeric(0),
      OUTPUT = "BIC"
    ))
  },
  finally = print(paste0(Method, ": complete")))
}

fit_rma_inf = function(X, Method, multilevel = F){
  
  out.RMA = rma(VALUE, 
                sei=STDERROR, 
                mods = ~MIN_TIME * IMM, 
                data = X
  )
  
  if(multilevel){
    if(Method == "intercept only: study and estimate"){
      #Single Random effect
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMM, 
                       data = X, random = list(~1|STUDY, 
                                               ~1|effect_id)
      )
    }
    
    if(Method == "slope-intercept correlation: study and estimate"){
      #Auto-correlated random effects following exp(-(time_1 - time_2)/rho)
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMM, struct = "GEN",
                       data = X, random = list(~MIN_TIME|STUDY,
                                               ~MIN_TIME|effect_id)
      )
    }
    
    if(Method == "slope-intercept-modulator correlation: study and estimate"){
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMM, struct = "GEN",
                       data = X, random = list(~IMM + MIN_TIME|STUDY,
                                               ~IMM + MIN_TIME|effect_id)
      )
    }
    
    if(Method == "intercept only: study only"){
      #Single Random effect
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMM, 
                       data = X, random = ~1|STUDY
      )
    }
    
    if(Method == "slope-intercept correlation: study only"){
      #Auto-correlated random effects following exp(-(time_1 - time_2)/rho)
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMM, struct = "GEN",
                       data = X, random = ~MIN_TIME|STUDY
      )
    }
    
    if(Method == "slope-intercept-modulator correlation: study only"){
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMM, struct = "GEN",
                       data = X, random = ~IMM + MIN_TIME|STUDY
      )
    }
    
    if(Method == "intercept only: estimate only"){
      #Single Random effect
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMM, 
                       data = X, random = ~1|effect_id
      )
    }
    
    if(Method == "slope-intercept correlation: estimate only"){
      #Auto-correlated random effects following exp(-(time_1 - time_2)/rho)
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMM, struct = "GEN",
                       data = X, random = ~MIN_TIME|effect_id
      )
    }
    
    if(Method == "slope-intercept-modulator correlation: estimate only"){
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMM, struct = "GEN",
                       data = X, random = ~IMM + MIN_TIME|effect_id
      )
    }
  }
  .out = as.data.table(predict.rma(out.RMA))
  
  .out$WEIGHT = weights(out.RMA)
  .out = cbind(.out[,.(pred, se, ci.lb, ci.ub, WEIGHT),], 
               X[,.(VALUE, LOWER, UPPER, IMM, MIN_TIME),])
  
  .out$OUTPUT = "Predicted"
  template_a = data.table(
    expand.grid(MIN_TIME = c(2, 30),
                #STRAIN = c("Omicron", "Pre-Omicron"),
                IMM = c("V", "I", "H")
    )
  )
  
  prediction.grid_a = model.matrix(~MIN_TIME * IMM, template_a)[,-1]
  prediction.grid_b = cbind(
    MIN_TIME = c(-28, -28, -28),
    `IMMH` = c(0,0,0),
    `IMMI` = c(0,0,0),
    `MIN_TIME:IMMH` = c(0, 0, -28),
    `MIN_TIME:IMMI` = c(0, -28, 0)
  )
  
  .out2a = as.data.table(predict.rma(out.RMA, newmods = prediction.grid_a))
  .out2b = as.data.table(predict.rma(out.RMA, newmods = prediction.grid_b, intercept = FALSE))
  
  .out2 = cbind(rbind(.out2a[,.(pred, se, ci.lb, ci.ub),], 
                      .out2b[,.(pred, se, ci.lb, ci.ub),]), 
                rbind(template_a, 
                      data.frame(
                        IMM = c("V",
                                "I", "H")), 
                      fill = T)
  )
  
  .out2$OUTPUT = "Standards"
  
  effect.grid = cbind(
    MIN_TIME = c(0, 0, 0, 0, 0, 0),
    `IMMH` = c(0, 0, 0, 1, 1, 0),
    `IMMI` = c(1, 1, 0, 0, 0, 0),
    `MIN_TIME:IMMH` = c(0,  0, 0, 2, 30, 1),
    `MIN_TIME:IMMI` = c(2, 30, 1, 0,  0, 0)
  )
  
  .out3 = as.data.table(predict.rma(out.RMA, newmods = effect.grid, intercept = FALSE))
  
  .out3 = cbind(.out3[,.(pred, se, ci.lb, ci.ub),], 
                MIN_TIME = c(2, 30, NA, 2, 30, NA), 
                IMM = rep(c("I vs V", "H vs V"), each = 3))
  .out3$OUTPUT = "Contrasts"
  
  .out4 = data.table(
    VALUE = out.RMA$fit.stats$REML,
    OUTPUT = row.names(out.RMA$fit.stats)
  )
  rbind(.out, .out2, .out3, .out4, fill = T)
  
}


fit_rma_inf_alt_wrap = function(X, multilevel = F){
  if(!multilevel){
    fit_rma_inf_alt(X, "single_random", multilevel)
  }
  
  .DT = data.table(method = c("intercept only: study and estimate",
                              "slope-intercept correlation: study and estimate",
                              "slope-intercept-modulator correlation: study and estimate",
                              "intercept only: study only",
                              "slope-intercept correlation: study only",
                              "slope-intercept-modulator correlation: study only",
                              "intercept only: estimate only",
                              "slope-intercept correlation: estimate only",
                              "slope-intercept-modulator correlation: estimate only"))
  .DT[,fit_rma_inf_alt_error_catch(X,method, multilevel = multilevel),method]
  
}

fit_rma_inf_alt_error_catch = function(X, Method, multilevel = F){
  result = tryCatch({
    fit_rma_inf_alt(X, Method, multilevel = multilevel)
  },
  error = function(e){
    print(e)
    return(data.table(
      pred = numeric(0),
      se = numeric(0),
      ci.lb = numeric(0),
      ci.ub = numeric(0),
      WEIGHT = numeric(0),
      VALUE = numeric(0),
      LOWER = numeric(0),
      UPPER = numeric(0),
      IMM = integer(0),
      MIN_TIME = numeric(0),
      OUTPUT = "BIC"
    ))
  },
  finally = print(paste0(Method, ": complete")))
}

## Special Case for D/Pre-Omicron ####
#Note that there is not enough data to include hybrid immunity in the case:
#TYPE = D and STRAIN = "Pre-Omicron" so we run this one separately
fit_rma_inf_alt = function(X, Method, multilevel = F){
  
  out.RMA = rma(VALUE, 
                sei=STDERROR, 
                mods = ~MIN_TIME * IMM, 
                data = X
  )
  
  if(multilevel){
    if(Method == "intercept only: study and estimate"){
      #Single Random effect
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMM, 
                       data = X, random = list(~1|STUDY, 
                                               ~1|effect_id)
      )
    }
    
    if(Method == "slope-intercept correlation: study and estimate"){
      #Auto-correlated random effects following exp(-(time_1 - time_2)/rho)
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMM, struct = "GEN",
                       data = X, random = list(~MIN_TIME|STUDY,
                                               ~MIN_TIME|effect_id)
      )
    }
    
    if(Method == "slope-intercept-modulator correlation: study and estimate"){
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMM, struct = "GEN",
                       data = X, random = list(~IMM + MIN_TIME|STUDY,
                                               ~IMM + MIN_TIME|effect_id)
      )
    }
    
    if(Method == "intercept only: study only"){
      #Single Random effect
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMM, 
                       data = X, random = ~1|STUDY
      )
    }
    
    if(Method == "slope-intercept correlation: study only"){
      #Auto-correlated random effects following exp(-(time_1 - time_2)/rho)
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMM, struct = "GEN",
                       data = X, random = ~MIN_TIME|STUDY
      )
    }
    
    if(Method == "slope-intercept-modulator correlation: study only"){
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMM, struct = "GEN",
                       data = X, random = ~IMM + MIN_TIME|STUDY
      )
    }
    
    if(Method == "intercept only: estimate only"){
      #Single Random effect
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMM, 
                       data = X, random = ~1|effect_id
      )
    }
    
    if(Method == "slope-intercept correlation: estimate only"){
      #Auto-correlated random effects following exp(-(time_1 - time_2)/rho)
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMM, struct = "GEN",
                       data = X, random = ~MIN_TIME|effect_id
      )
    }
    
    if(Method == "slope-intercept-modulator correlation: estimate only"){
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMM, struct = "GEN",
                       data = X, random = ~IMM + MIN_TIME|effect_id
      )
    }
  }
  .out = as.data.table(predict.rma(out.RMA))
  
  .out$WEIGHT = weights(out.RMA)
  .out = cbind(.out[,.(pred, se, ci.lb, ci.ub, WEIGHT),], 
               X[,.(VALUE, LOWER, UPPER, IMM, MIN_TIME),])
  
  .out$OUTPUT = "Predicted"
  template_a = data.table(
    expand.grid(MIN_TIME = c(2, 30),
                #STRAIN = c("Omicron", "Pre-Omicron"),
                IMM = c("V", "I")
    )
  )
  
  prediction.grid_a = model.matrix(~MIN_TIME * IMM, template_a)[,-1]
  prediction.grid_b = cbind(
    MIN_TIME = c(-28, -28),
    `IMMI` = c(0,0),
    `MIN_TIME:IMMI` = c(0, -28)
  )
  
  .out2a = as.data.table(predict.rma(out.RMA, newmods = prediction.grid_a))
  .out2b = as.data.table(predict.rma(out.RMA, newmods = prediction.grid_b, intercept = FALSE))
  
  .out2 = cbind(rbind(.out2a[,.(pred, se, ci.lb, ci.ub),], 
                      .out2b[,.(pred, se, ci.lb, ci.ub),]), 
                rbind(template_a, 
                      data.frame(
                        IMM = c("V","I")), 
                      fill = T)
  )
  
  .out2$OUTPUT = "Standards"
  
  effect.grid = cbind(
    MIN_TIME = c(0, 0, 0),
    `IMMI` = c(1,1,0),
    `MIN_TIME:IMMI` = c(2, 30, 1)
  )
  
  .out3 = as.data.table(predict.rma(out.RMA, newmods = effect.grid, intercept = FALSE))
  
  .out3 = cbind(.out3[,.(pred, se, ci.lb, ci.ub),], 
                MIN_TIME = c(2, 30, NA), IMM = "I vs V")
  .out3$OUTPUT = "Contrasts"
  
  .out4 = data.table(
    VALUE = out.RMA$fit.stats$REML,
    OUTPUT = row.names(out.RMA$fit.stats)
  )
  rbind(.out, .out2, .out3, .out4, fill = T)
  
}

# Model 4: Efficacy~TIME*Prior infecting strain ####

fit_rma_pstrain_wrap = function(X, multilevel = F){
  if(!multilevel){
    fit_rma_inf_alt(X, "single_random", multilevel)
  }
  
  .DT = data.table(method = c("intercept only: study and estimate",
                              "slope-intercept correlation: study and estimate",
                              "slope-intercept-modulator correlation: study and estimate",
                              "intercept only: study only",
                              "slope-intercept correlation: study only",
                              "slope-intercept-modulator correlation: study only",
                              "intercept only: estimate only",
                              "slope-intercept correlation: estimate only",
                              "slope-intercept-modulator correlation: estimate only"))
  .DT[,fit_rma_pstrain_error_catch(X,method, multilevel = multilevel),method]
  
}

fit_rma_pstrain_error_catch = function(X, Method, multilevel = F){
  result = tryCatch({
    fit_rma_pstrain(X, Method, multilevel = multilevel)
  },
  error = function(e){
    print(e)
    return(data.table(
      pred = numeric(0),
      se = numeric(0),
      ci.lb = numeric(0),
      ci.ub = numeric(0),
      WEIGHT = numeric(0),
      VALUE = numeric(0),
      LOWER = numeric(0),
      UPPER = numeric(0),
      IMM = integer(0),
      MIN_TIME = numeric(0),
      OUTPUT = "BIC"
    ))
  },
  finally = print(paste0(Method, ": complete")))
}

fit_rma_pstrain = function(X, Method, multilevel = F){
  out.RMA = rma(VALUE, 
                sei=STDERROR, 
                mods = ~MIN_TIME * IMMP, 
                data = X
  )
  
  if(multilevel){
    if(Method == "intercept only: study and estimate"){
      #Single Random effect
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMMP, 
                       data = X, random = list(~1|STUDY, 
                                               ~1|effect_id)
      )
    }
    
    if(Method == "slope-intercept correlation: study and estimate"){
      #Auto-correlated random effects following exp(-(time_1 - time_2)/rho)
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMMP, struct = "GEN",
                       data = X, random = list(~MIN_TIME|STUDY,
                                               ~MIN_TIME|effect_id)
      )
    }
    
    if(Method == "slope-intercept-modulator correlation: study and estimate"){
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMMP, struct = "GEN",
                       data = X, random = list(~IMMP + MIN_TIME|STUDY,
                                               ~IMMP + MIN_TIME|effect_id)
      )
    }
    
    if(Method == "intercept only: study only"){
      #Single Random effect
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMMP, 
                       data = X, random = ~1|STUDY
      )
    }
    
    if(Method == "slope-intercept correlation: study only"){
      #Auto-correlated random effects following exp(-(time_1 - time_2)/rho)
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMMP, struct = "GEN",
                       data = X, random = ~MIN_TIME|STUDY
      )
    }
    
    if(Method == "slope-intercept-modulator correlation: study only"){
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMMP, struct = "GEN",
                       data = X, random = ~IMMP + MIN_TIME|STUDY
      )
    }
    
    
    if(Method == "intercept only: estimate only"){
      #Single Random effect
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMMP, 
                       data = X, random = ~1|effect_id
      )
    }
    
    if(Method == "slope-intercept correlation: estimate only"){
      #Auto-correlated random effects following exp(-(time_1 - time_2)/rho)
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMMP, struct = "GEN",
                       data = X, random = ~MIN_TIME|effect_id
      )
    }
    
    if(Method == "slope-intercept-modulator correlation: estimate only"){
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * IMMP, struct = "GEN",
                       data = X, random = ~IMMP + MIN_TIME|effect_id
      )
    }
  }
  
  .out = as.data.table(predict.rma(out.RMA))
  
  .out$WEIGHT = weights(out.RMA)
  .out = cbind(.out[,.(pred, se, ci.lb, ci.ub, WEIGHT),], 
               X[,.(VALUE, LOWER, UPPER, IMMP, MIN_TIME),])
  
  .out$OUTPUT = "Predicted"
  template_a = data.table(
    expand.grid(MIN_TIME = c(2, 30),
                #STRAIN = c("Omicron", "Pre-Omicron"),
                IMMP = c("HPost - Omicron", "V", "IPre - Omicron", "IPost - Omicron", "HPre - Omicron")
    )
  )
  
  prediction.grid_a = model.matrix(~MIN_TIME * IMMP, template_a)[,-1]
  prediction.grid_b = cbind(
    MIN_TIME = c(-28, -28, -28, -28, -28),
    `IMMPV` = c(0,0,0,0,0),
    `IMMPIPre - Omicron` = c(0,0,0,0,0),
    `IMMPIPost - Omicron` = c(0,0,0,0,0),
    `IMMPHPre - Omicron` = c(0,0,0,0,0),
    `MIN_TIME:IMMPV` = c(0,-28,0,0,0),
    `MIN_TIME:IMMPIPre - Omicron` = c(0,0,-28,0,0),
    `MIN_TIME:IMMPIPost - Omicron` = c(0,0,0,-28,0),
    `MIN_TIME:IMMPHPre - Omicron` = c(0,0,0,0,-28)
  )
  
  .out2a = as.data.table(predict.rma(out.RMA, newmods = prediction.grid_a))
  .out2b = as.data.table(predict.rma(out.RMA, newmods = prediction.grid_b, intercept = FALSE))
  
  .out2 = cbind(rbind(.out2a[,.(pred, se, ci.lb, ci.ub),], 
                      .out2b[,.(pred, se, ci.lb, ci.ub),]), 
                rbind(template_a, 
                      data.frame(
                        IMMP = c("HPost - Omicron", "V", "IPre - Omicron", "IPost - Omicron", "HPre - Omicron")), 
                      fill = T)
  )
  
  .out2$OUTPUT = "Standards"
  
  
  effect.grid = cbind(
    MIN_TIME  = c( 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),
    `IMMPV` = c(-1, -1,  0, -1, -1,  0, -1, -1,  0, -1, -1,  0),
    `IMMPIPre - Omicron` = c( 1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0),
    `IMMPIPost - Omicron` = c( 0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0),
    `IMMPHPre - Omicron` = c( 0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0),
    `MIN_TIME:IMMPV` = c(-2,-30, -1, -2,-30, -1, -2,-30, -1, -2,-30, -1),
    `MIN_TIME:IMMPIPre - Omicron` = c( 2, 30,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0),
    `MIN_TIME:IMMPIPost - Omicron` = c( 0,  0,  0,  2, 30,  1,  0,  0,  0,  0,  0,  0),
    `MIN_TIME:IMMPHPre - Omicron` = c( 0,  0,  0,  0,  0,  0,  2, 30,  1,  0,  0,  0)
  )
  
  
  .out3 = as.data.table(predict.rma(out.RMA, newmods = effect.grid, intercept = FALSE))
  
  .out3 = cbind(.out3[,.(pred, se, ci.lb, ci.ub),], 
                MIN_TIME = c(2, 30, NA, 2, 30, NA, 2, 30, NA, 2, 30, NA), 
                IMMP = rep(c("I (Pre) vs V", "I (Post) vs V", "H (Pre) vs V", "H (Post) vs V"), each = 3))
  .out3$OUTPUT = "Contrasts"
  .out4 = data.table(
    VALUE = out.RMA$fit.stats$REML,
    OUTPUT = row.names(out.RMA$fit.stats)
  )
  rbind(.out, .out2, .out3, .out4, fill = T)
}

#Model 5: Efficacy~Time * Infection among those with booster doses ####

fit_rma_hboosted_wrap = function(X, multilevel = F){
  if(!multilevel){
    fit_rma_hboosted(X, "single_random", multilevel)
  }
  
  .DT = data.table(method = c("intercept only: study and estimate",
                              "slope-intercept correlation: study and estimate",
                              "slope-intercept-modulator correlation: study and estimate",
                              "intercept only: study only",
                              "slope-intercept correlation: study only",
                              "slope-intercept-modulator correlation: study only",
                              "intercept only: estimate only",
                              "slope-intercept correlation: estimate only",
                              "slope-intercept-modulator correlation: estimate only"))
  .DT[,fit_rma_hboosted(X,method, multilevel = multilevel),method]
  
}

fit_rma_hboosted = function(X, Method, multilevel = F){
  out.RMA = rma(VALUE, 
                sei=STDERROR, 
                mods = ~MIN_TIME * Infected, 
                data = X
  )
  
  if(multilevel){
    if(Method == "intercept only: study and estimate"){
      #Single Random effect
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * Infected, 
                       data = X, random = list(~1|STUDY, 
                                               ~1|effect_id)
      )
    }
    
    if(Method == "slope-intercept correlation: study and estimate"){
      #Auto-correlated random effects following exp(-(time_1 - time_2)/rho)
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * Infected, struct = "GEN",
                       data = X, random = list(~MIN_TIME|STUDY,
                                               ~MIN_TIME|effect_id)
      )
    }
    
    if(Method == "slope-intercept-modulator correlation: study and estimate"){
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * Infected, struct = "GEN",
                       data = X, random = list(~Infected + MIN_TIME|STUDY,
                                               ~Infected + MIN_TIME|effect_id)
      )
    }
    
    if(Method == "intercept only: study only"){
      #Single Random effect
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * Infected, 
                       data = X, random = ~1|STUDY
      )
    }
    
    if(Method == "slope-intercept correlation: study only"){
      #Auto-correlated random effects following exp(-(time_1 - time_2)/rho)
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * Infected, struct = "GEN",
                       data = X, random = ~MIN_TIME|STUDY
      )
    }
    
    if(Method == "slope-intercept-modulator correlation: study only"){
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * Infected, struct = "GEN",
                       data = X, random = ~Infected + MIN_TIME|STUDY
      )
    }
    
    if(Method == "intercept only: estimate only"){
      #Single Random effect
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * Infected, 
                       data = X, random = ~1|effect_id
      )
    }
    
    if(Method == "slope-intercept correlation: estimate only"){
      #Auto-correlated random effects following exp(-(time_1 - time_2)/rho)
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * Infected, struct = "GEN",
                       data = X, random = ~MIN_TIME|effect_id
      )
    }
    
    if(Method == "slope-intercept-modulator correlation: estimate only"){
      out.RMA = rma.mv(VALUE, 
                       V=STDERROR^2, 
                       mods = ~MIN_TIME * Infected, struct = "GEN",
                       data = X, random = ~Infected + MIN_TIME|effect_id
      )
    }
  }
  .out = as.data.table(predict.rma(out.RMA))
  
  .out$WEIGHT = weights(out.RMA)
  .out = cbind(.out[,.(pred, se, ci.lb, ci.ub, WEIGHT),], 
               X[,.(VALUE, LOWER, UPPER, Infected, MIN_TIME),])
  
  .out$OUTPUT = "Predicted"
  template_a = data.table(
    expand.grid(MIN_TIME = c(2, 30),
                #STRAIN = c("Omicron", "Pre-Omicron"),
                Infected = c(T, F)
    )
  )
  
  prediction.grid_a = model.matrix(~MIN_TIME * Infected, template_a)[,-1]
  prediction.grid_b = cbind(
    MIN_TIME = c(-28, -28),
    `InfectedTRUE` = c(0,0),
    `MIN_TIME:InfectedTRUE` = c(0, -28)
  )
  
  .out2a = as.data.table(predict.rma(out.RMA, newmods = prediction.grid_a))
  .out2b = as.data.table(predict.rma(out.RMA, newmods = prediction.grid_b, intercept = FALSE))
  
  .out2 = cbind(rbind(.out2a[,.(pred, se, ci.lb, ci.ub),], 
                      .out2b[,.(pred, se, ci.lb, ci.ub),]), 
                rbind(template_a, 
                      data.frame(
                        Infected = c(F, 
                                     T)), 
                      fill = T)
  )
  
  .out2$OUTPUT = "Standards"
  
  effect.grid = cbind(
    MIN_TIME = c(0, 0, 0),
    `InfectedTRUE` = c(1,1,0),
    `MIN_TIME:InfectedTRUE` = c(2, 30, 1)
  )
  
  .out3 = as.data.table(predict.rma(out.RMA, newmods = effect.grid, intercept = FALSE))
  
  .out3 = cbind(.out3[,.(pred, se, ci.lb, ci.ub),], MIN_TIME = c(2, 30, NA), Infected = "p-value")
  .out3$OUTPUT = "Contrasts"
  
  .out4 = data.table(
    VALUE = out.RMA$fit.stats$REML,
    OUTPUT = row.names(out.RMA$fit.stats)
  )
  rbind(.out, .out2, .out3, .out4, fill = T)
  
}

# Non-linear Model Efficacy~WANING*IMM ####
fit_rma_nonlinear_wrap = function(X, multilevel = F, decay_rate = 1/20){
  if(!multilevel){
    fit_rma_nonlinear(X, "single_random", multilevel)
  }
  
  .DT = data.table(method = c("intercept only: study and estimate",
                              "slope-intercept correlation: study and estimate",
                              #"slope-intercept-modulator correlation: study and estimate",
                              "intercept only: study only",
                              "slope-intercept correlation: study only",
                              #"slope-intercept-modulator correlation: study only",
                              "intercept only: estimate only",
                              "slope-intercept correlation: estimate only")
                              #"slope-intercept-modulator correlation: estimate only")
                   )
  .DT[,fit_rma_nonlinear_error_catch(X,method, 
                                     multilevel = multilevel, 
                                     decay_rate = decay_rate),
      method]
  
}

fit_rma_nonlinear_error_catch = function(X, Method, multilevel = F,
                                         decay_rate = 1/20){
  result = tryCatch({
    fit_rma_nonlinear(X, Method, multilevel = multilevel, 
                      decay_rate = decay_rate)
  },
  error = function(e){
    print(e)
    return(data.table(
      IMM = factor(character(0), levels = c("V", "I", "H", "B", "HB", "BvV", 
                                                     "BB", "BvB", "HvI", "HBvI", 
                                                     "Hvf", "HI", "HBvH")),
      WANING = numeric(0),
      VALUE = numeric(0),
      se = numeric(0),
      LOWER = numeric(0),
      UPPER = numeric(0),
      OUTPUT = "BIC"
    ))
  },
  finally = print(paste0(Method, ": complete")))
}

fit_rma_nonlinear = function(X, Method, multilevel = F, decay_rate = 1/20){
  out.RMA = rma(VALUE, 
                sei=STDERROR, 
                mods = ~WANING * IMM, 
                data = X[,.(VALUE, STDERROR, IMM, 
                            WANING = 1 - exp(-decay_rate * MIN_TIME))]
  )
  
  if(multilevel){
    if(Method == "intercept only: study and estimate"){
      #Single Random effect
      out.RMA = rma.mv(VALUE, 
                       V = STDERROR^2, 
                       mods = ~WANING * IMM, 
                       data = X[,.(VALUE, STDERROR, IMM, STUDY, effect_id, 
                                   WANING = 1 - exp(-decay_rate * MIN_TIME))], 
                       random = list(~1|STUDY, 
                                     ~1|effect_id)
      )
    }
    
    if(Method == "slope-intercept correlation: study and estimate"){
      #Auto-correlated random effects following exp(-(time_1 - time_2)/rho)
      out.RMA = rma.mv(VALUE, 
                       V = STDERROR^2, 
                       mods = ~WANING * IMM, 
                       struct = "GEN",
                       data = X[,.(VALUE, STDERROR, IMM, STUDY, effect_id, 
                                   WANING = 1 - exp(-decay_rate * MIN_TIME))],
                       random = list(~WANING|STUDY,
                                     ~WANING|effect_id)
      )
    }
    
    if(Method == "slope-intercept-modulator correlation: study and estimate"){
      out.RMA = rma.mv(VALUE, 
                       V = STDERROR^2, 
                       mods = ~WANING * IMM, 
                       struct = "GEN",
                       data = X[,.(VALUE, STDERROR, IMM, STUDY, effect_id, 
                                   WANING = 1 - exp(-decay_rate * MIN_TIME))],
                       random = list(~IMM + WANING|STUDY,
                                     ~IMM + WANING|effect_id)
      )
    }
    
    if(Method == "intercept only: study only"){
      #Single Random effect
      out.RMA = rma.mv(VALUE, 
                       V = STDERROR^2, 
                       mods = ~WANING * IMM,
                       data = X[,.(VALUE, STDERROR, IMM, STUDY, effect_id, 
                                   WANING = 1 - exp(-decay_rate * MIN_TIME))],
                       random = ~1|STUDY
      )
    }
    
    if(Method == "slope-intercept correlation: study only"){
      #Auto-correlated random effects following exp(-(time_1 - time_2)/rho)
      out.RMA = rma.mv(VALUE, 
                       V = STDERROR^2, 
                       mods = ~WANING * IMM, 
                       struct = "GEN",
                       data = X[,.(VALUE, STDERROR, IMM, STUDY, effect_id, 
                                   WANING = 1 - exp(-decay_rate * MIN_TIME))],
                       random = ~WANING|STUDY
      )
    }
    
    if(Method == "slope-intercept-modulator correlation: study only"){
      out.RMA = rma.mv(VALUE, 
                       V = STDERROR^2, 
                       mods = ~WANING * IMM, 
                       struct = "GEN",
                       data = X[,.(VALUE, STDERROR, IMM, STUDY, effect_id, 
                                   WANING = 1 - exp(-decay_rate * MIN_TIME))],
                       random = ~IMM + WANING|STUDY
      )
    }
    
    if(Method == "intercept only: estimate only"){
      #Single Random effect
      out.RMA = rma.mv(VALUE, 
                       V = STDERROR^2, 
                       mods = ~WANING * IMM, 
                       data = X[,.(VALUE, STDERROR, IMM, STUDY, effect_id, 
                                   WANING = 1 - exp(-decay_rate * MIN_TIME))],
                       random = ~1|effect_id
      )
    }
    
    if(Method == "slope-intercept correlation: estimate only"){
      #Auto-correlated random effects following exp(-(time_1 - time_2)/rho)
      out.RMA = rma.mv(VALUE, 
                       V = STDERROR^2, 
                       mods = ~WANING * IMM, 
                       struct = "GEN",
                       data = X[,.(VALUE, STDERROR, IMM, STUDY, effect_id, 
                                   WANING = 1 - exp(-decay_rate * MIN_TIME))],
                       random = ~WANING|effect_id
      )
    }
    
    if(Method == "slope-intercept-modulator correlation: estimate only"){
      out.RMA = rma.mv(VALUE, 
                       V = STDERROR^2, 
                       mods = ~WANING * IMM, 
                       struct = "GEN",
                       data = X[,.(VALUE, STDERROR, IMM, STUDY, effect_id, 
                                   WANING = 1 - exp(-decay_rate * MIN_TIME))],
                       random = ~IMM + WANING|effect_id
      )
    }
  }
  
  
  template = data.table(
    expand.grid(
      IMM = sort(unique(X$IMM)),
      WANING = c(0, 1)
    )
  )
  
  prediction_grid = 
    model.matrix( ~WANING * IMM, template)[,row.names(out.RMA$b)[-1]]
  
  model_prediction = as.data.table(predict.rma(out.RMA, newmods = prediction_grid))
  
  parameter_table = cbind(template, model_prediction[,.(VALUE = pred, se, 
                                                        LOWER = ci.lb, 
                                                        UPPER = ci.ub)])
  parameter_table$OUTPUT = "parameters"
  
  fit_summary = data.table(
    VALUE = out.RMA$fit.stats$REML,
    OUTPUT = row.names(out.RMA$fit.stats)
  )
  rbind(parameter_table, fit_summary, fill = T)
  
}