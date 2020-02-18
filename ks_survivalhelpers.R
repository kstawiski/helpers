suppressMessages(library(rmdformats))
suppressMessages(library(dplyr))
suppressMessages(library(plyr))
suppressMessages(library(tidyverse))
suppressMessages(library(reshape2))
suppressMessages(library(xlsx))
suppressMessages(library(psych))
suppressMessages(library(knitr))
suppressMessages(library(gplots))
suppressMessages(library(car))
suppressMessages(library(kableExtra))
suppressMessages(library(RVAideMemoire))
#suppressMessages(library(rms))
suppressMessages(library(survival))
suppressMessages(library(cmprsk))
# install.packages("biostatUZH", repos="http://R-Forge.R-project.org")
suppressMessages(library(biostatUZH))
suppressMessages(library(cr17))
suppressMessages(library(ggpubr))
suppressMessages(library(crrstep))
suppressMessages(library(mice))
suppressMessages(library(PairedData))
suppressMessages(library(Hmisc))
suppressMessages(library(corrplot))
suppressMessages(library(profileR))
suppressMessages(library(bestNormalize))
suppressMessages(library(mgcv))
suppressMessages(library(lm.beta))
suppressMessages(library(olsrr))
suppressMessages(library(inflection))
suppressMessages(library(data.table))
#suppressMessages(library(usdm))
suppressMessages(library(survminer))
# devtools::install_github("zabore/ezfun")
suppressMessages(library(ezfun))
library(tibble)
library(janitor)


# Funkcje generyczne


ks.competingrisk.univariate.nominal = function(daneOS, variable, variable_name, survival_name = "OS", time_name = "Months", our_cencode = "Alive", our_failcode = "Related", adjust_to = NA) {
  
  tempdaneOS = dplyr::select(daneOS, paste0(survival_name,c("","time")))
  colnames(tempdaneOS) = c("Event","Time") # czyli OS -> OS i OStime
  #tempdaneOS$Time = as.numeric(as.character(tempdaneOS$Time))
  tempdaneOS = cbind(tempdaneOS, data.frame(selected = as.factor(variable), adjust_to = adjust_to))
  tempdaneOS$selected[tempdaneOS$selected == ""] = NA
  tempdaneOS$selected = as.factor(as.character(tempdaneOS$selected))
  tempdaneOS = tempdaneOS[complete.cases(tempdaneOS),]
  resCumIncByDis <- cuminc(ftime = tempdaneOS$Time,  # failure time variable
                           fstatus = tempdaneOS$Event,  # variable with distinct codes for different causes of failure
                           group   = tempdaneOS$selected,  # estimates will calculated within groups
                           ## strata  = ,  # Tests will be stratified on this variable.
                           rho     = 0, # Power of the weight function used in the tests.
                           cencode = our_cencode,
                           ## subset = ,
                           ## na.action = na.omit
  )
  plot1 = ggcompetingrisks(resCumIncByDis, ggtheme = theme_bw(), ylab = paste0(survival_name," probability"), xlab = time_name, main = paste0(survival_name, " vs. ", variable_name), legend.title = "", multiple_panels = F)
  resCumIncByDis
  plot2 = ggcompetingrisks(resCumIncByDis, ggtheme = theme_bw(), ylab = paste0(survival_name," probability"), xlab = time_name, main = paste0(survival_name, " vs. ", variable_name), legend.title = "")
  
  if(!is.na(adjust_to)) {
    covs1 <- model.matrix(~ as.factor(selected) + as.numeric(adjust_to), data = tempdaneOS)[, -1]
  } else {
    covs1 <- model.matrix(~ as.factor(selected), data = tempdaneOS)[, -1] }
  shr_fit <- 
    crr(
      ftime = tempdaneOS$Time,
      fstatus = tempdaneOS$Event,
      cov1 = covs1,
      failcode = our_failcode,
      cencode = our_cencode
    )
  shr_fit
  
  tempdaneOS2 = as.data.frame(tempdaneOS)
  tempsurv = survfit(Surv(Time,Event == "Related") ~ as.factor(selected), data=tempdaneOS2, type = "kaplan-meier", error = "greenwood", conf.type = "log-log")
  names(tempsurv$strata) = levels(tempdaneOS$selected)
  plot_kaplan = ggsurvplot(
    fit = tempsurv,data=tempdaneOS2,
    xlab = time_name, 
    ylab = paste0(survival_name," probability"), 
    risk.table = "abs_pct",  # absolute number and percentage at risk.
    risk.table.y.text.col = T,# colour risk table text annotations.
    risk.table.y.text = T,# show bars instead of names in text annotations
    # in legend of risk table.
    ncensor.plot = F,      # plot the number of censored subjects at time t
    surv.median.line = "hv",  # add the median survival pointer.
  )
  plot_kaplan
  
  median_survival = surv_median(tempsurv)
  
  
  res = cbind("Parametr" = c(paste0(variable_name, " (", levels(tempdaneOS$selected)[2:length(levels(tempdaneOS$selected))], " vs. ", levels(tempdaneOS$selected)[1], ")"),"Adjustment")
              ,ezfun::mvcrrres(shr_fit), "Gray's test (p-value)" = resCumIncByDis$Tests[1, 2])
  list(result_table = res[-nrow(res),], model = shr_fit, cuminc = resCumIncByDis, plot1 = plot1, plot2 = plot2, plot_kaplan = plot_kaplan, median_surv = median_survival, used_data = tempdaneOS)
  
}

ks.competingrisk.univariate.continous = function(daneOS, variable, variable_name, survival_name = "OS", time_name = "Months", our_cencode = "Alive", our_failcode = "Related", adjust_to = daneOS$CzasDoRT) {
  tempdaneOS = dplyr::select(daneOS, paste0(survival_name,c("","time")))
  colnames(tempdaneOS) = c("Event","Time") # czyli OS -> OS i OStime
  tempdaneOS$Time = as.numeric(as.character(tempdaneOS$Time))
  tempdaneOS = cbind(tempdaneOS, data.frame(selected = as.character(variable), adjust_to2 = adjust_to))
  tempdaneOS$selected[tempdaneOS$selected == ""] = NA
  tempdaneOS$selected = as.numeric(tempdaneOS$selected)
  tempdaneOS = tempdaneOS[complete.cases(tempdaneOS),]
  
  if(!is.na(adjust_to)) {
    covs1 <- model.matrix(~ as.numeric(selected) + as.numeric(adjust_to2), data = tempdaneOS)[, -1]
  } else {
    covs1 <- model.matrix(~ as.numeric(selected), data = tempdaneOS)[, -1] }
  shr_fit <- 
    crr(
      ftime = tempdaneOS$Time,
      fstatus = tempdaneOS$Event,
      cov1 = covs1,
      cencode = our_cencode,
      failcode = our_failcode
    )
  shr_fit
  
  
  
  res = cbind("Parametr" = variable_name,ezfun::mvcrrres(shr_fit), "Gray's test (p-value)" = "NA (continous variable)")
  list(result_table = res[-nrow(res),], model = shr_fit, used_data = tempdaneOS)
}

ks.competingrisk.univariate = function(daneOS, nominal_variables = NA, continous_variables = NA, nominal_variables_names = nominal_variables, continous_variables_names = continous_variables, ...)
{
  univariate = list()
  if(!is.na(nominal_variables))
  {
    for(i in 1:length(nominal_variables))
    {
      tryCatch({
        univariate[[nominal_variables[i]]] = ks.competingrisk.univariate.nominal(daneOS, as.factor(unlist(dplyr::select(daneOS, nominal_variables[i])[,1])), nominal_variables_names[i], ...)})
    }
  }
  if(!is.na(continous_variables))
  {
    for(i in 1:length(continous_variables))
    {
      tryCatch({univariate[[continous_variables[i]]] = ks.competingrisk.univariate.continous(daneOS, as.numeric(unlist(dplyr::select(daneOS, continous_variables[i])[,1])), continous_variables_names[i], ...)})
    }
  }
  return(univariate)
}