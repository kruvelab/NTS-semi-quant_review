library(caret)
library(enviPat)
library(plotly)
library(RRF)
library(rcdk)
library(tidyverse)
library(webchem)

fn_isotope_distribution <- function(smiles){
  molecule <- parse.smiles(smiles)[[1]]
  formula <- get.mol2formula(molecule,charge=0)
  formula <- formula@string
  
  data(isotopes)
  pattern<-isopattern(isotopes,
                      formula,
                      threshold=0.1,
                      plotit=FALSE,
                      charge=FALSE,
                      algo=1)
  isotopes <- as.data.frame(pattern[[1]])
  isotope_dist <- as.numeric(sum(isotopes$abundance))
  return(isotope_dist)
}

fn_viscosity <- function(organic,organic_modifier){
  viscosity <- case_when(
    organic_modifier == "MeCN" ~ (-0.000103849885417527)*organic^2+0.00435719229180079*organic+0.884232851261593,
    organic_modifier == "MeOH" ~ (-0.00035908)*organic^2+0.031972067*organic+0.90273943)
  return(viscosity)
}

fn_surface_tension <- function(organic,organic_modifier){
  surface_tension <- case_when(
    organic_modifier == "MeCN" ~ 71.76-2.906*71.76*(organic/100)+(7.138*27.86+2.906*71.76-71.76)*(organic/100)^2+(27.86-7.138*27.86)*(organic/100)^3,
    organic_modifier == "MeOH" ~ 71.76-2.245*71.76*(organic/100)+(5.625*22.12+2.245*71.76-71.76)*(organic/100)^2+(22.12-5.625*22.12)*(organic/100)^3)
  return(surface_tension)
}

fn_polarity_index <- function(organic,organic_modifier){
  polarity_index <- case_when(
    organic_modifier == "MeCN" ~ (organic/100)*5.1+((100-organic)/100)*10.2,
    organic_modifier == "MeOH" ~ (organic/100)*5.1+((100-organic)/100)*10.2)
  return(polarity_index)
}

fn_organic_percentage <- function(eluent_parameters,ret_time){
  ApproxFun <- approxfun(x = eluent_parameters$time, y = eluent_parameters$B)
  organic <- ApproxFun(ret_time)
  return(organic)
}

assigning_closest_cal_comp_RF <- function(xRT, cal_compounds) {
  RF_cal <- cal_compounds %>% slice(which.min(abs(xRT - ret_time))) %>%select(RF)
  print(unlist(RF_cal))
}

assigning_closest_cal_comp <- function(xRT, cal_compounds) {
  Comp_cal <- cal_compounds %>% slice(which.min(abs(xRT - ret_time))) %>%select(Compound)
  print(unlist(Comp_cal))
}

assigning_closest_cal_comp_RT <- function(xRT, cal_compounds) {
  Comp_cal <- cal_compounds %>% slice(which.min(abs(xRT - ret_time))) %>%select(ret_time)
  print(unlist(Comp_cal))
}

SB1 <- read_delim('Quantem_SB1_w_compname.csv',
                  delim = ';',
                  col_names = TRUE)

SB1 <- SB1 %>%
  group_by(SMILES) %>%
  mutate(IC = fn_isotope_distribution(SMILES)) %>%
  ungroup()


SB1 <- SB1 %>%
  mutate(signal = signal*IC,
         
         RF = signal / concentration)

SA <- SB1 %>%
  filter(`SB/SA comp.` == "SA")

similarity_matrix <- read_delim('similarity_matrix.csv',
                                delim = ";",
                                col_names = TRUE)

similarity_matrix <- similarity_matrix %>%
 
  rename(SB_compound = Compound) %>%
  rename(Compound = Similar) %>%
  group_by(SB_compound) %>%
  slice(which.max(similarity_per_cent)) %>%
  ungroup() %>%
  left_join(SA) %>%
  select(SB_compound, Compound, RF, ret_time, similarity_per_cent, Maximum_expected_error) %>%
  rename(Similar = Compound) %>%
  rename(Compound = SB_compound,
         RF_similar_compound = RF,
         ret_time_similar_compound = ret_time)

SB1 <- SB1 %>%
  left_join(similarity_matrix) %>%
  mutate(c_pred_similar_compound = signal / RF_similar_compound)

cal_Comp_closeRT <- c()
for(i in 1:69){
 cal_Comp_closeRT <- c(cal_Comp_closeRT, assigning_closest_cal_comp(SB1[i,]$ret_time, SA))
}

cal_RF_closeRT <- c()
for(i in 1:69){
  cal_RF_closeRT <- c(cal_RF_closeRT, assigning_closest_cal_comp_RF(SB1[i,]$ret_time, SA))
}

cal_RT_closeRT <- c()
for(i in 1:69){
  cal_RT_closeRT <- c(cal_RT_closeRT, assigning_closest_cal_comp_RT(SB1[i,]$ret_time, SA))
}

SB1 <- data.frame(SB1, cal_Comp_closeRT, cal_RF_closeRT, cal_RT_closeRT)

SB1 <- SB1 %>%
  mutate(c_pred_close_ret_time = signal / cal_RF_closeRT)

regressor_pos <- read_rds("ESIpos_model_191116.rds")
descs_pos <-  read_rds("ESIpos_model_descs_191116.rds")

Padel_data <-  read_delim('descs200629.csv',
                          delim = ",",
                          col_names = TRUE,
                          trim_ws = TRUE) 

eluent_parameters <- read_delim('eluent.csv',
                                delim = ";",
                                col_names = TRUE)

organic_modifier <- "MeCN"
pH <- 2.7
NH4 <- 0

SB1 <- SB1 %>%
  mutate(
    organic_modifier = organic_modifier,
    organic = fn_organic_percentage(eluent_parameters, ret_time),
    pH.aq. = pH,
    NH4 = NH4,
    viscosity =  fn_viscosity(organic,organic_modifier),
    surface_tension = fn_surface_tension(organic,organic_modifier),
    polarity_index = fn_polarity_index(organic,organic_modifier)) 

prediction_set_model_pos <- SB1 %>%
  left_join(Padel_data) %>%
  select(Compound, SMILES, descs_pos) %>%
  na.omit() %>%
  mutate(logIE_pred = 0)
prediction <-  predict(regressor_pos, newdata = prediction_set_model_pos, predict.all = TRUE)
prediction <- prediction$aggregate
prediction_set_model_pos <- prediction_set_model_pos %>%
  mutate(logIE_pred = prediction) %>%
  select(Compound, SMILES,logIE_pred)

SB1 <- SB1 %>%
  left_join(prediction_set_model_pos) 

lin_fit_logRF <- lm(log(RF, 10) ~ logIE_pred, 
                    data = SB1 %>% filter(SB.SA.comp. == "SA") 
                    %>% filter(Compound != "Butocarboxim [M+NH4]+" &
                               Compound != "Clarithromycin [M+Na]+" &
                               Compound != "Ivermectin [M+NH4]+" &
                               Compound != "Nigericin [M+NH4]+" &
                               Compound != "Simvastatin [M+Na]+" &
                               Compound != "Simvastatin [M+NH4]+" &
                               Compound != "Sucralose [M+Na]+" &
                               Compound != "Sucralose [M+NH4]+"))

SB1 <- SB1 %>%
  mutate(logRF_pred = lin_fit_logRF$coefficients[2]*logIE_pred + lin_fit_logRF$coefficients[1],
         c_pred_IE = signal/(10^logRF_pred))

similarity_matrix <- read_delim('similarity_matrix.csv',
                                delim = ";",
                                col_names = TRUE)

similarity_matrix <- similarity_matrix %>%
  rename(SB_compound = Compound) %>%
  rename(Compound = Parent_compound) %>%
  group_by(SB_compound) %>%
  slice(which.max(similarity_per_cent)) %>%
  ungroup() %>%
  left_join(SA) %>%
  select(SB_compound, Compound, RF, ret_time, similarity_per_cent, Maximum_expected_error) %>%
  rename(Parent_compound = Compound) %>%
  rename(Compound = SB_compound,
         RF_parent_compound = RF,
         ret_time_parent_compound = ret_time)

SB1 <- SB1 %>%
  left_join(similarity_matrix) %>%
  mutate(c_pred_parent_compound = signal / RF_parent_compound)


SB1 <- SB1%>%
  mutate(similar_c_comparison = case_when(concentration < c_pred_similar_compound ~ c_pred_similar_compound / concentration,
                                          concentration > c_pred_similar_compound ~ concentration / c_pred_similar_compound),
         closeRT_c_comparison = case_when(concentration < c_pred_close_ret_time ~ c_pred_close_ret_time / concentration,
                                          concentration > c_pred_close_ret_time ~ concentration / c_pred_close_ret_time),
         IE_c_comparison = case_when(concentration < c_pred_IE ~ c_pred_IE / concentration,
                                     concentration > c_pred_IE ~ concentration / c_pred_IE),
         parent_c_comparison = case_when(concentration < c_pred_parent_compound ~ c_pred_parent_compound / concentration,
                                         concentration > c_pred_parent_compound ~ concentration / c_pred_parent_compound))

all_c_comparison_SB1 <- SB1%>%
  select(Compound, concentration, SB.SA.comp., c_pred_similar_compound,
         c_pred_close_ret_time, c_pred_IE, c_pred_parent_compound,
         similar_c_comparison, closeRT_c_comparison, IE_c_comparison, 
         parent_c_comparison, ret_time, RF, similarity_per_cent)

all_c_comparison_SB1_error <- gather(data = all_c_comparison_SB1, 
                                     key = "SQ_approach", 
                                     value = error, similar_c_comparison, 
                                     closeRT_c_comparison, IE_c_comparison, 
                                     parent_c_comparison) %>%
  select(Compound, concentration, SB.SA.comp., ret_time, RF, error, SQ_approach, similarity_per_cent)

all_c_comparison_SB1_pred_c <- gather(data = all_c_comparison_SB1, 
                                      key = "SQ_approach", 
                                      value = pred_conc, c_pred_similar_compound,
                                      c_pred_close_ret_time, c_pred_IE, c_pred_parent_compound, similarity_per_cent) %>%
  mutate(SQ_approach = case_when(
    SQ_approach == "c_pred_similar_compound" ~ "similar_c_comparison",
    SQ_approach == "c_pred_close_ret_time" ~ "closeRT_c_comparison",
    SQ_approach == "c_pred_IE" ~ "IE_c_comparison",
    SQ_approach == "c_pred_parent_compound" ~ "parent_c_comparison"
  )) %>%
  select(Compound, concentration, SB.SA.comp., ret_time, RF, pred_conc, SQ_approach)

all_c_comparison_SB1 <- all_c_comparison_SB1_error %>%
  left_join(all_c_comparison_SB1_pred_c)

all_c_comparison_SB1 <- all_c_comparison_SB1 %>%
  filter(SB.SA.comp. == "SB")%>%
  filter(error !="NA")

all_c_comparison_SB1 <- all_c_comparison_SB1 %>%
  mutate(sample = "SB1")

all_c_comparison_SB1 <- all_c_comparison_SB1 %>%
  select(Compound, concentration, ret_time, RF, similarity_per_cent, SQ_approach, error,
         pred_conc, sample)

error_stat_SB1 = all_c_comparison_SB1 %>%
  group_by(SQ_approach) %>%
  summarise(
    mean_error_c = mean(error),
    median_error_c = median(error),
    max_error_c = max(error),
    quantile_error_c = quantile(error, probs = c(0.95)),
    n_dp = length(error),
    n_dp_less_than_ten = length(error[error<10 & error>1]),
    percentage_less_than_ten = (length(error[error<10 & error>1]))/(length(error)))%>%
  mutate(sample = "SB1") %>%
  ungroup()
