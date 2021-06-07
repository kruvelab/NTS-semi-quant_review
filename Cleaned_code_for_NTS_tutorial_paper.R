library(tidyverse)
source("compound_eluent.R")


regressor_pos = readRDS("ESIpos_model_191116.rds")
regressor_neg = readRDS("regressor_neg_uus.rds")

mixM = read_delim("Mix_M_dataset.csv",
                  delim = ";",
                  col_names = TRUE)

SMILES = mixM %>%
  select(SMILES) %>%
  unique()

write_delim(SMILES,
            "SMILES_EmmaM.csv",
            delim = ",")

PaDEL = read_delim("PaDEL_EmmaM.csv",
                   delim = ",",
                   col_names = TRUE)

#reading the gradient program for reversed phase LC
eluent_parameters <- read_delim('eluent.csv',
                                delim = ";",
                                col_names = TRUE)

#organic modifier muutujana javast
organic_modifier = "MeCN"

mixM = mixM %>%
  rename(RT = ret_time) %>%
  na.omit() 

mixM = mixM %>%
  mutate(
    organic_modifier = organic_modifier,
    organic = organicpercentage(eluent_parameters, RT),
    pH.aq. = pH,
    viscosity =  viscosity(organic,organic_modifier),
    surface_tension = surfacetension(organic,organic_modifier),
    polarity_index = polarityindex(organic,organic_modifier)) %>%
  mutate(NH4 = case_when(
    pH == 10.0 ~ 1,
    pH == 8 ~ 1,
    TRUE ~ 0
  ))

mixM = mixM %>%
  left_join(PaDEL)


#------predicting IEs-------------
prediction_set_model = mixM %>%
  na.omit() %>%
  mutate(logIE_pred_pos = 0,
         logIE_pred_neg = 0)
prediction_pos =  predict(regressor_pos, newdata = prediction_set_model, predict.all = TRUE)
prediction_neg =  predict(regressor_neg, newdata = prediction_set_model, predict.all = TRUE)
prediction_pos = prediction_pos$aggregate
prediction_neg = prediction_neg$aggregate
prediction_set_model <- prediction_set_model %>%
  mutate(logIE_pred_pos = prediction_pos,
         logIE_pred_neg = prediction_neg) %>%
  #filter(RF > 2e18) %>%
  select(SMILES,logIE_pred_pos, logIE_pred_neg, everything())

#linear regression for the calibration compounds
lin_fit_logRF_pos <- lm(log10(slope) ~ logIE_pred_pos, 
                        data = prediction_set_model %>% filter(mode == "positive"))
lin_fit_logRF_neg <- lm(log10(slope) ~ logIE_pred_neg, 
                        data = prediction_set_model %>% filter(mode == "negative"))

prediction_set_model = prediction_set_model %>%
  mutate(logRF_pred_pos = lin_fit_logRF_pos$coefficients[2]*logIE_pred_pos + lin_fit_logRF_pos$coefficients[1],
         logRF_pred_neg = lin_fit_logRF_neg$coefficients[2]*logIE_pred_neg + lin_fit_logRF_neg$coefficients[1]) %>%
  mutate(c_pred_IE_pos = Peak_Area_IC/(10^logRF_pred_pos),
         c_pred_IE_neg = Peak_Area_IC/(10^logRF_pred_neg))

prediction_set_model = prediction_set_model %>%
  mutate(logIE_pred = case_when(
    mode == "positive" ~ logIE_pred_pos,
    TRUE ~ logIE_pred_neg
  )) %>%
  mutate(c_pred_IE = case_when(
    mode == "positive" ~ c_pred_IE_pos,
    TRUE ~ c_pred_IE_neg
  )) %>%
  mutate(error = case_when(
    conc_M > c_pred_IE ~ conc_M/c_pred_IE,
    TRUE ~ c_pred_IE/conc_M))

for (this_mode in levels(as.factor(prediction_set_model$mode))) {
  print(this_mode)
  for (this_pH in levels(as.factor(prediction_set_model$pH))) {
    data = prediction_set_model %>%
      filter(mode == this_mode & pH == this_pH)
    print(this_pH)
    print(mean(data$error))
  }
}

#prodictions




sum_pred = prediction_set_model %>%
  filter(mode == "positive" & pH == 2.7) %>%
  bind_rows(prediction_set_model %>%
              filter(mode == "negative" & pH == 10))

sum_pred = sum_pred %>%
  group_by(Compound_name, conc_M) %>%
  summarize(c_pred_IE = mean(c_pred_IE),
            count = n()) %>%
  ungroup() %>%
  mutate(mode = "mean") %>%
  filter(count > 1) %>%
  select(-count)

sum_pred = sum_pred %>%
  mutate(error = case_when(
    conc_M > c_pred_IE ~ conc_M/c_pred_IE,
    TRUE ~ c_pred_IE/conc_M))

mean(sum_pred$error)


choice_pred = prediction_set_model %>%
  filter(mode == "positive" & pH == 2.7) %>%
  bind_rows(prediction_set_model %>%
              filter(mode == "negative" & pH == 10))

choice_pred = choice_pred %>%
  group_by(Compound_name, conc_M) %>%
  mutate(c_pred_IE = case_when(
    Peak_Area == max(Peak_Area) ~ c_pred_IE,
    TRUE ~ 0),
    count = n()) %>%
  ungroup() %>%
  mutate(mode = "choice") %>%
  filter(c_pred_IE != 0) %>%
  filter(count > 1) %>%
  select(-count) 

choice_pred = choice_pred %>%
  mutate(error = case_when(
    conc_M > c_pred_IE ~ conc_M/c_pred_IE,
    TRUE ~ c_pred_IE/conc_M))


all_data = prediction_set_model %>%
  bind_rows(sum_pred) %>%
  bind_rows(choice_pred)




summary_df = all_data %>%
  filter((mode == "positive" & pH == 2.7) |
           (mode == "negative" & pH == 10) | 
           (mode == "mean") |
           (mode == "choice")) %>%
  group_by(Compound_name, conc_M) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  select(SMILES, count, everything())


summary_df %>%
  filter(count > 2) %>%
  group_by(mode) %>%
  summarize(mean_error = mean(error),
            median_error = median(error)) %>%
  ungroup()