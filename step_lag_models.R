library(tidyverse)
#library(Delgado)
library(lmerTest)
library(broom.mixed)
library(pscl);library(boot)

# Create the lag Variables

steps_active_df <- as.data.frame(steps_active_list)
steps_active_list_w_lag <-  #not using the lag of the steps - that was an error - 10-24-24
steps_active_df %>% 
  group_by(participant) %>% 
  mutate(
    step_lag_1 = lag(steps,1L),
    step_lag_1_scale = lag(steps,1L) %>% scale(.),
    step_lag_2_scale = lag(steps,2L) %>% scale(.),
    step_lag_3_scale = lag(steps,3L) %>% scale(.),
    step_lag_4_scale = lag(steps,4L) %>% scale(.),
    depressed_lag_1 = lag(depressed,1L)
  ) %>% 
  ungroup()


#### Within participant scaling ####

## function
within_lag.mod.list <- mapply(
  FUN = function(dv,lag_time){
    formulaz <- paste0("lag(",dv,",",lag_time,"L) ~ steps_scale_within + day + (1|participant)") %>% as.formula(.)
    
    fit <- lmer(formulaz,
                      data = steps_active_list_w_lag %>% 
                        filter(sub_group=="MA"),
                      REML = T)
    
    param.tbl <- tidy(fit,effects = "fixed",conf.int=T) %>% 
      mutate(
        stars = sapply(p.value,BiostatsALL::getstar)
      )
    
    results.list <- 
      list(
        'formula' = formulaz,
        'fit' = fit,
        'param.tbl' = param.tbl
      )
    
    return(results.list)
  },
  dv = c(
    rep("depressed",7),
    rep("distressed",7),
    rep("anxious",7)
    ),
  lag_time = rep(1:7,3),SIMPLIFY = F
) %>% 
  setNames(.,c(
    paste0("depressed_lag_",1:7),
    paste0("distressed_lag_",1:7),
    paste0("anxious_lag_",1:7)
  ))

all_lag_param.tbl <- lapply(within_lag.mod.list,function(x){x[["param.tbl"]]}) %>% 
  bind_rows(.,.id = "DV") %>% 
  mutate(
    term = factor(term, levels = c("(Intercept)","steps_scale_within","day"),
                  labels = c("(Intercept)","Within\nScaled Steps","Day")),
    lag_time = rep(rep(1:7,each=3),3) %>% factor(.,levels = 1:7,labels = paste0(1:7,"-period lag")),
    Outcome = rep(c("depressed","distressed","anxious"),each = 21) %>% 
      factor(.,
             levels = c("depressed","distressed","anxious"),
             labels = c("Depressed","Distressed","Anxious")
               )
  ) %>% 
  select(
    Outcome,lag_time,everything()
  )

all_lag_param.tbl

#write_csv(all_lag_param.tbl, '/Users/jackiebeltran/Documents/jmb_digital/results/within_scaling_lag_model.csv')

# reorder data table levels

all_lag_param.tbl$Outcome <- factor(all_lag_param.tbl$Outcome,
                                    levels = c("Anxious", "Distressed", "Depressed"))


## Plot ##
all_lag_MA_mod_parm_plot <- 
all_lag_param.tbl %>% 
  ggplot(aes(x = estimate, y = term,color = lag_time)) + 
  geom_point(size = 1.5) + 
  geom_line(aes(group = term)) +
  facet_wrap(~Outcome,ncol = 4) + 
  ggpubr::theme_pubclean() + 
  labs(
    x = "Estimate",
    y = "Regression Term",
    color = "Lag Time"
  ) +
  
  theme(
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text =  element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text = element_text(size = 14)
  )

all_lag_MA_mod_parm_plot

# ggsave(
#   plot = all_lag_MA_mod_parm_plot,
#   filename = "Plots/All Outcomes MA LMM Lag models.jpeg",
#   height = 4,width = 5.75
#   )

ema_labels <- c("Anxious" = "Anxiety", "Distressed" = "Distress", "Depressed" = "Depression")


# show only the step estimates
within_step_est_lag_mod_plot <- 
all_lag_param.tbl %>% 
  filter(term=="Within\nScaled Steps") %>% 
  ggplot(aes(x = estimate, y = lag_time)) + 
  geom_text(aes(label = stars,x = (conf.high*1.15),size = 20)) + 
  geom_errorbar(aes(xmin = conf.low,xmax = conf.high),width = 0.2) + 
  geom_point(aes(color= Outcome),size = 2.5,show.legend = F) + 
  geom_vline(xintercept = 0,linetype='dashed') + 
  facet_wrap(~Outcome,ncol = 4, labeller = labeller(Outcome = ema_labels)) + 
  ggpubr::theme_pubclean() + 
  labs(
    x = "Steps Estimate (within-participant scaling)",
    y = "Lag Period",
    color = "Lag Time"
  ) +
  scale_color_manual(values = c("#CC6633", "#669999", "#996699")) + 
  theme(
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),
    legend.position = "none"
  )


within_step_est_lag_mod_plot

# ggsave(
#   plot = step_est_lag_mod_plot,
#   filename = "Plots/All Outcomes MA LMM Lag models - Step Effect Only.jpeg",
#   height = 4,width = 8
# )

# Zip model# Zip model# Zip model
# zip_lag1 <-
#   zeroinfl(
#     depressed ~ steps_scale_within + day + sub_group | participant,
#     dist = "poisson",
#     link = "probit",
#     data = steps_active_list_w_lag
#   )

















































#### Between participant scaling ####

## function
bw_lag.mod.list <- mapply(
  FUN = function(dv,lag_time){
    formulaz <- paste0("lag(",dv,",",lag_time,"L) ~ steps_scale + day + (1|participant)") %>% as.formula(.)
    
    fit <- lmer(formulaz,
                data = steps_active_list_w_lag %>% 
                  filter(sub_group=="MA"),
                REML = T)
    
    param.tbl <- tidy(fit,effects = "fixed",conf.int=T) %>% 
      mutate(
        stars = sapply(p.value,BiostatsALL::getstar)
      )
    
    results.list <- 
      list(
        'formula' = formulaz,
        'fit' = fit,
        'param.tbl' = param.tbl
      )
    
    return(results.list)
  },
  dv = c(
    rep("depressed",7),
    rep("distressed",7),
    rep("anxious",7)
  ),
  lag_time = rep(1:7,3),SIMPLIFY = F
) %>% 
  setNames(.,c(
    paste0("depressed_lag_",1:7),
    paste0("distressed_lag_",1:7),
    paste0("anxious_lag_",1:7)
  ))

all_lag_param.tbl <- lapply(bw_lag.mod.list,function(x){x[["param.tbl"]]}) %>% 
  bind_rows(.,.id = "DV") %>% 
  mutate(
    term = factor(term, levels = c("(Intercept)","steps_scale","day"),
                  labels = c("(Intercept)","Between\nScaled Steps","Day")),
    lag_time = rep(rep(1:7,each=3),3) %>% factor(.,levels = 1:7,labels = paste0(1:7,"-period lag")),
    Outcome = rep(c("depressed","distressed","anxious"),each = 21) %>% 
      factor(.,
             levels = c("depressed","distressed","anxious"),
             labels = c("Depressed","Distressed","Anxious")
      )
  ) %>% 
  select(
    Outcome,lag_time,everything()
  )

all_lag_param.tbl

write_csv(all_lag_param.tbl, '/Users/jackiebeltran/Documents/jmb_digital/results/between_scaling_lag_model.csv')

# reorder data table levels

all_lag_param.tbl$Outcome <- factor(all_lag_param.tbl$Outcome,
                                    levels = c("Anxious", "Distressed", "Depressed"))


## Plot ##
all_lag_MA_mod_parm_plot <- 
  all_lag_param.tbl %>% 
  ggplot(aes(x = estimate, y = term,color = lag_time)) + 
  geom_point(size = 1.5) + 
  geom_line(aes(group = term)) +
  facet_wrap(~Outcome,ncol = 4) + 
  ggpubr::theme_pubclean() + 
  labs(
    x = "Estimate",
    y = "Regression Term",
    color = "Lag Time"
  ) +
  
  theme(
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text =  element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text = element_text(size = 14)
  )

all_lag_MA_mod_parm_plot

# ggsave(
#   plot = all_lag_MA_mod_parm_plot,
#   filename = "Plots/All Outcomes MA LMM Lag models.jpeg",
#   height = 4,width = 5.75
# )


# show only the step estimates
step_est_lag_mod_plot <- 
  all_lag_param.tbl %>% 
  filter(term=="Between\nScaled Steps") %>% 
  ggplot(aes(x = estimate, y = lag_time)) + 
  geom_text(aes(label = stars,x = (conf.high*1.15),size = 20)) + 
  geom_errorbar(aes(xmin = conf.low,xmax = conf.high),width = 0.2) + 
  geom_point(aes(color= Outcome),size = 2.5,show.legend = F) + 
  geom_vline(xintercept = 0,linetype='dashed') + 
  facet_wrap(~Outcome,ncol = 4) + 
  ggpubr::theme_pubclean() + 
  labs(
    x = "Steps Estimate (between-participant scaling)",
    y = "Lag Period",
    color = "Lag Time"
  ) +
  scale_color_manual(values = c("#CC6633", "#669999", "#996699")) + 
  theme(
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14),
    legend.position = "none"
  )


step_est_lag_mod_plot

# ggsave(
#   plot = step_est_lag_mod_plot,
#   filename = "Plots/All Outcomes MA LMM Lag models - Step Effect Only.jpeg",
#   height = 4,width = 8
# )

# Zip model# Zip model# Zip model
# zip_lag1 <-
#   zeroinfl(
#     depressed ~ steps_scale_within + day + sub_group | participant,
#     dist = "poisson",
#     link = "probit",
#     data = steps_active_list_w_lag
#   )















































