library(tidyverse)
library(lmerTest)
library(broom.mixed)
library(pscl);library(boot)

# Create the lag Variables

screen_active_df <- as.data.frame(screen_active_list)

#### Within participant scaling ####

## function
within_lag.mod.list <- mapply(
  FUN = function(dv,lag_time){
    formulaz <- paste0("lag(",dv,",",lag_time,"L) ~ screen_scale_within + day + (1|participant)") %>% as.formula(.)
    
    fit <- lmer(formulaz,
                data = screen_active_df %>% 
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

# get results table
all_lag_param.tbl <- lapply(within_lag.mod.list,function(x){x[["param.tbl"]]}) %>% 
  bind_rows(.,.id = "DV") %>% 
  mutate(
    term = factor(term, levels = c("(Intercept)","screen_scale_within","day"),
                  labels = c("(Intercept)","Within\nScaled Screentime","Day")),
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

write_csv(all_lag_param.tbl, '/Users/jackiebeltran/Documents/jmb_digital/results/within_scaling_screen_lag_model.csv')

# reorder data table levels

all_lag_param.tbl$Outcome <- factor(all_lag_param.tbl$Outcome,
                                    levels = c("Anxious", "Distressed", "Depressed"))


ema_labels <- c("Anxious" = "Anxiety", "Distressed" = "Distress", "Depressed" = "Depression")


## Plot ##
all_lag_MA_mod_parm_plot <- 
  all_lag_param.tbl %>% 
  ggplot(aes(x = estimate, y = term,color = lag_time)) + 
  geom_point(size = 1.5) + 
  geom_line(aes(group = term)) +
  facet_wrap(~Outcome,ncol = 4, labeller = labeller(Outcome = ema_labels)) + 
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

# show only the step estimates
screen_est_lag_mod_plot <- 
  all_lag_param.tbl %>% 
  filter(term=="Within\nScaled Screentime") %>% 
  ggplot(aes(x = estimate, y = lag_time)) + 
  geom_text(aes(label = stars,x = (conf.high*1.15),size = 20)) + 
  geom_errorbar(aes(xmin = conf.low,xmax = conf.high),width = 0.2) + 
  geom_point(aes(color= Outcome),size = 2.5,show.legend = F) + 
  geom_vline(xintercept = 0,linetype='dashed') + 
  facet_wrap(~Outcome,ncol = 4, labeller = labeller(Outcome = ema_labels)) + 
  ggpubr::theme_pubclean() + 
  labs(
    x = "Screentime Estimate (within-participant scaling)",
    y = "Lag Period",
    color = "Lag Time"
  ) +
  scale_color_manual(values = c("#CC6633", "#669999", "#996699")) + 
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 18),
    legend.position = "none"
  )


screen_est_lag_mod_plot


#### Between participant scaling ####

## function
between_lag.mod.list <- mapply(
  FUN = function(dv,lag_time){
    formulaz <- paste0("lag(",dv,",",lag_time,"L) ~ screen_scale + day + (1|participant)") %>% as.formula(.)
    
    fit <- lmer(formulaz,
                data = screen_active_df %>% 
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

# get results table
all_lag_param.tbl <- lapply(between_lag.mod.list,function(x){x[["param.tbl"]]}) %>% 
  bind_rows(.,.id = "DV") %>% 
  mutate(
    term = factor(term, levels = c("(Intercept)","screen_scale","day"),
                  labels = c("(Intercept)","Between\nScaled Screentime","Day")),
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

write_csv(all_lag_param.tbl, '/Users/jackiebeltran/Documents/jmb_digital/results/between_scaling_screen_lag_model.csv')

# reorder data table levels

all_lag_param.tbl$Outcome <- factor(all_lag_param.tbl$Outcome,
                                    levels = c("Anxious", "Distressed", "Depressed"))


ema_labels <- c("Anxious" = "Anxiety", "Distressed" = "Distress", "Depressed" = "Depression")


## Plot ##
all_lag_MA_mod_parm_plot <- 
  all_lag_param.tbl %>% 
  ggplot(aes(x = estimate, y = term,color = lag_time)) + 
  geom_point(size = 1.5) + 
  geom_line(aes(group = term)) +
  facet_wrap(~Outcome,ncol = 4, labeller = labeller(Outcome = ema_labels)) + 
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

# show only the screentime estimates
screen_est_lag_mod_plot <- 
  all_lag_param.tbl %>% 
  filter(term=="Between\nScaled Screentime") %>% 
  ggplot(aes(x = estimate, y = lag_time)) + 
  geom_text(aes(label = stars,x = (conf.high*1.15),size = 20)) + 
  geom_errorbar(aes(xmin = conf.low,xmax = conf.high),width = 0.2) + 
  geom_point(aes(color= Outcome),size = 2.5,show.legend = F) + 
  geom_vline(xintercept = 0,linetype='dashed') + 
  facet_wrap(~Outcome,ncol = 4, labeller = labeller(Outcome = ema_labels)) + 
  ggpubr::theme_pubclean() + 
  labs(
    x = "Screentime Estimate (between-participant scaling)",
    y = "Lag Period",
    color = "Lag Time"
  ) +
  scale_color_manual(values = c("#CC6633", "#669999", "#996699")) + 
  theme(
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 18),
    legend.position = "none"
  )


screen_est_lag_mod_plot
