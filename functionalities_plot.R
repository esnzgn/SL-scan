rm(list = ls())

library(reshape2)
library(ggplot2)
library(tidyr)
library(tidyverse)

# functionality checks
# load func data 
mn_SLscan_model = read.csv("../MetabolicSLOutput/core_g_446_imat_cons_fba_dim_ublb_ttestV.csv")
mn_SLscan_model = mn_SLscan_model %>% t() %>% data.frame()
mn_SLscan_model_dim = select(mn_SLscan_model, c("X13", "X14"))
mn_SLscan_model_dim = data.frame(Metabolites = (as.numeric(mn_SLscan_model_dim$X13)), Reactions = (as.numeric(mn_SLscan_model_dim$X14)))
mn_SLscan_model_dim = subset(mn_SLscan_model_dim, mn_SLscan_model_dim$Reactions < 5000)
png(paste0("../MetabolicSLOutput/stat_density2d_slscanMNM_dim.png"),width=400, height=400)
  ggplot(mn_SLscan_model_dim, aes(x = Metabolites, y = Reactions)) +
                 geom_point() +
                 stat_density2d() +
    labs(
      x = "Number of metabolites",
      y = "Number of reactions"
    ) +
    theme(
      axis.text = element_text(size = 2*rel(1)),
      axis.title = element_text(size = 2*rel(1))
    )
dev.off()

func_met = read.csv('../MetabolicSLOutput/Table S1 - func_test_ modelR204_cons_iMAT .csv')
# func_met$NoFunctionalities
df = read.csv('../MetabolicSLOutput/Table S2 - sequence of random funcs_modelR204_cons_iMAT.csv')

can_name = read.csv('../MetabolicSLOutput/R24_imat_cons_cantyp.csv')
rownames(df) = gsub('cons_R24_iMAT_', '', can_name[grep('cons_R24_iMAT_', can_name)])
df = rownames_to_column(df, var = "Cancer")

long = melt(df, id.vars = "Cancer")
long$value = log(long$value)
add_long = data.frame(cbind(long[1:14,-3] ,value = log(func_met$NoFunctionalities)))
pdf(paste0("../MetabolicSLOutput/Violin_plot_functionalities_check_14_MNM.pdf"),width=25, height=15)
# png(paste0("../MetabolicSLOutput/Violin_plot_functionalities_check_14_MNM.png"),width=1000, height=700)
# Plot the violin plots
plot = ggplot(long, aes(x = Cancer, y = value, color = value)) +
  geom_violin() +
  geom_jitter(width = 0.1, alpha = 0.5) +
  scale_y_continuous(trans = 'log', breaks = c(0, 1, 10, 100)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 2*rel(1)),
        axis.title.x = element_blank()
        ,plot.margin = unit(c(2, 2, 2, 2), "cm")) +
  labs(x = "Cancers", y = "Number of functionalities")

plot + geom_point(data = add_long, aes(x = Cancer, y = exp(value))
                 , color = "black", size = 2, alpha = 0.5)
dev.off()
  
# pdf(paste0("../MetabolicSLOutput/Violin_plot_functionalities_check_14_MNM.pdf"),width=25, height=15)


