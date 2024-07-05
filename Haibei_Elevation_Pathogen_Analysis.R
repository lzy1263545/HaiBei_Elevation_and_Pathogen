library(tidyverse)
library(scales)
library(readxl)
library(nlme)
library(lmerTest)
library(rsq)
library(sciplot)
library(ggrepel)
library(ggsci)
library(ggpubr)
library(MuMIn)
library(vegan)
library(ggpmisc)
library(MASS)
library(glmmTMB)
library(r2glmm)
library(jtools)
library(glmm.hp)
library(piecewiseSEM)
library(corrplot)
library(metafor)
library(ggnewscale)
library(car)

rm(list = ls())

##### ggplot2 theme #####

theme_set(theme_test())

mytheme = theme(plot.title = element_text(size = 10),
                plot.subtitle = element_text(size = 12),
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 12,
                                          face = "bold"),
                legend.title = element_text(size = 15,
                                            face = "bold"),
                legend.text = element_text(size = 12),
                legend.key.size = unit(0.5,"cm"))

##### Read data #####

dat = read_xlsx("HaibeiData.xlsx") %>% 
  mutate(PL = PL / 100,
         STE = STE / 100,
         SIV = SIV / 100,
         sfpRA = sfpRA / 100) %>%
  mutate(SIV = 0.5 * log((1 + SIV) / (1 - SIV)))
str(dat)
colnames(dat)

##### Principle component analysis for soil properties #####

soil = read_xlsx("soil.properties.xlsx")
soil.scaled = as.data.frame(scale(soil[,2:6]))
soil = cbind(soil[,1],
             soil.scaled)

soil.pca = rda(soil[,2:6])
soil.pca.result = summary(soil.pca)
pca.species = as.data.frame(soil.pca.result$species)
pca.sites = as.data.frame(soil.pca.result$sites) %>% 
  mutate(Elevation = soil$Elevation)

my.table = data.frame(Soil.properties = rownames(pca.species),
                      PC1 = round(pca.species$PC1,3),
                      PC2 = round(pca.species$PC2,3))

ggplot()+
  mytheme+
  geom_vline(xintercept = 0,
             linetype = 2,
             linewidth = 0.5,
             color = "gray25")+
  geom_hline(yintercept = 0,
             linetype = 2,
             linewidth = 0.5,
             color = "gray25")+
  geom_point(data = pca.sites,
             aes(x = PC1,
                 y = PC2,
                 fill = as.factor(Elevation)),
             shape = 21,
             size = 5,
             alpha = 0.75)+
  geom_segment(data = pca.species,
               aes(x = rep(0,5),
                   y = rep(0,5),
                   xend = PC1,
                   yend = PC2),
               linewidth = 1.25,
               arrow = arrow(angle = 25,
                             length = unit(0.25,"cm"),
                             type = "closed"),
               color = "#466983FF")+
  annotate(geom = "text",
           label = rownames(pca.species),
           x = c(1,1,-0.6,0.9,-0.75),
           y = c(-1,0.2,0.7,0.65,-1.1),
           size = 5,
           color = "#466983FF")+
  annotate(geom = "table",
           x = 3,
           y = 0,
           label = my.table,
           size = 4,
           color = "white",
           fill = "#466983FF",
           vjust = 1.4,
           hjust = 01)+
  labs(x = "Principle Component 1 (41.63%)",
       y = "Principle Component 2 (23.24%)")+
  scale_fill_frontiers()+
  theme(legend.position = "top",
        legend.title = element_blank())

ggsave("Soil.PCA.pdf",
       width = 8,
       height = 6)

##### Correlation matrix #####

corr = cor(dat %>% 
             dplyr::select(Elevation,
                           AGB,
                           BGB,
                           SR,
                           Evenness,
                           Soil_PCA1,
                           MDT,
                           MDH))

col = colorRampPalette(c("#809ECA","white","#953F56"))

corrplot.mixed(corr = corr,
               lower.col = col(50),
               upper.col = col(50),
               tl.pos = "lt",
               tl.col = "black")

##### Elevation effects on biotic & abiotic factors #####

col = c("AGB","BGB","SR","Evenness","Soil_PCA1","MDT","MDH")

plot.list = NULL
for (a in 1:length(col)) {
  
  foc.dat = dat %>% 
    select(Elevation,
           all_of(col[a]),
           lat,
           lon)
  colnames(foc.dat) = c("Elevation","Y","lon","lat")
  
  mod = lme(fixed = Y ~ Elevation,
            data = foc.dat, 
            correlation = corSpher(1,~ lat + lon),
            random = ~ 1 | Elevation,
            method = "REML")
  mod.r = summary(mod)
  mod.rsq = rsq.lmm(mod,
                    adj = TRUE)
  
  if(mod.r$tTable[2,5] < 0.05){
    
    reg.plot = ggplot(data = foc.dat,
                      aes(x = Elevation,
                          y = Y))+
      mytheme+
      geom_smooth(method = "lm",
                  linewidth = 1,
                  color = "#466983FF",
                  fill = "#466983FF")+
      geom_point(aes(fill = as.factor(Elevation)),
                 shape = 21,
                 alpha = 0.75,
                 size = 3)+
      labs(y = col[a],
           title = paste("p = ",
                         round(mod.r$tTable[2,5],3),
                         "   Marginal R^2 = ",
                         round(mod.rsq$fixed,3),
                         sep = ""))+
      scale_fill_frontiers()
  }
  
  if(mod.r$tTable[2,5] >= 0.05 & mod.r$tTable[2,5] < 0.1){
    
    reg.plot = ggplot(data = foc.dat,
                      aes(x = Elevation,
                          y = Y))+
      mytheme+
      geom_smooth(method = "lm",
                  linetype = 2,
                  se = FALSE,
                  linewidth = 1,
                  color = "#466983FF",
                  fill = "#466983FF")+
      geom_point(aes(fill = as.factor(Elevation)),
                 shape = 21,
                 alpha = 0.75,
                 size = 3)+
      labs(y = col[a],
           title = paste("p = ",
                         round(mod.r$tTable[2,5],3),
                         "   Marginal R^2 = ",
                         round(mod.rsq$fixed,3),
                         sep = ""))+
      scale_fill_frontiers()
  }
  
  if(mod.r$tTable[2,5] > 0.1){
    
    reg.plot = ggplot(data = foc.dat,
                      aes(x = Elevation,
                          y = Y))+
      mytheme+
      geom_point(aes(fill = as.factor(Elevation)),
                 shape = 21,
                 alpha = 0.75,
                 size = 3)+
      labs(y = col[a],
           title = paste("p = ",
                         round(mod.r$tTable[2,5],3),
                         "   Marginal R^2 = ",
                         round(mod.rsq$fixed,3),
                         sep = ""))+
      scale_fill_frontiers()
  }
  
  reg.plot = list(reg.plot)
  
  plot.list = c(plot.list,
                reg.plot)
}

ggarrange(plotlist = plot.list,
          ncol = 4,
          nrow = 2,
          align = "hv",
          common.legend = TRUE,
          legend = "right")

ggsave("Elevation.Effect.Factor.pdf",
       width = 12,
       height = 5)

##### Temperature effects on biotic & abiotic factors #####

col = c("AGB","BGB","SR","Evenness","Soil_PCA1")

plot.list = NULL
for (a in 1:length(col)) {
  
  foc.dat = dat %>% 
    dplyr::select(MDT,
                  Elevation,
                  all_of(col[a]),
                  lat,
                  lon)
  colnames(foc.dat) = c("Temperature","Elevation","Y","lon","lat")
  
  mod = lme(fixed = Y ~ Temperature,
            data = foc.dat, 
            correlation = corSpher(1,~ lat + lon),
            random = ~ 1 | Elevation,
            method = "REML")
  mod.r = summary(mod)
  mod.rsq = rsq.lmm(mod,
                    adj = TRUE)
  
  if(mod.r$tTable[2,5] < 0.05){
    
    reg.plot = ggplot(data = foc.dat,
                      aes(x = Temperature,
                          y = Y))+
      mytheme+
      geom_smooth(method = "lm",
                  linewidth = 1,
                  color = "#466983FF",
                  fill = "#466983FF")+
      geom_point(shape = 21,
                 alpha = 0.75,
                 size = 3)+
      labs(y = col[a],
           title = paste("p = ",
                         round(mod.r$tTable[2,5],3),
                         "   Marginal R^2 = ",
                         round(mod.rsq$fixed,3),
                         sep = ""))
  }
  
  if(mod.r$tTable[2,5] >= 0.05 & mod.r$tTable[2,5] < 0.1){
    
    reg.plot = ggplot(data = foc.dat,
                      aes(x = Temperature,
                          y = Y))+
      mytheme+
      geom_smooth(method = "lm",
                  linetype = 2,
                  se = FALSE,
                  linewidth = 1,
                  color = "#466983FF",
                  fill = "#466983FF")+
      geom_point(shape = 21,
                 alpha = 0.75,
                 size = 3)+
      labs(y = col[a],
           title = paste("p = ",
                         round(mod.r$tTable[2,5],3),
                         "   Marginal R^2 = ",
                         round(mod.rsq$fixed,3),
                         sep = ""))
  }
  
  if(mod.r$tTable[2,5] > 0.1){
    
    reg.plot = ggplot(data = foc.dat,
                      aes(x = Temperature,
                          y = Y))+
      mytheme+
      geom_point(shape = 21,
                 alpha = 0.75,
                 size = 3)+
      labs(y = col[a],
           title = paste("p = ",
                         round(mod.r$tTable[2,5],3),
                         "   Marginal R^2 = ",
                         round(mod.rsq$fixed,3),
                         sep = ""))
  }
  
  reg.plot = list(reg.plot)
  
  plot.list = c(plot.list,
                reg.plot)
}

ggarrange(plotlist = plot.list,
          ncol = 4,
          nrow = 2,
          align = "hv",
          common.legend = TRUE,
          legend = "right")

##### Elevation effects on pathogen #####

col = c("PL","STE","SIV","sfpOTUs","sfpRA")

plot.list = NULL
for (a in 1:length(col)) {
  
  foc.dat = dat %>% 
    dplyr::select(Elevation,
                  all_of(col[a]),
                  lat,
                  lon)
  colnames(foc.dat) = c("Elevation","Y","lon","lat")
  
  if(col[a] %in% c("PL","STE","sfpRA")){
    
    mod = glmmPQL(data = foc.dat,
                  fixed = Y ~ Elevation,
                  random = ~ 1 | Elevation,
                  correlation = corSpher(1,~ lat + lon),
                  family = beta_family(link = "logit"),
                  control = lmeControl(msMaxIter = 1000,
                                       msMaxEval = 1000,
                                       opt = "optim"))
    mod.r = summary(mod)
    mod.p = mod.r$tTable[2,5]
    mod.rsq = r2beta(mod)
    mod.rsq = mod.rsq$Rsq[2]
  }
  
  if(col[a] %in% c("SIV")){
    
    mod = lme(fixed = Y ~ Elevation,
              data = foc.dat, 
              correlation = corSpher(1,~ lat + lon),
              random = ~ 1 | Elevation,
              method = "REML",
              control = lmeControl(msMaxIter = 1000,
                                   msMaxEval = 1000,
                                   opt = "optim"))
    mod.r = summary(mod)
    mod.p = mod.r$tTable[2,5]
    mod.rsq = rsq.lmm(mod,
                      adj = TRUE)
    mod.rsq = mod.rsq$fixed
  }
  
  if(col[a] %in% c("sfpOTUs")){
    
    mod = glmmPQL(data = foc.dat,
                  fixed = Y ~ Elevation,
                  random = ~ 1 | Elevation,
                  correlation = corSpher(1,~ lat + lon),
                  family = poisson(link = "log"),
                  control = lmeControl(msMaxIter = 1000,
                                       msMaxEval = 1000,
                                       opt = "optim"))
    mod.r = summary(mod)
    mod.p = mod.r$tTable[2,5]
    mod.rsq = r2beta(mod)
    mod.rsq = mod.rsq$Rsq[2]
  }
  
  if(mod.p < 0.05){
    
    reg.plot = ggplot(data = foc.dat,
                      aes(x = Elevation,
                          y = Y))+
      mytheme+
      geom_smooth(method = "lm",
                  linewidth = 1,
                  color = "#466983FF",
                  fill = "#466983FF")+
      geom_point(aes(fill = as.factor(Elevation)),
                 shape = 21,
                 alpha = 0.75,
                 size = 3)+
      labs(y = col[a],
           title = paste("p = ",
                         round(mod.p,3),
                         "   Marginal R^2 = ",
                         round(mod.rsq,3),
                         sep = ""))+
      scale_fill_frontiers()
  }
  
  if(mod.p >= 0.05 & mod.p < 0.1){
    
    reg.plot = ggplot(data = foc.dat,
                      aes(x = Elevation,
                          y = Y))+
      mytheme+
      geom_smooth(method = "lm",
                  linetype = 2,
                  se = FALSE,
                  linewidth = 1,
                  color = "#466983FF",
                  fill = "#466983FF")+
      geom_point(aes(fill = as.factor(Elevation)),
                 shape = 21,
                 alpha = 0.75,
                 size = 3)+
      labs(y = col[a],
           title = paste("p = ",
                         round(mod.p,3),
                         "   Marginal R^2 = ",
                         round(mod.rsq,3),
                         sep = ""))+
      scale_fill_frontiers()
  }
  
  if(mod.p > 0.1){
    
    reg.plot = ggplot(data = foc.dat,
                      aes(x = Elevation,
                          y = Y))+
      mytheme+
      geom_point(aes(fill = as.factor(Elevation)),
                 shape = 21,
                 alpha = 0.75,
                 size = 3)+
      labs(y = col[a],
           title = paste("p = ",
                         round(mod.p,3),
                         "   Marginal R^2 = ",
                         round(mod.rsq,3),
                         sep = ""))+
      scale_fill_frontiers()
  }
  
  reg.plot = list(reg.plot)
  
  plot.list = c(plot.list,
                reg.plot)
}

ggarrange(plotlist = plot.list,
          ncol = 3,
          nrow = 2,
          align = "hv",
          common.legend = TRUE,
          legend = "right")

ggsave("Elevation.Effect.Pathogen.pdf",
       width = 10,
       height = 5)

##### Biotic & abiotic effects on pathogen #####

col = c("PL","STE","SIV","sfpOTUs","sfpRA")
X = c("AGB","BGB","SR","Evenness","Soil_PCA1","MDT","MDH")

plot.list = NULL
for (a in 1:length(col)){
  
  for (b in 1:length(X)) {
    
    foc.dat = dat %>% 
      dplyr::select(Elevation,
                    all_of(X[b]),
                    all_of(col[a]),
                    lat,
                    lon)
    colnames(foc.dat) = c("Elevation","X","Y","lon","lat")
    foc.dat$Elevation = as.factor(foc.dat$Elevation)
    
    if(col[a] %in% c("PL","STE","sfpRA")){
      
      mod = glmmPQL(data = foc.dat,
                    fixed = Y ~ X,
                    random = ~ 1 | Elevation,
                    correlation = corSpher(1,~ lat + lon),
                    family = beta_family(link = "logit"),
                    control = lmeControl(msMaxIter = 1000,
                                         msMaxEval = 1000,
                                         opt = "optim"))
      mod.r = summary(mod)
      mod.p = mod.r$tTable[2,5]
      mod.rsq = r2beta(mod)
      mod.rsq = mod.rsq$Rsq[2]
    }
    
    if(col[a] %in% c("SIV")){
      
      if(X[b] %nin% c("Soil_PCA1","MDT","MDH")){
        
        mod = lme(fixed = Y ~ X,
                  data = foc.dat, 
                  correlation = corSpher(1,~ lat + lon),
                  random = ~ 1 | Elevation,
                  method = "REML",
                  control = lmeControl(msMaxIter = 1000,
                                       msMaxEval = 1000,
                                       opt = "optim"))
        mod.r = summary(mod)
        mod.p = mod.r$tTable[2,5]
        mod.rsq = rsq.lmm(mod,
                          adj = TRUE)
        mod.rsq = mod.rsq$fixed
      }
      
      if(X[b] %in% c("Soil_PCA1","MDT","MDH")){
        
        mod = lme(fixed = Y ~ X,
                  data = foc.dat, 
                  correlation = corSpher(1,~ lat + lon),
                  random = ~ 1 | Elevation,
                  method = "REML",
                  control = lmeControl(msMaxIter = 1000,
                                       msMaxEval = 1000))
        mod.r = summary(mod)
        mod.p = mod.r$tTable[2,5]
        mod.rsq = rsq.lmm(mod,
                          adj = TRUE)
        mod.rsq = mod.rsq$fixed
      }
    }
    
    if(col[a] %in% c("sfpOTUs")){
      
      mod = glmmPQL(data = foc.dat,
                    fixed = Y ~ X,
                    random = ~ 1 | Elevation,
                    correlation = corSpher(1,~ lat + lon),
                    family = poisson(link = "log"),
                    control = lmeControl(msMaxIter = 1000,
                                         msMaxEval = 1000,
                                         opt = "optim"))
      mod.r = summary(mod)
      mod.p = mod.r$tTable[2,5]
      mod.rsq = r2beta(mod)
      mod.rsq = mod.rsq$Rsq[2]
    }
    
    if(mod.p < 0.05){
      
      reg.plot = ggplot(data = foc.dat,
                        aes(x = X,
                            y = Y))+
        mytheme+
        geom_smooth(method = "lm",
                    linewidth = 1,
                    color = "#466983FF",
                    fill = "#466983FF")+
        geom_point(aes(fill = as.factor(Elevation)),
                   shape = 21,
                   alpha = 0.75,
                   size = 3)+
        labs(y = col[a],
             x = X[b],
             title = paste("p = ",
                           round(mod.p,3),
                           "   Marginal R^2 = ",
                           round(mod.rsq,3),
                           sep = ""))+
        scale_fill_frontiers()
    }
    
    if(mod.p >= 0.05 & mod.p < 0.1){
      
      reg.plot = ggplot(data = foc.dat,
                        aes(x = X,
                            y = Y))+
        mytheme+
        geom_smooth(method = "lm",
                    linetype = 2,
                    se = FALSE,
                    linewidth = 1,
                    color = "#466983FF",
                    fill = "#466983FF")+
        geom_point(aes(fill = as.factor(Elevation)),
                   shape = 21,
                   alpha = 0.75,
                   size = 3)+
        labs(y = col[a],
             x = X[b],
             title = paste("p = ",
                           round(mod.p,3),
                           "   Marginal R^2 = ",
                           round(mod.rsq,3),
                           sep = ""))+
        scale_fill_frontiers()
    }
    
    if(mod.p > 0.1){
      
      reg.plot = ggplot(data = foc.dat,
                        aes(x = X,
                            y = Y))+
        mytheme+
        geom_point(aes(fill = as.factor(Elevation)),
                   shape = 21,
                   alpha = 0.75,
                   size = 3)+
        labs(y = col[a],
             x = X[b],
             title = paste("p = ",
                           round(mod.p,3),
                           "   Marginal R^2 = ",
                           round(mod.rsq,3),
                           sep = ""))+
        scale_fill_frontiers()
    }
    
   reg.plot = list(reg.plot)
  
  plot.list = c(plot.list,
                reg.plot) 
  }
}

all.plot = ggarrange(plotlist = plot.list,
                     ncol = 7,
                     nrow = 5,
                     align = "hv",
                     common.legend = TRUE,
                     legend = "top")

ggsave(plot = all.plot,
       filename = "Factor.Effect.Pathogen.pdf",
       width = 18,
       height = 10)

##### R-squared hierarchical partitioning & variables selection #####

###### [Pathogen load] ######

PL.model = glmmPQL(data = dat,
                   fixed = PL ~ Elevation + AGB + BGB + SR + Evenness + Soil_PCA1 + MDH, # removed MDT
                   random = ~ 1 | Elevation,
                   correlation = corSpher(1,~ lat + lon),
                   family = gaussian(),
                   control = lmeControl(msMaxIter = 1000,
                                        msMaxEval = 1000,
                                        opt = "optim"))
PL.model.result = summary(PL.model)
car::vif(PL.model)

PL.model.hp = glmm.hp(PL.model)
plot.glmmhp(PL.model.hp)

PL.model.hp$hierarchical.partitioning

###### [Species turnover effect] ######

STE.model = glmmPQL(data = dat,
                    fixed = STE ~ Elevation + AGB + BGB + SR + Evenness + Soil_PCA1 + MDH, # removed MDT
                    random = ~ 1 | Elevation,
                    correlation = corSpher(1,~ lat + lon),
                    family = gaussian(),
                    control = lmeControl(msMaxIter = 1000,
                                         msMaxEval = 1000))
STE.model.result = summary(STE.model)
car::vif(STE.model)

STE.model.hp = glmm.hp(STE.model)
plot.glmmhp(STE.model.hp)

STE.model.hp$hierarchical.partitioning

###### [Intraspecific variation] ######

SIV.model = lme(data = dat,
                fixed = SIV ~ Elevation + AGB + BGB + SR + Evenness + Soil_PCA1 + MDH, # removed MDT
                random = ~ 1 | Elevation,
                correlation = corSpher(1,~ lat + lon),
                control = lmeControl(msMaxIter = 1000,
                                     msMaxEval = 1000,
                                     opt = "optim"))
SIV.model.result = summary(SIV.model)
car::vif(SIV.model)

SIV.model.hp = glmm.hp(SIV.model)
plot.glmmhp(SIV.model.hp)

SIV.model.hp$hierarchical.partitioning

###### [Soil fungal pathogen OTUs] ######

sfpOTUs.model = glmmPQL(data = dat,
                        fixed = sfpOTUs ~ Elevation + AGB + BGB + SR + Evenness + Soil_PCA1 + MDH,  # removed MDT
                        random = ~ 1 | Elevation,
                        correlation = corSpher(1,~ lat + lon),
                        family = gaussian(),
                        control = lmeControl(msMaxIter = 1000,
                                             msMaxEval = 1000))
sfpOTUs.model.result = summary(sfpOTUs.model)
car::vif(sfpOTUs.model)

sfpOTUs.model.hp = glmm.hp(sfpOTUs.model)
plot.glmmhp(sfpOTUs.model.hp)

sfpOTUs.model.hp$hierarchical.partitioning

###### [Soil fungal pathogen relative abundance] ######

sfpRA.model = glmmPQL(data = dat,
                      fixed = sfpRA ~ Elevation + AGB + BGB + SR + Evenness + Soil_PCA1 + MDH,   # removed MDT
                      random = ~ 1 | Elevation,
                      correlation = cor.(1,~ lat + lon),
                      family = gaussian(),
                      control = lmeControl(msMaxIter = 1000,
                                           msMaxEval = 1000,
                                           opt = "optim"))
sfpRA.model.result = summary(sfpRA.model)
car::vif(sfpRA.model)

sfpRA.model.hp = glmm.hp(sfpRA.model)
plot.glmmhp(sfpRA.model.hp)

ad = as.data.frame(as.data.frame(PL.model.result$tTable) %>% 
                     mutate(Response = "Pathogen load",
                            Variable = rownames(PL.model.result$tTable),
                            CI = 1.96 * Std.Error) %>% 
                     filter(Variable != "(Intercept)"))

# ---- Plot

vsm.dat = rbind(as.data.frame(PL.model.result$tTable) %>% 
                  mutate(Response = "Pathogen load",
                         Variable = rownames(PL.model.result$tTable),
                         CI = 1.96 * Std.Error) %>% 
                  filter(Variable != "(Intercept)"),
                as.data.frame(STE.model.result$tTable) %>% 
                  mutate(Response = "Species turnover effect",
                         Variable = rownames(STE.model.result$tTable),
                         CI = 1.96 * Std.Error) %>% 
                  filter(Variable != "(Intercept)"),
                as.data.frame(SIV.model.result$tTable) %>% 
                  mutate(Response = "Intraspecific variation",
                         Variable = rownames(SIV.model.result$tTable),
                         CI = 1.96 * Std.Error) %>% 
                  filter(Variable != "(Intercept)"),
                as.data.frame(sfpOTUs.model.result$tTable) %>% 
                  mutate(Response = "Soil fungal pathogen OTUs",
                         Variable = rownames(sfpOTUs.model.result$tTable),
                         CI = 1.96 * Std.Error) %>% 
                  filter(Variable != "(Intercept)"),
                as.data.frame(sfpRA.model.result$tTable) %>% 
                  mutate(Response = "Soil fungal pathogen relative abundance",
                         Variable = rownames(sfpRA.model.result$tTable),
                         CI = 1.96 * Std.Error) %>% 
                  filter(Variable != "(Intercept)"))

colnames(vsm.dat) = c("mean","se","df","t","p","Response","Variable","CI")

vsm.dat$color = NA
vsm.dat$color[which(vsm.dat$mean < 0 & vsm.dat$p < 0.1)] = "neg.sig"
vsm.dat$color[which(vsm.dat$mean > 0 & vsm.dat$p < 0.1)] = "pos.sig"
vsm.dat$color[which(vsm.dat$p > 0.1)] = "no.sig"

hp.dat = rbind(as.data.frame(PL.model.hp$hierarchical.partitioning) %>% 
                 mutate(Response = "Pathogen load",
                        Variable = rownames(PL.model.hp$hierarchical.partitioning)),
               as.data.frame(STE.model.hp$hierarchical.partitioning) %>%
                 mutate(Response = "Species turnover effect",
                        Variable = rownames(STE.model.hp$hierarchical.partitioning)),
               as.data.frame(SIV.model.hp$hierarchical.partitioning) %>%
                 mutate(Response = "Intraspecific variation",
                        Variable = rownames(SIV.model.hp$hierarchical.partitioning)),
               as.data.frame(sfpOTUs.model.hp$hierarchical.partitioning) %>% 
                 mutate(Response = "Soil fungal pathogen OTUs",
                        Variable = rownames(sfpOTUs.model.hp$hierarchical.partitioning)),
               as.data.frame(sfpRA.model.hp$hierarchical.partitioning) %>% 
                 mutate(Response = "Soil fungal pathogen relative abundance",
                        Variable = rownames(sfpRA.model.hp$hierarchical.partitioning)))

colnames(hp.dat) = c("Unique","Average.share","Individual","Individual.effect","Response","Variable")

vsm.hp.dat = vsm.dat %>% 
  left_join(hp.dat,
            by = c("Response","Variable"))

(PL.es = ggplot(data = vsm.hp.dat %>% 
                   filter(Response == "Pathogen load"))+
    mytheme+
    geom_vline(xintercept = 0,
               linetype = 2,
               linewidth = 0.5,
               color = "lightblue4")+
    geom_point(aes(x = mean,
                   y = reorder(Variable,
                               Individual.effect),
                   color = color),
               shape = 16,
               size = 3,
               show.legend = FALSE)+
    geom_errorbar(aes(x = mean,
                      xmin = mean - CI,
                      xmax = mean + CI,
                      y = reorder(Variable,
                                  Individual.effect),
                      color = color),
                  linewidth = 1,
                  width = 0.25,
                  show.legend = FALSE)+
    labs(x = "Parameter estimates",
         title = "Pathogen load")+
    scale_color_manual(values = c("#B8B8B8FF","#7AA6DCFF","#CD534CFF"))+
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()))
PL.es.grob = ggplotGrob(PL.es)

(PL.imp = ggplot(data = vsm.hp.dat %>% 
                   filter(Response == "Pathogen load"),
                 aes(x = Individual.effect,
                     y = reorder(Variable,
                                 Individual.effect),
                     fill = Individual.effect))+
    mytheme+
    geom_bar(stat = "identity",
             linewidth = 1,
             alpha = 0.5,
             color = "black")+
    geom_text(aes(x = Individual.effect / 2,
                  y = reorder(Variable,
                              Individual.effect),
                  label = paste(Individual.effect,
                                "%",
                                sep = ""),
                  fontface = "bold"),
              size = 6)+
    labs(x = "Individual effect (%)")+
    scale_fill_gradientn(colours = c("steelblue","orange2","red3"))+
    theme(axis.title.y = element_blank(),
          axis.text = element_text(size = 16),
          axis.title.x = element_text(size = 18,
                                      face = "bold"))+
    annotation_custom(grob = PL.es.grob,
                      xmin = 14,xmax = 44,ymin = 0.5,ymax = 6.5))

(STE.es = ggplot(data = vsm.hp.dat %>% 
                  filter(Response == "Species turnover effect"))+
    mytheme+
    geom_vline(xintercept = 0,
               linetype = 2,
               linewidth = 0.5,
               color = "lightblue4")+
    geom_point(aes(x = mean,
                   y = reorder(Variable,
                               Individual.effect),
                   color = color),
               shape = 16,
               size = 3,
               show.legend = FALSE)+
    geom_errorbar(aes(x = mean,
                      xmin = mean - CI,
                      xmax = mean + CI,
                      y = reorder(Variable,
                                  Individual.effect),
                      color = color),
                  linewidth = 1,
                  width = 0.25,
                  show.legend = FALSE)+
    labs(x = "Parameter estimates",
         title = "Species turnover effect")+
    scale_color_manual(values = c("#B8B8B8FF","#CD534CFF","#7AA6DCFF"))+
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()))
STE.es.grob = ggplotGrob(STE.es)

(STE.imp = ggplot(data = vsm.hp.dat %>% 
                   filter(Response == "Species turnover effect"),
                 aes(x = Individual.effect,
                     y = reorder(Variable,
                                 Individual.effect),
                     fill = Individual.effect))+
    mytheme+
    annotation_custom(grob = STE.es.grob,
                      xmin = 13.8,xmax = 30,ymin = 0.5,ymax = 5.8)+
    geom_bar(stat = "identity",
             linewidth = 1,
             alpha = 0.5,
             color = "black")+
    geom_text(aes(x = Individual.effect / 2,
                  y = reorder(Variable,
                              Individual.effect),
                  label = paste(Individual.effect,
                                "%",
                                sep = ""),
                  fontface = "bold"),
              size = 6)+
    labs(x = "Individual effect (%)")+
    scale_fill_gradientn(colours = c("steelblue","orange2","red3"))+
    theme(axis.title.y = element_blank(),
          axis.text = element_text(size = 16),
          axis.title.x = element_text(size = 18,
                                      face = "bold")))

(SIV.es = ggplot(data = vsm.hp.dat %>% 
                   filter(Response == "Intraspecific variation"))+
    mytheme+
    geom_vline(xintercept = 0,
               linetype = 2,
               linewidth = 0.5,
               color = "lightblue4")+
    geom_point(aes(x = mean,
                   y = reorder(Variable,
                               Individual.effect),
                   color = color),
               shape = 16,
               size = 3,
               show.legend = FALSE)+
    geom_errorbar(aes(x = mean,
                      xmin = mean - CI,
                      xmax = mean + CI,
                      y = reorder(Variable,
                                  Individual.effect),
                      color = color),
                  linewidth = 1,
                  width = 0.25,
                  show.legend = FALSE)+
    labs(x = "Parameter estimates",
         title = "Intraspecific variation")+
    scale_color_manual(values = c("#B8B8B8FF","#CD534CFF","#7AA6DCFF"))+
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()))
SIV.es.grob = ggplotGrob(SIV.es)

(SIV.imp = ggplot(data = vsm.hp.dat %>% 
                    filter(Response == "Intraspecific variation"),
                  aes(x = Individual.effect,
                      y = reorder(Variable,
                                  Individual.effect),
                      fill = Individual.effect))+
    mytheme+
    annotation_custom(grob = SIV.es.grob,
                      xmin = 22,xmax = 64,ymin = 0.5,ymax = 6.75)+
    geom_bar(stat = "identity",
             linewidth = 1,
             alpha = 0.5,
             color = "black")+
    geom_text(aes(x = Individual.effect / 2,
                  y = reorder(Variable,
                              Individual.effect),
                  label = paste(Individual.effect,
                                "%",
                                sep = ""),
                  fontface = "bold"),
              size = 6)+
    labs(x = "Individual effect (%)")+
    scale_fill_gradientn(colours = c("steelblue","orange2","red3"))+
    theme(axis.title.y = element_blank(),
          axis.text = element_text(size = 16),
          axis.title.x = element_text(size = 18,
                                      face = "bold")))

(sfpOTUs.es = ggplot(data = vsm.hp.dat %>% 
                   filter(Response == "Soil fungal pathogen OTUs"))+
    mytheme+
    geom_vline(xintercept = 0,
               linetype = 2,
               linewidth = 0.5,
               color = "lightblue4")+
    geom_point(aes(x = mean,
                   y = reorder(Variable,
                               Individual.effect),
                   color = color),
               shape = 16,
               size = 3,
               show.legend = FALSE)+
    geom_errorbar(aes(x = mean,
                      xmin = mean - CI,
                      xmax = mean + CI,
                      y = reorder(Variable,
                                  Individual.effect),
                      color = color),
                  linewidth = 1,
                  width = 0.25,
                  show.legend = FALSE)+
    labs(x = "Parameter estimates",
         title = "Soil fungal pathogen OTUs")+
    scale_color_manual(values = c("#7AA6DCFF","#B8B8B8FF","#CD534CFF"))+
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()))
sfpOTUs.es.grob = ggplotGrob(sfpOTUs.es)

(sfpOTUs.imp = ggplot(data = vsm.hp.dat %>% 
                    filter(Response == "Soil fungal pathogen OTUs"),
                  aes(x = Individual.effect,
                      y = reorder(Variable,
                                  Individual.effect),
                      fill = Individual.effect))+
    mytheme+
    annotation_custom(grob = sfpOTUs.es.grob,
                      xmin = 28,xmax = 75,ymin = 0.5,ymax = 6.75)+
    geom_bar(stat = "identity",
             linewidth = 1,
             alpha = 0.5,
             color = "black")+
    geom_text(aes(x = Individual.effect / 2,
                  y = reorder(Variable,
                              Individual.effect),
                  label = paste(Individual.effect,
                                "%",
                                sep = ""),
                  fontface = "bold"),
              size = 6)+
    labs(x = "Individual effect (%)")+
    scale_fill_gradientn(colours = c("steelblue","orange2","red3"))+
    theme(axis.title.y = element_blank(),
          axis.text = element_text(size = 16),
          axis.title.x = element_text(size = 18,
                                      face = "bold")))

(sfpRA.es = ggplot(data = vsm.hp.dat %>% 
                       filter(Response == "Soil fungal pathogen relative abundance"))+
    mytheme+
    geom_vline(xintercept = 0,
               linetype = 2,
               linewidth = 0.5,
               color = "lightblue4")+
    geom_point(aes(x = mean,
                   y = reorder(Variable,
                               Individual.effect),
                   color = color),
               shape = 16,
               size = 5,
               show.legend = FALSE)+
    geom_errorbar(aes(x = mean,
                      xmin = mean - CI,
                      xmax = mean + CI,
                      y = reorder(Variable,
                                  Individual.effect),
                      color = color),
                  linewidth = 1,
                  width = 0.25,
                  show.legend = FALSE)+
    labs(x = "Parameter estimates",
         title = "Soil fungal pathogen relative abundance")+
    scale_color_manual(values = c("#B8B8B8FF","#7AA6DCFF","#CD534CFF"))+
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()))
sfpRA.es.grob = ggplotGrob(sfpRA.es)

(sfpRA.imp = ggplot(data = vsm.hp.dat %>% 
                        filter(Response == "Soil fungal pathogen relative abundance"),
                      aes(x = Individual.effect,
                          y = reorder(Variable,
                                      Individual.effect),
                          fill = Individual.effect))+
    mytheme+
    annotation_custom(grob = sfpRA.es.grob,
                      xmin = 14,xmax = 37,ymin = 0.5,ymax = 5.75)+
    geom_bar(stat = "identity",
             linewidth = 1,
             alpha = 0.5,
             color = "black")+
    geom_text(aes(x = Individual.effect / 2,
                  y = reorder(Variable,
                              Individual.effect),
                  label = paste(Individual.effect,
                                "%",
                                sep = ""),
                  fontface = "bold"),
              size = 6)+
    labs(x = "Individual effect (%)")+
    scale_fill_gradientn(colours = c("steelblue","orange2","red3"))+
    theme(axis.title.y = element_blank(),
          axis.text = element_text(size = 16),
          axis.title.x = element_text(size = 18,
                                      face = "bold")))

ggarrange(PL.imp,STE.imp,SIV.imp,sfpOTUs.imp,sfpRA.imp,
          nrow = 2,
          ncol = 3,
          common.legend = TRUE,
          legend = "right")

ggsave("Effect.Size.Importance.pdf",
       width = 18,
       height = 12)

##### Structural equation model #####

###### [Species turnover effect] ######

STE.sem.dat = dat %>% 
  dplyr::select(STE,Elevation,Evenness,Soil_PCA1,MDH) %>% 
  scale() %>% 
  as.data.frame() %>% 
  mutate(lat = dat$lat,
         lon = dat$lon)

STE.sem = psem(
  glmmPQL(data = STE.sem.dat,
          fixed = STE ~ Elevation + Evenness + Soil_PCA1 + MDH,
          random = ~ 1 | Elevation,
          correlation = corSpher(1,~ lat + lon),
          family = gaussian(),
          control = lmeControl(msMaxIter = 1000,
                               msMaxEval = 1000,
                               opt = "optim")),
  lme(data = STE.sem.dat,
          fixed = Evenness ~ Elevation,
          random = ~ 1 | Elevation,
          correlation = corSpher(1,~ lat + lon),
          control = lmeControl(msMaxIter = 1000,
                               msMaxEval = 1000,
                               opt = "optim")),
  lme(data = STE.sem.dat,
          fixed = Soil_PCA1 ~ Elevation,
          random = ~ 1 | Elevation,
          correlation = corSpher(1,~ lat + lon),
          control = lmeControl(msMaxIter = 1000,
                               msMaxEval = 1000,
                               opt = "optim")),
  lme(data = STE.sem.dat,
          fixed = MDH ~ Elevation,
          random = ~ 1 | Elevation,
          correlation = corSpher(1,~ lat + lon),
          control = lmeControl(msMaxIter = 1000,
                               msMaxEval = 1000,
                               opt = "optim"))
)

summary(STE.sem)

###### [Intraspecific variation] ######

SIV.sem.dat = dat %>% 
  dplyr::select(SIV,Elevation,Evenness,Soil_PCA1,MDT) %>% 
  scale() %>% 
  as.data.frame() %>% 
  mutate(lat = dat$lat,
         lon = dat$lon)

SIV.sem = psem(
  lme(data = SIV.sem.dat,
          fixed = SIV ~ Elevation +  Evenness + Soil_PCA1 + MDT,
          random = ~ 1 | Elevation,
      correlation = corSpher(1,~ lat + lon),
      control = lmeControl(msMaxIter = 1000,
                           msMaxEval = 1000)),
  lme(data = SIV.sem.dat,
      fixed = Evenness ~ Elevation,
      random = ~ 1 | Elevation,
      correlation = corSpher(1,~ lat + lon),
      control = lmeControl(msMaxIter = 1000,
                           msMaxEval = 1000)),
  lme(data = SIV.sem.dat,
      fixed = Soil_PCA1 ~ Elevation,
      random = ~ 1 | Elevation,
      correlation = corSpher(1,~ lat + lon),
      control = lmeControl(msMaxIter = 1000,
                           msMaxEval = 1000))
)
SIV.sem = update(SIV.sem,
                 Elevation %~~% MDT)

summary(SIV.sem)

###### [Soil fungal pathogen OTUs] ######

sfpOTUs.sem.dat = dat %>% 
  dplyr::select(sfpOTUs,Elevation,BGB,Soil_PCA1,MDT) %>% 
  scale() %>% 
  as.data.frame() %>% 
  mutate(lat = dat$lat,
         lon = dat$lon)

sfpOTUs.sem = psem(
  glmmPQL(data = sfpOTUs.sem.dat,
          fixed = sfpOTUs ~ Elevation + BGB + Soil_PCA1,
          random = ~ 1 | Elevation,
          correlation = corSpher(1,~ lat + lon),
          family = gaussian(),
          control = lmeControl(msMaxIter = 1000,
                               msMaxEval = 1000,
                               opt = "optim")),
  lme(data = sfpOTUs.sem.dat,
      fixed = BGB ~ Elevation,
      random = ~ 1 | Elevation,
      correlation = corSpher(1,~ lat + lon),
      control = lmeControl(msMaxIter = 1000,
                           msMaxEval = 1000,
                           opt = "optim")),
  lme(data = sfpOTUs.sem.dat,
      fixed = Soil_PCA1 ~ Elevation,
      random = ~ 1 | Elevation,
      correlation = corSpher(1,~ lat + lon),
      control = lmeControl(msMaxIter = 1000,
                           msMaxEval = 1000,
                           opt = "optim")))
sfpOTUs.sem = update(sfpOTUs.sem,
                     Elevation %~~% MDT)

summary(sfpOTUs.sem)

###### [Soil fungal pathogen relative abundance] ######

sfpRA.sem.dat = dat %>% 
  dplyr::select(sfpRA,Elevation,SR,Soil_PCA1,MDH) %>% 
  scale() %>% 
  as.data.frame() %>% 
  mutate(lat = dat$lat,
         lon = dat$lon)

sfpRA.sem = psem(
  glmmPQL(data = sfpRA.sem.dat,
          fixed = sfpRA ~ Elevation + SR + Soil_PCA1 + MDH,
          random = ~ 1 | Elevation,
          correlation = corSpher(1,~ lat + lon),
          family = gaussian(),
          control = lmeControl(msMaxIter = 1000,
                               msMaxEval = 1000,
                               opt = "optim")),
  lme(data = sfpRA.sem.dat,
      fixed = SR ~ Elevation,
      random = ~ 1 | Elevation,
      correlation = corSpher(1,~ lat + lon),
      control = lmeControl(msMaxIter = 1000,
                           msMaxEval = 1000,
                           opt = "optim")),
  lme(data = sfpRA.sem.dat,
      fixed = Soil_PCA1 ~ Elevation,
      random = ~ 1 | Elevation,
      correlation = corSpher(1,~ lat + lon),
      control = lmeControl(msMaxIter = 1000,
                           msMaxEval = 1000,
                           opt = "optim")),
  lme(data = sfpRA.sem.dat,
      fixed = MDH ~ Elevation,
      random = ~ 1 | Elevation,
      correlation = corSpher(1,~ lat + lon),
      control = lmeControl(msMaxIter = 1000,
                           msMaxEval = 1000,
                           opt = "optim"))
)

summary(sfpRA.sem)

##### Correlation among pathogen indices #####

mod = glmmPQL(data = dat,
              fixed = STE ~ SIV,
              random = ~ 1 | Elevation,
              correlation = corSpher(1,~ lat + lon),
              family = beta_family(),
              control = lmeControl(msMaxIter = 1000,
                                   msMaxEval = 1000,
                                   opt = "optim"))
summary(mod)
r2beta(mod)

mod = lme(data = dat,
          fixed = SIV ~ STE,
          random = ~ 1 | Elevation,
          correlation = corSpher(1,~ lat + lon),
          control = lmeControl(msMaxIter = 1000,
                               msMaxEval = 1000,
                               opt = "optim"))
summary(mod)
rsq(mod)


ggplot()+
  geom_point(data = dat,
             aes(x = STE,
                 y = SIV,
                 fill = as.factor(Elevation)),
             shape = 21,
             size = 3,
             alpha = 0.85)+
  theme_test()+
  scale_fill_frontiers()+
  mytheme+
  labs(x = "Species turnover effect",
       y = "Intraspecific effect",
       title = "p = 0.845   marginal R2 = 0.001")

ggsave("ReviewPlot.pdf",
       height = 4,
       width = 6)


(p1 = ggplot(data = dat %>% 
               gather("STE","SIV",
                      key = "X",value = "Value"),
             aes(x = Value,
                 y = PL,
                 group = X,
                 shape = X))+
    mytheme+
    geom_smooth(aes(color = X,
                    fill = X),
                method = "lm",
                se = TRUE,
                alpha = 0.25,
                linewidth = 1.5)+
    scale_color_manual(values = c("gray50","#6F99ADFF"))+
    scale_fill_manual(values = c("gray50","#6F99ADFF"))+
    new_scale_color()+
    new_scale_fill()+
    geom_point(aes(fill = as.factor(Elevation)),
               size = 3,
               alpha = 0.75)+
    labs(x = "Pathogen load components",
         y = "Pathogen load",
         title = paste("Species turnover effect: p = 0.005   Marginal R^2 = 0.283",
                       "\n",
                       "Intraspecific variation: p < 0.001   Marginal R^2 = 0.672",
                       sep = ""))+
    scale_fill_frontiers()+
    scale_shape_manual(values = c(21,22))+
    theme(plot.title = element_text(size = 12,
                                    face = "bold")))

(p2 = ggplot(data = dat %>% 
               gather(PL,STE,SIV,
                      key = "Response",value = "Value"))+
    mytheme+
    geom_density(aes(x = Value,
                     y = after_stat(density),
                     fill = Response,
                     color = Response),
                 alpha = 0.15,
                 linewidth = 1)+
    geom_segment(aes(x = Value,xend = Value,
                     y = Value - 0.5,yend = -2,
                     color = Response),
                 linewidth = 0.1)+
    scale_color_jama()+
    scale_fill_jama())

ggarrange(p1,p2,
          ncol = 1,
          align = "hv")

ggsave("Pathogen.Indices.Correlation.pdf",
       width = 6,
       height = 6)

##### Proportion #####

field.dat = read_xlsx("ProportionData.xlsx",
                      sheet = "field")
colnames(field.dat)

mod = lme(fixed = OTUs.Percentage ~ Elevation,
          data = field.dat, 
          correlation = corSpher(1,~ lat + lon),
          random = ~ 1 | Elevation,
          method = "REML")
summary(mod)
rsq.lmm(mod,
        adj = TRUE)

mod = glmmPQL(data = field.dat %>% 
                mutate(OTUs.Percentage = OTUs.Percentage / 100),
              fixed = Pathogenetic.OTUs ~ Elevation,
              random = ~ 1 | Elevation,
              correlation = corSpher(1,~ lat + lon),
              family = gaussian(),
              control = lmeControl(msMaxIter = 1000,
                                   msMaxEval = 1000,
                                   opt = "optim"))
summary(mod)
r2beta(mod)

(all.fungi = ggplot(data = field.dat)+
    geom_point(aes(x = Elevation,
                   y = Fungal.OTUs,
                   fill = as.factor(Elevation)),
               shape = 21,
               size = 5,
               alpha = 0.75)+
    labs(title = "p = 0.113   Marginal R2 = 0.154",
         y = "Soil fungal OTU richness")+
    scale_fill_frontiers()+
    theme(axis.title = element_text(size = 15,
                                    face = "bold"),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 12,
                                    face = "bold",
                                    hjust = 0.5)))

(p.fungi = ggplot(data = field.dat)+
    geom_smooth(aes(x = Elevation,
                    y = Pathogenetic.OTUs),
                method = "lm",
                linewidth = 1,
                color = "#6F99ADFF",
                fill = "#6F99ADFF",
                alpha = 0.25)+
    geom_point(aes(x = Elevation,
                   y = Pathogenetic.OTUs,
                   fill = as.factor(Elevation)),
               shape = 21,
               size = 5,
               alpha = 0.75)+
    scale_fill_frontiers()+
    labs(title = "p = 0.023   Marginal R2 = 0.429",
         y = "Soil fungal pathogenetic OTUs")+
    theme(axis.title = element_text(size = 15,
                                    face = "bold"),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 12,
                                    face = "bold",
                                    hjust = 0.5)))

(p.perc = ggplot(data = field.dat)+
    geom_smooth(aes(x = Elevation,
                    y = OTUs.Percentage),
                method = "lm",
                linewidth = 1,
                color = "#6F99ADFF",
                fill = "#6F99ADFF",
                alpha = 0.25)+
    geom_point(aes(x = Elevation,
                   y = OTUs.Percentage,
                   fill = as.factor(Elevation)),
               shape = 21,
               size = 5,
               alpha = 0.75)+
    labs(title = "p = 0.017   Marginal R2 = 0.462",
         y = "Percentage of\nfungal pathogenetic OTUs (%)")+
    scale_fill_frontiers()+
    theme(axis.title = element_text(size = 15,
                                    face = "bold"),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 12,
                                    face = "bold",
                                    hjust = 0.5)))

above.dat = read_xlsx("ProportionData.xlsx",
                      sheet = "field.aboveground")

vi.data = read_xlsx("ProportionData.xlsx",
                    sheet = "field.aboveground")
biomass.data = read_xlsx("ProportionData.xlsx",
                         sheet = "field.biomass")
pathogen.load = NULL
for (a in 1:30) {
  sum.n1 = NULL
  sum.n2 = NULL
  for (b in 1:73){
    foc.vi = as.numeric(vi.data[a,b + 1])
    foc.biomass = as.numeric(biomass.data[a,b + 1])
    if(is.na(foc.vi) == FALSE){
      calc.vi = foc.vi
      calc.biomass = foc.biomass
    }
    if(is.na(foc.vi) == TRUE){
      calc.vi = 0
      calc.biomass = 0
    }
    n1 = calc.vi * calc.biomass
    n2 = calc.biomass
    sum.n1 = sum(sum.n1,n1)
    sum.n2 = sum(sum.n2,n2)
  }
  PL = sum.n1 / sum.n2
  pathogen.load = c(pathogen.load,
                    PL)
}

turnover.effect = NULL
for (a in 1:30) {
  sum.n1 = NULL
  sum.n2 = NULL
  for (b in 1:73){
    calc.biomass = as.numeric(biomass.data[a,b + 1])
    if(is.na(calc.biomass) == FALSE){
      calc.biomass = calc.biomass
    }
    if(is.na(calc.biomass) == TRUE){
      calc.biomass = 0
    }
    pi = mean(c(vi.data[,b + 1])[[1]],na.rm = TRUE)
    n1 = pi * calc.biomass
    n2 = calc.biomass
    sum.n1 = sum(sum.n1,n1)
    sum.n2 = sum(sum.n2,n2)
  }
  turnover = sum.n1 / sum.n2
  turnover.effect = c(turnover.effect,
                      turnover)
}

intraspecific.variation = pathogen.load - turnover.effect

disease.dat = data.frame(PL = pathogen.load,
                         STE = turnover.effect,
                         SIV = intraspecific.variation)
write.csv(disease.dat,
          "disease.dat.csv")

proportion.dat = NULL
for (a in 1:30) {
  foc.ele = above.dat$Elevation[a]
  foc.vector = c(t(above.dat[a,]))[-1]
  n.all.sp = length(which(foc.vector != "NA"))
  n.uninfected.sp = length(which(foc.vector == 0))
  n.infected.sp = length(which(foc.vector > 0))
  proportion.df = data.frame(Elevation = foc.ele,
                             Uninfected = n.uninfected.sp / n.all.sp * 100,
                             Infected = n.infected.sp / n.all.sp * 100)
  proportion.dat = rbind(proportion.dat,
                         proportion.df)
}

proportion.dat = proportion.dat %>% 
  group_by(Elevation) %>% 
  summarise(Uninfected.mean = mean(Uninfected),
            Uninfected.se = se(Uninfected),
            Infected.mean = mean(Infected),
            Infected.se = se(Infected))

(proportion.plot = ggplot(data = proportion.dat)+
    geom_vline(xintercept = unique(proportion.dat$Elevation),
               linetype = 2,
               linewidth = 0.25,
               color = "lightblue4",
               alpha = 0.5)+
    geom_hline(yintercept = 50,
               linetype = 2,
               linewidth = 0.25,
               color = "lightblue4",
               alpha = 0.5)+
    geom_line(aes(x = Elevation,
                  y = Uninfected.mean),
              linewidth = 1,
              color = "#374E55FF")+
    geom_line(aes(x = Elevation,
                  y = Infected.mean),
              linewidth = 1,
              color = "#DF8F44FF")+
    geom_point(aes(x = Elevation,
                   y = Uninfected.mean),
               shape = 16,
               size = 5,
               color = "#374E55FF")+
    geom_errorbar(aes(x = Elevation,
                      ymin = Uninfected.mean - Uninfected.se,
                      ymax = Uninfected.mean + Uninfected.se),
                  width = 20,
                  linewidth = 1.5,
                  color = "#374E55FF")+
    geom_point(aes(x = Elevation,
                   y = Infected.mean),
               shape = 16,
               size = 5,
               color = "#DF8F44FF")+
    geom_errorbar(aes(x = Elevation,
                      ymin = Infected.mean - Infected.se,
                      ymax = Infected.mean + Infected.se),
                  width = 20,
                  linewidth = 1.5,
                  color = "#DF8F44FF")+
    geom_text(aes(x = Elevation,
                  y = 35,
                  label = round(Uninfected.mean,3)),
              color = "#374E55FF",
              size = 5)+
    geom_text(aes(x = Elevation,
                  y = 65,
                  label = round(Infected.mean,3)),
              color = "#DF8F44FF",
              size = 5)+
    annotate(geom = "text",
             label = c("Uninfected plant","Infected plant"),
             color = c("#374E55FF","#DF8F44FF"),
             x = c(3600,3600),
             y = c(45,55),
             size = 8)+
    labs(y = "Percentage of infected plant (%)")+
    theme(axis.title = element_text(size = 15,
                                    face = "bold"),
          axis.text = element_text(size = 12)))

meta.dat = read_xlsx("ProportionData.xlsx",
                     sheet = "meta")
colnames(meta.dat)

plot(meta.dat$Elevation,meta.dat$Fungal.OTUs)

calc.dat = meta.dat %>% 
  gather(Fungal.OTUs,Pathogenetic.OTUs,OTUs.Percentage,
         key = "Y",value = "value")

out.df = NULL
for (a in 1:length(unique(calc.dat$Citation))) {
  for (b in 1:length(unique(calc.dat$Y))) {
    foc.dat = calc.dat %>% 
      filter(Citation == unique(calc.dat$Citation)[a],
             Y == unique(calc.dat$Y)[b])
    N = nrow(foc.dat)
    corr = cor.test(foc.dat$Elevation,foc.dat$value,
                    method = "pearson")
    results = data.frame(Citation = unique(calc.dat$Citation)[a],
                         Y = unique(calc.dat$Y)[b],
                         r = corr$estimate,
                         N = N)
    out.df = rbind(out.df,
                   results)
  }
}
write.csv(out.df,
          "results.csv")

new.meta.dat = read_xlsx("ProportionData.xlsx",
                         sheet = "meta.es")
new.meta.dat$Y[which(new.meta.dat$Y == "Fungal.OTUs")] = "Soil fungal OTUs"
new.meta.dat$Y[which(new.meta.dat$Y == "Pathogenetic.OTUs")] = "Soil fungal pathogenetic OTUs"
new.meta.dat$Y[which(new.meta.dat$Y == "OTUs.Percentage")] = "Percentage of\n   fungal pathogenetic OTUs (%)"
unique(new.meta.dat$Y)

mod.meta = rma.uni(yi,vi,
                   subset = (new.meta.dat$Y == "Soil fungal OTUs"),
                   data = new.meta.dat,
                   method = "REML")
summary(mod.meta)
forest(mod.meta)

meta.plot.dat = data.frame(Y = c("Soil fungal OTUs",
                                 "Soil fungal pathogenetic OTUs",
                                 "Percentage of\n   fungal pathogenetic OTUs (%)"),
                           Mean = c(-0.2697,-0.2065,-0.006),
                           CI = c(0.2414,0.1986,0.2101))

(meta.plot = ggplot()+
    geom_hline(yintercept = 0,
               linetype = 2,
               linewidth = 0.5,
               color = "gray75")+
    geom_jitter(data = new.meta.dat,
                aes(x = Y,
                    y = yi,
                    size = 1 / vi),
                shape = 21,
                fill = "lightblue",
                alpha = 0.5,
                width = 0.1,
                show.legend = FALSE)+
    geom_pointrange(data = meta.plot.dat,
                    aes(x = Y,
                        y = Mean,
                        ymin = Mean - CI,
                        ymax = Mean + CI),
                    shape = 15,
                    size = 1,
                    linewidth = 2,
                    color = "#A73030FF")+
    annotate(geom = "text",
             label = c("*","*",'ns'),
             x = c(1,2,3),
             y = c(0.9,0.9,0.9),
             size = 6)+
    labs(y = "Effect size (Fisher's Z)")+
    scale_size(range = c(3,8))+
    scale_x_discrete(limits = c("Soil fungal OTUs",
                                "Soil fungal pathogenetic OTUs",
                                "Percentage of\n   fungal pathogenetic OTUs (%)"))+
    theme(axis.title = element_text(size = 15,
                                    face = "bold"),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12,
                                     angle = 10,
                                     hjust = 1),
          plot.title = element_text(size = 12,
                                    face = "bold",
                                    hjust = 0.5)))


ggpubr::ggarrange(ggarrange(all.fungi,p.fungi,p.perc,
                            nrow = 1,ncol = 3,
                            align = "hv",
                            legend = "none"),
                  ggarrange(proportion.plot,meta.plot,
                            nrow = 1,ncol = 2,
                            align = "hv"),
                  ncol = 1,nrow = 2,
                  align = "hv")

ggsave("Proportion.Plant.OTUs.pdf",
       height = 8,width = 12)

##### OTU rarecurve #####

oomy_zotutab = read_xlsx("HaibeiOTU.xlsx")

zotu = oomy_zotutab[,c(2:31)]
mean(colSums(zotu))
zotu_Flattening = as.data.frame(t(rrarefy(t(zotu),
                                          min(colSums(zotu)))))
colSums(zotu_Flattening)

rarecurve = rarecurve(t(zotu_Flattening),
                      step = 500,
                      label = TRUE,
                      col = "grey50")

rarecurve

Sample = attr(rarecurve[[1]],
              "Subsample")
plot_df_all = as.data.frame(NULL)
for (i in 1:30) {
  Sample_name = colnames(zotu_Flattening)[i]
  plot_df = data.frame(Sample, rarecurve[[i]],
                       rep(Sample_name,length(Sample)))
  plot_df_all = rbind(plot_df,plot_df_all)
}
colnames(plot_df_all) = c("Reads","Richness","Group")


plot_df_all$Elevation = NULL
plot_df_all$Elevation[which(plot_df_all$Group %in% paste("L0",1:6,sep = ""))] = "3200m"
plot_df_all$Elevation[which(plot_df_all$Group %in% c("L07","L08","L09","L10","L11","L12"))] = "3400m"
plot_df_all$Elevation[which(plot_df_all$Group %in% paste("L",13:18,sep = ""))] = "3600m"
plot_df_all$Elevation[which(plot_df_all$Group %in% paste("L",19:24,sep = ""))] = "3800m"
plot_df_all$Elevation[which(plot_df_all$Group %in% paste("L",25:30,sep = ""))] = "4000m"

max.dat = plot_df_all %>% 
  filter(Reads == 32001)

rarecurve_plot = ggplot(plot_df_all,aes(x = Reads,y = Richness,color = Elevation,group = Group))+
  theme_test()+
  mytheme+
  geom_line(linewidth = 0.5,alpha = 0.85)+
  geom_segment(data = max.dat,
               aes(x = rep(0,30),xend = rep(31979,30),
                   y = Richness,yend = Richness,
                   color = Elevation),
               linewidth = 0.25,
               linetype = 1,
               alpha = 0.85)+
  geom_vline(xintercept = 31979,linetype = "longdash",col = "steelblue4",linewidth = 1.0)+
  annotate('text',label = "Minimum sequencing depth = 31979 reads",x = 9000,y = 150,angle = 0,size = 4.5,color= "steelblue4",hjust = 0)+
  scale_y_continuous(expand = c(0,0),limit = c(0,1000),breaks = seq(0,1000,200))+
  scale_x_continuous(expand = c(0,0),limit = c(-1000,33001))+
  labs(x = "Number of sequences sampled",y = "Number of OTUs")+
  scale_color_frontiers()+
  theme(axis.title = element_text(size = 15),
        legend.key.size = unit(0.5,"cm"))
rarecurve_plot

ggsave("Rarecurve.plot.pdf",
       width = 8,
       height = 6)






