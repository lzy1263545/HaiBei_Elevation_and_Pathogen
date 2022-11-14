########## Analysis part ##########

library(vegan)
library(tidyverse)
library(readxl)
library(AICcmodavg)
library(piecewiseSEM)
library(nlme)

####:Note1 - Population-level disease severity (p), community-level pathogen load (PL), above- and belowground biomass (AB & BB) were directly calculated in Excel. Mean daily temperature and humidity (MDT & MDH) were extracted from the Temperature-Humidity Recorder Cos-03-0 (Renke Control Technology Co., Ltd., Jinan, Shandong, China). Soil fungal pathogen OTU richness and relative abundance (sfpOTUs & sfpRA) were directly calculated based on the OTU table in Excel; Note2 - Only community-level data of Haibei field study were open in this study

#####: Process the plant community survey data
Plant_community_data <- read.table("Plant_community.txt") #: The plant community original data is not open yet
Plant_community_data <- t(Plant_community_data)
###############: Calculate species richness (SR)
SR <- specnumber(Plant_community_data)
###############: Calculate Pielou's evenness index (Ev)
H <- diversity(Plant_community_data)
Ev <- H/log(SR)

#####: Process the soil properties data
Soil_properties_data <- read.table("Soil_properties.txt",header = T) #: A table contains 5 soil properties [i.e. pH, conductivity (C), moisture content (W), nitrate-nitrogen (NH4+) and ammonium-nitrogen (NO3-)] of each plot (not yet open)
Soil_properties_data <- Soil_properties_data %>% mutate(pH = scale(pH),
                                                        C = scale(C),
                                                        W = scale(W),
                                                        NH = scale(NH),
                                                        NO = scale(NO))
pca <- rda(Soil_properties_data)
pca.result <- summary(pca)
pca.result$species
###############: Extract the first principal component of soil properties (Soil PCA1)
Soil_PCA1 <- pca.result$sites[1:30]

#####: Read the data of Haibei field study
Field_study_data <- read_xlsx("HaibeiData.xlsx",sheet = "HaibeiData") #: We introduced the latitudes and longitudes for each sample plot to calculate the spatial matrix and introduce it into our linear mix effect models to resolve the spatial autocorrelation

Field_study_data <- Field_study_data %>% select(-Plot.ID) %>% mutate(PL = log(PL),
                                                                     p = log(p),
                                                                     sfpRA = log(sfpRA)) %>% #: Calculated log-transformed community pathogen load (PL), community proneness index (p) and soil fungal pathogen relative abundance (sfpRA) to achieve normal distribution
                                                              cbind.data.frame(Elevation = c("3200m","3200m","3200m","3200m","3200m","3200m","3400m","3400m","3400m","3400m","3400m","3400m","3600m","3600m","3600m","3600m","3600m","3600m","3800m","3800m","3800m","3800m","3800m","3800m","4000m","4000m","4000m","4000m","4000m","4000m")) #: We added "Elevation" (i.e. five elevational gradients) as the random effect term


#####: Linear regressions between multiple variables and AIC model selection (Note - the delta-AICc and wAICc value were directly calculated in Excel)
###############: Linear regressions and AIC model selection for community pathogen load (PL)
model.PL <- lme(fixed =  PL ~ X, data = Field_study_data, 
                correlation = corSpher(1,~lat+lon), #: We choose to use the corSpher (Spherical Correlation Structure) class for calculation of spatial matrix
                random = ~ 1 | Elevation, method = "ML")
#: X in linear mix effect model includes - Ele, AB, BB, SR, Ev, p, Soil_PCA1, MDT, MDH, ~1 (Null model)
summary(model.PL)
logLik(model.PL)
AICc(model.PL)
rsquared(model.PL)

###############: Linear regressions and AIC model selection for soil fungal pathogen OTU richness (sfpOTUs)
model.sfpOTUs <- lme(fixed =  sfpOTUs ~ X, data = Field_study_data, 
                correlation = corSpher(1,~lat+lon), #: We choose to use the corSpher (Spherical Correlation Structure) class for calculation of spatial matrix
                random = ~ 1 | Elevation, method = "ML")
#: X in linear mix effect model includes - Ele, AB, BB, SR, Ev, p, Soil_PCA1, MDT, MDH, ~1 (Null model)
summary(model.sfpOTUs)
logLik(model.sfpOTUs)
AICc(model.sfpOTUs)
rsquared(model.sfpOTUs)

###############: Linear regressions and AIC model selection for soil fungal pathogen OTU relative abundance (sfpRA)
model.sfpRA <- lme(fixed =  sfpRA ~ X, data = Field_study_data, 
                correlation = corSpher(1,~lat+lon), #: We choose to use the corSpher (Spherical Correlation Structure) class for calculation of spatial matrix
                random = ~ 1 | Elevation, method = "ML")
#: X in linear mix effect model includes - Ele, AB, BB, SR, Ev, p, Soil_PCA1, MDT, MDH, ~1 (Null model)
summary(model.sfpRA)
logLik(model.sfpRA)
AICc(model.sfpRA)
rsquared(model.sfpRA)

#####: PERMANOVA to test the compositional difference of soil fungal pathogens along the elevational gradient
OTU_site_matrix1 <- read.table("OTU_site_matrix1.txt",header = T) #: "OTU_site_matrix1" is the OTU table without elevation (not yet open)
OTU_site_matrix2 <- read.table("OTU_site_matrix2.txt",header = T)#: "OTU_site_matrix2" is the OTU table with elevation (not yet open)
disimilarity <- vegdist(OTU_site_matrix1,method = "bray")
adonis_result = adonis2(disimilarity~Elevation,OTU_site_matrix2, permutations = 999)
adonis_result

#####: Structural equation model (SEM) based on "piecewiseSEM" package
sem.model <- psem(
  lme(fixed =  sfpOTUs ~ Soil_PCA1+p+Ev+Ele, data = Field_study_data,correlation = corSpher(1,~lat+lon),
      random = ~ 1 | Elevation, method = "ML"),
  lme(fixed =  PL ~ Soil_PCA1+p+Ev+Ele, data = Field_study_data,correlation = corSpher(1,~lat+lon),
      random = ~ 1 | Elevation, method = "ML"),
  lme(fixed =  Soil_PCA1 ~ Ele, data = Field_study_data,correlation = corSpher(1,~lat+lon),
      random = ~ 1 | Elevation, method = "ML"),
  lme(fixed =  Ev ~ Ele, data = Field_study_data,correlation = corSpher(1,~lat+lon),
      random = ~ 1 | Elevation, method = "ML"),
  lme(fixed =  p ~ Ele, data = Field_study_data,correlation = corSpher(1,~lat+lon),
      random = ~ 1 | Elevation, method = "ML")
)

#: Report the results of final SEM
(sem.result <- summary(sem.model))
rsquared(sem.model)
write.csv(sem.result$coefficients,"sem.result.csv")


########## Figure part ##########

library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggthemes)
library(morris)
library(ochRe)
library(ggnewscale)
library(scales)
library(tidyverse)
library(ggcor)

#####: Plot biplot for principal component analysis (PCA) of soil properties
biplot.data <- as.data.frame(pca.result$sites[1:30,1:2]) %>% mutate(Ele = c("3200m","3200m","3200m","3200m","3200m","3200m","3400m","3400m","3400m","3400m","3400m","3400m","3600m","3600m","3600m","3600m","3600m","3600m","3800m","3800m","3800m","3800m","3800m","3800m","4000m","4000m","4000m","4000m","4000m","4000m")) #: Extract the first and second principle component of PCA and add the elevation for grouping
theme_set(theme_test())
biplot <- ggplot(data = biplot.data)+
  geom_point(aes(x = PC1,y = PC2,color = Ele,shape = Ele,fill = Ele),size = 4,alpha = 0.5,stroke = 1.5)+
  #stat_ellipse(aes(x = PC1,y = PC2,color = Ele,fill = Ele),geom ="polygon",level = 0.95,size = 1,alpha = 0.1)+
  scale_color_manual(values = c("#FDDBC7", "#F4A582",
                                "#D6604D", "#B2182B", "#67001F"))+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  scale_shape_manual(values = c(0,1,2,8,7))+
  geom_hline(yintercept = 0,size = 0.5, alpha = 0.4,linetype = 3)+
  geom_vline(xintercept = 0,size = 0.5, alpha = 0.4,linetype = 3)+
  geom_segment(aes(x = 0,y = 0,xend = 0.910,yend = -0.935),size = 1,color = "steelblue",arrow = arrow(length = unit(0.4,"cm")))+
  geom_segment(aes(x = 0,y = 0,xend = 0.927,yend = 0.335),size = 1,color = "steelblue",arrow = arrow(length = unit(0.4,"cm")))+
  geom_segment(aes(x = 0,y = 0,xend = -1.140,yend = 0.748),size = 1,color = "steelblue",arrow = arrow(length = unit(0.4,"cm")))+
  geom_segment(aes(x = 0,y = 0,xend = 1.196,yend = 0.529),size = 1,color = "steelblue",arrow = arrow(length = unit(0.4,"cm")))+
  geom_segment(aes(x = 0,y = 0,xend = -0.771,yend = -0.986),size = 1,color = "steelblue",arrow = arrow(length = unit(0.4,"cm")))+
  annotate(geom = "text", x = 1, y = -1,label = c("pH"),size = 10,family = "serif",color = "steelblue")+
  annotate(geom = "text", x = -0.75, y = -1.15,label = c("NH4+"),size = 10,family = "serif",color = "steelblue")+
  annotate(geom = "text", x = -1, y = 0.85,label = c("W"),size = 10,family = "serif",color = "steelblue")+
  annotate(geom = "text", x = 1.3, y = 0.7,label = c("NO3-"),size = 10,family = "serif",color = "steelblue")+
  annotate(geom = "text", x = 1, y = 0.25,label = c("C"),size = 10,family = "serif",color = "steelblue")
biplot
#ggsave("biplot.pdf",width = 8,height = 6)

#####: Plot correlation matrix among all variables and conduct Mantel test between soil fungal pathogen OTU richness (sfpOTUs) and multiple variables
Pathogen_variables <- Field_study_data %>% select(11:13) #: Three pathogen related variables (i.e. PL, sfpOTUs, sfpRA)
Other_variables <- Field_study_data %>% select(2:10) #: Other elevation, soil, plant community, environment variables
cor.coef <- read.csv("coefs.csv") #: Previously extract and process Pearson's correlation coefficients "r" (and corresponding P value) between "Pathogen_variables" and "Other_variables" (not provide here)
###############: Pearson correlation matrix
mantel1 <- mantel_test(spec = Other_variables,env = Pathogen_variables,mantel.fun = 'mantel.randtest',spec.dist.method = 'euclidean', env.dist.method = 'euclidean',
                       spec.select = list(PL = 1,
                                          sfpOTUs = 2,
                                          sfpRA = 3))%>%
  mutate(r_value = cut(cor.coef$rr, breaks = c(-Inf,0.1, 0.2, 0.3, Inf),
                       labels = c('<0.1','0.1-0.2','0.2-0.3','>=0.3'), right = FALSE),
         p_value = cut(corcoef$pp, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c('<0.001', '0.001-0.01', '0.01-0.05', '>=0.05'), right = FALSE),
         name = corcoef$d)
corplot1 <- quickcor(Pathogen_variables, type = "upper",grid.size = 0.5,grid.colour = "grey",show.diag = F) +
  geom_circle2() +
  anno_link(data = mantel1,aes(colour = p_value,size = r_value)) +
  scale_size_manual(values = c(0.3,1,1.5,2.5,4))+
  guides(size = guide_legend(title = "Pearson's r",title.theme = element_text(size = 12,face = "bold"),label.theme = element_text(size = 12,face = "bold"),order = 2),
         colour = guide_legend(title = "Pearson's p",title.theme = element_text(size = 12,face = "bold"),label.theme = element_text(size = 12,face = "bold"),order = 3),
         fill = guide_colorbar(title = "Pearson's r",title.theme = element_text(size = 12,face = "bold"),label.theme = element_text(size = 12,face = "bold"),order = 1))+
  theme_cor()+
  scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#CCCCCC"))+
  scale_fill_gradient2(midpoint = 0, low = "#053061", mid = "white",high = "#67001F", space = "Lab")+
  theme(
    legend.key.width = unit(0.8, "cm"),
    legend.key.height = unit(0.3, "cm")
  )
corplot1
###############: Mantel test (only for sfpOTUs)
mantel2 <- mantel_test(spec = Other_variables,env = Pathogen_variables,mantel.fun = 'mantel.randtest',spec.dist.method = 'bray', env.dist.method = 'euclidean',spec.select = list(sfpOTUs = 1))%>% mutate(r_value = cut(r, breaks = c(-Inf,0.1, 0.2, 0.3, Inf),labels = c('<0.1', '0.1-0.2','0.2-0.3','>=0.3'), right = FALSE))%>% mutate(p_value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),labels = c('<0.001', '0.001-0.01', '0.01-0.05', '>=0.05'), right = FALSE))
corplot2 <- quickcor(Pathogen_variables, type = "upper",grid.size = 0.5,grid.colour = "grey",show.diag = F) +
  geom_circle2() +
  anno_link(data = mantel2,aes(colour = p_value,size = r_value)) +
  scale_size_manual(values = c(0.3,1,1.5,2.5,4))+
  guides(size = guide_legend(title = "Mantel's r",title.theme = element_text(size = 12,face = "bold"),label.theme = element_text(size = 12,face = "bold"),order = 2),
         colour = guide_legend(title = "Mantel's p",title.theme = element_text(size = 12,face = "bold"),label.theme = element_text(size = 12,face = "bold"),order = 3),
         fill = guide_colorbar(title = "Pearson's r",title.theme = element_text(size = 12,face = "bold"),label.theme = element_text(size = 12,face = "bold"),order = 1))+
  theme_cor()+
  scale_color_manual(values=c("#E69F00", "#999999"))+
  scale_fill_gradient2(midpoint = 0, low = "#053061", mid = "white",high = "#67001F", space = "Lab")+
  theme(
    legend.key.width = unit(0.8, "cm"),
    legend.key.height = unit(0.4, "cm")
  )
corplot2
###############: Joint plot
ggpubr::ggarrange(corplot1,corplot2,ncol = 1,nrow = 2)
#ggsave("Correlation_matrix.pdf",width = 10,height = 10)

#: set the theme for plotting
theme_set(theme_light())
mytheme = theme(panel.border = element_rect(size = 1.5,color = "black"),
                plot.subtitle = element_text(size = 14,family = "serif"),
                axis.text.x = element_text(size = 12,face = "bold",family = "serif"),
                axis.text.y = element_text(size = 12,face = "bold",family = "serif"),
                axis.title.x = element_text(size = 16,face = "bold",family = "serif"),
                axis.title.y = element_text(size = 16,face = "bold",family = "serif"),
                legend.position = "none")
#####: Linear mix effect model regression plots between elevation and various (soil, plant community, environment) variables
###############: Aboveground biomass ~ Elevation
Elevation_plot1 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = Ele,y = AB,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  geom_smooth(aes(x = Ele,y = AB),method = "lm",se = T,size = 1.5,linetype = 1,alpha = 0.3,color = "#3C5488FF",fill = "#3C5488FF")+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                                "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.006   Marginal R2 = 0.733")+
  ylab(label = "Aboveground biomass")+
  xlab(label = "Elevation")
Elevation_plot1
###############: Belowground biomass ~ Elevation
Elevation_plot2 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = Ele,y = BB,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.112   Marginal R2 = 0.223")+
  ylab(label = "Belowground biomass")+
  xlab(label = "Elevation")
Elevation_plot2
###############: Species richness ~ Elevation
Elevation_plot3 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = Ele,y = SR,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  geom_smooth(aes(x = Ele,y = SR),method = "lm",se = T,size = 1.5,linetype = 1,alpha = 0.3,color = "#3C5488FF",fill = "#3C5488FF")+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.028   Marginal R2 = 0.482")+
  ylab(label = "Species  richness")+
  xlab(label = "Elevation")
Elevation_plot3
###############: Pielou's evenness index ~ Elevation
Elevation_plot4 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = Ele,y = Ev,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.798   Marginal R2 = 0.004")+
  ylab(label = "Pielou's evenness index")+
  xlab(label = "Elevation")
Elevation_plot4
###############: Community proneness ~ Elevation
Elevation_plot5 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = Ele,y = p,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.176   Marginal R2 = 0.103")+
  ylab(label = "Community proneness")+
  xlab(label = "Elevation")
Elevation_plot5
###############: Soil PCA1 ~ Elevation
Elevation_plot6 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = Ele,y = Soil_PCA1,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  geom_smooth(aes(x = Ele,y = Soil_PCA1),method = "lm",se = T,size = 1.5,linetype = 1,alpha = 0.3,color = "#3C5488FF",fill = "#3C5488FF")+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.038   Marginal R2 = 0.610")+
  ylab(label = "Soil PCA1")+
  xlab(label = "Elevation")
Elevation_plot6
###############: Mean daily temperature ~ Elevation
Elevation_plot7 <- ggplot(data = Field_study_data)+
  geom_jitter(aes(x = Ele,y = MDT,fill = Elevation),shape = 21,size = 7,alpha = 0.8,height = 0.1,width = 10)+
  geom_smooth(aes(x = Ele,y = MDT),method = "lm",se = T,size = 1.5,linetype = 1,alpha = 0.3,color = "#3C5488FF",fill = "#3C5488FF")+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.001   Marginal R2 = 0.971")+
  ylab(label = "Mean daily temperature")+
  xlab(label = "Elevation")
Elevation_plot7
###############: Mean daily humidity ~ Elevation
Elevation_plot8 <- ggplot(data = Field_study_data)+
  geom_jitter(aes(x = Ele,y = MDH,fill = Elevation),shape = 21,size = 7,alpha = 0.8,height = 0.05,width = 10)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.612   Marginal R2 = 0.066")+
  ylab(label = "Mean daily humidity")+
  xlab(label = "Elevation")
Elevation_plot8
###############: Joint plot
ggarrange(Elevation_plot1,Elevation_plot2,Elevation_plot3,Elevation_plot4,Elevation_plot5,Elevation_plot6,Elevation_plot7,Elevation_plot8,ncol = 4,nrow = 2, labels="AUTO",common.legend = T,legend = "bottom")
ggsave("Elevation_effect_on_variables.pdf",width = 16,height = 8)

#####: Linear mix effect model regression plots between elevation and pathogen related variables (PL, sfpOTUs, sfpRA)
###############: Community pathogen load ~ Elevation
Elevation_plot_p1 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = Ele,y = PL,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.403   Marginal R2 = 0.034")+
  ylab(label = "Pathogen load")+
  xlab(label = "Elevation")
Elevation_plot_p1
###############: Soil fungal pathogen OTU richness ~ Elevation
Elevation_plot_p2 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = Ele,y = sfpOTUs,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  geom_smooth(aes(x = Ele,y = sfpOTUs),method = "lm",se = T,size = 1.5,linetype = 1,alpha = 0.3,color = "#3C5488FF",fill = "#3C5488FF")+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.020   Marginal R2 = 0.433")+
  ylab(label = "Soil fungal pathogen OTU richness")+
  xlab(label = "Elevation")
Elevation_plot_p2
###############: Soil fungal pathogen relative abundance ~ Elevation
Elevation_plot_p3 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = Ele,y = sfpRA,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.973   Marginal R2 < 0.001")+
  ylab(label = "Soil fungal pathogen relative abundance")+
  xlab(label = "Elevation")
Elevation_plot_p3
###############: Joint plot
ggarrange(Elevation_plot_p1,Elevation_plot_p2,Elevation_plot_p3,ncol = 3,nrow = 1, labels="AUTO",common.legend = T,legend = "bottom")
ggsave("Elevation_effect_on_pathogen.pdf",width = 12,height = 4)

#####: Linear mix effect model regression plots between various (soil, plant community, environment) variables and community pathogen load (PL)
###############: Pathogen load ~ Aboveground biomass
PL_plot1 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = AB,y = PL,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.376   Marginal R2 = 0.029")+
  ylab(label = "Pathogen load")+
  xlab(label = "Aboveground biomass")
PL_plot1
###############: Pathogen load ~ Belowground biomass
PL_plot2 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = BB,y = PL,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.242   Marginal R2 = 0.050")+
  ylab(label = "Pathogen load")+
  xlab(label = "Belowground biomass")
PL_plot2
###############: Pathogen load ~ Species richness
PL_plot3 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = SR,y = PL,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.269   Marginal R2 = 0.045")+
  ylab(label = "Pathogen load")+
  xlab(label = "Species richness")
PL_plot3
###############: Pathogen load ~ Pielou's evenness index
PL_plot4 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = Ev,y = PL,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.514   Marginal R2 = 0.016")+
  ylab(label = "Pathogen load")+
  xlab(label = "Pielou's evenness index")
PL_plot4
###############: Pathogen load ~ Community proneness
PL_plot5 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = p,y = PL,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  geom_smooth(aes(x = p,y = PL),method = "lm",se = T,size = 1.5,linetype = 1,alpha = 0.3,color = "#3C5488FF",fill = "#3C5488FF")+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P < 0.001   Marginal R2 = 0.348")+
  ylab(label = "Pathogen load")+
  xlab(label = "Community proneness")
PL_plot5
###############: Pathogen load ~ Soil PCA1
PL_plot6 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = Soil_PCA1,y = PL,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.112   Marginal R2 = 0.091")+
  ylab(label = "Pathogen load")+
  xlab(label = "Soil PCA1")
PL_plot6
###############: Pathogen load ~ Mean daily temperature
PL_plot7 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = MDT,y = PL,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.329   Marginal R2 = 0.035")+
  ylab(label = "Pathogen load")+
  xlab(label = "Mean daily temperature")
PL_plot7
###############: Pathogen load ~ Mean daily humidity
PL_plot8 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = MDH,y = PL,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.259   Marginal R2 = 0.047")+
  ylab(label = "Pathogen load")+
  xlab(label = "Mean daily humidity")
PL_plot8
###############: Joint plot
ggarrange(PL_plot1,PL_plot2,PL_plot3,PL_plot4,PL_plot5,PL_plot6,PL_plot7,PL_plot8,ncol = 4,nrow = 2, labels="AUTO",common.legend = T,legend = "bottom")
ggsave("Variables_effect_on_PL.pdf",width = 16,height = 8)

#####: Linear mix effect model regression plots between various (soil, plant community, environment) variables and soil fungal pathogen OTU richness (sfpOTUs)
###############: sfpOTUs ~ Aboveground biomass
sfpOTUs_plot1 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = AB,y = sfpOTUs,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.318   Marginal R2 = 0.057")+
  ylab(label = "Soil fungal pathogen OTU richness")+
  xlab(label = "Aboveground biomass")
sfpOTUs_plot1
###############: sfpOTUs ~ Belowground biomass
sfpOTUs_plot2 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = BB,y = sfpOTUs,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.247   Marginal R2 = 0.036")+
  ylab(label = "Soil fungal pathogen OTU richness")+
  xlab(label = "Belowground biomass")
sfpOTUs_plot2
###############: sfpOTUs ~ Species richness
sfpOTUs_plot3 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = SR,y = sfpOTUs,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.851   Marginal R2 = 0.002")+
  ylab(label = "Soil fungal pathogen OTU richness")+
  xlab(label = "Species richness")
sfpOTUs_plot3
###############: sfpOTUs ~ Pielou's evenness index
sfpOTUs_plot4 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = Ev,y = sfpOTUs,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.920   Marginal R2 < 0.001")+
  ylab(label = "Soil fungal pathogen OTU richness")+
  xlab(label = "Pielou's evenness index")
sfpOTUs_plot4
###############: sfpOTUs ~ Community proneness
sfpOTUs_plot5 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = p,y = sfpOTUs,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.530   Marginal R2 = 0.009")+
  ylab(label = "Soil fungal pathogen OTU richness")+
  xlab(label = "Community proneness")
sfpOTUs_plot5
###############: sfpOTUs ~ Soil PCA1
sfpOTUs_plot6 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = Soil_PCA1,y = sfpOTUs,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.112   Marginal R2 = 0.091")+
  ylab(label = "Soil fungal pathogen OTU richness")+
  xlab(label = "Soil PCA1")
sfpOTUs_plot6
###############: sfpOTUs ~ Mean daily temperature
sfpOTUs_plot7 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = MDT,y = sfpOTUs,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  geom_smooth(aes(x = MDT,y = sfpOTUs),method = "lm",se = T,size = 1.5,linetype = 1,alpha = 0.3,color = "#3C5488FF",fill = "#3C5488FF")+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.001   Marginal R2 = 0.379")+
  ylab(label = "Soil fungal pathogen OTU richness")+
  xlab(label = "Mean daily temperature")
sfpOTUs_plot7
###############: sfpOTUs ~ Mean daily humidity
sfpOTUs_plot8 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = MDH,y = sfpOTUs,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.649   Marginal R2 = 0.023")+
  ylab(label = "Soil fungal pathogen OTU richness")+
  xlab(label = "Mean daily humidity")
sfpOTUs_plot8
###############: Joint plot
ggarrange(sfpOTUs_plot1,sfpOTUs_plot2,sfpOTUs_plot3,sfpOTUs_plot4,sfpOTUs_plot5,sfpOTUs_plot6,sfpOTUs_plot7,sfpOTUs_plot8,ncol = 4,nrow = 2, labels="AUTO",common.legend = T,legend = "bottom")
ggsave("Variables_effect_on_sfpOTUs.pdf",width = 16,height = 8)

#####: Linear mix effect model regression plots between various (soil, plant community, environment) variables and soil fungal pathogen relative abundance (sfpRA)
###############: sfpRA ~ Aboveground biomass
sfpRA_plot1 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = AB,y = sfpRA,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.418   Marginal R2 = 0.024")+
  ylab(label = "Soil fungal pathogen relative abundance")+
  xlab(label = "Aboveground biomass")
sfpRA_plot1
###############: sfpRA ~ Belowground biomass
sfpRA_plot2 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = BB,y = sfpRA,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.993   Marginal R2 < 0.001")+
  ylab(label = "Soil fungal pathogen relative abundance")+
  xlab(label = "Belowground biomass")
sfpRA_plot2
###############: sfpRA ~ Species richness
sfpRA_plot3 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = SR,y = sfpRA,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.111   Marginal R2 = 0.099")+
  ylab(label = "Soil fungal pathogen relative abundance")+
  xlab(label = "Species richness")
sfpRA_plot3
###############: sfpRA ~ Pielou's evenness index
sfpRA_plot4 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = Ev,y = sfpRA,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.360   Marginal R2 = 0.031")+
  ylab(label = "Soil fungal pathogen relative abundance")+
  xlab(label = "Pielou's evenness index")
sfpRA_plot4
###############: sfpRA ~ Community proneness
sfpRA_plot5 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = p,y = sfpRA,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.974   Marginal R2 < 0.001")+
  ylab(label = "Soil fungal pathogen relative abundance")+
  xlab(label = "Community proneness")
sfpRA_plot5
###############: sfpRA ~ Soil PCA1
sfpRA_plot6 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = Soil_PCA1,y = sfpRA,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.617   Marginal R2 = 0.009")+
  ylab(label = "Soil fungal pathogen relative abundance")+
  xlab(label = "Soil PCA1")
sfpRA_plot6
###############: sfpRA ~ Mean daily temperature
sfpRA_plot7 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = MDT,y = sfpRA,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.916   Marginal R2 < 0.001")+
  ylab(label = "Soil fungal pathogen relative abundance")+
  xlab(label = "Mean daily temperature")
sfpRA_plot7
###############: sfpRA ~ Mean daily humidity
sfpRA_plot8 <- ggplot(data = Field_study_data)+
  geom_point(aes(x = MDH,y = sfpRA,fill = Elevation),shape = 21,size = 7,alpha = 0.8)+
  scale_fill_manual(values = c("#FDDBC7", "#F4A582",
                               "#D6604D", "#B2182B", "#67001F"))+
  mytheme+
  labs(subtitle = "P = 0.893   Marginal R2 = 0.001")+
  ylab(label = "Soil fungal pathogen relative abundance")+
  xlab(label = "Mean daily humidity")
sfpRA_plot8
###############: Joint plot
ggarrange(sfpRA_plot1,sfpRA_plot2,sfpRA_plot3,sfpRA_plot4,sfpRA_plot5,sfpRA_plot6,sfpRA_plot7,sfpRA_plot8,ncol = 4,nrow = 2, labels="AUTO",common.legend = T,legend = "bottom")
ggsave("Variables_effect_on_sfpRA.pdf",width = 16,height = 8)
