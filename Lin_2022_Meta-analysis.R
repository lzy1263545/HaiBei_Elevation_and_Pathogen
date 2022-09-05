########## Analysis part ##########

library(metafor)
library(forestplot)
library(AICcmodavg)
library(lme4)
library(MuMIn)
library(nlme)
library(readxl)
library(tidyverse)

#####: Read the data of Meta-analysis
es.all <- read_xlsx("MetaData.xlsx")
#####: Set "Paper" and "Study" as factor to create ramdom effect
es.all$Paper <- as.factor(es.all$Paper)
es.all$Study <- as.factor(es.all$Study)
#####: Create subset for above- (foliar) and belowground (soil) data
es.above <- es.all %>% filter(Position == "foliage")
es.below <- es.all %>% filter(Position == "soil")

#####: Calculate overall effect of elevation on the foliar fungal disease based on random effect model
random1 <- rma.uni(yi,vi,subset = (es.above$Response.variable == "foliar.fungal.disease"),data = es.above,method = "REML")
summary(random1)
forest(random1)
funnel(random1)
radial(random1)
###############: Subset - Forest ecosystem
es.foliar.fungal.disease = es.all %>% filter(Response.variable == "foliar.fungal.disease")
random1_Forest <- rma.uni(yi,vi,subset = (es.foliar.fungal.disease$Ecosystem == "forest"),data = es.foliar.fungal.disease,method = "REML")
summary(random1_Forest)
###############: Subset - Grassland ecosystem
random1_Grassland <- rma.uni(yi,vi,subset = (es.foliar.fungal.disease$Ecosystem == "grassland"),data = es.foliar.fungal.disease,method = "REML")
summary(random1_Grassland)

#####: Calculate overall effect of elevation on the richness of foliar fungal pathogen (ffpOTUs) based on random effect model
random2 <- rma.uni(yi,vi,subset = (es.above$Response.variable == "ffpOTUs"),data = es.above,method = "REML")
summary(random2)
forest(random2)
funnel(random2)
radial(random2)

#####: Calculate overall effect of elevation on the richness of soil fungal pathogen (sfpOTUs) based on random effect model
random3 <- rma.uni(yi,vi,subset = (es.below$Response.variable == "sfpOTUs"),data = es.below,method = "REML")
summary(random3)
fsn(yi,vi,data = es.below %>% filter(Response.variable == "sfpOTUs"),type = "Rosenthal",alpha = 0.05, digits = 4)#: calculate the fail-safe number for the significant effect size
forest(random3)
funnel(random3)
radial(random3)

#####: Calculate overall effect of elevation on the relative abundance of soil fungal pathogen (sfpRA) based on random effect model
random4 <- rma.uni(yi,vi,subset = (es.below$Response.variable == "sfpRA"),data = es.below,method = "REML")
summary(random4)
forest(random4)
funnel(random4)
radial(random4)

#####: Test respective influence of multiple vriables on the elevational effect size on foliar fungal disease based on random effect model
###############: Mean Annual Temperature (MAT)
random1.1 <- rma.mv(yi,vi,subset = (es.above$Response.variable == "foliar.fungal.disease"), data = es.above, method = "REML", mods = ~MAT, random = ~1|Paper/Study)
summary(random1.1)
funnel(random1.1)
AICc(random1.1)
###############: Mean Annual Precipitation (MAP)
random1.2 <- rma.mv(yi,vi,subset = (es.above$Response.variable == "foliar.fungal.disease"), data = es.above, method = "REML", mods = ~MAP, random = ~1|Paper/Study)
summary(random1.2)
funnel(random1.2)
AICc(random1.2)
###############: Latitude
random1.3 <- rma.mv(yi,vi,subset = (es.above$Response.variable == "foliar.fungal.disease"), data = es.above, method = "REML", mods = ~abs(Latitude), random = ~1|Paper/Study)
summary(random1.3)
funnel(random1.3)
AICc(random1.3)
###############: Elevational Range of Sampling (ERS)
random1.4 <- rma.mv(yi,vi,subset = (es.above$Response.variable == "foliar.fungal.disease"), data = es.above, method = "REML", mods = ~ERS, random = ~1|Paper/Study)
summary(random1.4)
funnel(random1.4)
AICc(random1.4)
###############: Publish Year (Year)
random1.5 <- rma.mv(yi,vi,subset = (es.above$Response.variable == "foliar.fungal.disease"), data = es.above, method = "REML", mods = ~Year, random = ~1|Paper/Study)
summary(random1.5)
funnel(random1.5)
AICc(random1.5)
###############: Journal Impact Factor (IF)
random1.6 <- rma.mv(yi,vi,subset = (es.above$Response.variable == "foliar.fungal.disease"), data = es.above, method = "REML", mods = ~as.numeric(IF), random = ~1|Paper/Study)
summary(random1.6)
funnel(random1.6)
AICc(random1.6)
###############: Null model
random1.7 <- rma.mv(yi,vi,subset = (es.above$Response.variable == "foliar.fungal.disease"), data = es.above, method = "REML", mods = ~1, random = ~1|Paper/Study)
summary(random1.7)
funnel(random1.7)
AICc(random1.7)

#####: Test respective influence of multiple vriables on the elevational effect size on the richness of foliar fungal pathogen (ffpOTUs) based on random effect model
###############: Mean Annual Temperature (MAT)
random2.1 <- rma.mv(yi,vi,subset = (es.above$Response.variable == "ffpOTUs"), data = es.above, method = "REML", mods = ~MAT, random = ~1|Paper/Study)
summary(random2.1)
funnel(random2.1)
AICc(random2.1)
###############: Mean Annual Precipitation (MAP)
random2.2 <- rma.mv(yi,vi,subset = (es.above$Response.variable == "ffpOTUs"), data = es.above, method = "REML", mods = ~MAP, random = ~1|Paper/Study)
summary(random2.2)
funnel(random2.2)
AICc(random2.2)
###############: Latitude
random2.3 <- rma.mv(yi,vi,subset = (es.above$Response.variable == "ffpOTUs"), data = es.above, method = "REML", mods = ~abs(Latitude), random = ~1|Paper/Study)
summary(random2.3)
funnel(random2.3)
AICc(random2.3)
###############: Elevational Range of Sampling (ERS)
random2.4 <- rma.mv(yi,vi,subset = (es.above$Response.variable == "ffpOTUs"), data = es.above, method = "REML", mods = ~ERS, random = ~1|Paper/Study)
summary(random2.4)
funnel(random2.4)
AICc(random2.4)
###############: Publish Year (Year)
random2.5 <- rma.mv(yi,vi,subset = (es.above$Response.variable == "ffpOTUs"), data = es.above, method = "REML", mods = ~Year, random = ~1|Paper/Study)
summary(random2.5)
funnel(random2.5)
AICc(random2.5)
###############: Journal Impact Factor (IF)
random2.6 <- rma.mv(yi,vi,subset = (es.above$Response.variable == "ffpOTUs"), data = es.above, method = "REML", mods = ~as.numeric(IF), random = ~1|Paper/Study)
summary(random2.6)
funnel(random2.6)
AICc(random2.6)
###############: Null model
random2.7 <- rma.mv(yi,vi,subset = (es.above$Response.variable == "ffpOTUs"), data = es.above, method = "REML", mods = ~1, random = ~1|Paper/Study)
summary(random2.7)
funnel(random2.7)
AICc(random2.7)

#####: Test respective influence of multiple vriables on the elevational effect size on the richness of soil fungal pathogen (sfpOTUs) based on random effect model
###############: Mean Annual Temperature (MAT)
random3.1 <- rma.mv(yi,vi,subset = (es.below$Response.variable == "sfpOTUs"), data = es.below, method = "REML", mods = ~MAT, random = ~1|Paper/Study)
summary(random3.1)
funnel(random3.1)
AICc(random3.1)
###############: Mean Annual Precipitation (MAP)
random3.2 <- rma.mv(yi,vi,subset = (es.below$Response.variable == "sfpOTUs"), data = es.below, method = "REML", mods = ~MAP, random = ~1|Paper/Study)
summary(random3.2)
funnel(random3.2)
AICc(random3.2)
###############: Latitude
random3.3 <- rma.mv(yi,vi,subset = (es.below$Response.variable == "sfpOTUs"), data = es.below, method = "REML", mods = ~abs(Latitude), random = ~1|Paper/Study)
summary(random3.3)
funnel(random3.3)
AICc(random3.3)
###############: Elevational Range of Sampling (ERS)
random3.4 <- rma.mv(yi,vi,subset = (es.below$Response.variable == "sfpOTUs"), data = es.below, method = "REML", mods = ~ERS, random = ~1|Paper/Study)
summary(random3.4)
funnel(random3.4)
AICc(random3.4)
###############: Publish Year (Year)
random3.5 <- rma.mv(yi,vi,subset = (es.below$Response.variable == "sfpOTUs"), data = es.below, method = "REML", mods = ~Year, random = ~1|Paper/Study)
summary(random3.5)
funnel(random3.5)
AICc(random3.5)
###############: Journal Impact Factor (IF)
random3.6 <- rma.mv(yi,vi,subset = (es.below$Response.variable == "sfpOTUs"), data = es.below, method = "REML", mods = ~as.numeric(IF), random = ~1|Paper/Study)
summary(random3.6)
funnel(random3.6)
AICc(random3.6)
###############: Null model
random3.7 <- rma.mv(yi,vi,subset = (es.below$Response.variable == "sfpOTUs"), data = es.below, method = "REML", mods = ~1, random = ~1|Paper/Study)
summary(random3.7)
funnel(random3.7)
AICc(random3.7)

#####: Test respective influence of multiple vriables on the elevational effect size on the relative abundance of soil fungal pathogen (sfpRA) based on random effect model
###############: Mean Annual Temperature (MAT)
random4.1 <- rma.mv(yi,vi,subset = (es.below$Response.variable == "sfpRA"), data = es.below, method = "REML", mods = ~MAT, random = ~1|Paper/Study)
summary(random4.1)
funnel(random4.1)
AICc(random4.1)
###############: Mean Annual Precipitation (MAP)
random4.2 <- rma.mv(yi,vi,subset = (es.below$Response.variable == "sfpRA"), data = es.below, method = "REML", mods = ~MAP, random = ~1|Paper/Study)
summary(random4.2)
funnel(random4.2)
AICc(random4.2)
###############: Latitude
random4.3 <- rma.mv(yi,vi,subset = (es.below$Response.variable == "sfpRA"), data = es.below, method = "REML", mods = ~abs(Latitude), random = ~1|Paper/Study)
summary(random4.3)
funnel(random4.3)
AICc(random4.3)
###############: Elevational Range of Sampling (ERS)
random4.4 <- rma.mv(yi,vi,subset = (es.below$Response.variable == "sfpRA"), data = es.below, method = "REML", mods = ~ERS, random = ~1|Paper/Study)
summary(random4.4)
funnel(random4.4)
AICc(random4.4)
###############: Publish Year (Year)
random4.5 <- rma.mv(yi,vi,subset = (es.below$Response.variable == "sfpRA"), data = es.below, method = "REML", mods = ~Year, random = ~1|Paper/Study)
summary(random4.5)
funnel(random4.5)
AICc(random4.5)
###############: Journal Impact Factor (IF)
random4.6 <- rma.mv(yi,vi,subset = (es.below$Response.variable == "sfpRA"), data = es.below, method = "REML", mods = ~as.numeric(IF), random = ~1|Paper/Study)
summary(random4.6)
funnel(random4.6)
AICc(random4.6)
###############: Null model
random4.7 <- rma.mv(yi,vi,subset = (es.below$Response.variable == "sfpRA"), data = es.below, method = "REML", mods = ~1, random = ~1|Paper/Study)
summary(random4.7)
funnel(random4.7)
AICc(random4.7)

#####: Calculate Kendall's rank correlation for bias test (test for funnel plot asymmetry)
ranktest(random1)
ranktest(random2)
ranktest(random3)
ranktest(random4)
ranktest(random1_Forest)
ranktest(random1_Grassland)

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
library(plotbiomes)

#####: Create subsets based on different variables
es.foliar.fungal.disease = es.all %>% filter(Response.variable == "foliar.fungal.disease")
es.ffpOTUs = es.all %>% filter(Response.variable == "ffpOTUs")
es.sfpOTUs = es.all %>% filter(Response.variable == "sfpOTUs")
es.sfpRA = es.all %>% filter(Response.variable == "sfpRA")

#####: Plot the map and Whittaker's biomes to show the distribution of studies in Meta-analysis
###############: Map
theme_set(theme_light())
mymap <- NULL
mapworld <- borders("world",colour = "gray40",fill = "gray90",size = 0.01)
mymap <- ggplot(es.all)+
  mapworld+
  geom_point(aes(x = Longitude,y = Latitude,size = n,fill = Response.variable),alpha = 0.45,shape = 21)+
  coord_cartesian(xlim = c(-180,180))+
  coord_cartesian(ylim = c(-90,90))+
  scale_x_continuous(breaks = seq(-180,180,60))+
  scale_y_continuous(breaks = seq(-90,90,30))+
  scale_fill_manual(values = c("#0072B5FF","#20854EFF","#E18727FF","#BC3C29FF"))+
  scale_size(range=c(4,8))+
  theme(legend.position = "bottom")
mymap
#ggsave("map.pdf",width = 9.5,height = 6)
###############: Whittaker's biomes
theme_set(theme_test())
data(Whittaker_biomes)
as_tibble(Whittaker_biomes)
biome <- ggplot(data = Whittaker_biomes,aes(x = temp_c,y = precp_cm,fill = biome,alpha = 0.1))+
  geom_polygon(color = "white",size = 1.2)+
  scale_fill_ochre()+
  geom_point(aes(x = MAT,y = MAP/10,size = n),fill = "#20854EFF",alpha = 0.65,shape = 21,data = es.foliar.fungal.disease)+
  geom_point(aes(x = MAT,y = MAP/10,size = n),fill = "#0072B5FF",alpha = 0.65,shape = 21,data = es.ffpOTUs)+
  geom_point(aes(x = MAT,y = MAP/10,size = n),fill = "#E18727FF",alpha = 0.65,shape = 21,data = es.sfpOTUs)+
  geom_point(aes(x = MAT,y = MAP/10,size = n),fill = "#BC3C29FF",alpha = 0.65,shape = 21,data = es.sfpRA)+
  scale_size(range=c(3,8))+
  theme(legend.background = element_blank(),
        legend.position = c(0.25,0.49))
biome
#ggsave("biomes.pdf",width = 6,height = 6)

#####: Plot the overall effect of elevation on above- and belowground plant pathogens
theme_set(theme_test())
overall.es <- ggplot()+
  geom_pointrange(aes(x = c(0.7,2.3,5.3,7.3,3,3.7),
                      y = c(-0.157,-0.047,-0.257,-0.101,-0.127,0.023),ymin = c(-0.419,-0.289,-0.429,-0.283,-0.567,-0.130),ymax = c(0.106,0.195,-0.085,0.082,0.313,0.175)),
                  size = c(2,2,2,2,1,1),color = c("#0072B5FF","#20854EFF","#E18727FF","#BC3C29FF","#20854EFF","#20854EFF"),shape = c(16,16,16,16,15,15),alpha = 1)+
  coord_cartesian(ylim = c(-0.75,0.75),xlim = c(0,8))+
  scale_y_continuous(breaks = seq(-0.75,0.75,0.3))+
  geom_hline(yintercept = 0,color = "#3C5488FF",size = 0.75, alpha = 0.4,linetype = 2)+
  theme(panel.border = element_rect(fill = NA,color = "black",size = 0.5),
        axis.text.x = element_text(color="black",size = 13),
        axis.text.y = element_text(color="black",size = 13),
        axis.ticks.y = element_line(color="black",size=0.5,lineend = 2),
        axis.ticks.x = element_line(color="black",size=0.5,lineend = 2))+
  xlab("Position")+ylab("Effect size")+
  theme(plot.margin = unit(c(1,1,0.5,0.5),"cm"))
overall.es
#ggsave("Overall.es.pdf",width = 10,height = 7)

#####: Plot the forest plot of effect of elevation on above- and belowground plant pathogens
theme_set(theme_test())
###############: Create data frame of upper and lower boader and mean value of effect size
forest.foliar.fungal.disease <- data.frame(z.ub = c(0.17,0.71,-0.33,0.44,0.68,1.04,0.80,0.46,0.59,0.30,0.71,-0.04,0.26,-0.12,0.47,1.94,-0.42,0.60,0.36,-0.36,-1.05,0.06,0.70,3.46),z.lb = c(-0.06,0.43,-0.69,-0.29,-0.33,-0.56,-0.80,-0.00,-0.24,-0.29,-0.67,-0.31,-0.20,-0.77,-0.18,1.16,-1.32,-2.18,-3.56,-2.62,-4.97,-0.95,-0.54,-0.46),z = c(0.06,0.57,-0.51,0.08,0.18,0.24,0.00,0.23,0.18,0.01,0.02,-0.17,0.03,-0.45,0.15,1.55,-0.87,-0.79,-1.60,-1.49,-3.01,-0.44,0.08,1.50)) %>% arrange(z)
forest.ffpOTUs <- data.frame(z.ub = c(-0.26,1.06,0.07,0.04,0.14,0.41),z.lb = c(-1.02,-0.69,-0.30,-1.45,-0.17,-0.13),z = c(-0.64,0.19,-0.11,-0.70,-0.01,0.14)) %>% arrange(z)
forest.sfpOTUs <- data.frame(z.ub = c(0.15,0.54,-0.28,0.61,0.36,0.03,-0.49,-0.35,0.23,-0.11,0.54,-0.33,0.70,0.05,-0.65,0.06,0.15),z.lb = c(-0.65,-0.55,-0.87,0.04,-0.24,-0.72,-1.12,-0.90,-0.23,-0.63,-0.48,-0.66,-0.31,-0.90,-2.13,-0.22,-0.15),z = c(-0.25,-0.01,-0.57,0.33,0.06,-0.34,-0.80,-0.63,0.00,-0.37,0.03,-0.50,0.19,-0.43,-1.39,-0.08,-0.00)) %>% arrange(z)
forest.sfpRA <- data.frame(z.ub = c(0.33,1.01,-0.02,0.14,0.59,-0.10,0.47,-0.27,0.43,0.96,-0.29,0.77,0.17,-0.40,0.18),z.lb = c(-0.47,-0.07,-0.61,-0.43,-0.01,-0.86,-0.15,-0.82,-0.03,-0.05,-0.62,-0.24,-1.31,-0.67,-0.12),z = c(-0.07,0.47,-0.32,-0.14,0.29,-0.48,0.16,-0.54,0.20,0.45,-0.46,0.26,-0.57,-0.54,0.03)) %>% arrange(z)
forest.es <- ggplot()+
  geom_pointrange(aes(x = c(3.5,20.5,43,61),
                      y = c(-0.157,-0.047,-0.257,-0.101),ymin = c(-0.419,-0.289,-0.429,-0.283),ymax = c(0.106,0.195,-0.085,0.082)),
                  size = c(0.8,0.8,0.8,0.8),color = c("#0072B5FF","#20854EFF","#E18727FF","#BC3C29FF"),shape = c(16,16,16,16),alpha = 1)+
  geom_pointrange(data = forest.ffpOTUs,aes(x = c(1:6),
                      y = z,ymin = z.lb,ymax = z.ub),
                  size = 0.3,color = "#0072B5FF",alpha = 0.4)+
  geom_pointrange(data = forest.foliar.fungal.disease,aes(x = c(9:32),
                                            y = z,ymin = z.lb,ymax = z.ub),
                  size = 0.3,color = "#20854EFF",alpha = 0.4)+
  geom_pointrange(data = forest.sfpOTUs,aes(x = c(35:51),
                                            y = z,ymin = z.lb,ymax = z.ub),
                  size = 0.3,color = "#E18727FF",alpha = 0.4)+
  geom_pointrange(data = forest.sfpRA,aes(x = c(54:68),
                                            y = z,ymin = z.lb,ymax = z.ub),
                  size = 0.3,color = "#BC3C29FF",alpha = 0.4)+
  coord_cartesian(ylim = c(-5,4))+
  scale_y_continuous(breaks = seq(-5,4,1))+
  geom_hline(yintercept = 0,color = "gray40",size = 0.2, alpha = 0.4,linetype = 2)+
  geom_vline(xintercept = 7.5,color = "black",size = 0.3, alpha = 0.4,linetype = 1)+
  geom_vline(xintercept = 33.5,color = "black",size = 0.3, alpha = 0.4,linetype = 1)+
  geom_vline(xintercept = 52.5,color = "black",size = 0.3, alpha = 0.4,linetype = 1)+
  theme(panel.border = element_rect(fill = NA,color = "black",size = 0.5),
        axis.text.y = element_text(color="black",size = 13),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color="black",size=0.5,lineend = 2))+
  xlab("Position")+ylab("Effect size")+
  theme(plot.margin = unit(c(1,1,0.5,0.5),"cm"))
forest.es
#ggsave("Forest.pdf",width = 10,height = 10)

#####: Meta-regressions between multiple variables and effect size
theme_set(theme_light())
###############: Mean Annual Temperature (MAT)
MAT.yi <- ggplot()+
  geom_point(data = es.ffpOTUs,aes(x = MAT, y = yi,size = n),fill = "#0072B5FF",color = "black",shape = 21,alpha = 0.5)+
  geom_point(data = es.foliar.fungal.disease,aes(x = MAT, y = yi,size = n),fill = "#20854EFF",color = "black",shape = 21,alpha = 0.5)+
  geom_point(data = es.sfpOTUs,aes(x = MAT, y = yi,size = n),fill = "#E18727FF",color = "black",shape = 21,alpha = 0.5)+
  geom_point(data = es.sfpRA,aes(x = MAT, y = yi,size = n),fill = "#BC3C29FF",color = "black",shape = 21,alpha = 0.5)+
  scale_size_continuous(range = c(3,10))+
  xlab("MAT")+ylab("Effect size")+
  coord_cartesian(xlim = c(-1,26))+
  scale_x_continuous(breaks = seq(-1,26,4.5))+
  theme(panel.border = element_rect(fill = NA,color = "black",size = 0.5),
    plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
MAT.yi
###############: Mean Annual Precipitation (MAP)
MAP.yi <- ggplot()+
  geom_point(data = es.ffpOTUs,aes(x = MAP, y = yi,size = n),fill = "#0072B5FF",color = "black",shape = 21,alpha = 0.5)+
  geom_point(data = es.foliar.fungal.disease,aes(x = MAP, y = yi,size = n),fill = "#20854EFF",color = "black",shape = 21,alpha = 0.5)+
  geom_point(data = es.sfpOTUs,aes(x = MAP, y = yi,size = n),fill = "#E18727FF",color = "black",shape = 21,alpha = 0.5)+
  geom_point(data = es.sfpRA,aes(x = MAP, y = yi,size = n),fill = "#BC3C29FF",color = "black",shape = 21,alpha = 0.5)+
  scale_size_continuous(range = c(3,10))+
  xlab("MAP")+ylab("Effect size")+
  coord_cartesian(xlim = c(0,2700))+
  scale_x_continuous(breaks = seq(0,2700,450))+
  theme(panel.border = element_rect(fill = NA,color = "black",size = 0.5),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
MAP.yi
###############: Latitude
Latitude.yi <- ggplot()+
  geom_point(data = es.ffpOTUs,aes(x = abs(Latitude), y = yi,size = n),fill = "#0072B5FF",color = "black",shape = 21,alpha = 0.5)+
  geom_point(data = es.foliar.fungal.disease,aes(x = abs(Latitude), y = yi,size = n),fill = "#20854EFF",color = "black",shape = 21,alpha = 0.5)+
  geom_point(data = es.sfpOTUs,aes(x = abs(Latitude), y = yi,size = n),fill = "#E18727FF",color = "black",shape = 21,alpha = 0.5)+
  geom_point(data = es.sfpRA,aes(x = abs(Latitude), y = yi,size = n),fill = "#BC3C29FF",color = "black",shape = 21,alpha = 0.5)+
  geom_smooth(data = es.ffpOTUs,aes(x = abs(Latitude), y = yi,weight = 1/vi),fill = "#0072B5FF",color = "#0072B5FF",method = "lm",size = 0.5,alpha = 0.25)+ #: Add regression line for significant result
  scale_size_continuous(range = c(3,10))+
  xlab("Latitude")+ylab("Effect size")+
  coord_cartesian(xlim = c(0,70))+
  scale_x_continuous(breaks = seq(0,70,10))+
  theme(panel.border = element_rect(fill = NA,color = "black",size = 0.5),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
Latitude.yi
###############: Elevational Range of Sampling (ERS)
ERS.yi <- ggplot()+
  geom_point(data = es.ffpOTUs,aes(x = ERS, y = yi,size = n),fill = "#0072B5FF",color = "black",shape = 21,alpha = 0.5)+
  geom_point(data = es.foliar.fungal.disease,aes(x = ERS, y = yi,size = n),fill = "#20854EFF",color = "black",shape = 21,alpha = 0.5)+
  geom_point(data = es.sfpOTUs,aes(x = ERS, y = yi,size = n),fill = "#E18727FF",color = "black",shape = 21,alpha = 0.5)+
  geom_point(data = es.sfpRA,aes(x = ERS, y = yi,size = n),fill = "#BC3C29FF",color = "black",shape = 21,alpha = 0.5)+
  geom_smooth(data = es.sfpOTUs,aes(x = ERS, y = yi,weight = 1/vi),fill = "#E18727FF",color = "#E18727FF",method = "lm",size = 0.5,alpha = 0.25)+ #: Add regression line for significant result
  geom_smooth(data = es.sfpRA,aes(x = ERS, y = yi,weight = 1/vi),fill = "#BC3C29FF",color = "#BC3C29FF",method = "lm",size = 0.5,alpha = 0.25)+ #: Add regression line for significant result
  scale_size_continuous(range = c(3,10))+
  xlab("ERS")+ylab("Effect size")+
  coord_cartesian(xlim = c(0,3500))+
  scale_x_continuous(breaks = seq(0,3500,700))+
  theme(panel.border = element_rect(fill = NA,color = "black",size = 0.5),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
ERS.yi
###############: Joint plot
ggarrange(MAT.yi,MAP.yi,Latitude.yi,ERS.yi,nrow = 2,ncol = 2,legend = "none")
#ggsave("Meta-regressions.pdf",width = 10,height = 8)

#####: Meta-regressions between publish year and journal impact factor and effect size - publish bias test
###############: Publish year
Publish.year.yi <- ggplot()+
  geom_point(data = es.ffpOTUs,aes(x = Year, y = yi,size = n),fill = "#0072B5FF",color = "black",shape = 21,alpha = 0.5)+
  geom_point(data = es.foliar.fungal.disease,aes(x = Year, y = yi,size = n),fill = "#20854EFF",color = "black",shape = 21,alpha = 0.5)+
  geom_point(data = es.sfpOTUs,aes(x = Year, y = yi,size = n),fill = "#E18727FF",color = "black",shape = 21,alpha = 0.5)+
  geom_point(data = es.sfpRA,aes(x = Year, y = yi,size = n),fill = "#BC3C29FF",color = "black",shape = 21,alpha = 0.5)+
  scale_size_continuous(range = c(3,10))+
  xlab("Publish year")+ylab("Effect size")+
  coord_cartesian(xlim = c(1992,2022))+
  scale_x_continuous(breaks = seq(1992,2022,10))+
  theme(panel.border = element_rect(fill = NA,color = "black",size = 0.5),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
Publish.year.yi
###############: Journal Impact Factor (IF)
IF.yi <- ggplot()+
  geom_point(data = es.ffpOTUs,aes(x = as.numeric(IF), y = yi,size = n),fill = "#0072B5FF",color = "black",shape = 21,alpha = 0.5)+
  geom_point(data = es.foliar.fungal.disease,aes(x = as.numeric(IF), y = yi,size = n),fill = "#20854EFF",color = "black",shape = 21,alpha = 0.5)+
  geom_point(data = es.sfpOTUs,aes(x = as.numeric(IF), y = yi,size = n),fill = "#E18727FF",color = "black",shape = 21,alpha = 0.5)+
  geom_point(data = es.sfpRA,aes(x = as.numeric(IF), y = yi,size = n),fill = "#BC3C29FF",color = "black",shape = 21,alpha = 0.5)+
  scale_size_continuous(range = c(3,10))+
  xlab("IF")+ylab("Effect size")+
  coord_cartesian(xlim = c(0,12))+
  scale_x_continuous(breaks = seq(0,12,4))+
  theme(panel.border = element_rect(fill = NA,color = "black",size = 0.5),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
IF.yi
###############: Joint plot
ggarrange(Publish.year.yi,IF.yi,nrow = 2,ncol = 1,legend = "bottom")
#ggsave("Bias.test.pdf",width = 6,height = 9)

#####: Funnel plot
###############: For different response variables
funnel(random1,level = c(90,95,99),shade = c("#3C5488FF","#197EC0FF","#71D0F5FF"),legend = T,back = "gray95",lty = 2,col = "black",bg = "#E64B35FF",cex = 2.5,pch = 21)
funnel(random2,level = c(90,95,99),shade = c("#3C5488FF","#197EC0FF","#71D0F5FF"),legend = T,back = "gray95",lty = 2,col = "black",bg = "#E64B35FF",cex = 2.5,pch = 21)
funnel(random3,level = c(90,95,99),shade = c("#3C5488FF","#197EC0FF","#71D0F5FF"),legend = T,back = "gray95",lty = 2,col = "black",bg = "#E64B35FF",cex = 2.5,pch = 21)
funnel(random4,level = c(90,95,99),shade = c("#3C5488FF","#197EC0FF","#71D0F5FF"),legend = T,back = "gray95",lty = 2,col = "black",bg = "#E64B35FF",cex = 2.5,pch = 21)
###############: For all studies
random0 <- rma.uni(yi,vi,data = es.all,method = "REML")
funnel(random0,level = c(90,95,99),shade = c("#3C5488FF","#197EC0FF","#71D0F5FF"),legend = T,back = "gray95",lty = 2,col = "black",bg = "#E64B35FF",cex = 2.5,pch = 21)

