library(tidyverse)
library(readxl)
library(sciplot)
library(ggrepel)
library(ggsci)
library(ggpubr)
library(metafor)
library(plotbiomes)
library(ochRe)
library(ggnewscale)

rm(list = ls())

##### ggplot2 theme #####

theme_set(theme_test())

mytheme = theme(plot.title = element_text(size = 10),
                plot.subtitle = element_text(size = 12),
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 15,
                                          face = "bold"),
                legend.title = element_text(size = 15,
                                            face = "bold"),
                legend.text = element_text(size = 12),
                legend.key.size = unit(0.5,"cm"))

##### Read data #####

dat = read_xlsx("MetaData.xlsx") %>% 
  mutate(Paper = as.factor(Paper),
         Study = as.factor(Study))
str(dat)
colnames(dat)
with(dat %>% mutate(type = paste(Response.variable,
                                 Ecosystem,
                                 sep = "|")),
     table(type))

##### Modelling #####

mod.all = rma.uni(yi,vi,
                      data = dat,
                      method = "REML")

mod.ffpOTUs = rma.uni(yi,vi,
                      data = dat %>% 
                        filter(Response.variable == "ffpOTUs"),
                      method = "REML")
mod.ffpOTUs.result = summary(mod.ffpOTUs)

mod.disease = rma.uni(yi,vi,
               data = dat %>% 
                 filter(Response.variable == "foliar.fungal.disease"),
               method = "REML")
mod.disease.result = summary(mod.disease)

mod.disease.forest = rma.uni(yi,vi,
                             data = dat %>% 
                               filter(Response.variable == "foliar.fungal.disease",
                                      Ecosystem == "forest"),
                             method = "REML")
mod.disease.forest.result = summary(mod.disease.forest)

mod.disease.grassland = rma.uni(yi,vi,
                                data = dat %>% 
                                  filter(Response.variable == "foliar.fungal.disease",
                                         Ecosystem == "grassland"),
                                method = "REML")
mod.disease.grassland.result = summary(mod.disease.grassland)

mod.sfpOTUs = rma.uni(yi,vi,
                      data = dat %>% 
                        filter(Response.variable == "sfpOTUs"),
                      method = "REML")
mod.sfpOTUs.result = summary(mod.sfpOTUs)

mod.sfpRA = rma.uni(yi,vi,
                    data = dat %>% 
                      filter(Response.variable == "sfpRA"),
                    method = "REML")
mod.sfpRA.result = summary(mod.sfpRA)

##### Plot (Overall and category) #####

overall.dat = data.frame(response = c("Foliar fungal disease",
                                      "Foliar fungal disease forest",
                                      "Foliar fungal disease grassland",
                                      "Foliar fungal pathogen OTU richness",
                                      "Soil fungal pathogen relative abundance",
                                      "Soil fungal pathogen OTU richness"),
                         mean = c(mod.disease.result$b,
                                  mod.disease.forest.result$b,
                                  mod.disease.grassland.result$b,
                                  mod.ffpOTUs.result$b,
                                  mod.sfpRA$b,
                                  mod.sfpOTUs.result$b),
                         ci.lb = c(mod.disease.result$ci.lb,
                                   mod.disease.forest.result$ci.lb,
                                   mod.disease.grassland.result$ci.lb,
                                   mod.ffpOTUs.result$ci.lb,
                                   mod.sfpRA$ci.lb,
                                   mod.sfpOTUs.result$ci.lb),
                         ci.ub = c(mod.disease.result$ci.ub,
                                   mod.disease.forest.result$ci.ub,
                                   mod.disease.grassland.result$ci.ub,
                                   mod.ffpOTUs.result$ci.ub,
                                   mod.sfpRA$ci.ub,
                                   mod.sfpOTUs.result$ci.ub),
                         p = c(mod.disease.result$pval,
                               mod.disease.forest.result$pval,
                               mod.disease.grassland.result$pval,
                               mod.ffpOTUs.result$pval,
                               mod.sfpRA$pval,
                               mod.sfpOTUs.result$pval))

overall.dat$response = factor(overall.dat$response,
                              levels = c("Foliar fungal disease",
                                         "Foliar fungal disease forest",
                                         "Foliar fungal disease grassland",
                                         "Foliar fungal pathogen OTU richness",
                                         "Soil fungal pathogen relative abundance",
                                         "Soil fungal pathogen OTU richness"))

forest.dat = rbind(data.frame(z.ub = c(0.17,0.71,-0.33,0.44,0.68,1.04,0.80,0.46,0.59,0.30,0.71,-0.04,0.26,-0.12,0.47,1.94,-0.42,0.60,0.36,-0.36,-1.05,0.06,0.70,3.46),
                              z.lb = c(-0.06,0.43,-0.69,-0.29,-0.33,-0.56,-0.80,-0.00,-0.24,-0.29,-0.67,-0.31,-0.20,-0.77,-0.18,1.16,-1.32,-2.18,-3.56,-2.62,-4.97,-0.95,-0.54,-0.46),
                              z = c(0.06,0.57,-0.51,0.08,0.18,0.24,0.00,0.23,0.18,0.01,0.02,-0.17,0.03,-0.45,0.15,1.55,-0.87,-0.79,-1.60,-1.49,-3.01,-0.44,0.08,1.50),
                              response = "Foliar fungal disease"),
                   data.frame(z.ub = c(-0.26,1.06,0.07,0.04,0.14,0.41),
                              z.lb = c(-1.02,-0.69,-0.30,-1.45,-0.17,-0.13),
                              z = c(-0.64,0.19,-0.11,-0.70,-0.01,0.14),
                              response = "Foliar fungal pathogen OTU richness"),
                   data.frame(z.ub = c(0.33,1.01,-0.02,0.14,0.59,-0.10,0.47,-0.27,0.43,0.96,-0.29,0.77,0.17,-0.40,0.18),
                              z.lb = c(-0.47,-0.07,-0.61,-0.43,-0.01,-0.86,-0.15,-0.82,-0.03,-0.05,-0.62,-0.24,-1.31,-0.67,-0.12),
                              z = c(-0.07,0.47,-0.32,-0.14,0.29,-0.48,0.16,-0.54,0.20,0.45,-0.46,0.26,-0.57,-0.54,0.03),
                              response = "Soil fungal pathogen relative abundance"),
                   data.frame(z.ub = c(0.15,0.54,-0.28,0.61,0.36,0.03,-0.49,-0.35,0.23,-0.11,0.54,-0.33,0.70,0.05,-0.65,0.06,0.15),
                              z.lb = c(-0.65,-0.55,-0.87,0.04,-0.24,-0.72,-1.12,-0.90,-0.23,-0.63,-0.48,-0.66,-0.31,-0.90,-2.13,-0.22,-0.15),
                              z = c(-0.25,-0.01,-0.57,0.33,0.06,-0.34,-0.80,-0.63,0.00,-0.37,0.03,-0.50,0.19,-0.43,-1.39,-0.08,-0.00),
                              response = "Soil fungal pathogen OTU richness")) %>% 
  arrange(response,
          z) %>% 
  mutate(ID = c(1:nrow(forest.dat)))

forest.dat$response = factor(forest.dat$response,
                             levels = c("Foliar fungal disease",
                                        "Foliar fungal pathogen OTU richness",
                                        "Soil fungal pathogen relative abundance",
                                        "Soil fungal pathogen OTU richness"))

forest.dat$sig = NA
forest.dat$sig[which(forest.dat$z.ub < 0 & forest.dat$z.lb < 0)] = "sig"
forest.dat$sig[which(forest.dat$z.ub > 0 & forest.dat$z.lb > 0)] = "sig"
forest.dat$sig[which(forest.dat$z.ub > 0 & forest.dat$z.lb < 0)] = "no.sig"

(overall.plot = ggplot(data = overall.dat,
                       aes(x = c(1,1.5,2,3,4,5),
                           y = mean,
                           ymin = ci.lb,
                           ymax = ci.ub,
                           color = response))+
    mytheme+
    geom_hline(yintercept = 0,
               linetype = 2,
               linewidth = 0.5,
               color = "gray25")+
    geom_point(shape = 16,
               size = 5)+
    geom_errorbar(linewidth = 1,
                  width = 0.1)+
    geom_text(aes(y = c(0.25,0.47,0.24,0.17,0.15,-0.02),
                  label = paste("p = ",
                                round(p,3),
                                sep = ""),
                  fontface = "bold"),
              size = 5)+
    geom_text(aes(y = c(0.33,0.55,0.32,0.25,0.23,0.06),
                  label = c("N = 24","N = 14","N = 10","N = 6","N = 15","N = 17"),
                  fontface = "bold"),
              size = 5)+
    geom_text(aes(y = c(-0.32,-0.6,-0.16,-0.45,-0.32,-0.47),
                  label = c("(Overall)","(Forest)","(Grassland)","(Overall)","(Overall)","(Overall)")),
              size = 5)+
    labs(y = "Effect size (Z)")+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank(),
          panel.border = element_rect(linewidth = 2))+
    guides(color = guide_legend(nrow = 2))+
    scale_x_continuous(expand = c(0.1,0))+
    scale_color_manual(values = c(rep("#20854EFF",3),"#0072B5FF","#E18727FF","#BC3C29FF")))

ggsave("Meta.Overall.pdf",
       width = 8,
       height = 6)

(forest.plot = ggplot()+
    mytheme+
    geom_hline(yintercept = 0,
               linetype = 2,
               linewidth = 0.5,
               color = "gray25")+
    geom_point(data = forest.dat,
               aes(x = ID,
                   y = z,
                   color = response,
                   alpha = sig),
               shape = 16,
               size = 3)+
    geom_errorbar(data = forest.dat,
                  aes(x = ID,
                      y = z,
                      ymin = z.lb,
                      ymax = z.ub,
                      color = response,
                      alpha = sig),
                  linewidth = 1,
                  width = 0)+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.title = element_blank(),
          legend.position = "bottom",
          panel.border = element_blank(),
          axis.line = element_line(size = 1))+
    guides(color = guide_legend(nrow = 2))+
    labs(y = "Effect size (Z)")+
    scale_color_manual(values = c("#20854EFF","#0072B5FF","#E18727FF","#BC3C29FF"))+
    scale_alpha_manual(values = c(0.25,1)))

ggsave("Meta.Forest.pdf",
       width = 8,
       height = 6)

##### Map #####

meta.plot.data = dat %>% 
  dplyr::select(Longitude,
                Latitude,
                Response.variable,
                n,
                yi,
                vi,
                MAT,
                MAP,
                Year,
                IF,
                ERS)
unique(meta.plot.data$Response.variable)

meta.plot.data$Response.variable[which(meta.plot.data$Response.variable == "foliar.fungal.disease")] = "Foliar fungal disease"
meta.plot.data$Response.variable[which(meta.plot.data$Response.variable == "ffpOTUs")] = "Foliar fungal pathogen OTU richness"
meta.plot.data$Response.variable[which(meta.plot.data$Response.variable == "sfpRA")] = "Soil fungal pathogen relative abundance"
meta.plot.data$Response.variable[which(meta.plot.data$Response.variable == "sfpOTUs")] = "Soil fungal pathogen OTU richness"

meta.plot.data$Response.variable = factor(meta.plot.data$Response.variable,
                                          levels = c("Foliar fungal disease",
                                                     "Foliar fungal pathogen OTU richness",
                                                     "Soil fungal pathogen relative abundance",
                                                     "Soil fungal pathogen OTU richness"))

mapworld = borders("world",
                   colour = "gray25",
                   fill = "gray90",
                   size = 0.01)

(meta.map = ggplot(data = meta.plot.data)+
    mapworld+
    mytheme+
    geom_point(aes(x = Longitude,
                   y = Latitude,
                   size = n,
                   fill = Response.variable),
               alpha = 0.75,
               shape = 21)+
    coord_cartesian(xlim = c(-180,180))+
    coord_cartesian(ylim = c(-90,90))+
    scale_x_continuous(breaks = seq(-180,180,60))+
    scale_y_continuous(breaks = seq(-90,90,30))+
    scale_fill_manual(values = c("#20854EFF","#0072B5FF","#E18727FF","#BC3C29FF"))+
    scale_size(range = c(3,6))+
    theme(legend.position = "bottom")+
    labs(x = "Longitude",
         y = "Latitude"))

data(Whittaker_biomes)
as_tibble(Whittaker_biomes)

(meta.biome = ggplot()+
    mytheme+
    geom_polygon(data = Whittaker_biomes,
                     aes(x = temp_c,
                         y = precp_cm,
                         fill = biome),
                 color = "white",
                 linewidth = 2,
                 alpha = 0.6)+
    guides(fill = guide_legend(title = "Biomes",
                               ncol = 1))+
    scale_fill_ochre()+
    new_scale_fill()+
    geom_point(data = meta.plot.data,
               aes(x = MAT,
                   y = MAP/10,
                   size = n,
                   fill = Response.variable),
               alpha = 0.65,
               shape = 21,
               show.legend = FALSE)+
    scale_fill_manual(values = c("#20854EFF","#0072B5FF","#E18727FF","#BC3C29FF"))+
    scale_size(range=c(3,8))+
    labs(x = "Mean annual temperature",
         y = "Mean annual precipitation")+
    theme(legend.background = element_blank(),
          legend.position = c(0.22,0.75)))

ggarrange(meta.map,meta.biome,
          widths = c(1.5,1),
          ncol = 2,
          align = "hv")

ggsave("Meta.Map.pdf",
       width = 15,
       height = 6)

##### Funnel plot #####

funnel(mod.disease,
       level = c(90,95,99),
       shade = c("gray45","gray60","gray75"),
       legend = TRUE,
       back = "gray90",
       lty = 2,
       col = "black",
       bg = "#20854EFF",
       cex = 2.5,
       pch = 21)

funnel(mod.ffpOTUs,
       level = c(90,95,99),
       shade = c("gray45","gray60","gray75"),
       legend = TRUE,
       back = "gray90",
       lty = 2,
       col = "black",
       bg = "#0072B5FF",
       cex = 2.5,
       pch = 21)

funnel(mod.sfpRA,
       level = c(90,95,99),
       shade = c("gray45","gray60","gray75"),
       legend = TRUE,
       back = "gray90",
       lty = 2,
       col = "black",
       bg = "#E18727FF",
       cex = 2.5,
       pch = 21)

funnel(mod.sfpOTUs,
       level = c(90,95,99),
       shade = c("gray45","gray60","gray75"),
       legend = TRUE,
       back = "gray90",
       lty = 2,
       col = "black",
       bg = "#BC3C29FF",
       cex = 2.5,
       pch = 21)

##### Bias plot #####

(meta.year = ggplot(data = meta.plot.data,
                    aes(x = Year,
                        y = yi,
                        size = n,
                        fill = Response.variable))+
   mytheme+
   geom_point(shape = 21,
              alpha = 0.5,
              show.legend = FALSE)+
   scale_fill_manual(values = c("#20854EFF","#0072B5FF","#E18727FF","#BC3C29FF"))+
   scale_size_continuous(range = c(3,10))+
   xlab("Publish year")+
   ylab("Effect size (Z)"))

(meta.IF = ggplot(data = meta.plot.data,
                    aes(x = as.numeric(IF),
                        y = yi,
                        size = n,
                        fill = Response.variable))+
    mytheme+
    geom_point(shape = 21,
               alpha = 0.5,
               show.legend = FALSE)+
    scale_fill_manual(values = c("#20854EFF","#0072B5FF","#E18727FF","#BC3C29FF"))+
    scale_size_continuous(range = c(3,10))+
    xlab("Journal impact factor")+
    ylab("Effect size (Z)"))

ggarrange(meta.year,meta.IF,
          ncol = 2,
          align = "hv")

ggsave("Meta.Bias.pdf",
       width = 8,
       height = 4)

##### Meta-regression #####

mod.null = rma.mv(yi,vi,
                  subset = (dat$Response.variable == "sfpOTUs"),
                  data = dat,
                  method = "REML",
                  mods = ~ 1,
                  random = ~1|Paper/Study)
summary(mod.null)
AICc(mod.null)

mod = rma.mv(yi,vi,
             subset = (dat$Response.variable == "sfpOTUs"),
             data = dat,
             method = "REML",
             mods = ~ MAT,
             random = ~1|Paper/Study)
summary(mod)
AICc(mod)

(meta.MAT = ggplot()+
    mytheme+
    geom_point(data = meta.plot.data,
               aes(x = MAT,
                   y = yi,
                   size = n,
                   fill = Response.variable),
               shape = 21,
               alpha = 0.5,
               show.legend = FALSE)+
    scale_size_continuous(range = c(3,10))+
    scale_fill_manual(values = c("#20854EFF","#0072B5FF","#E18727FF","#BC3C29FF"))+
    xlab("Mean annual temperature")+
    ylab("Effect size (Z)"))

(meta.MAP = ggplot()+
    mytheme+
    geom_point(data = meta.plot.data,
               aes(x = MAP,
                   y = yi,
                   size = n,
                   fill = Response.variable),
               shape = 21,
               alpha = 0.5,
               show.legend = FALSE)+
    scale_size_continuous(range = c(3,10))+
    scale_fill_manual(values = c("#20854EFF","#0072B5FF","#E18727FF","#BC3C29FF"))+
    xlab("Mean annual precipitation")+
    ylab("Effect size (Z)"))


(meta.Latitude = ggplot()+
    mytheme+
    geom_point(data = meta.plot.data,
               aes(x = abs(Latitude),
                   y = yi,
                   size = n,
                   fill = Response.variable),
               shape = 21,
               alpha = 0.5,
               show.legend = FALSE)+
    geom_smooth(data = meta.plot.data,
                aes(x = abs(Latitude),
                    y = yi,
                    color = Response.variable,
                    fill = Response.variable),
                alpha = 0.25,
                method = "lm",
                se = TRUE,
                show.legend = FALSE)+
    scale_size_continuous(range = c(3,10))+
    scale_color_manual(values = c("#20854EFF","#0072B5FF","#E18727FF","#BC3C29FF"))+
    scale_fill_manual(values = c("#20854EFF","#0072B5FF","#E18727FF","#BC3C29FF"))+
    xlab("Absolute latitude")+
    ylab("Effect size (Z)"))

(meta.ERS = ggplot()+
    mytheme+
    geom_point(data = meta.plot.data,
               aes(x = ERS,
                   y = yi,
                   size = n,
                   fill = Response.variable),
               shape = 21,
               alpha = 0.5,
               show.legend = FALSE)+
    geom_smooth(data = meta.plot.data,
                aes(x = ERS,
                    y = yi,
                    color = Response.variable,
                    fill = Response.variable),
                alpha = 0.25,
                method = "lm",
                se = TRUE,
                show.legend = FALSE)+
    scale_size_continuous(range = c(3,10))+
    scale_color_manual(values = c("#20854EFF","#0072B5FF","#E18727FF","#BC3C29FF"))+
    scale_fill_manual(values = c("#20854EFF","#0072B5FF","#E18727FF","#BC3C29FF"))+
    xlab("Elevation range of sampling")+
    ylab("Effect size (Z)"))

ggarrange(meta.MAT,meta.MAP,meta.Latitude,meta.ERS,
          nrow = 2,
          ncol = 2,
          align = "hv")

ggsave("Meta.Regression.pdf",
       height = 8,
       width = 8)







