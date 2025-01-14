# Codes related to the mixed effect modeling analyses, Figure 4 - 6 and Figure S7-S8
# creates and saves model output and figures

library(lme4)
library(sjPlot)
library(ggplot2)
library(cowplot)
library(lmerTest)
library(Hmisc)
library(dplyr)
set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'Helvetica',   #To change the font type
          axis.title.size = 1.2,  #To change axis title size
          axis.textsize.x = 1.2,  #To change x axis text size
          axis.textsize.y = 1.2)  #To change y axis text size
# reading ROI
roi <- 'DG'
# read hipp subfields data
dsall <- read.csv(paste('/schizo/results/alltrial_', roi,'_imaging.csv', sep =""))

## scale variables
dsall$bvalues <- dsall$bvalues/100

## remove NaNs
dsna <- dsall%>% 
  filter(resp!='NaN')
  
## remove > 3.5 SD for b
ds <- dsna %>% 
  group_by(condition, subid) %>% 
  filter(abs(bvalues - mean(bvalues)) < (sd(bvalues) * 3.5))
ungroup(ds)

## check the distribution
hist(ds['bvalues'],plot = TRUE)

############################### across two sessions #################################
options(warn = 2)   
conditions <- list("s","o","n")
for (val in conditions)
{
df <- ds[ds$session == "ses-01",]
df_s1 <- df[df$condition == val,]
df_s1$session<-factor(df_s1$group, levels=c('patient','control'))
smodel1<-glmer(behavior~bvalues*group+(1|subid), data=df_s1, family=binomial)
summary(smodel1)

s1 <- summary(smodel1)
capture.output(s1, file = paste('/schizo_pattern/results/model_b_',val,'_ses01_',roi,'_1.txt'))
p1 <- plot_model(smodel1, type="pred" ,terms=c("bvalues[all]","group"),colors = c("#8AAD67","#755298"),line.size =1.25, axis.title = c("Beta values (centered)","Predicted performance"), title = "", show.p = TRUE, ci.lvl = .95)
#p1 + font_size(axis_title.y = 20)
df <- ds[ds$session == "ses-02",]
df_s2 <- df[df$condition == val,]
smodel2<-glmer(behavior~bvalues*group+(1|subid), data=df_s2, family=binomial)
s2 <- summary(smodel2)
capture.output(s2, file = paste('/schizo_pattern/results/model_b_',val,'_ses01_',roi,'_2.txt'))
p2 <- plot_model(smodel2, type="pred" ,terms=c("bvalues[all]","group"),colors = c("#8AAD67","#755298"),line.size =1.25, axis.title = c("Beta values (centered)","Predicted behavioral performance"),  title = "", show.p = TRUE, ci.lvl = .95)
p4 <- plot_grid(p1, p2, labels = c('similar trails session 1','similar trails session 2'), label_size = 12)
save_plot(path="/schizo_pattern/results/",paste('pc_b_',val,'_1_',roi,'.pdf', sep =""), p4, base_aspect_ratio = 2.6)
}

##################################### across two groups #####################################
conditions <- list("s","o")
for (val in conditions)
{
df <- bothSubset[bothSubset$group == "patient",]
df_s1 <- df[df$condition == val,]
df_s1$session<-factor(df_s1$session, levels=c('ses-02','ses-01'))
smodel1<-glmer(behavior~bvalues*session+(1|subid), data=df_s1, family=binomial)
s1 <- summary(smodel1)
capture.output(s1, file = paste('/schizo_pattern/results/model_b_',val,'_ses012_',roi,'_1.txt'))
p1 <- plot_model(smodel1, type="pred" ,terms=c("bvalues[all]","session"),colors = c("#D3C2E5","#755298"),line.size =1.25, axis.title = c("Beta values (centered)","Predicted behavioral performance"), title = "", show.p = TRUE, ci.lvl = .95)

df <- bothSubset[bothSubset$group == "control",]
df_s2 <- df[df$condition == val,]
df_s2$session<-factor(df_s2$session, levels=c('ses-02','ses-01'))
smodel2<-glmer(behavior~bvalues*session+(1|subid), data=df_s2, family=binomial)
s2 <- summary(smodel2)
capture.output(s2, file = paste('/schizo_pattern/results/model_b_',val,'_ses012_',roi,'_2.txt'))
p2 <- plot_model(smodel2, type="pred" ,terms=c("bvalues[all]","session"),colors = c("#DEEFCC","#8AAD67"),line.size =1.25, axis.title = c("Beta values (centered)","Predicted behavioral performance"), title = "", show.p = TRUE, ci.lvl = .95)

p4 <- plot_grid(p1, p2, labels = c('similar trails patients','similar trails controls'), label_size = 12)
save_plot(path="/schizo_pattern/results/",paste('pc_b_',val,'_12_', roi,'.pdf', sep =""), p4, base_aspect_ratio = 2.6)
}


