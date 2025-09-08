#load all packages
require(survival)
require(dplyr)
require(survminer)
require(ggplot2)
require(MASS)
require(viridis)
require(extrafont)
require(scales)
require(cowplot)
require(gridGraphics)
loadfonts(device = "win")

#load datasets
visit_dates <- read.csv("C:/Users/charl/Downloads/School/visit_status_impute_deid.csv")
tobacco_wide <- read_sas("C:/Users/charl/Downloads/School/cenic1p2_wide_2019_01_07.sas7bdat", NULL)
tobacco_long <- read.csv("C:/Users/charl/Downloads/School/ivr_tlfb_deidentified-2024-08-21.csv")


#create survival analysis datasets
tobacco_long[-c(10:61)]->tbclongsub

#defining abstinence variables. indicator for if previous (3 or 7) days have a total of 0 cpd. 
for (i in 7:nrow(tbclongsub)) {
  ifelse(tbclongsub$days_from_rand[i]>=7&&0==sum(tbclongsub[(i-6):i,11]),1,0)->tbclongsub$abst_7day[i]
}

tbclongsub$abst_3day<-NA

for (i in 3:nrow(tbclongsub)) {
  ifelse(tbclongsub$days_from_rand[i]>=3&&0==sum(tbclongsub[(i-2):i,11]),1,0)->tbclongsub$abst_3day[i]
}

tbclongsub$abst_7day->tbclongsub$abst_7daynoNA
tbclongsub[is.na(tbclongsub$abst_7day),17]<-0

tbclongsub$abst_3day->tbclongsub$abst_3daynoNA
tbclongsub[is.na(tbclongsub$abst_7day),18]<-0

#the trial ends at 20 weeks, so we truncate the data based on each patient's week 20 visit day relative to randomization
tbcRC<-data.frame()

for (i in 1:1250) {
  check<- tbclongsub %>% 
    filter(tbclongsub$studyid==visit_dates[i,1]&tbclongsub$days_from_rand<visit_dates[i,12]) 
  rbind(check,tbcRC)->tbcRC
}

tbcRC$site<-substr(tbcRC$studyid,3,3)

#defining each dataset. in order is 7 day abstinence complete case, 3 day complete case, 7 day imputed and 3 day imputed.

event<-tbcRC %>% 
  filter(abst_7day==1) %>% 
  group_by(studyid) %>% 
  slice_min(days_from_rand,n=1,with_ties = FALSE) %>% 
  ungroup()
noevent<-tbcRC %>% 
  filter(!abst_7day==1&!studyid %in% event$studyid) %>% 
  group_by(studyid) %>% 
  slice_max(days_from_rand,n=1,with_ties = FALSE) %>% 
  ungroup()
rbind(event,noevent)->tbc7dayCC

event<-tbcRC %>% 
  filter(tbcRC$abst_3day==1) %>% 
  group_by(studyid) %>% 
  slice(1) %>% 
  ungroup()

noevent<-tbcRC %>% 
  filter(!abst_3day==1&!studyid %in% event$studyid) %>% 
  group_by(studyid) %>% 
  slice_max(days_from_rand,n=1,with_ties = FALSE) %>% 
  ungroup()

rbind(event,noevent)->tbc3dayCC

event<-tbcRC %>% 
  filter(tbcRC$abst_3daynoNA==1) %>% 
  group_by(studyid) %>% 
  slice(1) %>% 
  ungroup()

noevent<-tbcRC %>% 
  filter(!abst_3daynoNA==1&!studyid %in% event$studyid) %>% 
  group_by(studyid) %>% 
  slice_max(days_from_rand,n=1,with_ties = FALSE) %>% 
  ungroup()
rbind(event,noevent)->tbc3daynoNA

event<-tbcRC %>% 
  filter(tbcRC$abst_7daynoNA==1) %>% 
  group_by(studyid) %>% 
  slice(1) %>% 
  ungroup()

noevent<-tbcRC %>% 
  filter(!abst_7daynoNA==1&!studyid %in% event$studyid) %>% 
  group_by(studyid) %>% 
  slice_max(days_from_rand,n=1,with_ties = FALSE) %>% 
  ungroup()
rbind(event,noevent)->tbc7daynoNA

#site variable, not included in long data
merge(tobacco_wide[,c(2,7)],tbc3dayCC, by = 'studyid')->tbc3dayCC
merge(tobacco_wide[,c(2,7)],tbc7dayCC, by = 'studyid')->tbc7dayCC
merge(tobacco_wide[,c(2,7)],tbc3daynoNA, by = 'studyid')->tbc3daynoNA
merge(tobacco_wide[,c(2,7)],tbc7daynoNA, by = 'studyid')->tbc7daynoNA

list<-list(tbc3daynoNA=tbc3daynoNA,tbc3dayCC=tbc3dayCC,tbc7daynoNA=tbc7daynoNA,tbc7dayCC=tbc7dayCC)

mapply(write.csv,list,file=paste0(names(list), '.csv'))

tbc3daynoNA<-tbc3daynoNA%>% 
  mutate(Treatment = factor(Treatment, levels = c("B", "A", "C")))
tbc3dayCC<-tbc3dayCC%>% 
  mutate(Treatment = factor(Treatment, levels = c("B", "A", "C")))
tbc7daynoNA<-tbc7daynoNA%>% 
  mutate(Treatment = factor(Treatment, levels = c("B", "A", "C")))
tbc7dayCC<-tbc7dayCC%>% 
  mutate(Treatment = factor(Treatment, levels = c("B", "A", "C")))

#fitting survival curves
survfit(Surv(days_from_rand,abst_7day)~Treatment,data=tbc7dayCC)->km7dayCC
survfit(Surv(days_from_rand,abst_3day)~Treatment,data=tbc3dayCC)->km3dayCC
survfit(Surv(days_from_rand,abst_7daynoNA)~Treatment,data=tbc7daynoNA)->km7daynoNA
survfit(Surv(days_from_rand,abst_3daynoNA)~Treatment,data=tbc3daynoNA)->km3daynoNA

#plots
ggsurvplot(km7daynoNA, 
           palette = viridis(4)[1:3], 
           legend.labs=c("Immediate","Gradual","Control"),
           ggtheme = theme_classic2(base_size=12, base_family = "Verdana"),
           font.family="Verdana",
           linetype = "strata",
           font.title=(c(16,"bold","black")),
           font.x = c(16, "bold", "black"),
           font.y = c(16, "bold", "black"),
           font.tickslab = c(13, "plain", "black"),
           ylab="7-day abstinence probability",
           xlab="Weeks from Randomization",
           risk.table = TRUE,
           break.time.by=28,
           xscale=7,
           fun='event',
           ylim=c(0,0.3),
           xlim=c(0,140),
           risk.table.height = 0.2
)->imputed7dayplot
ggsurvplot(km3daynoNA, 
           palette = viridis(4)[1:3], 
           legend.labs=c("Immediate","Gradual","Control"),
           ggtheme = theme_classic2(base_size=12, base_family = "Verdana"),
           font.family="Verdana",
           linetype = "strata",
           font.title=(c(16,"bold","black")),
           font.x = c(16, "bold", "black"),
           font.y = c(16, "bold", "black"),
           font.tickslab = c(13, "plain", "black"),
           ylab="3-day abstinence probability",
           xlab="Weeks from Randomization",
           risk.table = TRUE,
           break.time.by=28,
           xscale=7,
           fun='event',
           ylim=c(0,0.3),
           xlim=c(0,140),
           risk.table.height = 0.2
)->imputed3dayplot
ggsurvplot(km3dayCC, 
           palette = viridis(4)[1:3], 
           legend.labs=c("Immediate","Gradual","Control"),
           ggtheme = theme_classic2(base_size=12, base_family = "Verdana"),
           font.family="Verdana",
           linetype = "strata",
           font.title=(c(16,"bold","black")),
           font.x = c(16, "bold", "black"),
           font.y = c(16, "bold", "black"),
           font.tickslab = c(13, "plain", "black"),
           ylab="3-day abstinence probability",
           xlab="Weeks from Randomization",
           risk.table = TRUE,
           break.time.by=28,
           xscale=7,
           fun='event',
           ylim=c(0,0.3),
           xlim=c(0,140),
           risk.table.height = 0.2
)->CC3dayplot
ggsurvplot(km7dayCC, 
           palette = viridis(4)[1:3], 
           legend.labs=c("Immediate","Gradual","Control"),
           ggtheme = theme_classic2(base_size=12, base_family = "Verdana"),
           font.family="Verdana",
           linetype = "strata",
           font.title=(c(16,"bold","black")),
           font.x = c(16, "bold", "black"),
           font.y = c(16, "bold", "black"),
           font.tickslab = c(13, "plain", "black"),
           ylab="7-day abstinence probability",
           xlab="Weeks from Randomization",
           risk.table = TRUE,
           break.time.by=28,
           xscale=7,
           fun='event',
           ylim=c(0,0.3),
           xlim=c(0,140),
           risk.table.height = 0.2
)->CC7dayplot 

plots[[1]]<-imputed7dayplot
plots[[2]]<-CC7dayplot
plots[[3]]<-imputed3dayplot
plots[[4]]<-CC3dayplot
arrange_ggsurvplots(plots,nrow=1,ncol=2)

#results from stratified log-rank and stratified cox regression
survdiff(Surv(days_from_rand,abst_7day)~Treatment+strata(site),data=tbc7dayCC[tbc7dayCC$Treatment!='C',])
survdiff(Surv(days_from_rand,abst_3day)~Treatment+strata(site),data=tbc3dayCC[tbc3dayCC$Treatment!='C',])
survdiff(Surv(days_from_rand,abst_7daynoNA)~Treatment+strata(site),data=tbc7daynoNA[tbc7daynoNA$Treatment!='A',])
survdiff(Surv(days_from_rand,abst_3daynoNA)~Treatment+strata(site),data=tbc3daynoNA[tbc3daynoNA$Treatment!='C',])

coxph(Surv(days_from_rand,abst_7day)~Treatment+strata(site),data=tbc7dayCC)->cox7dayCC
coxph(Surv(days_from_rand,abst_3day)~Treatment+strata(site),data=tbc3dayCC)->cox3dayCC
coxph(Surv(days_from_rand,abst_7daynoNA)~Treatment+strata(site),data=tbc7daynoNA)->cox7daynoNA
coxph(Surv(days_from_rand,abst_3daynoNA)~Treatment+strata(site),data=tbc3daynoNA)->cox3daynoNA

#andersen-gills data definition
AG7dayImp <- tbcRC %>%
  group_by(studyid) %>%
  filter(abst_7daynoNA == 1 | row_number() == n()) %>%
  mutate(
    start = days_from_rand,
    stop = days_from_rand + 1,
    event = as.integer(abst_7daynoNA)
  ) %>%
  ungroup()

coxph(Surv(start, stop, abst_7daynoNA) ~ Treatment + factor(site) + cluster(studyid), data = AG7dayImp)
#repeat for 3-day and missing method 2.

#defining dataset for GLMs 
ifelse(tbcRC$TotalCPD_ivr!=0|is.na(tbcRC$TotalCPD_ivr),0,1)->tbcRC$cfd

nb7dayImp <- tbcRC %>%
  group_by(studyid) %>%
  mutate(total_abstinent_days = sum(cfd, na.rm = TRUE)) %>%
  mutate(any1day = ifelse(total_abstinent_days>0,1,0)) %>% 
  slice_max(days_from_rand,with_ties = FALSE) %>% 
  ungroup()

glm(abst_7daynoNA~Treatment+site,data=nb7dayImp,family = binomial(link='logit'))
glm(any1day~Treatment+site,data=nb7dayImp,family = binomial(link='logit'))
glm.nb(total_abstinent_days ~ Treatment + factor(site), data = nb7dayImp)


#power calculations
nb7dayImpA<-nb7dayImp[nb7dayImp$Treatment=='A',]
nb7dayImpB<-nb7dayImp[nb7dayImp$Treatment=='B',]
nb7dayImpC<-nb7dayImp[nb7dayImp$Treatment=='C',]
unvpwra<-data.frame()
for (j in c(0:10*50+150,675)) {
  pvals <- numeric(1000)
  
  for (i in 1:1000) {
    indicesB <- sample(1:nrow(nb7dayImpB), floor(2*j/3), replace = TRUE)
    indicesC <- sample(1:nrow(nb7dayImpC), floor(j/3), replace = TRUE)
    data<-rbind(nb7dayImpB[indicesB,],nb7dayImpC[indicesC,])
    glm(abst_7daynoNA~Treatment+site,data=data,family = binomial(link='logit'))->unvglm
    unvsum <- summary(unvglm)
    pvals[i] <- coef(unvsum)[2, 4]
  }
mean(pvals<0.05)->unvpwra[nrow(unvpwra)+1,1]
}

for (j in c(0:10*50+150,675)) {
  pvals <- numeric(1000)
  
  for (i in 1:1000) {
    indicesB <- sample(1:nrow(nb7dayImpB), floor(2*j/3), replace = TRUE)
    indicesC <- sample(1:nrow(nb7dayImpC), floor(j/3), replace = TRUE)
    data<-rbind(nb7dayImpB[indicesB,],nb7dayImpC[indicesC,])
    glm(abst_7daynoNA~Treatment+site,data=data,family = binomial(link='logit'))->unvglm
    unvsum <- summary(unvglm)
    pvals[i] <- coef(unvsum)[2, 4]
  }
  mean(pvals<0.05)->unvpwr[nrow(unvpwr)+1,1]
}
tbcdataA<-tbc7daynoNA[tbc7daynoNA$Treatment=='A',]
tbcdataB<-tbc7daynoNA[tbc7daynoNA$Treatment=='B',]
tbcdataC<-tbc7daynoNA[tbc7daynoNA$Treatment=='C',]
coxp<-data_frame()
coxpwr<-data.frame()
for (j in floor(c(0:10*50+150,675))) {
  for (i in 1:1000) {
    sample(1:nrow(tbcdataB),floor(3*j/3),replace=TRUE)->indicesB
    sample(1:nrow(tbcdataC),floor(j/2),replace=TRUE)->indicesC
    data<-rbind(tbcdataB[indicesB,],tbcdataC[indicesC,])
    summary(coxph(Surv(days_from_rand,abst_7daynoNA)~Treatment+strata(site),data=data))->coxsimsum
    coxsimsum$logtest["pvalue"]->coxp[i,1]
  }
  
  mean(coxp<0.05)->coxpwr[ceiling(j/50)-2,1]
}
AGA_split<-AG7dayImp[AG7dayImp$Treatment=='A',]
AGB_split<-AG7dayImp[AG7dayImp$Treatment=='B',]
AGC_split<-AG7dayImp[AG7dayImp$Treatment=='C',]
AGpwrBC<-data.frame()
for (j in floor(c(0:10*50+150,675)/.75)) {
  pvals <- numeric(1000)
  
  for (i in seq_len(1000)) {
    studyidB <- sample(names(AGB_split), floor(j / 2), replace = TRUE)
    studyidC <- sample(names(AGC_split), floor(j / 2), replace = TRUE)
    
    sampledB <- relabel_subjects(AGB_split, studyidB, "B", i)
    sampledC <- relabel_subjects(AGC_split, studyidC, "C", i)
    
    data <- rbind(sampledB, sampledC)
    
    fit <- coxph(Surv(start, stop, abst_3daynoNA) ~ Treatment + strata(site) + cluster(studyid), data = data)
    pvals[i] <- summary(fit)$logtest["pvalue"]
  }
  
  AGpwrBC[ceiling(j / 50)-2,1] <- mean(pvals < 0.05)
}


cbind(glmpwr,coxpwr,unvpwr,nbpwr,AGpwrBC)->curvedataBC
curvedataBC$SampleSize<-floor(c(0:10*50+150,675)/.75)
colnames(curvedataBC)<-c('glmpwr','coxpwr','kmpwr','any1glmpwr','nbpwr','AGpwrBC','SampleSize')

longdata <- curvedataBC %>%
  pivot_longer(
    cols = -SampleSize,
    names_to = "MethodCode",
    values_to = "Power"
  )

method_labels <- c(
  glmpwr = "Logistic Regression, Abstinence at Week 20",
  nbpwr = "Negative Binomial, Number of Cigarette Free Days",
  any1glmpwr = "Logistic Regression, Any Cigarette Free Days",
  coxpwr = "Stratified Cox Regression",
  AGpwrAB = "Andersen–Gill Model, Recurrent Time-to-Event"
)
longdata$Method <- method_labels[longdata$MethodCode]
ggplot(longdata, aes(x = SampleSize, y = Power, color = Method)) +
  geom_point(size = 4) +
  geom_line(size = 1) +
  xlim(150,675)+
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2),
                     labels = label_number(accuracy = 0.01))+ 
  labs(
    x = "Combined Sample Size",
    y = "Power"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_color_viridis_d(option = "D")->powerplotBC

#repeat for AC and AB. For AB, we sample from the two groups at a 1:1 ratio, meaning the combined sample size is larger than that of BC and AC. 

legend<-get_legend(powerplotAC)

p1_clean <- powerplotBC + theme(legend.position = "none")
p2_clean <- powerplotAC + theme(legend.position = "none")
p3_clean <- powerplotAB + theme(legend.position = "none")

top_row <- plot_grid(p1_clean, p2_clean, ncol = 2, align = "hv")
bottom_row <- plot_grid(p3_clean, legend, ncol = 2,rel_widths = c(1, 1))

final_plot <- plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(1, 1))
