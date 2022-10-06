library(tidyverse)
library(mirt)
library(haven)
library(rstudioapi)
library(psych)
library(MplusAutomation)
library(BifactorIndicesCalculator)
library(lordif)
library(ggthemr)
library(Hmisc)

setwd(dirname(getActiveDocumentContext()$path))  

data <- read_sav("reduction_new.sav")
ids <- read_csv("csc-ids.csv")
ids2 <- ids %>%
  mutate(ITT=1) %>%
  select(-ContactID)

join<- full_join(data, ids2, by="contactid")

join2<- join %>%
  mutate(study=case_when(contactid>1 ~ 1,
                         contactid<1 ~ 0)) %>%
  mutate(itt2=case_when(study==0 & ITT==1 ~ 1,
                        study==1 ~ 1),
         k6cat=case_when(m0k6>=13 ~ 1,
                                m0k6>=0 & m0k6<=12 ~ 0))


data2<- join2 %>%
  filter(itt2==1)

table(data2$k6cat)

##splitting dataset into calibrationa and test data 80/20
sample_size <- floor(0.8*nrow(data2))
set.seed(12345)

flagged <- sample(seq_len(nrow(data2)), size=sample_size)
calibration <- data2[flagged,]
validation <- data2[-flagged,]



###cronbach's alpha and mcdonalds omega for scale totals
chu9<-calibration%>%
  select(m0chu1, m0chu2,m0chu3,m0chu4,m0chu5,m0chu6,m0chu7,m0chu8,m0chu9)
omega(chu9)
gad7<-calibration%>%
  select(m0gad71,m0gad72,m0gad73,m0gad74,m0gad75,m0gad76,m0gad77)
omega(gad7)
bsi<-calibration%>%
  select(m0h1suic,m0h2lonely,m0h3sad,m0h4anhedonia,m0h5hopeless,m0h6worthless)
omega(bsi)
k6<-calibration%>%
  select(m0k61,m0k62,m0k63,m0k64,m0k65,m0k66)
omega(k6)
phq8<-calibration%>%
  select(m0phq1,m0phq2,m0phq3,m0phq4,m0phq5,m0phq6,m0phq7,m0phq8)
omega(phq8)
sdq<-calibration%>%
  select(m0sd13,m0sd16,m0sd24,m0sd3,m0sd8)
omega(sdq, nfactors=2)
alpha(sdq)

comb1<-calibration%>%
  select(m0chu1, m0chu2,m0chu3,m0chu4,m0chu5,m0chu6,m0chu7,m0chu8,m0chu9,
         m0sd13,m0sd16,m0sd24,m0sd3,m0sd8,m0k61,m0k62,m0k63,m0k64,m0k65,m0k66)
omega(comb1)
comb2<-calibration%>%
  select(m0gad71,m0gad72,m0gad73,m0gad74,m0gad75,m0gad76,m0gad77,
         m0sd13,m0sd16,m0sd24,m0sd3,m0sd8,m0k61,m0k62,m0k63,m0k64,m0k65,m0k66)
omega(comb2)
comb3<-calibration%>%
  select(m0h1suic,m0h2lonely,m0h3sad,m0h4anhedonia,m0h5hopeless,m0h6worthless,
         m0sd13,m0sd16,m0sd24,m0sd3,m0sd8,m0k61,m0k62,m0k63,m0k64,m0k65,m0k66)
omega(comb3)
comb4<-calibration%>%
  select(m0phq1,m0phq2,m0phq3,m0phq4,m0phq5,m0phq6,m0phq7,m0phq8,
         m0sd13,m0sd16,m0sd24,m0sd3,m0sd8,m0k61,m0k62,m0k63,m0k64,m0k65,m0k66)
omega(comb4)
comb5<-calibration%>%
  select(m0minisp1,m0minisp2,m0minisp3,
         m0sd13,m0sd16,m0sd24,m0sd3,m0sd8,m0k61,m0k62,m0k63,m0k64,m0k65,m0k66)
omega(comb5)


##setting up dataset that just contains items for IRT model. Reverse scoring the K6 items so higher scores = higher distress

data0<-calibration %>%
  select(m0chu1,m0chu2,m0chu3,m0chu4,m0chu5,m0chu6,m0chu7,m0chu8,m0chu9,m0gad71,m0gad72,m0gad73,m0gad74,m0gad75,m0gad76,m0gad77,
         m0h1suic,m0h2lonely,m0h3sad,m0h4anhedonia,m0h5hopeless,m0h6worthless,
         m0k61,m0k62,m0k63,m0k64,m0k65,m0k66,
         m0phq1,m0phq2,m0phq3,m0phq4,m0phq5,m0phq6,m0phq7,m0phq8,
         m0sd13,m0sd16,m0sd24,m0sd3,m0sd8) %>%
  mutate(m0k61=4-m0k61,
         m0k62=4-m0k62,
         m0k63=4-m0k63,
         m0k64=4-m0k64,
         m0k65=4-m0k65,
         m0k66=4-m0k66)

datamiss<-data0 %>%
  mutate(miss= rowSums(!is.na(data0)))%>%
  filter(miss>0)

mean(datamiss$miss)
sd(datamiss$miss)
summary(datamiss)
summary(data0)

# Load function
source("http://pcwww.liv.ac.uk/~william/R/crosstab.r")

crosstab(calibration, row.vars="study", col.vars="m0chu1", type="f")
crosstab(calibration, row.vars="study", col.vars="m0chu2", type="f")
crosstab(calibration, row.vars="study", col.vars="m0chu3", type="f")
crosstab(calibration, row.vars="study", col.vars="m0chu4", type="f")
crosstab(calibration, row.vars="study", col.vars="m0chu5", type="f")
crosstab(calibration, row.vars="study", col.vars="m0chu6", type="f")
crosstab(calibration, row.vars="study", col.vars="m0chu7", type="f")
crosstab(calibration, row.vars="study", col.vars="m0chu8", type="f")
crosstab(calibration, row.vars="study", col.vars="m0chu9", type="f")
crosstab(calibration, row.vars="study", col.vars="m0gad71", type="f")
crosstab(calibration, row.vars="study", col.vars="m0gad72", type="f")
crosstab(calibration, row.vars="study", col.vars="m0gad73", type="f")
crosstab(calibration, row.vars="study", col.vars="m0gad74", type="f")
crosstab(calibration, row.vars="study", col.vars="m0gad75", type="f")
crosstab(calibration, row.vars="study", col.vars="m0gad76", type="f")
crosstab(calibration, row.vars="study", col.vars="m0gad77", type="f")
crosstab(calibration, row.vars="study", col.vars="m0h1suic", type="f")
crosstab(calibration, row.vars="study", col.vars="m0h2lonely", type="f")
crosstab(calibration, row.vars="study", col.vars="m0h3sad", type="f")
crosstab(calibration, row.vars="study", col.vars="m0h4anhedonia", type="f")
crosstab(calibration, row.vars="study", col.vars="m0h5hopeless", type="f")
crosstab(calibration, row.vars="study", col.vars="m0h6worthless", type="f")
crosstab(calibration, row.vars="study", col.vars="m0k61", type="f")
crosstab(calibration, row.vars="study", col.vars="m0k62", type="f")
crosstab(calibration, row.vars="study", col.vars="m0k63", type="f")
crosstab(calibration, row.vars="study", col.vars="m0k64", type="f")
crosstab(calibration, row.vars="study", col.vars="m0k65", type="f")
crosstab(calibration, row.vars="study", col.vars="m0k66", type="f")
crosstab(calibration, row.vars="study", col.vars="m0phq1", type="f")
crosstab(calibration, row.vars="study", col.vars="m0phq2", type="f")
crosstab(calibration, row.vars="study", col.vars="m0phq3", type="f")
crosstab(calibration, row.vars="study", col.vars="m0phq4", type="f")
crosstab(calibration, row.vars="study", col.vars="m0phq5", type="f")
crosstab(calibration, row.vars="study", col.vars="m0phq6", type="f")
crosstab(calibration, row.vars="study", col.vars="m0phq7", type="f")
crosstab(calibration, row.vars="study", col.vars="m0phq8", type="f")
crosstab(calibration, row.vars="study", col.vars="m0sd13", type="f")
crosstab(calibration, row.vars="study", col.vars="m0sd16", type="f")
crosstab(calibration, row.vars="study", col.vars="m0sd24", type="f")
crosstab(calibration, row.vars="study", col.vars="m0sd3", type="f")
crosstab(calibration, row.vars="study", col.vars="m0sd8", type="f")


##descriptive statistics and correlation between items
datamat<-as.matrix(data0)
describe(datamat)


cors<-cor.ci(data0, poly=TRUE)
corPlot(cors, xlas=2, main="Correlation Plot", numbers=FALSE)


## fitting graded response IRT model to all items first round
data0uni<-mirt(data0,1, itemtype="graded", SE=FALSE)
print(data0uni)
summary(data0uni)
coef(data0uni,printSE=TRUE,IRTpars=TRUE)
coef(data0uni, simplify=TRUE,IRTpars=TRUE)
distresseap<-fscores(data0uni, method="EAP", full.scores.SE=TRUE)

totinfo<-plot(data0uni, type="info")
totinfo
infotot<-as.data.frame(cbind(totinfo$panel.args[[1]][[1]],totinfo$panel.args[[1]][[2]]))

ggthemr('flat')
ggplot(infotot, aes(x=V1, y=V2))+
  geom_line(size=1.0) +
  xlab("Theta score")+
  ylab("Information")+
  xlim(-6,6)+
  scale_y_continuous(breaks = c(0,10,20,30,40,50,60))+
  scale_x_continuous(breaks=c(-6,-3,0,3,6),
                     sec.axis=sec_axis(~ . *10+50, name="T-score", breaks=c(20,50,80)))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))+
  scale_colour_ggthemr_d()



plot(data0uni, type="trace")
plot(data0uni, type="infotrace")
plot(data0uni, type="score")
totse<-plot(data0uni, type="SE")
setot<-as.data.frame(cbind(totse$panel.args[[1]][[1]],totse$panel.args[[1]][[2]]))
ggthemr('flat')
ggplot(setot, aes(x=V1, y=V2))+
  geom_line(size=0.7) +
  xlab("Theta score")+
  ylab("SE")+
  xlim(-6,6)+
  geom_hline(yintercept=0.44, linetype="dashed",color="green", size=0.7)+
  geom_hline(yintercept=0.54, linetype="dashed",color="orange", size=0.7)+
  geom_hline(yintercept=0.64, linetype="dashed",color="red", size=0.7)+
  #scale_y_continuous(breaks = c(0,10,20,30,40,50,60))+
  scale_x_continuous(breaks=c(-6,-3,0,3,6),
                     sec.axis=sec_axis(~ . *10+50, name="T-score", breaks=c(20,50,80)))+
  scale_colour_ggthemr_d()

####Item probability map for CHU9D
Theta <-as.numeric(matrix(seq(-2,4, by=.01)))
extr.1 <- extract.item(data0uni, 1)
extr.2 <- extract.item(data0uni, 2)
extr.3 <- extract.item(data0uni, 3)
extr.4 <- extract.item(data0uni, 4)
extr.5 <- extract.item(data0uni, 5)
extr.6 <- extract.item(data0uni, 6)
extr.7 <- extract.item(data0uni, 7)
extr.8 <- extract.item(data0uni, 8)
extr.9 <- extract.item(data0uni, 9)
item1<-probtrace(extr.1, Theta)
item2<-probtrace(extr.2, Theta)
item3<-probtrace(extr.3, Theta)
item4<-probtrace(extr.4, Theta)
item5<-probtrace(extr.5, Theta)
item6<-probtrace(extr.6, Theta)
item7<-probtrace(extr.7, Theta)
item8<-probtrace(extr.8, Theta)
item9<-probtrace(extr.9, Theta)
tracelines1<-as.data.frame(cbind(Theta, colnames(item1)[max.col(item1, ties.method = "first")]))
tracelines2<-as.data.frame(cbind(Theta, colnames(item2)[max.col(item2, ties.method = "first")]))
tracelines3<-as.data.frame(cbind(Theta, colnames(item3)[max.col(item3, ties.method = "first")]))
tracelines4<-as.data.frame(cbind(Theta, colnames(item4)[max.col(item4, ties.method = "first")]))
tracelines5<-as.data.frame(cbind(Theta, colnames(item5)[max.col(item5, ties.method = "first")]))
tracelines6<-as.data.frame(cbind(Theta, colnames(item6)[max.col(item6, ties.method = "first")]))
tracelines7<-as.data.frame(cbind(Theta, colnames(item7)[max.col(item7, ties.method = "first")]))
tracelines8<-as.data.frame(cbind(Theta, colnames(item8)[max.col(item8, ties.method = "first")]))
tracelines9<-as.data.frame(cbind(Theta, colnames(item9)[max.col(item9, ties.method = "first")]))

combtrace<-as.data.frame(cbind(tracelines1$Theta,tracelines1$V2,tracelines2$V2,tracelines3$V2,tracelines4$V2,
                               tracelines5$V2,tracelines6$V2,tracelines7$V2,tracelines8$V2,tracelines9$V2))
combtrace$V1<-as.numeric(combtrace$V1)
combtrace <- combtrace %>%
  mutate(Theta=V1,
         Item1=V2,
         Item2=V3,
         Item3=V4,
         Item4=V5,
         Item5=V6,
         Item6=V7,
         Item7=V8,
         Item8=V9,
         Item9=V10) %>%
  select(-V1:-V10)

combtrace<-combtrace %>%
  pivot_longer(cols=starts_with("I"), names_to = "Item", values_to = "Category")

Mattplot_df<-group_by(combtrace, Item, Category) %>% mutate(start=min(Theta)) %>% mutate(end=max(Theta))
Mattplot_df <-Mattplot_df%>% distinct(start,end, .keep_ll=TRUE)

Mattplot_df$Category<-factor(Mattplot_df$Category,
                             levels=c("P.1","P.2","P.3","P.4","P.5"),
                             label=c("Response option 1","Response option 2","Response option 3","Response option 4","Response option 5"))
Mattplot_df$Item<-as_factor(Mattplot_df$Item)

ggthemr('flat')
ggplot(Mattplot_df) +
  geom_segment(aes(x=start, xend=end, y=Item, yend=Item, colour=Category), size=10) +
  labs(x="Theta", y="", title="CHU-9D Item probability map") +
  #scale_color_manual("Category",values=c("#34495e", "#3498db", "#1abc9c","#f39c12", "#d35400"))+
  scale_x_continuous(breaks=c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4),
                     sec.axis=sec_axis(~ . *10+50, name="T-score", breaks=c(10,20,30,40,50,60,70,80,90)))+
  scale_y_discrete(limits=rev, labels=c("Item1"="chu1","Item2"="chu2","Item3"="chu3","Item4"="chu4",
                                        "Item5"="chu5","Item6"="chu6", "Item7"="chu7", "Item8"="chu8",
                                        "Item9"="chu9"))+
  geom_vline(xintercept=0, linetype="dashed", color="tomato2", size=1.0)+
  theme(legend.position="bottom", legend.title=element_blank(),legend.text = element_text(size=8),
        panel.border = element_rect(colour = "grey", fill=NA, size=0.5))+
  scale_colour_ggthemr_d()


### item probability map for GAD7##############################
Theta <-as.numeric(matrix(seq(-1,3, by=.01)))
extr.10 <- extract.item(data0uni, 10)
extr.11 <- extract.item(data0uni, 11)
extr.12 <- extract.item(data0uni, 12)
extr.13 <- extract.item(data0uni, 13)
extr.14 <- extract.item(data0uni, 14)
extr.15 <- extract.item(data0uni, 15)
extr.16 <- extract.item(data0uni, 16)
item10<-probtrace(extr.10, Theta)
item11<-probtrace(extr.11, Theta)
item12<-probtrace(extr.12, Theta)
item13<-probtrace(extr.13, Theta)
item14<-probtrace(extr.14, Theta)
item15<-probtrace(extr.15, Theta)
item16<-probtrace(extr.16, Theta)
tracelines10<-as.data.frame(cbind(Theta, colnames(item10)[max.col(item10, ties.method = "first")]))
tracelines11<-as.data.frame(cbind(Theta, colnames(item11)[max.col(item11, ties.method = "first")]))
tracelines12<-as.data.frame(cbind(Theta, colnames(item12)[max.col(item12, ties.method = "first")]))
tracelines13<-as.data.frame(cbind(Theta, colnames(item13)[max.col(item13, ties.method = "first")]))
tracelines14<-as.data.frame(cbind(Theta, colnames(item14)[max.col(item14, ties.method = "first")]))
tracelines15<-as.data.frame(cbind(Theta, colnames(item15)[max.col(item15, ties.method = "first")]))
tracelines16<-as.data.frame(cbind(Theta, colnames(item16)[max.col(item16, ties.method = "first")]))

combtracegad<-as.data.frame(cbind(tracelines10$Theta,tracelines10$V2,tracelines11$V2,tracelines12$V2,tracelines13$V2,
                                  tracelines14$V2,tracelines15$V2,tracelines16$V2))
combtracegad$V1<-as.numeric(combtracegad$V1)
combtracegad <- combtracegad %>%
  mutate(Theta=V1,
         Item1=V2,
         Item2=V3,
         Item3=V4,
         Item4=V5,
         Item5=V6,
         Item6=V7,
         Item7=V8) %>%
  select(-V1:-V8)

combtracegad<-combtracegad %>%
  pivot_longer(cols=starts_with("I"), names_to = "Item", values_to = "Category")

Mattplot_dfgad<-group_by(combtracegad, Item, Category) %>% mutate(start=min(Theta)) %>% mutate(end=max(Theta))
Mattplot_dfgad <-Mattplot_dfgad%>% distinct(start,end, .keep_ll=TRUE)

Mattplot_dfgad$Category<-factor(Mattplot_dfgad$Category,
                                levels=c("P.1","P.2","P.3","P.4"),
                                label=c("Not at all","Several days","More than half the days","Nearly every day"))
Mattplot_dfgad$Item<-as_factor(Mattplot_dfgad$Item)

ggthemr('flat')
ggplot(Mattplot_dfgad) +
  geom_segment(aes(x=start, xend=end, y=Item, yend=Item, colour=Category), size=10) +
  labs(x="Theta", y="", title="GAD7 Item probability map") +
  #scale_color_manual("Category",values=c("#34495e", "#3498db", "#1abc9c","#f39c12", "#d35400"))+
  scale_x_continuous(breaks=c(-1,-0.5,0,0.5,1,1.5,2,2.5,3),
                     sec.axis=sec_axis(~ . *10+50, name="T-score", breaks=c(10,20,30,40,50,60,70,80,90)))+
  scale_y_discrete(limits=rev)+
  geom_vline(xintercept=0, linetype="dashed", color="tomato2", size=1.0)+
  theme(legend.position="bottom", legend.title=element_blank(),legend.text = element_text(size=8),
        panel.border = element_rect(colour = "grey", fill=NA, size=0.5))+
  scale_colour_ggthemr_d()


### item probability map for BSI-depression##############################
Theta <-as.numeric(matrix(seq(-1,3, by=.01)))
extr.17 <- extract.item(data0uni, 17)
extr.18 <- extract.item(data0uni, 18)
extr.19 <- extract.item(data0uni, 19)
extr.20 <- extract.item(data0uni, 20)
extr.21 <- extract.item(data0uni, 21)
extr.22 <- extract.item(data0uni, 22)
item17<-probtrace(extr.17, Theta)
item18<-probtrace(extr.18, Theta)
item19<-probtrace(extr.19, Theta)
item20<-probtrace(extr.20, Theta)
item21<-probtrace(extr.21, Theta)
item22<-probtrace(extr.22, Theta)
tracelines17<-as.data.frame(cbind(Theta, colnames(item17)[max.col(item17, ties.method = "first")]))
tracelines18<-as.data.frame(cbind(Theta, colnames(item18)[max.col(item18, ties.method = "first")]))
tracelines19<-as.data.frame(cbind(Theta, colnames(item19)[max.col(item19, ties.method = "first")]))
tracelines20<-as.data.frame(cbind(Theta, colnames(item20)[max.col(item20, ties.method = "first")]))
tracelines21<-as.data.frame(cbind(Theta, colnames(item21)[max.col(item21, ties.method = "first")]))
tracelines22<-as.data.frame(cbind(Theta, colnames(item22)[max.col(item22, ties.method = "first")]))


combtracebsi<-as.data.frame(cbind(tracelines17$Theta,tracelines17$V2,tracelines18$V2,tracelines19$V2,tracelines20$V2,
                                  tracelines21$V2,tracelines22$V2))
combtracebsi$V1<-as.numeric(combtracebsi$V1)
combtracebsi <- combtracebsi %>%
  mutate(Theta=V1,
         Item1=V2,
         Item2=V3,
         Item3=V4,
         Item4=V5,
         Item5=V6,
         Item6=V7) %>%
  select(-V1:-V7)

combtracebsi<-combtracebsi %>%
  pivot_longer(cols=starts_with("I"), names_to = "Item", values_to = "Category")

Mattplot_dfbsi<-group_by(combtracebsi, Item, Category) %>% mutate(start=min(Theta)) %>% mutate(end=max(Theta))
Mattplot_dfbsi <-Mattplot_dfbsi%>% distinct(start,end, .keep_ll=TRUE)

Mattplot_dfbsi$Category<-factor(Mattplot_dfbsi$Category,
                                levels=c("P.1","P.2","P.3","P.4","P.5"),
                                label=c("Not at all","A little bit","Moderately","Quite a bit", "Often"))
Mattplot_dfbsi$Item<-as_factor(Mattplot_dfbsi$Item)

ggthemr('flat')
ggplot(Mattplot_dfbsi) +
  geom_segment(aes(x=start, xend=end, y=Item, yend=Item, colour=Category), size=10) +
  labs(x="Theta", y="", title="BSI Depression Item probability map") +
  #scale_color_manual("Category",values=c("#34495e", "#3498db", "#1abc9c","#f39c12", "#d35400"))+
  scale_x_continuous(breaks=c(-1,-0.5,0,0.5,1,1.5,2,2.5,3),
                     sec.axis=sec_axis(~ . *10+50, name="T-score", breaks=c(10,20,30,40,50,60,70,80,90)))+
  scale_y_discrete(limits=rev, labels=c("Item1"="h1suic","Item2"="h2lonely","Item3"="h3sad","Item4"="h4anhedonia",
                                        "Item5"="h5hopeless","Item6"="h6worthless"))+
  geom_vline(xintercept=0, linetype="dashed", color="tomato2", size=1.0)+
  theme(legend.position="bottom", legend.title=element_blank(),legend.text = element_text(size=8),
        panel.border = element_rect(colour = "grey", fill=NA, size=0.5))+
  scale_colour_ggthemr_d()



### item probability map for K6##############################
Theta <-as.numeric(matrix(seq(-2,3, by=.01)))
extr.23 <- extract.item(data0uni, 23)
extr.24 <- extract.item(data0uni, 24)
extr.25 <- extract.item(data0uni, 25)
extr.26 <- extract.item(data0uni, 26)
extr.27 <- extract.item(data0uni, 27)
extr.28 <- extract.item(data0uni, 28)
item23<-probtrace(extr.23, Theta)
item24<-probtrace(extr.24, Theta)
item25<-probtrace(extr.25, Theta)
item26<-probtrace(extr.26, Theta)
item27<-probtrace(extr.27, Theta)
item28<-probtrace(extr.28, Theta)
tracelines23<-as.data.frame(cbind(Theta, colnames(item23)[max.col(item23, ties.method = "first")]))
tracelines24<-as.data.frame(cbind(Theta, colnames(item24)[max.col(item24, ties.method = "first")]))
tracelines25<-as.data.frame(cbind(Theta, colnames(item25)[max.col(item25, ties.method = "first")]))
tracelines26<-as.data.frame(cbind(Theta, colnames(item26)[max.col(item26, ties.method = "first")]))
tracelines27<-as.data.frame(cbind(Theta, colnames(item27)[max.col(item27, ties.method = "first")]))
tracelines28<-as.data.frame(cbind(Theta, colnames(item28)[max.col(item28, ties.method = "first")]))


combtracek6<-as.data.frame(cbind(tracelines23$Theta,tracelines23$V2,tracelines24$V2,tracelines25$V2,tracelines26$V2,
                                 tracelines27$V2,tracelines28$V2))
combtracek6$V1<-as.numeric(combtracek6$V1)
combtracek6 <- combtracek6 %>%
  mutate(Theta=V1,
         Item1=V2,
         Item2=V3,
         Item3=V4,
         Item4=V5,
         Item5=V6,
         Item6=V7) %>%
  select(-V1:-V7)

combtracek6<-combtracek6 %>%
  pivot_longer(cols=starts_with("I"), names_to = "Item", values_to = "Category")

Mattplot_dfk6<-group_by(combtracek6, Item, Category) %>% mutate(start=min(Theta)) %>% mutate(end=max(Theta))
Mattplot_dfk6 <-Mattplot_dfk6%>% distinct(start,end, .keep_ll=TRUE)

Mattplot_dfk6$Category<-factor(Mattplot_dfk6$Category,
                               levels=c("P.1","P.2","P.3","P.4","P.5"),
                               label=c("None of the time","A little of the time","Some of the time","Most of the time", "All of the time"))
Mattplot_dfk6$Item<-as_factor(Mattplot_dfk6$Item)

ggthemr('flat')
ggplot(Mattplot_dfk6) +
  geom_segment(aes(x=start, xend=end, y=Item, yend=Item, colour=Category), size=10) +
  labs(x="Theta", y="", title="K6 Item probability map") +
  #scale_color_manual("Category",values=c("#34495e", "#3498db", "#1abc9c","#f39c12", "#d35400"))+
  scale_x_continuous(breaks=c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3),
                     sec.axis=sec_axis(~ . *10+50, name="T-score", breaks=c(10,20,30,40,50,60,70,80,90)))+
  scale_y_discrete(limits=rev, labels=c("Item1"="K61","Item2"="K62","Item3"="K63","Item4"="K64",
                                        "Item5"="K65","Item6"="K66"))+
  geom_vline(xintercept=0, linetype="dashed", color="tomato2", size=1.0)+
  theme(legend.position="bottom", legend.title=element_blank(),legend.text = element_text(size=8),
        panel.border = element_rect(colour = "grey", fill=NA, size=0.5))+
  scale_colour_ggthemr_d()


### item probability map for phq##############################
Theta <-as.numeric(matrix(seq(-1,3, by=.01)))
extr.29 <- extract.item(data0uni, 29)
extr.30 <- extract.item(data0uni, 30)
extr.31 <- extract.item(data0uni, 31)
extr.32 <- extract.item(data0uni, 32)
extr.33 <- extract.item(data0uni, 33)
extr.34 <- extract.item(data0uni, 34)
extr.35 <- extract.item(data0uni, 35)
extr.36 <- extract.item(data0uni, 36)
item29<-probtrace(extr.29, Theta)
item30<-probtrace(extr.30, Theta)
item31<-probtrace(extr.31, Theta)
item32<-probtrace(extr.32, Theta)
item33<-probtrace(extr.33, Theta)
item34<-probtrace(extr.34, Theta)
item35<-probtrace(extr.35, Theta)
item36<-probtrace(extr.36, Theta)
tracelines29<-as.data.frame(cbind(Theta, colnames(item29)[max.col(item29, ties.method = "first")]))
tracelines30<-as.data.frame(cbind(Theta, colnames(item30)[max.col(item30, ties.method = "first")]))
tracelines31<-as.data.frame(cbind(Theta, colnames(item31)[max.col(item31, ties.method = "first")]))
tracelines32<-as.data.frame(cbind(Theta, colnames(item32)[max.col(item32, ties.method = "first")]))
tracelines33<-as.data.frame(cbind(Theta, colnames(item33)[max.col(item33, ties.method = "first")]))
tracelines34<-as.data.frame(cbind(Theta, colnames(item34)[max.col(item34, ties.method = "first")]))
tracelines35<-as.data.frame(cbind(Theta, colnames(item35)[max.col(item35, ties.method = "first")]))
tracelines36<-as.data.frame(cbind(Theta, colnames(item36)[max.col(item36, ties.method = "first")]))

combtracephq<-as.data.frame(cbind(tracelines29$Theta,tracelines29$V2,tracelines30$V2,tracelines31$V2,tracelines32$V2,
                                  tracelines33$V2,tracelines34$V2,tracelines35$V2,tracelines36$V2))
combtracephq$V1<-as.numeric(combtracephq$V1)
combtracephq <- combtracephq %>%
  mutate(Theta=V1,
         Item1=V2,
         Item2=V3,
         Item3=V4,
         Item4=V5,
         Item5=V6,
         Item6=V7,
         Item7=V8,
         Item8=V9) %>%
  select(-V1:-V9)

combtracephq<-combtracephq %>%
  pivot_longer(cols=starts_with("I"), names_to = "Item", values_to = "Category")

Mattplot_dfphq<-group_by(combtracephq, Item, Category) %>% mutate(start=min(Theta)) %>% mutate(end=max(Theta))
Mattplot_dfphq <-Mattplot_dfphq%>% distinct(start,end, .keep_ll=TRUE)

Mattplot_dfphq$Category<-factor(Mattplot_dfphq$Category,
                                levels=c("P.1","P.2","P.3","P.4"),
                                label=c("Not at all","Several days","More than half the days","Nearly every day"))
Mattplot_dfphq$Item<-as_factor(Mattplot_dfphq$Item)

ggthemr('flat')
ggplot(Mattplot_dfphq) +
  geom_segment(aes(x=start, xend=end, y=Item, yend=Item, colour=Category), size=10) +
  labs(x="Theta", y="", title="PHQ Item probability map") +
  #scale_color_manual("Category",values=c("#34495e", "#3498db", "#1abc9c","#f39c12", "#d35400"))+
  scale_x_continuous(breaks=c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3),
                     sec.axis=sec_axis(~ . *10+50, name="T-score", breaks=c(10,20,30,40,50,60,70,80,90)))+
  scale_y_discrete(limits=rev, labels=c("Item1"="PHQ1","Item2"="PHQ2","Item3"="PHQ3","Item4"="PHQ4",
                                        "Item5"="PHQ5","Item6"="PHQ6","Item7"="PHQ7","Item8"="PHQ8"))+
  geom_vline(xintercept=0, linetype="dashed", color="tomato2", size=1.0)+
  theme(legend.position="bottom", legend.title=element_blank(),legend.text = element_text(size=8),
        panel.border = element_rect(colour = "grey", fill=NA, size=0.5))+
  scale_colour_ggthemr_d()


### item probability map for sdq##############################
Theta <-as.numeric(matrix(seq(-1,3, by=.01)))
extr.37<- extract.item(data0uni, 37)
extr.38 <- extract.item(data0uni, 38)
extr.39 <- extract.item(data0uni, 39)
extr.40 <- extract.item(data0uni, 40)
extr.41 <- extract.item(data0uni, 41)
item37<-probtrace(extr.37, Theta)
item38<-probtrace(extr.38, Theta)
item39<-probtrace(extr.39, Theta)
item40<-probtrace(extr.40, Theta)
item41<-probtrace(extr.41, Theta)
tracelines37<-as.data.frame(cbind(Theta, colnames(item37)[max.col(item37, ties.method = "first")]))
tracelines38<-as.data.frame(cbind(Theta, colnames(item38)[max.col(item38, ties.method = "first")]))
tracelines39<-as.data.frame(cbind(Theta, colnames(item39)[max.col(item39, ties.method = "first")]))
tracelines40<-as.data.frame(cbind(Theta, colnames(item40)[max.col(item40, ties.method = "first")]))
tracelines41<-as.data.frame(cbind(Theta, colnames(item41)[max.col(item41, ties.method = "first")]))

combtracesdq<-as.data.frame(cbind(tracelines37$Theta,tracelines37$V2,tracelines38$V2,tracelines39$V2,tracelines40$V2,
                                  tracelines41$V2))
combtracesdq$V1<-as.numeric(combtracesdq$V1)
combtracesdq <- combtracesdq %>%
  mutate(Theta=V1,
         Item1=V2,
         Item2=V3,
         Item3=V4,
         Item4=V5,
         Item5=V6) %>%
  select(-V1:-V6)

combtracesdq<-combtracesdq %>%
  pivot_longer(cols=starts_with("I"), names_to = "Item", values_to = "Category")

Mattplot_dfsdq<-group_by(combtracesdq, Item, Category) %>% mutate(start=min(Theta)) %>% mutate(end=max(Theta))
Mattplot_dfsdq <-Mattplot_dfsdq%>% distinct(start,end, .keep_ll=TRUE)

Mattplot_dfsdq$Category<-factor(Mattplot_dfsdq$Category,
                                levels=c("P.1","P.2","P.3"),
                                label=c("Not true","Somewhat true","Certainly true"))
Mattplot_dfsdq$Item<-as_factor(Mattplot_dfsdq$Item)

ggthemr('flat')
ggplot(Mattplot_dfsdq) +
  geom_segment(aes(x=start, xend=end, y=Item, yend=Item, colour=Category), size=10) +
  labs(x="Theta", y="", title="SDQ-Emotional Item probability map") +
  #scale_color_manual("Category",values=c("#34495e", "#3498db", "#1abc9c","#f39c12", "#d35400"))+
  scale_x_continuous(breaks=c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3),
                     sec.axis=sec_axis(~ . *10+50, name="T-score", breaks=c(10,20,30,40,50,60,70,80,90)))+
  scale_y_discrete(limits=rev, labels=c("Item1"="sd13","Item2"="sd16","Item3"="sd24","Item4"="sd3",
                                        "Item5"="sd8"))+
  geom_vline(xintercept=0, linetype="dashed", color="tomato2", size=1.0)+
  theme(legend.position="bottom", legend.title=element_blank(),legend.text = element_text(size=8),
        panel.border = element_rect(colour = "grey", fill=NA, size=0.5))+
  scale_colour_ggthemr_d()


###generating datafile for Mplus
datampls<-calibration %>%
  select(contactid,m0chu1, m0chu2,m0chu3,m0chu4,m0chu5,m0chu6,m0chu7,m0chu8,m0chu9,m0gad71,m0gad72,m0gad73,m0gad74,m0gad75,m0gad76,m0gad77,
         m0h1suic,m0h2lonely,m0h3sad,m0h4anhedonia,m0h5hopeless,m0h6worthless,m0k61,m0k62,m0k63,m0k64,m0k65,m0k66,
         m0minisp1,m0minisp2,m0minisp3,
         m0phq1,m0phq2,m0phq3,m0phq4,m0phq5,m0phq6,m0phq7,m0phq8,
         m0sd13,m0sd16,m0sd24,m0sd3,m0sd8) %>%
  mutate(m0k61=4-m0k61,
         m0k62=4-m0k62,
         m0k63=4-m0k63,
         m0k64=4-m0k64,
         m0k65=4-m0k65,
         m0k66=4-m0k66) %>%
  mutate_all(funs(replace_na(.,-9)))

write.csv(datampls, "datampls.csv", row.names=FALSE)
##### USE MPLUS TO FIT EXPLORATORY BIFACTOR MODELS AND COMPARE TO UNIDIMENSIONAL




##inspecting residual correlation in subpopulation 1
data1<-calibration %>%
  select(m0chu1, m0chu2,m0chu3,m0chu4,m0chu5,m0chu6,m0chu7,m0chu8,m0chu9,m0gad71,m0gad72,m0gad73,m0gad74,m0gad75,m0gad76,m0gad77,
         m0k61,m0k62,m0k63,m0k64,m0k65,m0k66,m0phq1,m0phq2,m0phq3,m0phq4,m0phq5,m0phq6,m0phq7,m0phq8,
         m0sd13,m0sd16,m0sd24,m0sd3,m0sd8) %>%
  mutate(m0k61=4-m0k61,
         m0k62=4-m0k62,
         m0k63=4-m0k63,
         m0k64=4-m0k64,
         m0k65=4-m0k65,
         m0k66=4-m0k66)
data1uni<-mirt(data1,1, itemtype="graded", technical=list(removeEmptyRows=TRUE))
summary(data1uni)
coef(data1uni, simplify=TRUE, printSE=TRUE)
residmat1<-M2(data1uni,residmat=TRUE, na.rm=TRUE, suppress=0.3)


##inspecting residual correlation in subpopulation 2
data22<-calibration %>%
  select(m0h1suic,m0h2lonely,m0h3sad,m0h4anhedonia,m0h5hopeless,m0h6worthless,m0k61,m0k62,m0k63,m0k64,m0k65,m0k66,
         m0sd13,m0sd16,m0sd24,m0sd3,m0sd8) %>%
  mutate(m0k61=4-m0k61,
         m0k62=4-m0k62,
         m0k63=4-m0k63,
         m0k64=4-m0k64,
         m0k65=4-m0k65,
         m0k66=4-m0k66)
data2uni<-mirt(data22,1, itemtype="graded", technical=list(removeEmptyRows=TRUE))
summary(data2uni)
coef(data2uni, simplify=TRUE)
residmat2<-M2(data2uni,residmat=TRUE, na.rm=TRUE, suppress=0.3)




##dif testing common items sample
options(scipen=999)

datadif<-calibration %>%
  mutate(m0k61=4-m0k61,
         m0k62=4-m0k62,
         m0k63=4-m0k63,
         m0k64=4-m0k64,
         m0k65=4-m0k65,
         m0k66=4-m0k66,
         sample=ifelse(contactid>=1,1,0))%>%
  select(sample,
         m0k61,m0k62,m0k63,m0k64,m0k65,m0k66,
         m0sd13,m0sd16,m0sd24,m0sd3,m0sd8) 

datadif2<-datadif%>%
  select(-sample)

datadif3<-as.matrix(datadif%>%
                      select(-sample))
groupdif<-datadif$sample
sex<-as.numeric(calibration$m0sex)

groupdif2<-factor(datadif$sample)

difftab<-lordif(datadif3, groupdif,model="GRM",criterion="R2", R2.change=0.02)
difftab


datak6<-as.matrix(calibration %>%
                    mutate(m0k61=4-m0k61,
                           m0k62=4-m0k62,
                           m0k63=4-m0k63,
                           m0k64=4-m0k64,
                           m0k65=4-m0k65,
                           m0k66=4-m0k66)%>%
                    select(
                      m0k61,m0k62,m0k63,m0k64,m0k65,m0k66))
difk6grp<-lordif(datak6,groupdif,model="GRM", criterion="R2", R2.change=0.02)
summary(difk6grp)
difk6sex<-lordif(datak6,sex,model="GRM", criterion="R2", R2.change=0.02)
summary(difk6sex)

datasd<-as.matrix(calibration %>%
                    select(m0sd13,m0sd16,m0sd24,m0sd3,m0sd8))
difsdgrp<-lordif(datasd,groupdif,model="GRM", criterion="R2", R2.change=0.02)
summary(difsdgrp)
difsdsex<-lordif(datasd,sex,model="GRM", criterion="R2", R2.change=0.02)
summary(difsdsex)




####Creating crosswalks using common metric

##setting up individual scale datasets
chudat<-calibration %>%
  select(m0chu1, m0chu2,m0chu3,m0chu4,m0chu5,m0chu6,m0chu7,m0chu8,m0chu9)
gaddat<-calibration %>%
  select(m0gad71,m0gad72,m0gad73,m0gad74,m0gad75,m0gad76,m0gad77)
bsidepdat<-calibration %>%
  select(m0h1suic,m0h2lonely,m0h3sad,m0h4anhedonia,m0h5hopeless,m0h6worthless)
k6dat<-calibration %>%
  select(m0k61,m0k62,m0k63,m0k64,m0k65,m0k66)%>%
  mutate(m0k61=4-m0k61,
         m0k62=4-m0k62,
         m0k63=4-m0k63,
         m0k64=4-m0k64,
         m0k65=4-m0k65,
         m0k66=4-m0k66)
minispdat<-calibration %>%
  select(m0minisp1,m0minisp2,m0minisp3)
phqdat<-calibration %>%
  select(m0phq1,m0phq2,m0phq3,m0phq4,m0phq5,m0phq6,m0phq7,m0phq8)
sdqdat<-calibration %>%
  select(m0sd13,m0sd16,m0sd24,m0sd3,m0sd8)


coef0uni<-coef(data0uni, as.data.frame=TRUE)

#crosswalk and figures for the CHU9D
chusv1<-mirt(chudat,1, itemtype="graded", pars="values", technical=list(removeEmptyRows=TRUE))
chusv2<-chusv1
chusv2$value<-c(coef0uni[1:45,], 0, 1)

chusv2$est<- FALSE

chuwalk<-mirt(chudat,1, itemtype="graded", pars=chusv2, technical=list(removeEmptyRows=TRUE))
coef(chuwalk)

infochufixed<-plot(chuwalk, type="info")
infochufixed
infofixedchu<-as.data.frame(cbind(infochufixed$panel.args[[1]][[1]],infochufixed$panel.args[[1]][[2]]))

sechufixed<-plot(chuwalk, type="SE")
sechufixed
sefixedchu<-as.data.frame(cbind(sechufixed$panel.args[[1]][[1]],sechufixed$panel.args[[1]][[2]]))

tccchufixed<-plot(chuwalk, type="score")
tccchufixed
scorefixedchu<-as.data.frame(cbind(tccchufixed$panel.args[[1]][[1]],tccchufixed$panel.args[[1]][[2]]))

chueap<-fscores(chuwalk, method="EAP", full.scores.SE=TRUE)
chueapsum<-fscores(chuwalk, method="EAPsum", full.scores = FALSE,na.rm=TRUE)

#cross walk and figures for the GAD7
gadsv1<-mirt(gaddat,1, itemtype="graded", pars="values", technical=list(removeEmptyRows=TRUE))
gadsv2<-gadsv1
gadsv2$value<-c(coef0uni[46:73], 0, 1)

gadsv2$est<- FALSE

gadwalk<-mirt(gaddat,1, itemtype="graded", pars=gadsv2, technical=list(removeEmptyRows=TRUE))
coef(gadwalk)

infogadfixed<-plot(gadwalk, type="info")
infogadfixed
infofixedgad<-as.data.frame(cbind(infogadfixed$panel.args[[1]][[1]],infogadfixed$panel.args[[1]][[2]]))

segadfixed<-plot(gadwalk, type="SE")
segadfixed
sefixedgad<-as.data.frame(cbind(segadfixed$panel.args[[1]][[1]],segadfixed$panel.args[[1]][[2]]))

tccgadfixed<-plot(gadwalk, type="score")
tccgadfixed
scorefixedgad<-as.data.frame(cbind(tccgadfixed$panel.args[[1]][[1]],tccgadfixed$panel.args[[1]][[2]]))

gadeap<-fscores(gadwalk, method="EAP", full.scores.SE=TRUE)
gadeapsum<-fscores(gadwalk, method="EAPsum", full.scores = FALSE,na.rm=TRUE)

#cross walk for the BSI depression scale
bsisv1<-mirt(bsidepdat,1, itemtype="graded", pars="values", technical=list(removeEmptyRows=TRUE))
bsisv2<-bsisv1
bsisv2$value<-c(coef0uni[74:103], 0, 1)

bsisv2$est<- FALSE

bsiwalk<-mirt(bsidepdat,1, itemtype="graded", pars=bsisv2, technical=list(removeEmptyRows=TRUE))
coef(bsiwalk)

infobsifixed<-plot(bsiwalk, type="info")
infobsifixed
infofixedbsi<-as.data.frame(cbind(infobsifixed$panel.args[[1]][[1]],infobsifixed$panel.args[[1]][[2]]))

sebsifixed<-plot(bsiwalk, type="SE")
sebsifixed
sefixedbsi<-as.data.frame(cbind(sebsifixed$panel.args[[1]][[1]],sebsifixed$panel.args[[1]][[2]]))

tccbsifixed<-plot(bsiwalk, type="score")
tccbsifixed
scorefixedbsi<-as.data.frame(cbind(tccbsifixed$panel.args[[1]][[1]],tccbsifixed$panel.args[[1]][[2]]))

bsieap<-fscores(bsiwalk, method="EAP", full.scores.SE=TRUE)
bsieapsum<-fscores(bsiwalk, method="EAPsum", full.scores = FALSE,na.rm=TRUE)

#cross walk for the K6 scale
k6sv1<-mirt(k6dat,1, itemtype="graded", pars="values", technical=list(removeEmptyRows=TRUE))
k6sv2<-k6sv1
k6sv2$value<-c(coef0uni[104:133], 0, 1)

k6sv2$est<- FALSE

k6walk<-mirt(k6dat,1, itemtype="graded", pars=k6sv2, technical=list(removeEmptyRows=TRUE))
coef(k6walk)

infok6fixed<-plot(k6walk, type="info")
infok6fixed
infofixedk6<-as.data.frame(cbind(infok6fixed$panel.args[[1]][[1]],infok6fixed$panel.args[[1]][[2]]))

sek6fixed<-plot(k6walk, type="SE")
sek6fixed
sefixedk6<-as.data.frame(cbind(sek6fixed$panel.args[[1]][[1]],sek6fixed$panel.args[[1]][[2]]))

tcck6fixed<-plot(k6walk, type="score")
tcck6fixed
scorefixedk6<-as.data.frame(cbind(tcck6fixed$panel.args[[1]][[1]],tcck6fixed$panel.args[[1]][[2]]))

k6eap<-fscores(k6walk, method="EAP", full.scores.SE=TRUE)
k6eapsum<-fscores(k6walk, method="EAPsum", full.scores = FALSE,na.rm=TRUE)

#cross walk for the Mini spin scale
#minispsv1<-mirt(minispdat,1, itemtype="graded", pars="values", technical=list(removeEmptyRows=TRUE))
#minispsv2<-minispsv1
#minispsv2$value<-c(coef0uni[134:148], 0, 1)

#minispsv2$est<- FALSE

#minispwalk<-mirt(minispdat,1, itemtype="graded", pars=minispsv2, technical=list(removeEmptyRows=TRUE))
#coef(minispwalk)

#infominispfixed<-plot(minispwalk, type="info")
#infominispfixed
#infofixedminisp<-as.data.frame(cbind(infominispfixed$panel.args[[1]][[1]],infominispfixed$panel.args[[1]][[2]]))

#seminispfixed<-plot(minispwalk, type="SE")
#seminispfixed
#sefixedminisp<-as.data.frame(cbind(seminispfixed$panel.args[[1]][[1]],seminispfixed$panel.args[[1]][[2]]))

#tccminispfixed<-plot(minispwalk, type="score")
#tccminispfixed
#scorefixedminisp<-as.data.frame(cbind(tccminispfixed$panel.args[[1]][[1]],tccminispfixed$panel.args[[1]][[2]]))

#minispeap<-fscores(minispwalk, method="EAP", full.scores.SE=TRUE)
#minispeapsum<-fscores(minispwalk, method="EAPsum", full.scores = FALSE,na.rm=TRUE)

#cross walk for the PHQ8 scale
phqsv1<-mirt(phqdat,1, itemtype="graded", pars="values", technical=list(removeEmptyRows=TRUE))
phqsv2<-phqsv1
phqsv2$value<-c(coef0uni[134:165], 0, 1)

phqsv2$est<- FALSE

phqwalk<-mirt(phqdat,1, itemtype="graded", pars=phqsv2, technical=list(removeEmptyRows=TRUE))
coef(phqwalk)

infophqfixed<-plot(phqwalk, type="info")
infophqfixed
infofixedphq<-as.data.frame(cbind(infophqfixed$panel.args[[1]][[1]],infophqfixed$panel.args[[1]][[2]]))

sephqfixed<-plot(phqwalk, type="SE")
sephqfixed
sefixedphq<-as.data.frame(cbind(sephqfixed$panel.args[[1]][[1]],sephqfixed$panel.args[[1]][[2]]))

tccphqfixed<-plot(phqwalk, type="score")
tccphqfixed
scorefixedphq<-as.data.frame(cbind(tccphqfixed$panel.args[[1]][[1]],tccphqfixed$panel.args[[1]][[2]]))

phqeap<-fscores(phqwalk, method="EAP", full.scores.SE=TRUE)
phqeapsum<-fscores(phqwalk, method="EAPsum", full.scores = FALSE,na.rm=TRUE)

#cross walk for the sdq emo scale
sdqsv1<-mirt(sdqdat,1, itemtype="graded", pars="values", technical=list(removeEmptyRows=TRUE))
sdqsv2<-sdqsv1
sdqsv2$value<-c(coef0uni[166:182])

sdqsv2$est<- FALSE

sdqwalk<-mirt(sdqdat,1, itemtype="graded", pars=sdqsv2, technical=list(removeEmptyRows=TRUE))
coef(sdqwalk)

infosdqfixed<-plot(sdqwalk, type="info")
infosdqfixed
infofixedsdq<-as.data.frame(cbind(infosdqfixed$panel.args[[1]][[1]],infosdqfixed$panel.args[[1]][[2]]))

sesdqfixed<-plot(sdqwalk, type="SE")
sesdqfixed
sefixedsdq<-as.data.frame(cbind(sesdqfixed$panel.args[[1]][[1]],sesdqfixed$panel.args[[1]][[2]]))

tccsdqfixed<-plot(sdqwalk, type="score")
tccsdqfixed
scorefixedsdq<-as.data.frame(cbind(tccsdqfixed$panel.args[[1]][[1]],tccsdqfixed$panel.args[[1]][[2]]))

sdqeap<-fscores(sdqwalk, method="EAP", full.scores.SE=TRUE)
sdqeapsum<-fscores(sdqwalk, method="EAPsum", full.scores = FALSE,na.rm=TRUE)

bsieapsum2<-bsieapsum %>%
  mutate(Scale="BSI") %>%
  select(c(Sum.Scores, F1, Scale))
chueapsum2<-chueapsum %>%
  mutate(Scale="CHU") %>%
  select(c(Sum.Scores, F1, Scale))
gadeapsum2<-gadeapsum %>%
  mutate(Scale="GAD") %>%
  select(c(Sum.Scores, F1, Scale))
k6eapsum2<-k6eapsum %>%
  mutate(Scale="K6") %>%
  select(c(Sum.Scores, F1, Scale))
phqeapsum2<-phqeapsum %>%
  mutate(Scale="PHQ") %>%
  select(c(Sum.Scores, F1, Scale))
sdqeapsum2<-sdqeapsum %>%
  mutate(Scale="SDQ") %>%
  select(c(Sum.Scores, F1, Scale))

walktable<-bind_rows(bsieapsum2,chueapsum2,gadeapsum2,k6eapsum2,phqeapsum2,sdqeapsum2)
walktable$Scale<-factor(walktable$Scale)
walktable$Sum.Scores<-as.numeric(walktable$Sum.Scores)

library(ggthemr)

ggthemr('flat')
ggplot(walktable, aes(Scale, F1, label = Sum.Scores, color = Scale)) +
  geom_point(shape = "-", size = 7) +
  geom_text(nudge_x=0.1, size = 3, alpha = 1, check_overlap =TRUE) +
  scale_x_discrete(position = "top", name = "Scale Sum Score") +
  guides(color = "none") +
  scale_y_continuous(name="Theta Score", limits=c(-2,4), n.breaks=15)+
  theme(panel.border = element_rect(colour = "grey", fill=NA, size=0.5))+
  scale_colour_ggthemr_d()

datatest

###creating figures of all scales
library(ggthemr)

datfig11<-as_tibble(cbind(sefixedchu$V1,sefixedchu$V2,sefixedgad$V2,sefixedbsi$V2,sefixedk6$V2,sefixedphq$V2,sefixedsdq$V2))
datfig1<-datfig11%>%
  mutate(Theta=V1,
         CHU9D=V2,
         GAD7=V3,
         BSIDEP=V4,
         K6=V5,
         PHQ8=V6,
         SDQEMO=V7)%>%
  select(-V1,-V2,-V3,-V4,-V5,-V6,-V7)

datfig1long<-gather(datfig1, scale, se, CHU9D:SDQEMO)

datfig33<-as_tibble(cbind(infofixedchu$V1,infofixedchu$V2,infofixedgad$V2,infofixedbsi$V2,infofixedk6$V2,infofixedphq$V2,infofixedsdq$V2))
datfig3<-datfig33%>%
  mutate(Theta=V1,
         CHU9D=V2,
         GAD7=V3,
         BSIDEP=V4,
         K6=V5,
         PHQ8=V6,
         SDQEMO=V7)%>%
  select(-V1,-V2,-V3,-V4,-V5,-V6,-V7)

datfig3long<-gather(datfig3, scale, info, CHU9D:SDQEMO)

ggthemr('flat')

ggplot(datfig1long, aes(x=Theta, y=se, color=scale))+
  geom_line(size=1) +
  xlab("Theta score")+
  ylab("SE")+
  geom_hline(yintercept=0.44, linetype="dashed",color="green")+
  geom_hline(yintercept=0.54, linetype="dashed",color="orange")+
  geom_hline(yintercept=0.64, linetype="dashed",color="red")+
  xlim(-6,6)+
  ylim(0,5)+
  scale_x_continuous(breaks=c(-6,-3,0,3,6),
                     sec.axis=sec_axis(~ . *10+50, name="T-score", breaks=c(20,50,80)))+
  theme(panel.border = element_rect(colour = "grey", fill=NA, size=0.5))+
  scale_colour_ggthemr_d()

ggplot(datfig3long, aes(x=Theta, y=info, color=scale))+
  geom_line(size=1) +
  xlab("Theta score")+
  ylab("Standard Error")+
  geom_hline(yintercept=11, linetype="dashed",color="green")+
  xlim(-6,6)+
  ylim(0,5)+
  scale_x_continuous(breaks=c(-6,-3,0,3,6),
                     sec.axis=sec_axis(~ . *10+50, name="T-score", breaks=c(20,50,80)))+
  theme(panel.border = element_rect(colour = "grey", fill=NA, size=0.5))+
  scale_colour_ggthemr_d()

datfig22<-as_tibble(cbind(scorefixedchu$V1,scorefixedchu$V2,scorefixedgad$V2,scorefixedbsi$V2,scorefixedk6$V2,scorefixedphq$V2,scorefixedsdq$V2))
datfig2<-datfig22%>%
  mutate(Theta=V1,
         CHU9D=V2,
         GAD7=V3,
         BSIDEP=V4,
         K6=V5,
         PHQ8=V6,
         SDQEMO=V7)%>%
  select(-V1,-V2,-V3,-V4,-V5,-V6,-V7)

datfig2long<-gather(datfig2, scale, sum, CHU9D:SDQEMO)


ggplot(datfig2long, aes(x=Theta, y=sum, color=scale))+
  geom_line(size=1.0) +
  xlim(-6,6)+
  xlab("Theta score")+
  ylab("Total score")+
  #geom_vline(xintercept=1.35678392, color="#aaa488", linetype=2, size=1.0)+
 # geom_vline(xintercept=1.47738693, color="#b2432f", linetype=2, size=1.0)+
 # geom_vline(xintercept=1.17587940, color="#3a6589", linetype=2, size=1.0)+
 # geom_vline(xintercept=2.02010050, color="#9b5672", linetype=2, size=1.0)+
  scale_x_continuous(breaks=c(-6,-3,0,3,6),
                     sec.axis=sec_axis(~ . *10+50, name="T-score", breaks=c(20,50,80)))+
  theme(panel.border = element_rect(colour = "grey", fill=NA, size=0.5))+
  scale_colour_ggthemr_d()

swatch()


### agreement between EAP score in calibration sample
# CHU9D and Full items
ICC(as.matrix(cbind(chueap[,1],distresseap[,1])))
cor(as.matrix(cbind(chueap[,1],distresseap[,1])),use="complete.obs")
mean(chueap[,1]-distresseap[,1], na.rm=TRUE)
SD(chueap[,1]-distresseap[,1])

Avg <- (chueap[,1] + distresseap[,1]) / 2
Dif <- chueap[,1] - distresseap[,1]
sample_df_chu<-as.data.frame(cbind(Avg, Dif, 1))

ggplot(sample_df_chu, aes(x = Avg, y = Dif)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(sample_df_chu$Dif, na.rm=TRUE), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_chu$Dif, na.rm=TRUE) - (1.96 * SD(sample_df_chu$Dif)), colour = "red", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_chu$Dif, na.rm=TRUE) + (1.96 * SD(sample_df_chu$Dif)), colour = "red", size = 0.5) +
  ylab("Diff. Between Measures") +
  xlab("Average Measure")

# GAD and Full items
ICC(as.matrix(cbind(gadeap[,1],distresseap[,1])))
mean(gadeap[,1]-distresseap[,1], na.rm=TRUE)
SD(gadeap[,1]-distresseap[,1])

Avg <- (gadeap[,1] + distresseap[,1]) / 2
Dif <- gadeap[,1] - distresseap[,1]
sample_df_gad<-as.data.frame(cbind(Avg, Dif, 2))

ggplot(sample_df_gad, aes(x = Avg, y = Dif)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(sample_df_gad$Dif, na.rm=TRUE), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_gad$Dif, na.rm=TRUE) - (1.96 * SD(sample_df_gad$Dif)), colour = "red", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_gad$Dif, na.rm=TRUE) + (1.96 * SD(sample_df_gad$Dif)), colour = "red", size = 0.5) +
  ylab("Diff. Between Measures") +
  xlab("Average Measure")

# BSI and Full items
ICC(as.matrix(cbind(bsieap[,1],distresseap[,1])))
mean(bsieap[,1]-distresseap[,1], na.rm=TRUE)
SD(bsieap[,1]-distresseap[,1])

Avg <- (bsieap[,1] + distresseap[,1]) / 2
Dif <- bsieap[,1] - distresseap[,1]
sample_df_bsi<-as.data.frame(cbind(Avg, Dif, 3))

ggplot(sample_df_bsi, aes(x = Avg, y = Dif)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(sample_df_bsi$Dif, na.rm=TRUE), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_bsi$Dif, na.rm=TRUE) - (1.96 * SD(sample_df_bsi$Dif)), colour = "red", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_bsi$Dif, na.rm=TRUE) + (1.96 * SD(sample_df_bsi$Dif)), colour = "red", size = 0.5) +
  ylab("Diff. Between Measures") +
  xlab("Average Measure")

# k6 and Full items
ICC(as.matrix(cbind(k6eap[,1],distresseap[,1])))
mean(k6eap[,1]-distresseap[,1], na.rm=TRUE)
SD(k6eap[,1]-distresseap[,1])

Avg <- (k6eap[,1] + distresseap[,1]) / 2
Dif <- k6eap[,1] - distresseap[,1]
sample_df_k6<-as.data.frame(cbind(Avg, Dif, 4))

ggplot(sample_df_k6, aes(x = Avg, y = Dif)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(sample_df_k6$Dif, na.rm=TRUE), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_k6$Dif, na.rm=TRUE) - (1.96 * SD(sample_df_k6$Dif)), colour = "red", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_k6$Dif, na.rm=TRUE) + (1.96 * SD(sample_df_k6$Dif)), colour = "red", size = 0.5) +
  ylab("Diff. Between Measures") +
  xlab("Average Measure")

# phq and Full items
ICC(as.matrix(cbind(phqeap[,1],distresseap[,1])))
mean(phqeap[,1]-distresseap[,1], na.rm=TRUE)
SD(phqeap[,1]-distresseap[,1])

Avg <- (phqeap[,1] + distresseap[,1]) / 2
Dif <- phqeap[,1] - distresseap[,1]
sample_df_phq<-as.data.frame(cbind(Avg, Dif, 5))

ggplot(sample_df_phq, aes(x = Avg, y = Dif)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(sample_df_phq$Dif, na.rm=TRUE), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_phq$Dif, na.rm=TRUE) - (1.96 * SD(sample_df_phq$Dif)), colour = "red", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_phq$Dif, na.rm=TRUE) + (1.96 * SD(sample_df_phq$Dif)), colour = "red", size = 0.5) +
  ylab("Diff. Between Measures") +
  xlab("Average Measure")

# sdq and Full items
ICC(as.matrix(cbind(sdqeap[,1],distresseap[,1])))
mean(sdqeap[,1]-distresseap[,1], na.rm=TRUE)
SD(sdqeap[,1]-distresseap[,1])

Avg <- (sdqeap[,1] + distresseap[,1]) / 2
Dif <- sdqeap[,1] - distresseap[,1]
sample_df_sdq<-as.data.frame(cbind(Avg, Dif, 6))

ggplot(sample_df_sdq, aes(x = Avg, y = Dif)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(sample_df_sdq$Dif, na.rm=TRUE), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_sdq$Dif, na.rm=TRUE) - (1.96 * SD(sample_df_sdq$Dif)), colour = "red", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_sdq$Dif, na.rm=TRUE) + (1.96 * SD(sample_df_sdq$Dif)), colour = "red", size = 0.5) +
  ylab("Diff. Between Measures") +
  xlab("Average Measure")

# producing all bland altman plots in single graph
blandplotc <- rbind(sample_df_chu, sample_df_gad, sample_df_bsi, sample_df_k6, sample_df_phq, sample_df_sdq)
blandplotc$V3 <- factor(blandplotc$V3, labels = c("CHU v Full","GAD v Full","BSI v Full","K6 v Full","PHQ v Full","SDQ v Full"))

baMean <- blandplotc %>%
  group_by(V3) %>%
  summarise(MN=mean(Dif, na.rm=TRUE),
            Std=SD(Dif))

ggplot(blandplotc) +
  geom_point(aes(x = Avg, y = Dif), alpha = 0.5) +
  geom_hline(data=baMean, aes(yintercept = MN), colour = "blue", size = 0.5) +
  geom_hline(data=baMean, aes(yintercept = MN - (1.96 * Std)), colour = "red", size = 0.5) +
  geom_hline(data=baMean, aes(yintercept = MN + (1.96 * Std)), colour = "red", size = 0.5) +
  ylab("Diff. Between Measures") +
  xlab("Average Measure")+
  facet_wrap(~V3)

### testing crosswalks in validation sample using calibration sample parameters

####Creating crosswalks using common metric
data0v<-validation %>%
  select(m0chu1,m0chu2,m0chu3,m0chu4,m0chu5,m0chu6,m0chu7,m0chu8,m0chu9,m0gad71,m0gad72,m0gad73,m0gad74,m0gad75,m0gad76,m0gad77,
         m0h1suic,m0h2lonely,m0h3sad,m0h4anhedonia,m0h5hopeless,m0h6worthless,
         m0k61,m0k62,m0k63,m0k64,m0k65,m0k66,
         m0phq1,m0phq2,m0phq3,m0phq4,m0phq5,m0phq6,m0phq7,m0phq8,
         m0sd13,m0sd16,m0sd24,m0sd3,m0sd8) %>%
  mutate(m0k61=4-m0k61,
         m0k62=4-m0k62,
         m0k63=4-m0k63,
         m0k64=4-m0k64,
         m0k65=4-m0k65,
         m0k66=4-m0k66)

datamiss<-data0v %>%
  mutate(miss= rowSums(!is.na(data0v)))%>%
  filter(miss>0)

##setting up individual scale datasets
chudatv<-validation %>%
  select(m0chu1, m0chu2,m0chu3,m0chu4,m0chu5,m0chu6,m0chu7,m0chu8,m0chu9)
gaddatv<-validation %>%
  select(m0gad71,m0gad72,m0gad73,m0gad74,m0gad75,m0gad76,m0gad77)
bsidepdatv<-validation %>%
  select(m0h1suic,m0h2lonely,m0h3sad,m0h4anhedonia,m0h5hopeless,m0h6worthless)
k6datv<-validation %>%
  select(m0k61,m0k62,m0k63,m0k64,m0k65,m0k66)%>%
  mutate(m0k61=4-m0k61,
         m0k62=4-m0k62,
         m0k63=4-m0k63,
         m0k64=4-m0k64,
         m0k65=4-m0k65,
         m0k66=4-m0k66)
minispdatv<-validation %>%
  select(m0minisp1,m0minisp2,m0minisp3)
phqdatv<-validation %>%
  select(m0phq1,m0phq2,m0phq3,m0phq4,m0phq5,m0phq6,m0phq7,m0phq8)
sdqdatv<-validation %>%
  select(m0sd13,m0sd16,m0sd24,m0sd3,m0sd8)

##crosswalk using all items
distressval1<-mirt(data0v,1, itemtype="graded", pars="values")
distressval2<-distressval1
distressval2$value<-c(coef0uni)

distressval2$est<- FALSE

distresswalkv<-mirt(data0v,1, itemtype="graded", pars=distressval2)
distresseapv<-fscores(distresswalkv, method="EAP", full.scores.SE=TRUE)

##crosswalk using chu9d
chusv1v<-mirt(chudatv,1, itemtype="graded", pars="values")
chusv2v<-chusv1v
chusv2v$value<-c(coef0uni[1:45,], 0, 1)

chusv2v$est<- FALSE

chuwalkv<-mirt(chudatv,1, itemtype="graded", pars=chusv2v)
chueapv<-fscores(chuwalkv, method="EAP", full.scores.SE=TRUE)

##crosswalk using gad
gadsv1v<-mirt(gaddatv,1, itemtype="graded", pars="values")
gadsv2v<-gadsv1
gadsv2v$value<-c(coef0uni[46:73], 0, 1)

gadsv2v$est<- FALSE

gadwalkv<-mirt(gaddatv,1, itemtype="graded", pars=gadsv2v)

gadeapv<-fscores(gadwalkv, method="EAP", full.scores.SE=TRUE)

##crosswalk using bsi
bsisv1v<-mirt(bsidepdatv,1, itemtype="graded", pars="values")
bsisv2v<-bsisv1v
bsisv2v$value<-c(coef0uni[74:103], 0, 1)

bsisv2v$est<- FALSE

bsiwalkv<-mirt(bsidepdatv,1, itemtype="graded", pars=bsisv2v)

bsieapv<-fscores(bsiwalkv, method="EAP", full.scores.SE=TRUE)

##crosswalk using k6
k6sv1v<-mirt(k6datv,1, itemtype="graded", pars="values")
k6sv2v<-k6sv1v
k6sv2v$value<-c(coef0uni[104:133], 0, 1)

k6sv2v$est<- FALSE

k6walkv<-mirt(k6datv,1, itemtype="graded", pars=k6sv2v)
k6eapv<-fscores(k6walkv, method="EAP", full.scores.SE=TRUE)

##crosswalk using phq
phqsv1v<-mirt(phqdatv,1, itemtype="graded", pars="values")
phqsv2v<-phqsv1v
phqsv2v$value<-c(coef0uni[134:165], 0, 1)

phqsv2v$est<- FALSE

phqwalkv<-mirt(phqdatv,1, itemtype="graded", pars=phqsv2v)
phqeapv<-fscores(phqwalkv, method="EAP", full.scores.SE=TRUE)

##crosswalk using sdq
sdqsv1v<-mirt(sdqdatv,1, itemtype="graded", pars="values")
sdqsv2v<-sdqsv1v
sdqsv2v$value<-c(coef0uni[166:182])

sdqsv2v$est<- FALSE

sdqwalkv<-mirt(sdqdatv,1, itemtype="graded", pars=sdqsv2v)
sdqeapv<-fscores(sdqwalkv, method="EAP", full.scores.SE=TRUE)

### agreement between EAP score in validation sample
# chu and Full items
ICC(as.matrix(cbind(chueapv[,1],distresseapv[,1])))
mean(chueapv[,1]-distresseapv[,1], na.rm=TRUE)
SD(chueapv[,1]-distresseapv[,1])

Avg <- (chueapv[,1] + distresseapv[,1]) / 2
Dif <- chueapv[,1] - distresseapv[,1]
sample_df_chuv<-as.data.frame(cbind(Avg, Dif, 1))

ggplot(sample_df_chuv, aes(x = Avg, y = Dif)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(sample_df_chuv$Dif, na.rm=TRUE), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_chuv$Dif, na.rm=TRUE) - (1.96 * SD(sample_df_chuv$Dif)), colour = "red", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_chuv$Dif, na.rm=TRUE) + (1.96 * SD(sample_df_chuv$Dif)), colour = "red", size = 0.5) +
  ylab("Diff. Between Measures") +
  xlab("Average Measure")

# gad and Full items
ICC(as.matrix(cbind(gadeapv[,1],distresseapv[,1])))
mean(gadeapv[,1]-distresseapv[,1], na.rm=TRUE)
SD(gadeapv[,1]-distresseapv[,1])

Avg <- (gadeapv[,1] + distresseapv[,1]) / 2
Dif <- gadeapv[,1] - distresseapv[,1]
sample_df_gadv<-as.data.frame(cbind(Avg, Dif, 2))

ggplot(sample_df_gadv, aes(x = Avg, y = Dif)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(sample_df_gadv$Dif, na.rm=TRUE), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_gadv$Dif, na.rm=TRUE) - (1.96 * SD(sample_df_gadv$Dif)), colour = "red", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_gadv$Dif, na.rm=TRUE) + (1.96 * SD(sample_df_gadv$Dif)), colour = "red", size = 0.5) +
  ylab("Diff. Between Measures") +
  xlab("Average Measure")

# bsi and Full items
ICC(as.matrix(cbind(bsieapv[,1],distresseapv[,1])))
mean(bsieapv[,1]-distresseapv[,1], na.rm=TRUE)
SD(bsieapv[,1]-distresseapv[,1])

Avg <- (bsieapv[,1] + distresseapv[,1]) / 2
Dif <- bsieapv[,1] - distresseapv[,1]
sample_df_bsiv<-as.data.frame(cbind(Avg, Dif, 3))

ggplot(sample_df_bsiv, aes(x = Avg, y = Dif)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(sample_df_bsiv$Dif, na.rm=TRUE), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_bsiv$Dif, na.rm=TRUE) - (1.96 * SD(sample_df_bsiv$Dif)), colour = "red", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_bsiv$Dif, na.rm=TRUE) + (1.96 * SD(sample_df_bsiv$Dif)), colour = "red", size = 0.5) +
  ylab("Diff. Between Measures") +
  xlab("Average Measure")

# k6 and Full items
ICC(as.matrix(cbind(k6eapv[,1],distresseapv[,1])))
mean(k6eapv[,1]-distresseapv[,1], na.rm=TRUE)
SD(k6eapv[,1]-distresseapv[,1])

Avg <- (k6eapv[,1] + distresseapv[,1]) / 2
Dif <- k6eapv[,1] - distresseapv[,1]
sample_df_k6v<-as.data.frame(cbind(Avg, Dif, 4))

ggplot(sample_df_k6v, aes(x = Avg, y = Dif)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(sample_df_k6v$Dif, na.rm=TRUE), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_k6v$Dif, na.rm=TRUE) - (1.96 * SD(sample_df_k6v$Dif)), colour = "red", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_k6v$Dif, na.rm=TRUE) + (1.96 * SD(sample_df_k6v$Dif)), colour = "red", size = 0.5) +
  ylab("Diff. Between Measures") +
  xlab("Average Measure")

# phq and Full items
ICC(as.matrix(cbind(phqeapv[,1],distresseapv[,1])))
mean(phqeapv[,1]-distresseapv[,1], na.rm=TRUE)
SD(phqeapv[,1]-distresseapv[,1])

Avg <- (phqeapv[,1] + distresseapv[,1]) / 2
Dif <- phqeapv[,1] - distresseapv[,1]
sample_df_phqv<-as.data.frame(cbind(Avg, Dif, 5))

ggplot(sample_df_phqv, aes(x = Avg, y = Dif)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(sample_df_phqv$Dif, na.rm=TRUE), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_phqv$Dif, na.rm=TRUE) - (1.96 * SD(sample_df_phqv$Dif)), colour = "red", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_phqv$Dif, na.rm=TRUE) + (1.96 * SD(sample_df_phqv$Dif)), colour = "red", size = 0.5) +
  ylab("Diff. Between Measures") +
  xlab("Average Measure")


# sdq and Full items
ICC(as.matrix(cbind(sdqeapv[,1],distresseapv[,1])))
mean(sdqeapv[,1]-distresseapv[,1], na.rm=TRUE)
SD(sdqeapv[,1]-distresseapv[,1])

Avg <- (sdqeapv[,1] + distresseapv[,1]) / 2
Dif <- sdqeapv[,1] - distresseapv[,1]
sample_df_sdqv<-as.data.frame(cbind(Avg, Dif, 6))

ggplot(sample_df_sdqv, aes(x = Avg, y = Dif)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = mean(sample_df_sdqv$Dif, na.rm=TRUE), colour = "blue", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_sdqv$Dif, na.rm=TRUE) - (1.96 * SD(sample_df_sdqv$Dif)), colour = "red", size = 0.5) +
  geom_hline(yintercept = mean(sample_df_sdqv$Dif, na.rm=TRUE) + (1.96 * SD(sample_df_sdqv$Dif)), colour = "red", size = 0.5) +
  ylab("Diff. Between Measures") +
  xlab("Average Measure")

# producing all bland altman plots in single graph
blandplotv <- rbind(sample_df_chuv, sample_df_gadv, sample_df_bsiv, sample_df_k6v, sample_df_phqv, sample_df_sdqv)
blandplotv$V3 <- factor(blandplotv$V3, labels = c("CHU v Full","GAD v Full","BSI v Full","K6 v Full","PHQ v Full","SDQ v Full"))

baMeanv <- blandplotv %>%
  group_by(V3) %>%
  summarise(MN=mean(Dif, na.rm=TRUE),
            Std=SD(Dif))

ggplot(blandplotv) +
  geom_point(aes(x = Avg, y = Dif), alpha = 0.5) +
  geom_hline(data=baMeanv, aes(yintercept = MN), colour = "blue", size = 0.5) +
  geom_hline(data=baMeanv, aes(yintercept = MN - (1.96 * Std)), colour = "red", size = 0.5) +
  geom_hline(data=baMeanv, aes(yintercept = MN + (1.96 * Std)), colour = "red", size = 0.5) +
  ylab("Diff. Between Measures") +
  xlab("Average Measure")+
  facet_wrap(~V3)




