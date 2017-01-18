### CDA Final 
### Authour: Issac Li
### Created: Dec 2016
###  
require(ResourceSelection)
require(ggplot2)
require(reshape2)
require(data.table)
source("http://peterhaschke.com/Code/multiplot.R")

## Read dataset 
ml <- fread("~/Downloads/final dataset.txt")
dataset<- subset(ml, select = -c(Treatment,Radiation))

## Add the response outcome variable
## Cases: outcome = 0 Controls: outcome = 1
dataset$outcome=-1
dataset$outcome[dataset$Vital_status == 0 & dataset$Survival_mo<=48] <- 0
dataset$outcome[dataset$Survival_mo>48] <- 1

nrow(dataset)
nrow(dataset[dataset$outcome==1,])
nrow(dataset[dataset$outcome==-1,])

## Remove data points with outcome = -1

sample = dataset[outcome>=0]
sample=data.table(data.frame(unclass(sample),stringsAsFactors=TRUE))

## Make weeks: 7-day interval (roughly)
## Adjusted boundaries for balance of observations in each "week"
sample$Chemo_weeks = ifelse(sample$Chemo_days<=28,4,NA)
sample[Chemo_days>28]$Chemo_weeks = floor(sample[Chemo_days>28,Chemo_days]/7)+1
sample[Chemo_weeks==11]$Chemo_weeks = 10
sample[Chemo_weeks>=12]$Chemo_weeks = 11

head(sample[,list(Chemo_days,Chemo_weeks)],10)
par(mfrow=c(1,2))
hist(log(sample$Chemo_days),nclass=20,xlab = "Log Days before Chemotherapy",main = "A                                                ")
hist(sample$Chemo_weeks,breaks = c(3,4,5,6,7,8,9,10,11),xlab = "Weeks before Chemotherapy",main="B                                                ")

## Graphically Explore the difference between control group and case group
tmp <- melt(data.frame(sample[Tumor_size<=200])[, c("outcome", "Tumor_size", "Chemo_days", "Age", "Chemo_weeks")],
            id.vars="outcome")
  

p1<-ggplot(tmp, aes(factor(outcome), y = value, fill=factor(outcome))) +
  geom_boxplot() +
  theme(axis.text.x = element_text(size=12))
  facet_wrap(~variable, scales="free_y")

### No apparent association between chemo_days/weeks and survival outcome
p1


## There is violation of linearity of predictors in the data 
### Histology and Sex
### Make are more prone to have squamous cell carcinoma
p2<- ggplot(data=sample[Grade!="Unknown"],aes(x=outcome,y=Chemo_days,alpha=Age,size=Tumor_size,color=Histology))+
  geom_jitter()+
  geom_boxplot(aes(group=c(outcome)),alpha=0.3)+
  guides(alpha="none")+
  facet_grid(Sex~Grade)

### Pneumonectomy leads to longer Chemo_days 
p3<-ggplot(data=sample[Grade!="Unknown"],aes(x=Surgery,y=Chemo_days))+geom_jitter()+geom_boxplot(alpha=0.3)+facet_wrap(~Sex)


### Histology also influence tumor size
p4<-ggplot(data=sample[Grade!="Unknown"&Tumor_size<250],aes(x=Histology,y=Tumor_size,fill=Histology))+
  geom_boxplot(alpha=0.3)+facet_wrap(~Sex)

### Larger ages lead to more chemo days
p5<-ggplot(data=sample[Grade!="Unknown"&Tumor_size<250],aes(x=Age,y=Chemo_days),fill=Histology)+geom_point(alpha=0.3)+
  stat_smooth()+
  facet_wrap(CD_score~Sex)

multiplot(p2,p3,p4,p5,cols = 2)

## There are sufficient number of unknowns in grade which should be imputed first
## Explore the data with no unknown grade
impute=droplevels(subset(sample,Grade=='Unknown'))
impute_train=droplevels(sample[Grade!='Unknown'])

table1=impute_train[,.N, by=Grade]
table2=sample[,.N,by=Age]
table3=sample[,mean(Chemo_days),by=c("Sex","Histology","Path_stage","CD_score","outcome")][order(Sex,Histology,Path_stage)]

## Tumor size greater than 200 are definitely outliers and should be removed.
hist(sample$Tumor_size[sample$Tumor_size<=200],breaks = seq(0,200,5))

## Impute unknown grade by 
library(MASS)
require(scales)
formula_grade=Grade ~ Path_stage + Tumor_size+ Histology + CD_score+ Primary_site+Surgery+Age+Sex+(Facility_type+Facility_location)^2
train_ind <- sample(seq_len(nrow(impute_train)), size = floor(nrow(impute_train)*0.75))
impute_test=impute_train[-train_ind,]
impute_train=impute_train[train_ind,]


olg=polr(formula=formula_grade, data=impute_train,Hess=TRUE)
summary(olg)

olg.bwd=step(olg,direction = "backward")

olg.null=polr(formula = Grade ~ Path_stage + Tumor_size, data = impute_train)
summary(olg.null)

ctable <- coef(summary(olg.bwd))
### calculate and store p values
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2

### combined table
(ctable <- cbind(ctable, "p value" = p))

### odds ratio
exp(coef(olg.bwd))

pred_g=predict(olg.bwd,impute_test,type="class")
pred_g=predict(olg.bwd,impute_train,type="class")

confMat = table(impute_train$Grade,pred_g)
confMat
rowSums(confMat)
colSums(confMat)
### Gives prediction accuracy
sum(diag(confMat))/sum(confMat)

anova(olg.null,olg.bwd,test = "Chisq")
anova(olg.bwd,olg,test = "Chisq")

pred_gu=predict(olg.bwd,impute,type="class")
roc_pred <- prediction(pred_gu, impute_train$Grade)
plot(performance(roc_pred, measure="tpr", x.measure="fpr"), colorize=TRUE)
abline(0,1,lty=2)
## Add back imputed grades

imputed=impute
imputed$Grade=pred_gu

sample.g=droplevels(rbind(sample[Grade!="Unknown"],imputed))

# Socioeconomic Status ----------------------

## It is interesting to see that in most of the cases people using a private insurance have better 
## outcome
ggplot(sample.g,aes(x=Insurance,y=..count.., fill=as.factor(outcome)))+ 
  geom_bar(alpha=c(0.6),position='dodge') +
  facet_grid(Path_stage ~.,margins = FALSE) +
  labs(fill="Outcome")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))

## Check out for education 

ggplot(sample.g,aes(x=Education,y=..count.., fill=as.factor(outcome)))+ 
  geom_bar(alpha=c(0.6),position='dodge') +
  facet_grid(Path_stage ~ Chemo_weeks,margins = FALSE) +
  labs(fill="Outcome")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))

####  Diagnose of cancer peaked at 2008?
ggplot(data=tmp2,aes(Year_diag,N,color=Insurance))+geom_line()

ggplot(sample.g[Year_diag<=2010],aes(x=Year_diag,y=..count.., fill=as.factor(outcome)))+ 
  geom_bar(alpha=c(0.6),position='dodge') +
  facet_grid(Facility_type ~ Facility_location,margins = FALSE,scales="free_y") +
  labs(fill="Outcome")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))

ggplot(sample.g[Year_diag<=2010],aes(x=Facility_location,y=Survival_mo, fill=Facility_location))+
         geom_boxplot()+
  facet_grid(Facility_type ~ .,margins = FALSE,scales="free_y")+
  coord_flip()
  

socioecon=glm(formula = outcome ~ Education+Insurance+Income+Facility_location+Facility_type,data = sample.g[Education!="Unknown"],family = "binomial")

### Get aod for wald test
require(aod)
wald.test(b = coef(socioecon), Sigma = vcov(socioecon), Terms = 2:4)

wald.test(b = coef(socioecon), Sigma = vcov(socioecon), Terms = 5:9)

wald.test(b = coef(socioecon), Sigma = vcov(socioecon), Terms = 10:11)

wald.test(b = coef(socioecon), Sigma = vcov(socioecon), Terms = 12:15)

### Facility Type
wald.test(b = coef(socioecon), Sigma = vcov(socioecon), Terms = 16)

### ROC curve
roc_pred <- prediction(socioecon$fitted.values, sample.g[Education!="Unknown"]$outcome)
plot(performance(roc_pred, measure="tpr", x.measure="fpr"), colorize=TRUE)
abline(0,1,lty=2)


# Mixed Effects Model ------------------------

tmp1 <- sample.g[Tumor_size<=200,.N,by=c("Facility_type,Facility_location,Insurance")][order(Facility_type,Facility_location,Insurance)]
tmp2<-sample.g[Tumor_size<=200,.N,by=c("Year_diag,Insurance")][order(Year_diag)]

require(lme4)
# estimate the model and store results in mixed effect model
sample.g$sAge = scale(sample.g$Age)
sample.g$sTumor_size = scale(sample.g$Tumor_size)

mmod <- glmer(outcome ~ Path_stage + sTumor_size + Age+Chemo_weeks+Histology+Grade+
                (1|Facility_location)+(1|Year_diag)+( 1|Facility_type)+(1|Income), data = sample.g[Tumor_size<=250&Year_diag<2011], family = binomial, control = glmerControl(optimizer = "bobyqa"),
           nAGQ = 1)

# print the model results without correlations among fixed effects
print(mmod, corr = FALSE)

se <- sqrt(diag(vcov(mmod)))
# table of estimates with 95% CI
(tab <- cbind(Est = fixef(mmod), LL = fixef(mmod) - 1.96 * se, UL = fixef(mmod) + 1.96 *
                se))
# Odds Ratios
exp(tab)
pred_m=predict(mmod,sample.g[Tumor_size<=250&Year_diag<2011],type="response")  
predAccu(data = sample.g[Tumor_size<=250&Year_diag<2011],prob = pred_m,p_i = 0.6)

roc_pred <- prediction(pred_m,sample.g[Tumor_size<=250&Year_diag<2011]$outcome)
plot(performance(roc_pred, measure="tpr", x.measure="fpr"), colorize=TRUE)

df=data.frame(cbind(sample.g[Tumor_size<=250&Year_diag<2011],pred_m))
df$pred = df$pred_m
df$outcome = as.numeric(df$outcome)-1

plot_pred_type_distribution <- function(df,threshold) {
  v <- rep(NA, nrow(df))
  v <- ifelse(df$pred >= threshold & df$outcome == 1, "TP", v)
  v <- ifelse(df$pred >= threshold & df$outcome == 0, "FP", v)
  v <- ifelse(df$pred < threshold & df$outcome == 1, "FN", v)
  v <- ifelse(df$pred < threshold & df$outcome == 0, "TN", v)
  
  df$pred_type <- v
  
  ggplot(data=df, aes(x=as.factor(outcome), y=pred)) + 
    geom_jitter(aes(color=pred_type), alpha=0.6) +
    geom_violin(fill=rgb(1,1,1,alpha=0.3), color=NA) + 
    geom_hline(yintercept=threshold, color="red", alpha=0.6) +
    scale_color_discrete(name = "type") +
    labs(title=sprintf("Threshold at %.2f", threshold))
}

plot_pred_type_distribution(df,0.5)


sample.g$Chemo_weeks=factor(sample.g$Chemo_weeks)

## Ordinary logistic regressions on two sub models ---------

formula_log = outcome ~ (Age+Sex+CD_score+Path_stage+Grade+Tumor_size+Surgery+Chemo_weeks+Histology+Primary_site)^2

c1=data.table(sample.g[Facility_type!="Academic"])[,c("Facility_type","Facility_location"):=NULL]
c2=data.table(sample.g[Facility_type=="Academic"])[,c("Facility_type","Facility_location"):=NULL]

logfit_c1f<-glm(formula_log,data=c1,family="binomial")
logfit_c2f<-glm(formula_log,data=c2,family="binomial")
logfit_null <- glm(outcome~1,data=c2,family = "binomial")

bwd.c1<-step(logfit_c1f,direction = "backward")
bwd.c2<-step(logfit_null,direction = "forward",scope = formula(bwd.c1))


formula_new = outcome ~ Age + Sex + CD_score + Path_stage + Grade + Tumor_size + 
  Surgery + Chemo_weeks + Histology + Primary_site + Age:CD_score + 
  Sex:CD_score + Sex:Surgery  + Tumor_size:Histology 



logfit_c1<-glm(formula=formula_new,data=c1,family="binomial")
logfit_c2<-glm(formula_new,data=c2,family="binomial")


### Hosmer-Lemeshow Test
### p-value >0.05 meaning no significant lack of fit
hoslem.test(as.numeric(c1$outcome)-1, logfit_c1$fitted.values, g=10)
hoslem.test(as.numeric(c2$outcome)-1, logfit_c2$fitted.values, g=10)


anova(logfit_c1,logfit_c1f,test="Chisq")
anova(logfit_c2,logfit_c2f,test="Chisq")

pred_c1=predict(logfit_c1f,c1,type="response")  
predAccu(data = c1,prob = pred_c1,p_i = 0.5)

pred_c2=predict(logfit_c2,c2,type="response")
predAccu(data = c2,prob = pred_c2,p_i = 0.5)

par(mfrow=c(2,1),margin = 1)
roc_pred <- prediction(logfit_c1$fitted.values, c1$outcome)
plot(performance(roc_pred, measure="tpr", x.measure="fpr"), colorize=TRUE)
abline(0,1,lty=2)
roc_pred <- prediction(logfit_c2$fitted.values, c2$outcome)
plot(performance(roc_pred, measure="tpr", x.measure="fpr"), colorize=TRUE)
abline(0,1,lty=2)


# Graphically Explore Mixed Effect model -----------------
## Predict Probobilities in various cases
tmpdat <- sample.g[Tumor_size<=250&Year_diag<2011,list(Path_stage, sTumor_size, Age, Histology, Chemo_weeks,Facility_type,Facility_location,Income,Year_diag,Education,Grade)]
jvalues <- with(sample.g, seq(from = min(as.numeric(Chemo_weeks)), to = max(as.numeric(Chemo_weeks)), length.out = 6))

pp <- lapply(jvalues, function(j) {
  tmpdat$Chemo_weeks <- as.integer(j)
  predict(mmod, newdata = tmpdat, type = "response")
})

pp=predict(mmod, newdata = tmpdat, type = "response")
tmptb=cbind(tmpdat,unlist(pp))

sapply(pp[c(1, 2, 3, 4, 5, 6)], mean)

plotdat=tmptb[,list(Mean=mean(V2),LB=quantile(V2,0.25),UB=quantile(V2,0.75)),by=list(Path_stage,Chemo_weeks)]

p6<-ggplot(plotdat, aes(x = as.numeric(Chemo_weeks), y = Mean)) +
  geom_line(aes(colour = Path_stage), size = 1)+
  geom_ribbon(aes(ymin = LB, ymax = UB, fill = Path_stage), alpha = .35) +
  labs(y="Predicted Prob. of Survival",x="Weeks After Resection")+
  scale_x_discrete(limits=c("",5,6,7,8,9,10))+
  theme(axis.text.x = element_text(hjust=0,size=12),
        axis.title.x = element_text(size=13),
        axis.title.y = element_text(size=13))+
  ylim(c(0.3, 0.8)) + facet_wrap(~Path_stage)

p7<-ggplot(plotdat1, aes(x = as.numeric(Chemo_weeks), y = Mean)) +
  geom_line(aes(colour = Facility_type), size = 1)+
  geom_ribbon(aes(ymin = LB, ymax = UB, fill = Facility_type), alpha = .35) +
  labs(y="Predicted Prob. of Survival",x="Weeks After Resection")+
  scale_x_discrete(limits=c("",5,6,7,8,9,10))+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=13),
        axis.title.y = element_text(size=13))+
  ylim(c(0.3, 0.8)) + facet_wrap(~Facility_type)
plotdat1=tmptb[,list(Mean=mean(V2),LB=quantile(V2,0.25),UB=quantile(V2,0.75)),by=list(Facility_type,Chemo_weeks)]

plotdat2=tmptb[,list(Mean=mean(V2),LB=quantile(V2,0.25),UB=quantile(V2,0.75)),by=list(Year_diag,Chemo_weeks)]

plotdat3=tmptb[,list(Mean=mean(V2),LB=quantile(V2,0.25),UB=quantile(V2,0.75)),by=list(Facility_location,Chemo_weeks)]

sample_proportion=sample.g[Tumor_size<=250&Year_diag<2011,list(Outcome=sum(as.numeric(outcome)-1),Num=.N),by=list(Chemo_weeks,Facility_type,Facility_location,Year_diag)] 

p8<-ggplot(plotdat2[Year_diag>2004&Year_diag<2011], aes(x = as.numeric(Chemo_weeks), y = Mean)) +
  geom_line(aes(colour = as.character(Year_diag)), size = 1)+
  #geom_line(data=sample_proportion[Year<2011,prop=sum(V1)/.N,by=Year_diag],aes(x=Chemo_weeks,y=prop))+
  geom_ribbon(aes(ymin = LB, ymax = UB),fill="grey", alpha = .35) +
  labs(y="Predicted Prob. of Survival",x="Weeks After Resection")+
  scale_x_discrete(limits=c("",5,6,7,8,9,10))+
  scale_color_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=13),
        axis.title.y = element_text(size=13))+
  facet_wrap(~Year_diag)+
  ylim(c(0.2, 0.75))

p9<-ggplot(plotdat3, aes(x = as.numeric(Chemo_weeks), y = Mean)) +
  geom_line(aes(colour = Facility_location), size = 1)+
  geom_ribbon(aes(ymin = LB, ymax = UB),fill="grey", alpha = .35) +
  labs(y="Predicted Prob. of Survival",x="Weeks After Resection")+
  scale_x_discrete(limits=c("",5,6,7,8,9,10))+
  scale_color_brewer(palette = "Set1")+
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=13),
        axis.title.y = element_text(size=13))+
  facet_wrap(~Facility_location)+
  ylim(c(0.2, 0.75))

multiplot(p6,p7,p8,p9,cols = 2)

# Results Characteristic Table --------------
require(tableone)
listVars = c("Age","Survival_mo","Tumor_size","Chemo_days","Sex","CD_score","Path_stage","Chemo_weeks","Grade","Histology",
             "Primary_site","Surgery","Facility_type","Facility_location","Insurance","Income","Education","outcome")

catVars = c("Sex","CD_score","Path_stage","Chemo_weeks","Grade","Histology",
            "Primary_site","Surgery","Facility_type","Facility_location","Insurance","Income","Education","outcome")



group1=sample.g[as.numeric(Chemo_weeks)<4]
group2=sample.g[as.numeric(Chemo_weeks)==4]
group3=sample.g[as.numeric(Chemo_weeks)>4]

group1$Interval="Early Start"
group2$`Interval`="Optimal"
group3$`Interval`="Delayed Start"

groups = rbind(group1,group2,group3)

table_1 <- CreateTableOne(listVars, groups, catVars, strata = c("Interval"))

