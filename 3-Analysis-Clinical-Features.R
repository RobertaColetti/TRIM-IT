### Step 3 of TRIM-IT: ANALYSIS OF CLUSTER CLINICAL CHARACTERISTICS AND SURVIVAL 

# Library section
library(ggsurvfit) 
library(survival)
library(survminer) 

# Load data and prior results.
load("~/Results/Kmeans-GBM-results.RData") #K-means results
load("~/Data/Glioma-clinic-TCGA.RData") # clinical TCGA data

# Prepare the data.
#create a dataframe with name of the patients and clusters 
df.km.MOMIP <- as.data.frame(cbind(names=rownames(as.data.frame(res.km.MOMIP$cluster)), cluster = res.km.MOMIP$cluster))
#reduce the clinical data to the patients of interest
reduced.clinic=Glioma_clinic[which(Glioma_clinic$names %in% df.km.MOMIP$names),]
reduced.clinic=reduced.clinic[order(reduced.clinic$names),] #sort in alphabetical order
all(reduced.clinic$names==df.km.MOMIP$names) #it MUST be TRUE
reduced.clinic$cluster=df.km.MOMIP$cluster #add cluster column
#create column with simplified histological information (LGG vs GBM)
reduced.clinic$simplified.hist=NA
reduced.clinic$simplified.hist[reduced.clinic$histological_typ=="treated primary gbm"] = "gbm"
reduced.clinic$simplified.hist[reduced.clinic$histological_typ=="untreated primary (de novo) gbm"] = "gbm"
reduced.clinic$simplified.hist[reduced.clinic$histological_typ=="astrocytoma"] = "lgg"
reduced.clinic$simplified.hist[reduced.clinic$histological_typ=="oligoastrocytoma"] = "lgg"
reduced.clinic$simplified.hist[reduced.clinic$histological_typ=="oligodendroglioma"] = "lgg"
#check cluster assignement (both conditions MUST be TRUE)
all(reduced.clinic$names==df.km.MOMIP$names)
all(reduced.clinic$cluster==df.km.MOMIP$cluster)

## .............Correlation analysis...................
# Clinical features to consider: histology, gender.
#visualization
table(reduced.clinic$cluster,reduced.clinic$gender)
table(reduced.clinic$cluster,reduced.clinic$simplified.hist)
#correlation with histology
tbl.hist = as.matrix(table(reduced.clinic$cluster,reduced.clinic$simplified.hist))
dimnames(tbl.hist) = list(Clust=c(1, 2,3), Histology=c('GBM', 'LGG'))
chi2.hist = chisq.test(tbl.hist, correct=F)
c(chi2.hist$statistic, chi2.hist$p.value)
sqrt(chi2.hist$statistic / sum(tbl.hist)) #cramer V (the smaller v, the lower the correlation)
#correlation with gender
tbl.gen = as.matrix(table(reduced.clinic$cluster,reduced.clinic$gender))
dimnames(tbl.gen) = list(Clust=c(1, 2,3), Sex=c('F','M'))
chi2.gen = chisq.test(tbl.gen, correct=F)
c(chi2.gen$statistic, chi2.gen$p.value)
sqrt(chi2.gen$statistic / sum(tbl.gen)) #cramer V (the smaller v, the lower the correlation)


## .............Survival Analysis...................
#
# Prepare survival status to be used by the R functions.
#
# The current clinical information retrieved by TCGA contains the time of survival divided into two columns:
# if the event has occurred (survival_status = 1), the time is recorded in the "days_to_death" column;
# if the event has not been identified (vital_status = 0) the information is reported in the "days_to_last_followup" column.
# To ensure compatibility with the survfit2 function in R, which we are using for survival analysis, we need to consolidate this information into a single column, "time".
# 
reduced.clinic$days_to_death[is.na(reduced.clinic$days_to_death)]=0 #change NA to 0
reduced.clinic$days_to_last_followup[is.na(reduced.clinic$days_to_last_followup)]=0  #change NA to 0
#We overwrite the 1st column of the clinical data (not used in this analysis) with the time of survival by summing the values into the two columns "days_to_death" and "days_to_last_followup"
reduced.clinic[,1]=as.numeric(reduced.clinic$days_to_death)+as.numeric(reduced.clinic$days_to_last_followup)
colnames(reduced.clinic)[1]="time" 
reduced.clinic$vital_status=as.numeric(reduced.clinic$vital_status) 

#Filter the data for the survival analysis to retain only the necessary information
data_surv=reduced.clinic[,c(1,3,15,14)] #days followup (time), survival status, clusters, sample name
#Kaplan Meier curves of the GBM clusters
fit2 <- survfit2(Surv(time, vital_status) ~ cluster, data = data_surv)
#plot curves
ggsurvfit(fit2)+ add_pvalue(caption = "Log-rank {p.value}") 


#.........Comparison of cluster survival curves........
#
# Perform pairwise log-rank tests
pairwise_survdiff(Surv(time, vital_status) ~ cluster, data = data_surv, p.adjust.method = "BH") #BH= Benjamini-Hochberg adjustment 

# Cox regression to further study cluster 2 vs cluster 3 
data_surv$cluster <- as.factor(data_surv$cluster)  #convert to factor
data_surv$cluster <- relevel(data_surv$cluster, ref = "2")  # set cluster 2 as the reference
# Fit the Cox model
cox_model2 <- coxph(Surv(time, vital_status) ~ cluster, data = data_surv)
summary(cox_model2) #results

