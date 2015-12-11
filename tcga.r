#Project TCGA Breast Cancer
#INCLUDES
require(xgboost)
library(glasso)
library(lars)
library(survival)
library(kernlab)
library(ROCR)
library(e1071)
library(glmnet)

#Surv()
#LARS -> CLUSTER -> CLASSIFICATION
#READ PREPROCESSED DATA #rm(list=ls())

mrna_dir = "C:/Users/son/Desktop/2015CCA/mrna_ver1.txt"
micro_rna_dir = "C:/Users/son/Desktop/2015CCA/cnv_ver1.txt"
death_dir = "C:/Users/son/Desktop/2015CCA/patient_ver1.txt"
clinic_dir = "C:/Users/son/Desktop/2015CCA/summary_ver2.txt"

temp_m <- read.table(file=mrna_dir, header = F, nrows = 5,sep="\t")
classes <- sapply(temp_m, class)
mrna <- read.table(file=mrna_dir, header = F, colClasses = classes,sep="\t")

#temp_mi <- read.table(file=micro_rna_dir, header = F, nrows = 5,sep="\t")
#classes2 <- sapply(temp_mi, class)
#micro_rna <- read.table(file=micro_rna_dir, header = F, colClasses = classes2,sep="\t")
micro_rna <- read.table(file=micro_rna_dir, header = F,sep="\t")


death <- read.table(file=death_dir, header = F, sep="\t")
clinic <- read.table(file=clinic_dir, header = F, sep="\t")
death <- rbind(replace(death[1,],1,"Patient"),death[2:nrow(death),] )
clinic <- rbind(replace(clinic[1,],1,"Patient"),clinic[2:nrow(clinic),] )


#SO FAR WE HAVE
# MRNA, micro_RNA, DEATH, CLINIC / TABLEs

#
mrna_t = t(mrna)
micro_rna_t = t(micro_rna)
# PREPROCESS INPUTS

mrna_pat = mrna_t[,1]
micro_rna_pat = micro_rna_t[,1]
death_pat = death[,1]

make01 <- function (x) paste(x,"-01",sep="")
death_name = sapply(death_pat,make01)
death_name = replace(death_name,1,"Patient")
fin_death = cbind(death_name, death[,2:ncol(death)])
fin_death = fin_death[2:nrow(fin_death),]
colnames(fin_death) <- c("Patient","vital_status","days_to_last_followup","days_to_death","radiation_therapy")

### MADE DEATH

namefunc <- function(x) substr(x,1,15)
mrna_pat_name = sapply(mrna_pat[2:length(mrna_pat)], namefunc)
micro_rna_pat_name = sapply(micro_rna_pat[2:length(micro_rna_pat)], namefunc)

fin_mrna = cbind(mrna_pat_name, mrna_t[2:nrow(mrna_t),2:ncol(mrna_t)])
fin_mirna = cbind(micro_rna_pat_name, 
		micro_rna_t[2:nrow(micro_rna_t),2:ncol(micro_rna_t)])

colnames(fin_mrna) <- replace(as.vector(mrna_t[1,]),1,"Patient")
colnames(fin_mirna) <- replace(as.vector(micro_rna_t[1,]),1,"Patient")


clinic <- clinic[,c(1,3:ncol(clinic))]
temp <- as.vector(clinic[1,])
clinic <- clinic[2:nrow(clinic),]
colnames(clinic) <- c("Patient","ER_Status","PR_Status",
			"HER2_Status","PAM50_TCGA_mRNASeq_AWG",
			"PAM50_TCGA_2012","PAM50_BHI_RNASeq_Raw",
			"PAM50_BHI_RNASeq_Log2","Age","Age_Diagnosis",
			"Race","sample_type","Menopause_Status")

##COMPLETE fin_death, fin_mirna, fin_mrna, clinic

##################### todo
m=fin_mrna[,1]
mi=fin_mirna[,1]
dea=fin_death[,1]
cli = clinic[,1]


tp1=intersect(m, mi)
tp2=intersect(dea,cli)
tp3=intersect(tp1,tp2)

mrow = match(tp3,m)
mirow = match(tp3,mi)
dearow = match(tp3,dea)
clirow = match(tp3,cli)
#is.unsorted?!

##NORMALIZE with common
fin_mrna = fin_mrna[mrow,]
fin_mirna = fin_mirna[mirow,]
fin_death = fin_death[dearow,]
clinic = clinic[clirow,]
#1077 patients -> 198

#############################################################
#MAKING OBJECTIVE VECTOR

#death's col 3 & 4 is ..
year = 365
days_alive = sapply(sapply(fin_death[,3],as.character),as.numeric)
days_death = sapply(sapply(fin_death[,4],as.character),as.numeric)
alive_col = which(days_alive > year*5)
death_col = which(!is.na(days_death))

#totcol = 5 year as longest
totcol = sort(c(alive_col,death_col))

## RENEWING DATA
fin_mrna = fin_mrna[totcol,]
fin_mirna = fin_mirna[totcol,]
fin_death = fin_death[totcol,]
clinic = clinic[totcol,]


#EXAMPLE CODE
#make all attributes numeric
#apply(extra4, 2, as.numeric)
#sapply(extra4, 2, as.numeric)
#class(extra4) <- "as.numeric"
#storage.mode(extra4) <- "numeric"

num_mrna <- apply(fin_mrna[,2:ncol(fin_mrna)], 2, as.numeric)
num_mirna <- apply(fin_mirna[,2:ncol(fin_mirna)], 2, as.numeric)

## NUM_MRNA , NUM_MIRNA as numerical values

days_alive = sapply(sapply(fin_death[,3],as.character),as.numeric)
days_death = sapply(sapply(fin_death[,4],as.character),as.numeric)

days_death[is.na(days_death)] = days_alive[is.na(days_death)]

#View "days_death" as objective vector
# "goal"

#(a,b]#
death_break2 = c(0,1*year,5*year,50*year)
death_label2 = c("Early","Medium","Long")
death_break = c(0,5*year,50*year)
death_label = c("Early","Long")
death_label_num = c("0","1")
goal = cut(days_death, labels= death_label, breaks = death_break)
goal_num = cut(days_death, labels= death_label_num, breaks = death_break)
goal_num =as.numeric(as.character(goal_num))
goal2 = cut(days_death, labels= death_label2, breaks = death_break2)


parked <- as.matrix(cbind(num_mrna,num_mirna))
bad_parked <- as.matrix(num_mrna)

# to estimate
dlasso <- lars(parked,as.matrix(days_death),type="lasso",use.Gram=FALSE)
mlasso <- lars(as.matrix(num_mrna),as.matrix(days_death),type="lasso",use.Gram=FALSE)
milasso <- glmnet(as.matrix(num_mirna),as.matrix(days_death),alpha=1)

#dlasso
plot(dlasso, cex=0.5)
plot(dlasso, plottype="Cp")
#Lasso bit hard to understand
#dl_coef <- coef(dlasso)
#lasso_glm <- glmnet(parked[1:2,],as.matrix(days_death[1:2]),alpha=1)
plot(milasso)

#gLasso Fail 15.3G
#glasso_parked <- as.matrix(cbind(days_death,parked)) # age,test_time, motor ...
#cov_matrix <- cov(glasso_parked)
#glasso(cov_matrix,1) # L1 regularization

svp
attributes(svp)
alpha(svp)
b(svp)
plot(svp,data= parked)
print(cross(svp))

# try gLasso
# k-means
# SVM + XGBOOST
# clinic (9=age)/(10=age of diagnosis)
good_clinic <- (cbind(as.numeric(as.character(clinic[,9])),
			as.numeric(as.character(clinic[,10]))))

new_parked <- as.matrix(cbind(num_mrna,num_mirna,good_clinic))

svp <- ksvm(parked,as.matrix(goal),type="C-svc",kernel="rbfdot",C=100,scaled=c(),cross=20)
#for demo
svp_two <- ksvm(parked,as.matrix(goal),type="C-svc",kernel="rbfdot",C=20,scaled=c(),cross=20)

#best so far
svp_poly <- ksvm(parked,as.matrix(goal),type="C-svc",kernel="rbfdot",C=10,scaled=c(),cross=20)
svp_cone <- ksvm(parked,as.matrix(goal),type="C-svc",kernel="rbfdot",C=1,scaled=c(),cross=20)
############################################################
SVM with clinic

svp_clinic <- ksvm(new_parked,as.matrix(goal),type="C-svc",kernel="rbfdot",C=100,scaled=c(),cross=20)

##Compare SVP and svp-clinic

svp_bad <- ksvm(bad_parked,as.matrix(goal),type="C-svc",kernel="rbfdot",C=20,scaled=c(),cross=20)


#############################################################

#NEW HORIZON
#XGBOOST

dtest_sub <- xgb.DMatrix(data = parked[11:nrow(parked),],
				 label=goal_num[11:length(goal_num)])
dtrain_sub <- xgb.DMatrix(data = parked[1:10,],
				 label=goal_num[1:10])
dtrain <- xgb.DMatrix(data = parked, label=goal_num)
watchlist <- list(train=dtrain_sub, test=dtest_sub)

bst <- xgboost(data = parked[11:nrow(parked),],label=goal_num[11:length(goal_num)], 
		max.depth=3, eta=1, nthread = 5, nround=20,
	watchlist=watchlist, eval.metric = "error",objective = "binary:logistic")

#eval.metric = "logloss", ommited



history2 <- xgb.cv(data = dtrain, nround=50, nthread = 2, nfold = 20, 
			metrics=list("error","rmse","auc"),max.depth =8, eta = 1, objective = "binary:logistic")

best <- xgb.cv(data = dtrain, nround=20, nthread = 2, nfold = 20, 
			metrics=list("error","rmse","auc"),max.depth =5, eta = 1, objective = "binary:logistic")

history <- xgb.cv(data = dtrain, nround=20, nthread = 2, nfold = 20, 
			metrics=list("error","rmse","auc"),max.depth =3, eta = 1, objective = "binary:logistic")
print(history)

#history as best model
history <- xgb.cv(data = dtrain, nround=20, nthread = 2, nfold = 20, 
			metrics=list("error","rmse","auc"),max.depth =3, eta = 1, objective = "binary:logistic")


#with clinic data
dtrain_cli <- xgb.DMatrix(data = new_parked, label=goal_num)

history_cli <- xgb.cv(data = dtrain_cli, nround=20, nthread = 5, nfold = 30, 
			metrics=list("error"),max.depth =5, 
			eta = 1, objective = "binary:logistic")




names <- c(mrna_t[1,2:ncol(mrna_t)],micro_rna_t[1,2:ncol(micro_rna_t)])
xgb.plot.tree(c(mrna_t[1,2:ncol(mrna_t)],micro_rna_t[1,2:ncol(micro_rna_t)]),model = bst,n_first_tree=1)

# Compute feature importance matrix
importance_matrix <- xgb.importance(names, model = bst)

# Nice graph
xgb.plot.importance(importance_matrix[1:10,])

importance_matrix <- xgb.importance(model = bst)
print(importance_matrix)
xgb.plot.importance(importance_matrix = importance_matrix)

########
hc <-hclust(dist(parked)) 
#hc2 <-hclust(dist(parked),"ave")
#hc_cen <-hclust(dist(parked),"cent")
#plot(hc)
#plot(hc2)
#plot(hc_cen)
#plot(hc,hang=-1)
alable<-cutree(hc, k = 10)
#colMeans(rbind(c(1,2),c(2,6)), na.rm=TRUE) works
#dist(rbind(c(0.5,5),colMeans(rbind(c(1,2),c(2,6)), na.rm=TRUE)))[1]
outs <- matrix(nrow=198, ncol=10)
sq_outs <- matrix(nrow=198, ncol=10)

for(i in 1:10)
{
	tempor <- which(alable==i)
	if(length(tempor)!=1)
	{
		tcol <- colMeans(parked[tempor,], na.rm=TRUE)
	}
	if(length(tempor)==1)
	{
		tcol <- parked[tempor,]
	}
	#tcol as mean vector of cluster
	for(j in 1:198)
	{
		outs[j,i]<-dist(rbind(parked[j,],tcol))[1]
		sq_outs[j,i]<-sqrt(dist(rbind(parked[j,],tcol))[1])
	}	
}

############################################################
#only with new vect
fucked <- as.matrix(cbind(num_mrna,num_mirna,good_clinic,outs))
dtrain_fuck <- xgb.DMatrix(data = fucked, label=goal_num)
recent <- xgb.cv(data = dtrain_fuck, nround=50, nthread = 4, nfold = 20, 
			metrics=list("error"),max.depth =40, eta = 1, 
			objective = "binary:logistic")

#only with mrna
damm <- as.matrix(cbind(num_mrna,good_clinic,outs))
dtrain_damm <- xgb.DMatrix(data = damm, label=goal_num)
recent_dam <- xgb.cv(data = dtrain_damm, nround=20, nthread = 5, nfold = 20, 
			metrics=list("error"),max.depth =40, eta = 1, 
			objective = "binary:logistic")

#only with mirna
damm2 <- as.matrix(cbind(num_mirna,good_clinic,outs))
dtrain_damm2 <- xgb.DMatrix(data = damm2, label=goal_num)
recent_dam2 <- xgb.cv(data = dtrain_damm2, nround=20, nthread = 5, nfold = 20, 
			metrics=list("error"),max.depth =40, eta = 1, 
			objective = "binary:logistic")

#without patch
damm3 <- as.matrix(cbind(num_mrna))
dtrain_damm3 <- xgb.DMatrix(data = damm3, label=goal_num)
recent_dam3 <- xgb.cv(data = dtrain_damm3, nround=20, nthread = 5, nfold = 20, 
			metrics=list("error"),max.depth =40, eta = 1, 
			objective = "binary:logistic")

#with clinic data



#############################################################

# USEFUL TACTICS
# which function
# paste-> as.character
# 
# factor term := mrna_t[1,]
mrnafile = rbind( mrna_t[1,],fin_mrna)
mirnafile = rbind( micro_rna_t[1,],fin_mirna)

write.csv(mrnafile, "C:/Users/son/Desktop/2015CCA/mrna_fin.csv",sep="\t")

######

data(agaricus.train, package='xgboost')

#Both dataset are list with two items, a sparse matrix and labels
#(labels = outcome column which will be learned).
#Each column of the sparse Matrix is a feature in one hot encoding format.
trainz <- agaricus.train

bstz <- xgboost(data = trainz$data, label = trainz$label, max.depth = 2,
               eta = 1, nthread = 5, nround = 10,objective = "binary:logistic")

#agaricus.test$data@Dimnames[[2]] represents the column names of the sparse matrix.
xgb.plot.tree(agaricus.train$data@Dimnames[[2]], model = bstz)

