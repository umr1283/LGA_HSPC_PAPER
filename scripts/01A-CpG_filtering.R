#deter threshold for the different quality metrics [TO update]

source("scripts/utils/methyl_utils.R")

data_all<-fread("datasets/CD34_angle_119_noEmptyLocis_withConfScore_withoutChrXY.txt",header = T)

data_all<-data.frame(data_all)
rownames(data_all)<-data_all$id
head(data_all)

dim(data_all) #1709224     132
samples<-names(data_all)[str_detect(names(data_all),"CBP")]


batch<-read.csv2("datasets/batch_CD34_library_date_090420.csv",header=T,row.names = 1)

#d'abord en fct msp1c et nbNA
deterSeuilQC(data_all,metrique = "msp1c",qTestes = 1:9/40) #exclu locis < q0.125
quantile(data_all$msp1c,0.125) #6.532433e-08

deterSeuilQC(data_all,metrique = "pctNA",qTestes = 0:5,test = "brut") #seuil exclu ; no NA

mat<-as.matrix(data_all[,samples])
data_F<-data_all[data_all$msp1c>quantile(data_all$msp1c,0.125)&
                   rowSums(is.na(mat))==0,]

#puis on retire les locis full methylated
data_F$nbNonFullMethyl<-rowSums(data_F[,samples]>10,na.rm = T)

deterSeuilQC(data_F,metrique = "nbNonFullMethyl",qTestes = 0:10,test = "brut",lowerThan = F) #exclu locis avec nbNonFUllmethyl<5

data_F<-data_F[rowSums(data_F[,samples]>10)>4,]
nrow(data_F) #989522
#plus conf Score, nbMethylNonzeros dans pct0 elevÃ© :

names(data_F)
deterSeuilQC(data_F,metrique = "confidenceScore",qTestes = 1:9/10) #exclu locis < q0.2
data_F<-data_F[data_F$confidenceScore>quantile(data_F$confidenceScore,0.2),]
nrow(data_F) #791613

deterSeuilQC(data_F[data_F$pct0>0.7,],metrique = "nbMethylNonZeros",qTestes = 0:5,test = "brut") #exclu locis avec pct0>0.7 et nbMethylVraizers==0


data_F<-data_F[!(data_F$pct0>0.7&data_F$nbMethylNonZeros==0),]
nrow(data_F) #fre
locisF<-rownames(data_F)
#gain en qualitÃ©
#avant
deterQual(mat) #41% des locis avec Vrais zeros
deterQual(mat[locisF,]) #54%>64% de locis avec vrais zeros

deterQual2(mat,batch) #"PC 1  ( 16.9 % de la variance) a R2 avec Group_Complexity = 0.75 et pval = 10^ -36.1671884642798"


deterQual2(mat[locisF,],batch) #"PC 1  ( 19.4 % de la variance a R2 avec Group_Complexity = 0.72 et pval = 10^ -33.9025713104284"

deterQual2(mat[!(rownames(mat)%in%locisF),],batch)
# "PC 1  ( 10.6 % de la variance a R2 avec Group_Complexity = 0.73 et pval = 10^ -34.8942186739751"

deterQual2(mat[sample(rownames(mat),length(locisF)),],batch)
#"PC 1  ( 16.9 % de la variance a R2 avec Group_Complexity = 0.75 et pval = 10^ -36.2830174521662"
deterQual2(mat[sample(rownames(mat),length(locisF)),],batch)
#PC 1  ( 16.9 % de la variance a R2 avec Group_Complexity = 0.75 et pval = 10^ -36.0886977688336"

#visualisation de nouvelle distribution

dim(data_F) #1029401     132

plot(density(as.matrix(data_F[,samples])))
plot(density(data_all$mean))
lines(density(data_F$mean),col=2)
plot(density(data_all$sd))
lines(density(data_F$sd),col=2)
plot(density(data_all$pct0))
lines(density(data_F$pct0),col=2)

hist(data_F$RankConfidenceScore,breaks = 100)

#comparÃ© aux hasard : 
hist(data_all[sample(rownames(mat),length(locisF)),"RankConfidenceScore"],breaks = 100)
hist(data_all[sample(rownames(mat),length(locisF)),"RankConfidenceScore"],breaks = 100)
