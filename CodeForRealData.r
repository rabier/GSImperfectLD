# Code regarding Table 7 of manuscript

### This code will fill a Table, called summary Table
## and analyze rice real data from Spindel et al.
### these real data are available from dryad at https://doi.org/10.5061/dryad.7369p 
# you need to consider the folder random subsets once you have downloaded Spindel' s data

## in the code below, we can  choose between different density of markers
# either 448, or 781 , 1553 , 3076 
### each column gives an estimated accuracy according to a given method
### each row is one data set ( a different draw of training and test individuals)
## finally, all results will be in a generated file called  MyTableREALDATAInperfectLD...

##########################
#### EmpiricalAccuracyWithoutQTL : empirical Accuracy 
## TheoreticalAccuracyWithoutQTLLASSO  :  rhohat(X*hat, BetaLASSO*hat)
## TheoreticalAccuracyWithoutQTLADLASSO : rhohat(X*hat, BetaADLASSO*hat)
## TheoreticalAccuracyScandinavianTRN : rhohatPerfectLD(BetaADLASSO*hat)
## TheoreticalAccuracyScandinavianTEST : rhocheckPerfectLD(BetaADLASSO*hat)
## BernardoProxyH2Based : Proxy Lian and Bernardo
## OurProxyH2Based: Proxy Lian and Bernardo adapted by Rabier Mangin
#########################


### on this tutorial, we will consider 448 markers

library(rrBLUP)
library(gglasso)
library(parcor)

### To be replaced by the location of your current directory
currentDirectory<-"/Users/charley/Desktop/PARTIALLD/ForOurUsers/"

setwd(currentDirectory)
source("functionsForAccuracy.R")

### choose the number of markers you want to use for your Test set ( on line 172, you ll have also to 
### do the change appropriately)
MyTab448<-read.csv("Spindel/448/geno_random_sub_448_SNPs_1.csv", header=TRUE)
#dim(MyTab448)
# 332 449

### we can also choose another set with 448 markers, cf the folder Spindel/448/
#MyTab448<-read.csv("Spindel/448/geno_random_sub_448_SNPs_3.csv", header=TRUE)


#MyTab781<-read.csv("Spindel/781/geno_random_sub_781_SNPs_3.csv", header=TRUE)
#dim(MyTab781)
# 332 782

#MyTab1553<-read.csv("Spindel/1553/geno_random_sub_1553_SNPs_3.csv", header=TRUE)

#dim(MyTab1553)
#332 1554

#MyTab3076<-read.csv("Spindel/3076/geno_random_sub_3076_SNPs_10.csv", header=TRUE)

#dim(MyTab3076)
#332 3077

 
 
 

### we will use this data set on TRN set, for finding QTLs. In other words, it is our dense map
MyTabAll<-read.csv("Spindel/MET_crfilt_.90_outliers_removed_for_RRBlup_line_corrected.csv", header=TRUE)
#dim(MyTabAll)
# 332 73148
 
Tab2012<-read.csv("Spindel/corrected_RYT2012DS_plotdata_by_GHID.csv", header=TRUE,sep=",", colClasses = "character", na.strings = "NA")
#dim(Tab2012)
# 712  12

####  Remove column 1 and 2 since they represent season and year 
TabPheno<-Tab2012[,-(1:2)]

Ind<-which(TabPheno[,1] %in% MyTabAll[,1])

TabPheno2<-TabPheno[Ind,]

IndLDG<- which(TabPheno2$LDG == 0)

TabPheno2LDG0<-TabPheno2[IndLDG,]

IndOrdered<-order(TabPheno2LDG0[,1])

TabPheno3<-TabPheno2LDG0[IndOrdered,]

TabPheno3$PLTHGT <- as.numeric(TabPheno3$PLTHGT)
TabPheno3$FL <- as.numeric(TabPheno3$FL)
TabPheno3$Plot_Yld <- as.numeric(TabPheno3$Plot_Yld)
TabPheno3$MAT <- as.numeric(TabPheno3$MAT)
TabPheno3$TILLER <- as.numeric(TabPheno3$TILLER)
TabPheno3$LDG<- as.numeric(TabPheno3$LDG)
TabPheno3$Plot_Yld<- as.numeric(TabPheno3$Plot_Yld)
TabPheno3$Kg_ha<- as.numeric(TabPheno3$Kg_ha)

TabPheno4<-data.frame(GHID=unique(TabPheno3$GHID), PLTHGT=tapply(TabPheno3$PLTHGT, TabPheno3$GHID, mean),
FL=tapply(TabPheno3$FL, TabPheno3$GHID, mean), Plot_Yld=tapply(TabPheno3$Plot_Yld, TabPheno3$GHID, mean),
MAT=tapply(TabPheno3$MAT, TabPheno3$GHID, mean), TILLER=tapply(TabPheno3$TILLER, TabPheno3$GHID, mean),
LDG=tapply(TabPheno3$LDG, TabPheno3$GHID, mean), Plot_Yld=tapply(TabPheno3$Plot_Yld, TabPheno3$GHID, mean),
Kg_ha=tapply(TabPheno3$Kg_ha, TabPheno3$GHID, mean))

### we compute the average on replicates for each indivduals  
####################################################################

Ind2<-which(MyTabAll[,1] %in% TabPheno4[,1])
TabSNP2<-MyTabAll[Ind2,]
IndOrderedSNP<-order(TabSNP2[,1])
TabSNP3<-TabSNP2[IndOrderedSNP,]

### TABSNP4 is the same table but without ID of individuals
TabSNP4<-TabSNP3[,-1]
TabSNP5<-data.matrix(TabSNP4)
nbSNPS<-dim(TabSNP5)[2]
###########################################################

### Remove SNP not unique
### i.e. identical SNPs along the chromosome are filtered out, keeping only the first occurence of that
##### SNP on the chromosome
indSNPDup<-which(duplicated(t(TabSNP5)))
allGenomesWithoutDup<-removeSNPNotUnique(indSNPDup, TabSNP5)

### keep only SNPS with polymorphisms
sumCol<-apply(allGenomesWithoutDup[,],2,sum)
markNotPoly<-which(sumCol==dim(allGenomesWithoutDup)[1] | sumCol==0)
allGenomesWithoutDupPoly<-keepSNPSPolymorphic(markNotPoly, allGenomesWithoutDup)
NbRemainingMarkers<-dim(allGenomesWithoutDupPoly)[2]


#####################Compute the Lambda, ie the tuning parameter by REML , with QTL

#nbTRN is the numer of training inidviduals
nbTRN<-round(dim(allGenomesWithoutDupPoly)[1]*0.8)
#nbTestis the numer of test inidviduals
nbTest<-dim(allGenomesWithoutDupPoly)[1]-nbTRN

nbSim<-100
summaryTable = data.frame(IDsim=1:nbSim, EmpiricalAccuracyWithoutQTL=rep(0,nbSim), 
OurProxyH2Based=rep(0,nbSim), BernardoProxyH2Based=rep(0,nbSim), 
TheoreticalAccuracyWithoutQTLLASSO=rep(0,nbSim), 
TheoreticalAccuracyWithoutQTLADLASSO=rep(0,nbSim), TheoreticalAccuracyScandinavianTRN=rep(0,nbSim), TheoreticalAccuracyScandinavianTEST=rep(0,nbSim))

IndIndiv<-seq(1,nbTRN+nbTest,1)
 

set.seed(4)

for (idSample in 1:nbSim) {


IndShuff<-sample(IndIndiv, nbTRN+nbTest, replace = FALSE)
IndTRN<-IndShuff[1:nbTRN]
IndTest<-IndShuff[(nbTRN+1):length(IndShuff)]

GenomeTRNCentered<-calcCenteredGenomes(allGenomesWithoutDupPoly[IndTRN,])

data <- data.frame(phenoTrain=TabPheno4$FL[IndTRN], gid=seq(1,nbTRN))
rownames(GenomeTRNCentered)<-seq(1,nbTRN)
myK<- GenomeTRNCentered %*% t(GenomeTRNCentered)
#myK is my Kinship matrix

##### we will USE REML 
### use kin.blup to obtain variance components
u<-kin.blup(data=data, geno="gid", pheno="phenoTrain", K=myK)
### lambda REML
lambda<-u$Ve/u$Vg
 
##################################################################################################
############################################# Without QTL  

Ind2<-which(MyTabAll[,1] %in% TabPheno4[,1])

#You need to choose the density of markers here
MyTabWithoutQTL<-MyTab448
#MyTabWithoutQTL<-MyTab3076
#MyTabWithoutQTL<-MyTab1553
#MyTabWithoutQTL<-MyTab781



Ind2WithoutQTL<-which(MyTabWithoutQTL[,1] %in% TabPheno4[,1]) 

TabSNP2WithoutQTL<-MyTabWithoutQTL[Ind2WithoutQTL,]

IndOrderedSNPWithoutQTL<-order(TabSNP2WithoutQTL[,1])

TabSNP3WithoutQTL<-TabSNP2WithoutQTL[IndOrderedSNPWithoutQTL,]

### TABSNP4 est la table de SNP sans les ID de mes individus
TabSNP4WithoutQTL<-TabSNP3WithoutQTL[,-1]

TabSNP5WithoutQTL<-data.matrix(TabSNP4WithoutQTL)
 
### It is time for filtering !!!
### Remove SNP not unique
### i.e. identical SNPs along the chromosome are filtered out, keeping only the first occurence of that
##### SNP on the chromosome
indSNPDupWithoutQTL<-which(duplicated(t(TabSNP5WithoutQTL)))
allGenomesWithoutDupWithoutQTL<-removeSNPNotUnique(indSNPDupWithoutQTL, TabSNP5WithoutQTL)

### keep only SNPS with polymorphisms

sumColWithoutQTL<-apply(allGenomesWithoutDupWithoutQTL[,],2,sum)
markNotPolyWithoutQTL<-which(sumColWithoutQTL==dim(allGenomesWithoutDupWithoutQTL)[1] | sumColWithoutQTL==0)
allGenomesWithoutDupPolyWithoutQTL<-keepSNPSPolymorphic(markNotPolyWithoutQTL, allGenomesWithoutDupWithoutQTL)
NbRemainingMarkersWithoutQTL<-dim(allGenomesWithoutDupPolyWithoutQTL)[2]

#### Let us compute the average r2. It will be useful for computing bernardo proxy

R2Bernardo<-0
consecutiveR<-0

for (i in 1:(NbRemainingMarkersWithoutQTL-1)){

consecutiveR<-cor(allGenomesWithoutDupPolyWithoutQTL[,i],allGenomesWithoutDupPolyWithoutQTL[,(i+1)])

R2Bernardo<- R2Bernardo + consecutiveR^2

}

R2Bernardo<- R2Bernardo/(NbRemainingMarkersWithoutQTL-1)

## since we are analyzing doubled haploid
R2Bernardo<-sqrt(R2Bernardo)

### end  r2 Bernardo

###########################################################
## handle phenotypes

YTRNCentered<-TabPheno4$FL[IndTRN]-mean(TabPheno4$FL[IndTRN])
YTestCentered<-TabPheno4$FL[IndTest]-mean(TabPheno4$FL[IndTest])

#####################Compute the Lambda, ie the tuning parameter  by REML 
 
GenomeTRNCenteredWithoutQTL<-calcCenteredGenomes(allGenomesWithoutDupPolyWithoutQTL[IndTRN,])
GenomeTestCenteredWithoutQTL<-calcCenteredGenomes(allGenomesWithoutDupPolyWithoutQTL[IndTest,])

dataWithoutQTL <- data.frame(phenoTrain=TabPheno4$FL[IndTRN], gid=seq(1,nbTRN))
rownames(GenomeTRNCenteredWithoutQTL)<-seq(1,nbTRN)
myKWithoutQTL<- GenomeTRNCenteredWithoutQTL %*% t(GenomeTRNCenteredWithoutQTL)
#myKWithoutQTL is my Kinship matrix

##### we  USE REML 
### use kin.blup to obtain variance components
uWithoutQTL<-kin.blup(data=dataWithoutQTL, geno="gid", pheno="phenoTrain", K=myKWithoutQTL)
### lambda REML
lambdaWithoutQTL<-uWithoutQTL$Ve/uWithoutQTL$Vg

##### compute Ridge regression
matrixTWithoutQTL <- solve(myKWithoutQTL + lambdaWithoutQTL*diag(nbTRN))

#################################### Empirical Accuracy Without QTL
#### I will use REML here 
### Prediction for Test individual
TestPredictionWithoutQTL<-GenomeTestCenteredWithoutQTL[,] %*% t(GenomeTRNCenteredWithoutQTL[,])  %*% matrixTWithoutQTL %*% YTRNCentered

#Empirical Accuracy
EmpiricalAccuracyWithoutQTL<-cor(YTestCentered , TestPredictionWithoutQTL)

summaryTable$EmpiricalAccuracyWithoutQTL[idSample]<-EmpiricalAccuracyWithoutQTL



#################### Imperfect LD proxies #############################################
 

######### Computing Lian/Bernardo proxy 

NbCut <- floor( dim(allGenomesWithoutDupPolyWithoutQTL)[2]/dim(allGenomesWithoutDupPolyWithoutQTL)[1] ) + 1

MeLiJi<-calcMeLiJi(NbCut,NbRemainingMarkersWithoutQTL,allGenomesWithoutDupPolyWithoutQTL)
 

## value according to Table 1 from Spindel
EstimatedH2<-0.4378

DenomProxyBernardo <- MeLiJi/nbTRN  + R2Bernardo*EstimatedH2/(1-EstimatedH2)

NumProxyBernardo <- sqrt(EstimatedH2/(1-EstimatedH2))


summaryTable$BernardoProxyH2Based[idSample] <- R2Bernardo*sqrt(EstimatedH2)* NumProxyBernardo / ( sqrt(DenomProxyBernardo))


### end Lian/Bernardo Proxy


#############################################################
 
LDCorrectedWithQTLTRNBis<-t(GenomeTRNCenteredWithoutQTL[,]) %*% matrixTWithoutQTL[,]  %*% GenomeTRNCentered

########## Focus on  LASSO
###### So, Let us compute rhohat(X*hat, BetaLASSO*hat), that we call TheoreticalAccuracyWithoutQTLLASSO 

MyLassoCV<-cv.glmnet(GenomeTRNCentered, t(YTRNCentered), standardize =FALSE, intercept=FALSE)


MyLassoCoeff<-coef(MyLassoCV, MyLassoCV$lambda.1se)


RRWithoutQTLTRN<-GenomeTRNCenteredWithoutQTL[,] %*% t(GenomeTRNCenteredWithoutQTL[,]) %*% matrixTWithoutQTL[,]
DesignTermWithoutQTLTRN<-sum(RRWithoutQTLTRN^2)/nbTRN

TheoreticalAccuracyWithoutQTLLASSO<-computeTheoreticalAccuracy(GenomeTRNCenteredWithoutQTL, GenomeTRNCentered, LDCorrectedWithQTLTRNBis, MyLassoCoeff[2:(NbRemainingMarkers+1)]
, DesignTermWithoutQTLTRN,u$Ve)

summaryTable$TheoreticalAccuracyWithoutQTLLASSO[idSample]<-TheoreticalAccuracyWithoutQTLLASSO


########## Focus on Adaptive LASSO
###### So, Let us compute rhohat(X*hat, BetaADLASSO*hat), that we call TheoreticalAccuracyWithoutQTLADLASSO


MyVec<-seq(1,(NbRemainingMarkers+1),1)

for (i in 1:(NbRemainingMarkers+1)){

if (MyLassoCoeff[i]!=0){
 MyVec[i] <- 1/MyLassoCoeff[i]}
else{ 
  MyVec[i]<-5000
}

}

adlasso <- cv.glmnet(GenomeTRNCentered, t(YTRNCentered), standardize =FALSE, intercept=FALSE, 
                  penalty.factor = MyVec )

MyAdaptLassoCoeff<-coef(adlasso , s=adlasso$lambda.1se)


TheoreticalAccuracyWithoutQTLADLASSO<-computeTheoreticalAccuracy(GenomeTRNCenteredWithoutQTL, GenomeTRNCentered, LDCorrectedWithQTLTRNBis,
 MyAdaptLassoCoeff[2:(NbRemainingMarkers+1)], DesignTermWithoutQTLTRN,u$Ve)


summaryTable$TheoreticalAccuracyWithoutQTLADLASSO[idSample]<-TheoreticalAccuracyWithoutQTLADLASSO


#################### Perfect LD proxies #############################################


########## Focus on  Adaptive LASSO
###### So, Let us compute  rhohatPerfectLD(BetaADLASSO*hat), that we call  TheoreticalAccuracyScandinavianTRN
 
LDCorrectedWithoutQTLTRN<-t(GenomeTRNCenteredWithoutQTL[,]) %*% matrixTWithoutQTL[,]  %*% GenomeTRNCenteredWithoutQTL[,]

MyAdaptTRNWithout<-adalasso(GenomeTRNCenteredWithoutQTL,t(YTRNCentered),k=10)
MyAdaptLassoCoeffTRNWithout<-MyAdaptTRNWithout$coefficients.adalasso

TheoreticalAccuracyScandinavianTRN<-computeTheoreticalAccuracy(GenomeTRNCenteredWithoutQTL, GenomeTRNCenteredWithoutQTL, LDCorrectedWithoutQTLTRN, MyAdaptLassoCoeffTRNWithout, DesignTermWithoutQTLTRN,uWithoutQTL$Ve)

summaryTable$TheoreticalAccuracyScandinavianTRN[idSample]<-TheoreticalAccuracyScandinavianTRN


###### So, Let us compute  rhocheckPerfectLD(BetaADLASSO*hat), that we call TheoreticalAccuracyScandinavianTEST

RRWithoutQTL<-GenomeTestCenteredWithoutQTL[,] %*% t(GenomeTRNCenteredWithoutQTL[,]) %*% matrixTWithoutQTL[,]
DesignTermWithoutQTL<-sum(RRWithoutQTL^2)/nbTest

TheoreticalAccuracyScandinavianTEST<-computeTheoreticalAccuracy(GenomeTestCenteredWithoutQTL, GenomeTestCenteredWithoutQTL, LDCorrectedWithoutQTLTRN, MyAdaptLassoCoeffTRNWithout, DesignTermWithoutQTL,uWithoutQTL$Ve)


summaryTable$TheoreticalAccuracyScandinavianTEST[idSample]<-TheoreticalAccuracyScandinavianTEST


############## Compute our simple proxy  based on Lian/Bernardo and on Rabier Mangin
 
DenomOurProxy <- DesignTermWithoutQTL + R2Bernardo*EstimatedH2/(1-EstimatedH2)
 
summaryTable$OurProxyH2Based[idSample]<-R2Bernardo*sqrt(EstimatedH2)* NumProxyBernardo / ( sqrt(DenomOurProxy))


##################################

print(idSample)

######### Write Table in File
write.table(summaryTable, file="MyTableREALDATAInperfectLD448seed1", quote=F, sep="\t", row.names=F, col.names=T)


}










