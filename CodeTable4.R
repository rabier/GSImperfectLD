### code regarding Table 4 of the manuscript

### This code will fill a Table, called summary Table
### each column gives an estimated accuracy according to a given method
### each row is one simulated data set
## finally, all results will be in a geenrated file called MyTableResults

#####################################################

#### EmpiricalAccuracyWithoutQTL : empirical Accuracy
## TheoreticalAccuracyWithoutQTL :  rhocheck(X*,X*new,Beta*)
## TheoreticalAccuracyWithoutQTLGroupLASSO :  rhohat(X*hat, BetaGPLASSO*hat)
## TheoreticalAccuracyWithoutQTLLASSO  :  rhohat(X*hat, BetaLASSO*hat)
## TheoreticalAccuracyWithoutQTLADLASSO : rhohat(X*hat, BetaADLASSO*hat)
## TheoreticalAccuracyScandinavianTRNADLASSO : rhohatPerfectLD(BetaADLASSO*hat)
## TheoreticalAccuracyScandinavianTESTADLASSO : rhocheckPerfectLD(BetaADLASSO*hat)
## TheoreticalAccuracyScandinavianTRNLASSO : rhohatPerfectLD(BetaLASSO*hat)
## TheoreticalAccuracyScandinavianTESTLASSO :  rhocheckPerfectLD(BetaLASSO*hat)
## EmpiricalAccuracyWithoutQTLCompressedLASSO : empirical Accuracy using model
### based on the improved predictor associated with LASSO
## EmpiricalAccuracyWithoutQTLCompressedADLASSO : empirical Accuracy using model
### based on the improved predictor associated with AdaptiveLASSO
 

##############################################################
 
 
library(hypred)
library(rrBLUP)
library(glmnet)
library(gglasso)
library(parcor)
library(Rlab)
 

### To be replaced by the location of your current directory
currentDirectory<-"/Users/charley/Desktop/PARTIALLD/ForOurUsers"
setwd(currentDirectory)

source("getAccuracy.R") 
source("functionsForAccuracy.R")
source("functionsForMultiFounders.R")  
source("getSmartAxesPartial.R")


chroLength<-4
nbSNPS<-4000
#nbSNPS<-8000

#chroLength<-6
#nbSNPS<-6000
#nbSNPS<-12000
  

## we place QTLs on the first Morgan , some of them do not have any effects on the trait,
## but it is useful for creating loci with hypred package
QTLlocations<-seq(2*chroLength/nbSNPS,1,2*chroLength/nbSNPS)
length(QTLlocations)
 
QTLeffects<-rep(0,length(QTLlocations))
 
# 2 large QTLs will be located at 3 and 80cM

if (nbSNPS/chroLength==1000) { 
	QTLeffects[15]<-2
	QTLeffects[400]<--4	

} else if (nbSNPS/chroLength ==2000) { 

	QTLeffects[30]<-2
 	QTLeffects[800]<--4  
}


 
### n is the number of Training individuals that are not full sibs
n<-500
#n<-1000

## VarEnv is the environmental variance
VarEnv<-4

### nbFullSibs is the number of FullSibs considered in the Training set (at least >=2)
nbFullSibs<-2 

### nbTRN is the number of Training individuals
nbTRN<-n+nbFullSibs


## nbTest is the number of Test individuals
nbTest<-100

nbgeneration<-100

set.seed(289)
nbSim<-100 

summaryTable = data.frame(IDsim=1:nbSim, TheoreticalAccuracyWithoutQTL=rep(0,nbSim), EmpiricalAccuracyWithoutQTL=rep(0,nbSim), TheoreticalAccuracyWithoutQTLGroupLASSO=rep(0,nbSim),
 TheoreticalAccuracyWithoutQTLLASSO=rep(0,nbSim), TheoreticalAccuracyWithoutQTLADLASSO=rep(0,nbSim), TheoreticalAccuracyScandinavianTRNADLASSO=rep(0,nbSim), TheoreticalAccuracyScandinavianTESTADLASSO=rep(0,nbSim), TheoreticalAccuracyScandinavianTRNLASSO=rep(0,nbSim), TheoreticalAccuracyScandinavianTESTLASSO=rep(0,nbSim),
 EmpiricalAccuracyWithoutQTLCompressedLASSO=rep(0,nbSim), EmpiricalAccuracyWithoutQTLCompressedADLASSO=rep(0,nbSim))
 
 
## 1 chromosome of length chroLength and nbSNPS SNPs
genomeDef <- hypredGenome(1, chroLength, nbSNPS)

###create a genetic map with markers at fixed location
step<-chroLength/nbSNPS 

change_map <- seq(step,chroLength,step)
genomeDefFixedMap <- hypredNewMap(genomeDef,change_map)

###find the QTL index which corresponds to QTL locations
QTLindex<-QTLlocations/step
### QTL has to be on markers, they will be defined explicitely with hypredNewQTL

genomeDefFixedMapWithQTL <- hypredNewQTL(genomeDefFixedMap,
new.id.add = round(QTLindex),
new.id.dom = NULL, ## default
new.eff.add = QTLeffects,
new.eff.dom = NULL ## default
)
   


## count the number of loci
nloci<-genomeDefFixedMapWithQTL@num.snp.chr + genomeDefFixedMapWithQTL@num.add.qtl.chr
genomesOvertime<-array(0,dim=c(nbgeneration,n,nloci))
##genomesOvertime[g,i,j] denotes the genome of individu i at locus j at generation g


### we consider 8 founders for the population
nbfounder<-8
genomesOvertime<-WrightFisherPopulationMultiFounders(n,nbgeneration,nloci,genomeDefFixedMapWithQTL,nbfounder)
#### generate genomes of the n individuals according to Wright Fisher
## Be careful, at the first generation , it is a population of n individuals who are descendents of nbfounder founders 


## generate full sibs
genomeFullSibs<-array(0,dim=c(nbFullSibs,nloci))
#genomeFullSibs[i,j] will denote genome of sib i at loci j

#Choose 2 parents randomly among the n individuals
indparentsFullSibs<-sample(seq(1,n), 2, replace = FALSE)

genomeFullSibs<-generateFullSibs(nbFullSibs,nloci,genomesOvertime[nbgeneration,indparentsFullSibs[1],], genomesOvertime[nbgeneration,indparentsFullSibs[2],],genomeDefFixedMapWithQTL)

####### Generate Test individuals, ie phenotypes, and genome matrix,

for (idSample in 1:nbSim) {
 

#### genomesTest is the genome of Test individuals
genomesTest<-array(0,dim=c(nbTest,nloci))

### Test individuals are nbTest extra individual from the random mating, so go back to generation nbgeneration-1
genomesTest<-generateTestIndividuals(genomesOvertime[nbgeneration-1,,],n,nbTest,nloci,genomeDefFixedMapWithQTL)
 
#################################
### generate phenotypes

QTLeffectsMatrix<-matrix(genomeDefFixedMapWithQTL@add.and.dom.eff$add,ncol=1)
Y<-rep(0,n)
Y[1:n]<-generatePhenotypes(genomesOvertime[nbgeneration, , genomeDefFixedMapWithQTL@pos.add.qtl$ID],QTLeffectsMatrix,VarEnv)
 
#generate phenotypes of full Sibs
YFullSibs<-rep(0,nbFullSibs)
YFullSibs<-generatePhenotypes(genomeFullSibs[,genomeDefFixedMapWithQTL@pos.add.qtl$ID],QTLeffectsMatrix,VarEnv)

YTRN<-rep(0,nbTRN)
YTRN<-c(Y,YFullSibs)
## YTRN contains phenotypes of the TRN
## we center the phenotypes of TRN
YTRNCentered<-YTRN-mean(YTRN)

#generate phenotypes of Tests
YTest<-rep(0,nbTest)
YTest<-generatePhenotypes(genomesTest[,genomeDefFixedMapWithQTL@pos.add.qtl$ID],QTLeffectsMatrix,VarEnv)
## we center the phenotypes of Test
YTestCentered<-YTest -  mean(YTest)

################################################
allGenomes<-array(0,dim=c(n+nbFullSibs+nbTest,nloci))
### Fill all genomes, by respectively n indiv, Full Sibs, and Tests
#recall that we focus only on last generation
allGenomes<-mergeGenomes(genomesOvertime[nbgeneration, ,],genomeFullSibs,genomesTest)
 
allGenomesPreFiltered<-allGenomes


##########################################

### Remove SNP not unique
### i.e. identical SNPs along the chromosome are filtered out, keeping only the first occurence of that
##### SNP on the chromosome
indSNPDup<-which(duplicated(t(allGenomesPreFiltered)))
allGenomesWithoutDup<-removeSNPNotUnique(indSNPDup, allGenomesPreFiltered)

### keep only SNPS with polymorphisms (since former SNPs are only markers, I should have called them markers)

sumCol<-apply(allGenomesWithoutDup[,],2,sum)
markNotPoly<-which(sumCol==dim(allGenomesWithoutDup)[1] | sumCol==0)
allGenomesWithoutDupPoly<-keepSNPSPolymorphic(markNotPoly, allGenomesWithoutDup)
NbRemainingMarkers<-dim(allGenomesWithoutDupPoly)[2]

##########Center TRN genomes and Test genomes and Extra TRN

GenomeTRNCentered<-array(0,dim=c(nbTRN,dim(allGenomesWithoutDupPoly)[2]))
#genomeTRNCentered[i,j] denotes the genome of indiv i of TRN at locus j , but centered
## Recall that in TRN we consider that we have the n indiv and Full Sibs
GenomeTRNCentered<-calcCenteredGenomes(allGenomesWithoutDupPoly[1:nbTRN,])

indFirstTest<-n+nbFullSibs+1
indLastTest<-n+nbFullSibs+nbTest
GenomeTestCentered<-calcCenteredGenomes(allGenomesWithoutDupPoly[indFirstTest:indLastTest,])
#genomeTestCentered[i,j] denotes the genome of indiv i of Test at locus j , but centered
 
 
#####################Compute the Lambda, ie the tuning parameter . We will use REML.

data <- data.frame(phenoTrain=YTRNCentered, gid=seq(1,nbTRN))
rownames(GenomeTRNCentered)<-seq(1,nbTRN)
myK<- GenomeTRNCentered %*% t(GenomeTRNCentered)
#myK is my Kinship matrix
 
 
##### we  use REML 
### use kin.blup to obtain variance components
u<-kin.blup(data=data, geno="gid", pheno="phenoTrain", K=myK)
### lambda REML
lambda<-u$Ve/u$Vg
 

################################### Evaluation of the prediction model on Test individuals
##### compute Ridge regression
matrixT <- solve(myK + lambda*diag(nbTRN))

 
################## We need to obtain the SNPs matrix without the putative QTLs

allGenomesPreFilteredWithoutQTL<-allGenomes[ , -genomeDefFixedMapWithQTL@pos.add.qtl$ID]

### Remove SNP not unique
### i.e. identical SNPs along the chromosome are filtered out, keeping only the first occurence of that
##### SNP on the chromosome
indSNPDupWithoutQTL<-which(duplicated(t(allGenomesPreFilteredWithoutQTL)))
allGenomesWithoutDupWithoutQTL<-removeSNPNotUnique(indSNPDupWithoutQTL, allGenomesPreFilteredWithoutQTL)

### keep only SNPS with polymorphisms

sumColWithoutQTL<-apply(allGenomesWithoutDupWithoutQTL[,],2,sum)
markNotPolyWithoutQTL<-which(sumColWithoutQTL==dim(allGenomesWithoutDupWithoutQTL)[1] | sumColWithoutQTL==0)
allGenomesWithoutDupPolyWithoutQTL<-keepSNPSPolymorphic(markNotPolyWithoutQTL, allGenomesWithoutDupWithoutQTL)
NbRemainingMarkersWithoutQTL<-dim(allGenomesWithoutDupPolyWithoutQTL)[2]

########## center TRN genomes WithoutQTL and Test genomes WithoutQTL

GenomeTRNCenteredWithoutQTL<-array(0,dim=c(nbTRN,dim(allGenomesWithoutDupPolyWithoutQTL)[2]))
#genomeTRNCentered[i,j] denotes the genome of indiv i of TRN at locus j , but centered
## Recall that in TRN we consider that we have the n indiv and Full Sibs
GenomeTRNCenteredWithoutQTL<-calcCenteredGenomes(allGenomesWithoutDupPolyWithoutQTL[1:nbTRN,])

indFirstTest<-n+nbFullSibs+1
indLastTest<-n+nbFullSibs+nbTest
GenomeTestCenteredWithoutQTL<-calcCenteredGenomes(allGenomesWithoutDupPolyWithoutQTL[indFirstTest:indLastTest,])


#####################Compute the Lambda REML Based Without QTL

dataWithoutQTL <- data.frame(phenoTrain=YTRNCentered, gid=seq(1,nbTRN))
rownames(GenomeTRNCenteredWithoutQTL)<-seq(1,nbTRN)
myKWithoutQTL<- GenomeTRNCenteredWithoutQTL %*% t(GenomeTRNCenteredWithoutQTL)

##### we use REML 
### use kin.blup to obtain variance components
u<-kin.blup(data=dataWithoutQTL, geno="gid", pheno="phenoTrain", K=myKWithoutQTL)
### lambda REML
lambdaWithoutQTL<-u$Ve/u$Vg

##### compute Ridge regression
matrixTWithoutQTL <- solve(myKWithoutQTL + lambdaWithoutQTL*diag(nbTRN))

#################################### Empirical Accuracy Without QTL, ie in imperfect LD

### Prediction for Test individual
TestPredictionWithoutQTL<-GenomeTestCenteredWithoutQTL[,] %*% t(GenomeTRNCenteredWithoutQTL[,])  %*% matrixTWithoutQTL %*% YTRNCentered

#Empirical Accuracy
EmpiricalAccuracyWithoutQTL<-cor(YTestCentered , TestPredictionWithoutQTL)


summaryTable$EmpiricalAccuracyWithoutQTL[idSample] <- EmpiricalAccuracyWithoutQTL
 

#################### Imperfect LD proxies #############################################

###### So, Let us compute rhocheck(X*,X*new,Beta*), that we call TheoreticalAccuracyWithoutQTL
 

RRWithoutQTL<-GenomeTestCenteredWithoutQTL[,] %*% t(GenomeTRNCenteredWithoutQTL[,]) %*% matrixTWithoutQTL[,]
DesignTermWithoutQTL<-sum(RRWithoutQTL^2)/nbTest

## Linkage disequilibrium corrected for relatedness , between SNPS and QTLs, for the Training
genomeAtQTLTRNCentered<-calcCenteredGenomes(allGenomes[1:nbTRN, genomeDefFixedMapWithQTL@pos.add.qtl$ID])
LDCorrectedWithQTLTRN<-t(GenomeTRNCenteredWithoutQTL[,]) %*% matrixTWithoutQTL[,]  %*% genomeAtQTLTRNCentered

### Center the genome at QTL of Test 
genomeAtQTLTestCentered<-calcCenteredGenomes(genomesTest[,genomeDefFixedMapWithQTL@pos.add.qtl$ID])
##note that genomeTest is the genome of Test before filtering

###Compute Theoretical Accuracy
TheoreticalAccuracyWithoutQTL<-computeTheoreticalAccuracy(GenomeTestCenteredWithoutQTL, genomeAtQTLTestCentered, LDCorrectedWithQTLTRN, QTLeffectsMatrix, DesignTermWithoutQTL,VarEnv)

summaryTable$TheoreticalAccuracyWithoutQTL[idSample]<-TheoreticalAccuracyWithoutQTL

 

########## Focus on Adaptive LASSO
###### So, Let us compute rhohat(X*hat, BetaADLASSO*hat), that we call TheoreticalAccuracyWithoutQTLADLASSO

LDCorrectedWithQTLTRNBis<-t(GenomeTRNCenteredWithoutQTL[,]) %*% matrixTWithoutQTL[,]  %*% GenomeTRNCentered

MyAdapt<-adalasso(GenomeTRNCentered,YTRNCentered,k=10)
MyAdaptLassoCoeff<-MyAdapt$coefficients.adalasso

RRWithoutQTLTRN<-GenomeTRNCenteredWithoutQTL[,] %*% t(GenomeTRNCenteredWithoutQTL[,]) %*% matrixTWithoutQTL[,]
DesignTermWithoutQTLTRN<-sum(RRWithoutQTLTRN^2)/nbTRN

TheoreticalAccuracyWithoutQTLADLASSO<-computeTheoreticalAccuracy(GenomeTRNCenteredWithoutQTL, GenomeTRNCentered, LDCorrectedWithQTLTRNBis, MyAdaptLassoCoeff, DesignTermWithoutQTLTRN,VarEnv)

summaryTable$TheoreticalAccuracyWithoutQTLADLASSO[idSample]<-TheoreticalAccuracyWithoutQTLADLASSO

########## Focus on  LASSO
###### So, Let us compute rhohat(X*hat, BetaLASSO*hat), that we call TheoreticalAccuracyWithoutQTLLASSO 
 

MyLassoCV<-cv.glmnet(GenomeTRNCentered, YTRNCentered, standardize =FALSE, intercept=FALSE, alpha=1)

MyLassoCoeff<-coef(MyLassoCV, MyLassoCV$lambda.1se)

TheoreticalAccuracyWithoutQTLLASSO<-computeTheoreticalAccuracy(GenomeTRNCenteredWithoutQTL, GenomeTRNCentered, LDCorrectedWithQTLTRNBis, MyLassoCoeff[2:(NbRemainingMarkers+1)]
, DesignTermWithoutQTLTRN,VarEnv)

summaryTable$TheoreticalAccuracyWithoutQTLLASSO[idSample]<-TheoreticalAccuracyWithoutQTLLASSO

########## Focus on  Group LASSO
###### So, Let us compute  rhohat(X*hat, BetaGPLASSO*hat), that we call TheoreticalAccuracyWithoutQTLGroupLASSO
 
NbGroupsForGroupLASSO<-floor(NbRemainingMarkers/10)

MyGroups= rep(1:NbGroupsForGroupLASSO,each=10)

if  ( (NbRemainingMarkers%%10) != 0 ) {

MyGroupToAdd<-rep((NbGroupsForGroupLASSO+1), NbRemainingMarkers%%10)

MyGroups<-c(MyGroups,MyGroupToAdd)
}

cvGroupLASSO <- cv.gglasso(GenomeTRNCentered, YTRNCentered, MyGroups, loss="ls",
pred.loss="L2", lambda.factor=0.05, nfolds=5)

MyGroupLassoCoeff<-coef(cvGroupLASSO,s=cvGroupLASSO$lambda.1se)


TheoreticalAccuracyWithoutQTLGroupLASSO<-computeTheoreticalAccuracy(GenomeTRNCenteredWithoutQTL, GenomeTRNCentered, LDCorrectedWithQTLTRNBis, MyGroupLassoCoeff[2:(NbRemainingMarkers+1)]
, DesignTermWithoutQTLTRN,VarEnv)

summaryTable$TheoreticalAccuracyWithoutQTLGroupLASSO[idSample]<-TheoreticalAccuracyWithoutQTLGroupLASSO

#################### Perfect LD proxies #############################################

########## Focus on  Adaptive LASSO
###### So, Let us compute  rhohatPerfectLD(BetaADLASSO*hat), that we call  TheoreticalAccuracyScandinavianTRNADLASSO 
 
LDCorrectedWithoutQTLTRN<-t(GenomeTRNCenteredWithoutQTL[,]) %*% matrixTWithoutQTL[,]  %*% GenomeTRNCenteredWithoutQTL[,]

MyAdaptTRNWithout<-adalasso(GenomeTRNCenteredWithoutQTL,YTRNCentered,k=10)
MyAdaptLassoCoeffTRNWithout<-MyAdaptTRNWithout$coefficients.adalasso

TheoreticalAccuracyScandinavianTRNADLASSO<-computeTheoreticalAccuracy(GenomeTRNCenteredWithoutQTL, GenomeTRNCenteredWithoutQTL, LDCorrectedWithoutQTLTRN, MyAdaptLassoCoeffTRNWithout, DesignTermWithoutQTLTRN,VarEnv)


summaryTable$TheoreticalAccuracyScandinavianTRNADLASSO[idSample]<-TheoreticalAccuracyScandinavianTRNADLASSO

###### So, Let us compute  rhocheckPerfectLD(BetaADLASSO*hat), that we call TheoreticalAccuracyScandinavianTESTADLASSO 

 
TheoreticalAccuracyScandinavianTESTADLASSO<-computeTheoreticalAccuracy(GenomeTestCenteredWithoutQTL, GenomeTestCenteredWithoutQTL, LDCorrectedWithoutQTLTRN, MyAdaptLassoCoeffTRNWithout, DesignTermWithoutQTL,VarEnv)

summaryTable$TheoreticalAccuracyScandinavianTESTADLASSO[idSample]<-TheoreticalAccuracyScandinavianTESTADLASSO

################ Focus on LASSO
###### So, Let us compute  rhohatPerfectLD(BetaLASSO*hat), that we call TheoreticalAccuracyScandinavianTRNLASSO 
  
MyLassoCVForScand<-cv.glmnet(GenomeTRNCenteredWithoutQTL, YTRNCentered, standardize =FALSE, intercept=FALSE, alpha=1)

MyLassoCoeffForScand<-coef(MyLassoCVForScand, MyLassoCVForScand$lambda.1se)


TheoreticalAccuracyScandinavianTRNLASSO<-computeTheoreticalAccuracy(GenomeTRNCenteredWithoutQTL, GenomeTRNCenteredWithoutQTL, LDCorrectedWithoutQTLTRN, MyLassoCoeffForScand[2:( dim(LDCorrectedWithoutQTLTRN)[2]+1)], DesignTermWithoutQTLTRN,VarEnv)


summaryTable$TheoreticalAccuracyScandinavianTRNLASSO[idSample]<-TheoreticalAccuracyScandinavianTRNLASSO

######  So, Let us compute  rhocheckPerfectLD(BetaLASSO*hat), that we call TheoreticalAccuracyScandinavianTESTLASSO
 
TheoreticalAccuracyScandinavianTESTLASSO<-computeTheoreticalAccuracy(GenomeTestCenteredWithoutQTL, GenomeTestCenteredWithoutQTL, LDCorrectedWithoutQTLTRN,  MyLassoCoeffForScand[2:( dim(LDCorrectedWithoutQTLTRN)[2]+1)], DesignTermWithoutQTL,VarEnv)

summaryTable$TheoreticalAccuracyScandinavianTESTLASSO[idSample]<-TheoreticalAccuracyScandinavianTESTLASSO


#########################################################################################
############# From here, we focus on our modified predictor


#####   SVD now

svdMatTrainingsWithoutQTL<-svd(GenomeTRNCenteredWithoutQTL)

RankGenomeTRNCenteredWithoutQTL <- rankMatrix(GenomeTRNCenteredWithoutQTL,method = "qrLINPACK", tol=1e-20)


#######################  Modified predictor relying on Adaptive LASSO  

### we choose the axes smartly 

MyP1<-c(0.7,0.8,0.9, 0.925, 0.95, 0.975, 0.99)

MyCompressedAccuracyAdapt<-0
IDAxesAdapt<-0

for (compteur in 1:length(MyP1)){

IDAxesAdaptNew=getSmartAxesPartial(RankGenomeTRNCenteredWithoutQTL[1],MyAdaptLassoCoeff,svdMatTrainingsWithoutQTL,MyP1[compteur])

InterestingTermsAdapt<-ComputeCompressedAccuracySmartAxes(IDAxesAdaptNew,nbTRN,MyAdaptLassoCoeff,svdMatTrainingsWithoutQTL,GenomeTRNCenteredWithoutQTL,matrixTWithoutQTL,lambda,GenomeTRNCentered)

#VarGenetTRN<-var(genomeAtQTLTRNCentered %*% QTLeffectsMatrix)

VarGenetTRN<-var(GenomeTRNCentered %*% MyAdaptLassoCoeff)

DenomCompressedAdapt<-sqrt( VarGenetTRN + VarEnv) *  sqrt(VarEnv*InterestingTermsAdapt[[2]] + InterestingTermsAdapt[[3]])

MyCompressedAccuracyAdaptNew<-InterestingTermsAdapt[[1]]/DenomCompressedAdapt

 

	if ( MyCompressedAccuracyAdaptNew > MyCompressedAccuracyAdapt ) {
 
		
		IDAxesAdapt<-IDAxesAdaptNew	
		## on conserve les nouveaux axes pour l emp acc reduced calculee apres	

		MyCompressedAccuracyAdapt <- MyCompressedAccuracyAdaptNew

		summaryTable$MyCompressedAccuracyAdapt[idSample]<-MyCompressedAccuracyAdapt
	 }


}

  
### empirical Accuracy using model
### based on the improved predictor associated with AdaptiveLASSO

# Prediction for Test individual
TestPredictionWithoutQTLCompressedAdaptLASSO<-GenomeTestCenteredWithoutQTL[,] %*% t(GenomeTRNCenteredWithoutQTL[,])  %*% matrixTWithoutQTL %*% 
svdMatTrainingsWithoutQTL$u[,IDAxesAdapt] %*% t(svdMatTrainingsWithoutQTL$u[,IDAxesAdapt]) %*% YTRNCentered

#Empirical Accuracy
EmpiricalAccuracyWithoutQTLCompressedADLASSO<-cor(YTestCentered , TestPredictionWithoutQTLCompressedAdaptLASSO)


summaryTable$EmpiricalAccuracyWithoutQTLCompressedADLASSO[idSample]<-EmpiricalAccuracyWithoutQTLCompressedADLASSO
 

#########  Modified predictor relying on  LASSO now

MyP1<-c(0.7,0.8,0.9, 0.925, 0.95, 0.975, 0.99)
MyCompressedAccuracyLASSO<-0
IDAxesLASSO<-0

for (compteur in 1:length(MyP1)){

IDAxesLASSONew=getSmartAxesPartial(RankGenomeTRNCenteredWithoutQTL[1],MyLassoCoeff[2:(NbRemainingMarkers+1)],svdMatTrainingsWithoutQTL,MyP1[compteur])

InterestingTermsLASSO<-ComputeCompressedAccuracySmartAxes(IDAxesLASSONew,nbTRN,MyLassoCoeff[2:(NbRemainingMarkers+1)],svdMatTrainingsWithoutQTL,GenomeTRNCenteredWithoutQTL,matrixTWithoutQTL,lambda,GenomeTRNCentered)


VarGenetTRN<-var(GenomeTRNCentered %*% MyLassoCoeff[2:(NbRemainingMarkers+1)])

DenomCompressedLASSO<-sqrt( VarGenetTRN + VarEnv) *  sqrt(VarEnv*InterestingTermsLASSO[[2]] + InterestingTermsLASSO[[3]])

MyCompressedAccuracyLASSONew <- InterestingTermsLASSO[[1]]/DenomCompressedLASSO
 

	if ( MyCompressedAccuracyLASSONew > MyCompressedAccuracyLASSO ) {

		
		IDAxesLASSO<-IDAxesLASSONew	 

		MyCompressedAccuracyLASSO <- MyCompressedAccuracyLASSONew

		summaryTable$MyCompressedAccuracyLASSO[idSample]<-MyCompressedAccuracyLASSO
	 }


}


### empirical Accuracy using model
### based on the improved predictor associated with  LASSO

# Prediction for Test individual
TestPredictionWithoutQTLCompressedLASSO<-GenomeTestCenteredWithoutQTL[,] %*% t(GenomeTRNCenteredWithoutQTL[,])  %*% matrixTWithoutQTL %*% 
svdMatTrainingsWithoutQTL$u[,IDAxesLASSO] %*% t(svdMatTrainingsWithoutQTL$u[,IDAxesLASSO]) %*% YTRNCentered

#Empirical Accuracy
EmpiricalAccuracyWithoutQTLCompressedLASSO<-cor(YTestCentered , TestPredictionWithoutQTLCompressedLASSO)

summaryTable$EmpiricalAccuracyWithoutQTLCompressedLASSO[idSample]<-EmpiricalAccuracyWithoutQTLCompressedLASSO


######### Write Table in File
write.table(summaryTable, file="MyTableResults", quote=F, sep="\t", row.names=F, col.names=T)


}



