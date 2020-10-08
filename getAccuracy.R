
getAccuracy<-
function (REML=FALSE, REMOVEQTL=FALSE, nbSNPS, chroLength, nbQTLs, QTLlocations, QTLeffects, heritability, n, nbFullSibs, nbTest, nbgeneration, VarEnv) {

#Function that gives the accuracy realted to only one sample

#INPUT

# REML : Boolean, TRUE if we want to compute the tuning parameter by REML 
# vs FALSE if we want a tuning parameter heritability based
# REMOVEQTL : Boolean, TRUE if we want to remove the QTL from the design matrix
# vs FALSE if we do not want to remove QTLs
# nbSNPS : number of SNPS
# chroLength : length of the chromosome in Morgan
# nbQTLs : number of QTLs
# QTLlocations : vector containing location of the QTLs, they have to be on markers
# QTLeffects : vector containing QTL effects
# heritability : heritability of the considered Trait
# n : number of individuals present in the TRN and that are not full sibs
# nbFullSibs : number of Full Sibs that will be considered in the TRN set
# nbTest : number of Test individuals
# nbgeneration : number of generation during which the  population will evolve by random mating

#OUTPUT

# List containing in the following order
# TheoreticalAccuracy : Theoretical Accuracy according to our formula
# EmpiricalAccuracy : Empirical accuracy
# EmpiricalBias : Empirical bias
# TheoreticalQuadraticError : Theoretical Quadratic Error
# EmpiricalQuadraticError : Empirical Quadratic Error 
# LDCorrectedWithQTLTRN : Corrected LD between markers and QTLs, for the Training
# NewProxy : Our New Proxy
# proxyMe1 : Proxy based on Me1 (Ne based on Hill and Weirr)
# proxyMe2 : Proxy based on Me2 (Ne based on Hill and Weirr)
# proxyMe3 : Proxy based on Me3 (Ne based on Hill and Weirr)
# ProxyLiJi : Li and Ji's Proxy
# NewMe : Me corresponding to our method
# Me1 : Me1 (Ne based on Hill and Weirr)
# Me2 : Me2 (Ne based on Hill and Weirr)
# Me3 :  Me2 (Ne based on Hill and Weirr)
# MeLiJi : MeLiJi corresponding to Li and Ji

## 1 chromosome of length chroLength and nbSNPS SNPs
genomeDef <- hypredGenome(1, chroLength, nbSNPS)

###create a genetic map with markers at fixed location
step<-chroLength/nbSNPS
change_map <- seq(step,chroLength,step)
genomeDefFixedMap <- hypredNewMap(genomeDef,change_map)

###find the QTL index which corresponds to QTL locations
QTLindex<-QTLlocations/step
### QTL has to be on markers, they will be defined explicitely with hypredNewQTL

##we assume no dominance
###warning:  have to use round since bug de R QTLindex are not considered as integer by R
genomeDefFixedMapWithQTL <- hypredNewQTL(genomeDefFixedMap,
new.id.add = round(QTLindex),
new.id.dom = NULL, ## default
new.eff.add = QTLeffects,
new.eff.dom = NULL ## default
)

## produce two haploid founder line genomes
founder <- hypredFounder(genomeDefFixedMapWithQTL,1)
### these 2 lines are completely gentically different

## count the number of loci
nloci<-genomeDefFixedMapWithQTL@num.snp.chr + genomeDefFixedMapWithQTL@num.add.qtl.chr
genomesOvertime<-array(0,dim=c(nbgeneration,n,nloci))
##genomesOvertime[g,i,j] denotes the genome of individu i at locus j at generation g
#### generate genomes of the n individuals according to Wright Fisher
## Be careful, at the first generation , it a population of n individuals who are descendents of the 2 founders
genomesOvertime<-WrightFisherPopulation(n,nbgeneration,nloci,founder,genomeDefFixedMapWithQTL)

#######################################################################
## generate full sibs

genomeFullSibs<-array(0,dim=c(nbFullSibs,nloci))
#genomeFullSibs[i,j] will denote genome of sib i at loci j

#Choose 2 parents randomly among the n individuals
indparentsFullSibs<-sample(seq(1,n), 2, replace = FALSE)
genomeFullSibs<-array(0,dim=c(nbFullSibs,nloci))
genomeFullSibs<-generateFullSibs(nbFullSibs,nloci,genomesOvertime[nbgeneration,indparentsFullSibs[1],], genomesOvertime[nbgeneration,indparentsFullSibs[2],],genomeDefFixedMapWithQTL)


####### Generate Test individuals

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


nbTRN<-n+nbFullSibs
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


#######################################
## check if the user asked to remove the QTLs from the SNP matrix

if  (REMOVEQTL==TRUE) {
### Remove QTLs
allGenomesPreFiltered<-allGenomes[ , -genomeDefFixedMapWithQTL@pos.add.qtl$ID]
} else {
###Do not remove QTLs
allGenomesPreFiltered<-allGenomes
}


##########################################

### Remove SNP not unique
### i.e. identical SNPs along the chromosome are filtered out, keeping only the first occurence of that
##### SNP on the chromosome
indSNPDup<-which(duplicated(t(allGenomesPreFiltered)))
allGenomesWithoutDup<-removeSNPNotUnique(indSNPDup, allGenomesPreFiltered)

### keep only SNPS with polymorphisms

sumCol<-apply(allGenomesWithoutDup[,],2,sum)
markNotPoly<-which(sumCol==dim(allGenomesWithoutDup)[1] | sumCol==0)
allGenomesWithoutDupPoly<-keepSNPSPolymorphic(markNotPoly, allGenomesWithoutDup)
NbRemainingMarkers<-dim(allGenomesWithoutDupPoly)[2]

##########Center TRN genomes and Test genomes

GenomeTRNCentered<-array(0,dim=c(nbTRN,dim(allGenomesWithoutDupPoly)[2]))
#genomeTRNCentered[i,j] denotes the genome of indiv i of TRN at locus j , but centered
## Recall that in TRN we consider that we have the n indiv and Full Sibs
GenomeTRNCentered<-calcCenteredGenomes(allGenomesWithoutDupPoly[1:nbTRN,])


indFirstTest<-n+nbFullSibs+1
indLastTest<-n+nbFullSibs+nbTest
GenomeTestCentered<-calcCenteredGenomes(allGenomesWithoutDupPoly[indFirstTest:indLastTest,])
#genomeTestCentered[i,j] denotes the genome of indiv i of Test at locus j , but centered

#####################Compute the Lambda, ie the tuning parameter either by REML or heritability based

data <- data.frame(phenoTrain=YTRNCentered, gid=seq(1,nbTRN))
rownames(GenomeTRNCentered)<-seq(1,nbTRN)
myK<- GenomeTRNCentered %*% t(GenomeTRNCentered)
#myK is my Kinship matrix

if (REML==TRUE){
##### we want to USE REML 
### use kin.blup to obtain variance components
u<-kin.blup(data=data, geno="gid", pheno="phenoTrain", K=myK)
### lambda REML
lambda<-u$Ve/u$Vg
} else {
### The tuning parameter is heritability based
lambda <- ( ( 1 - heritability) / (heritability) ) * dim(allGenomesWithoutDupPoly)[2] * mean( diag(var(allGenomesWithoutDupPoly[1:nbTRN,1:dim(allGenomesWithoutDupPoly)[2]]) ) )
}

################################### Evaluation of the prediction model on Test individuals
##### compute Ridge regression
matrixT <- solve(myK + lambda*diag(nbTRN))


######################################### Empirical Bias, Empirical Accuracy and Empirical Quadratic Error
#### I will use REML here , fix it later
### Prediction for Test individual
TestPrediction<-GenomeTestCentered[,] %*% t(GenomeTRNCentered[,])  %*% matrixT %*% YTRNCentered

##Empirical Bias
EmpiricalBias= mean(YTestCentered - TestPrediction)

#Empirical Accuracy
EmpiricalAccuracy<-cor(YTestCentered , TestPrediction)

### Empirical Quadratic Error
EmpiricalQuadraticError= mean( (YTestCentered- TestPrediction)^2 ) 

###################################### Theoretical Bias, Accuracy and Quadratic Error

RR<-GenomeTestCentered[,] %*% t(GenomeTRNCentered[,]) %*% matrixT[,]
DesignTerm<-sum(RR^2)/nbTest

## Linkage disequilibrium corrected for relatedness , between SNPS and QTLs, for the Training
genomeAtQTLTRNCentered<-calcCenteredGenomes(allGenomes[1:nbTRN, genomeDefFixedMapWithQTL@pos.add.qtl$ID])
LDCorrectedWithQTLTRN<-t(GenomeTRNCentered[,]) %*% matrixT[,]  %*% genomeAtQTLTRNCentered

### Center the genome at QTL of Test 
genomeAtQTLTestCentered<-calcCenteredGenomes(genomesTest[,genomeDefFixedMapWithQTL@pos.add.qtl$ID])

#### Compute Theoretical Bias
TheoreticalBias<-mean(GenomeTestCentered %*% LDCorrectedWithQTLTRN %*% QTLeffectsMatrix - genomeAtQTLTestCentered %*% QTLeffectsMatrix)

#### Compute Theoretical Quadratic Error
### Be careful , we assume that environmental variance is equal to one
TheoreticalQuadraticError<-computeTheoreticalQuadraticError(GenomeTestCentered,LDCorrectedWithQTLTRN,QTLeffectsMatrix,genomeAtQTLTestCentered,DesignTerm)

###Compute Theoretical Accuracy
TheoreticalAccuracy<-computeTheoreticalAccuracy(GenomeTestCentered, genomeAtQTLTestCentered, LDCorrectedWithQTLTRN, QTLeffectsMatrix, DesignTerm,VarEnv)

##################################################    Proxies

### Our New Proxy
NewMe<-DesignTerm*nbTRN
NewProxy<-computeProxy(heritability, NewMe, nbTRN)

### Proxy using Li and Ji for estimating Me

NbCut <- floor( dim(allGenomesWithoutDupPoly)[2]/dim(allGenomesWithoutDupPoly)[1] ) + 1
MeLiJi<-calcMeLiJi(NbCut,NbRemainingMarkers,allGenomesWithoutDupPoly)
ProxyLiJi<-computeProxy(heritability, MeLiJi, nbTRN)

##############################################################
###### Hill and Weirr

#### Compute r2 among all markers

LDUsualTRN<-array(0,dim=c(NbRemainingMarkers,NbRemainingMarkers))
LDUsualTRN<-cor(GenomeTRNCentered[,], GenomeTRNCentered[,])

VectLDUsualTRN<-as.vector(LDUsualTRN)
VectLDUsualTRNSquare<-VectLDUsualTRN^2
VectDist<-distanceBetweenRemainingSNPS(nbSNPS, NbRemainingMarkers, indSNPDup, markNotPoly, step)
NeHillWeir<-HillWeir.est(VectLDUsualTRNSquare,VectDist, nbTRN)
###NeHillWeirr is the effective size


### Compare 3 differents ways of estimating Me using NeHillWeirr

Me1<-2*NeHillWeir*chroLength/(log(4*NeHillWeir*chroLength))
Me2<-2*NeHillWeir*chroLength/(log(2*NeHillWeir*chroLength))
Me3<-2*NeHillWeir*chroLength/(log(NeHillWeir*chroLength))

proxyHillWeirMe1<-computeProxy(heritability, Me1, nbTRN)
proxyHillWeirMe2<-computeProxy(heritability, Me2, nbTRN)
proxyHillWeirMe3<-computeProxy(heritability, Me3, nbTRN)

return(list(TheoreticalAccuracy=TheoreticalAccuracy, EmpiricalAccuracy=EmpiricalAccuracy, TheoreticalBias=TheoreticalBias,
 EmpiricalBias=EmpiricalBias, TheoreticalQuadraticError=TheoreticalQuadraticError,
 EmpiricalQuadraticError=EmpiricalQuadraticError, LDCorrectedWithQTLTRN=LDCorrectedWithQTLTRN, NewProxy=NewProxy,
proxyMe1=proxyHillWeirMe1, proxyMe2=proxyHillWeirMe2, proxyMe3=proxyHillWeirMe3, ProxyLiJi=ProxyLiJi, 
	 NewMe=NewMe, Me1=Me1, Me2=Me2, Me3=Me3, MeLiJi=MeLiJi))

}

























