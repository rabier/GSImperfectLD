
##### Useful functions for computing the accuracy
#############################################################
#############################################################

WrightFisherPopulation<-function(n,nbgeneration,nloci,founder,genomeDefFixedMapWithQTL){
##  Generate individuals according to Wright Fisher
## Be careful at the early beginning we only have two founders
## input : n number of individuals to generate
### nbgeneration : number of generations
#### nloci: number of loci
## founder : matrix of size 2 X nloci, containing genomes of the two founders


genomesOvertime<-array(0,dim=c(nbgeneration,n,nloci))
#generate genomes of n individuals that are descendents of the 2 founders
# i.e. it is the first generation of descendents
# genomesOvertime[1,i,] genome of individual i , for generation 1
#recombination according to Haldane
for (i in 1:n) { 
	gamete <- hypredRecombine(genomeDefFixedMapWithQTL,
		genomeA = founder[1,],
		genomeB = founder[2,],
		mutate = FALSE,
		block = FALSE)
	genomesOvertime[1,i,]<-as.matrix(gamete)

}


# generate generations from 2 to nbgeneration by random mating, and recombination according to Haldane
#genomesOvertime[g,i,j] denotes the genome of individu i at locus j at generation g
for (g in 2:nbgeneration){

	for (i in 1:n) { 

		indparents<-sample(seq(1,n), 2, replace = FALSE)

		gamete <- hypredRecombine(genomeDefFixedMapWithQTL,
			genomeA = genomesOvertime[g-1,indparents[1],],
			genomeB = genomesOvertime[g-1,indparents[2],],
			mutate = FALSE,
			block = FALSE)

		genomesOvertime[g,i,]<-as.matrix(gamete)

	}
} 


return(genomesOvertime)


}



##################################################

generateFullSibs<-function(nbFullSibs,nloci,genomeParent1, genomeParent2,genomeDefFixedMapWithQTL){
### generates genomes of  Full Sibs, descendents of Parent1 and 2
#### input : 
## genomeParent1 matrix of size 1 X nbloci, representing genome of Parent 1
## genomeParent2 matrix of size 1 X nbloci, representing genome of Parent 2


genomeFullSibs<-array(0,dim=c(nbFullSibs,nloci))
#genomeFullSibs[i,j] will denote genome of sib i at loci j

for (i in 1:nbFullSibs) { 
gamete1 <- hypredRecombine(genomeDefFixedMapWithQTL,
genomeA = genomeParent1,
genomeB = genomeParent2,
mutate = FALSE,
block = FALSE)
genomeFullSibs[i,]<-as.matrix(gamete1)
}

return(genomeFullSibs)

}

##################################################################################

generateTestIndividuals<-function(genomeLastButOneGeneration,n,nbTest,nloci,genomeDefFixedMapWithQTL) {

## function that generates Test individuals, they are extra individual from the random mating
## input :  
## genomeLastButOneGeneration : matrix of size n X nloci
## containing the genomes of the n potential parents
##nbTest: number of Test individuals
### n : number of regular indiv obtained by random mating. We will pick up parents of Tests from those individuals at last but one generation

#### output: genomesTest is the genome of Test individuals
genomesTest<-array(0,dim=c(nbTest,nloci))

### Test individuals are nbTest extra individual from the random mating, so go back to generation nbgeneration-1
for (i in 1:nbTest) { 

indparents<-sample(seq(1,n), 2, replace = FALSE)

gamete <- hypredRecombine(genomeDefFixedMapWithQTL,
genomeA = genomeLastButOneGeneration[indparents[1],],
genomeB = genomeLastButOneGeneration[indparents[2],],
mutate = FALSE,
block = FALSE)

genomesTest[i,]<-as.matrix(gamete)

}


return(genomesTest)

}













########################################################################################

generatePhenotypes<-function(genomeAtQTL,QTLeffectsMatrix,VarEnv){
### output : generate phenotypes
### input : 
### genomeAtQTL , matrix of size n X nbQTLS
### QTLeffectsMatrix : matrix of size nbQTLs x 1, containing the QTLs effects
### VarEnv : environmental variance

n<-dim(genomeAtQTL)[1]
Y<-rep(0,n)
Y[1:n]<-genomeAtQTL %*% QTLeffectsMatrix + sqrt(VarEnv)*rnorm(n)

return(Y)

}


##################################################################################

mergeGenomes<-function(genomeWF, genomeFullSibs, genomesTest){

##### input
#### genomeWF: matrix of size n X nbmarkers (it includes the QTL)
#### genomesTest: matrix of size nbTest X nbmarkers (it includes the QTL)
#### genomeFullSibs: matrix of size nbFullSibs X nbmarkers (it includes the QTL)

### output
### allGenomes, matrix containing all the genomes

n<-dim(genomeWF)[1]
nloci<-dim(genomeWF)[2]
nbFullSibs<-dim(genomeFullSibs)[1]
nbTest<-dim(genomesTest)[1]
allGenomes<-array(0,dim=c(n+nbFullSibs+nbTest,nloci))

allGenomes[1:n,1:nloci] =  genomeWF[,]

#add full Sibs
begin<-n+1
end<-n+nbFullSibs
allGenomes[begin:end, 1:nloci] =  genomeFullSibs[,]

## add TEST
beginBis<-end+1
endBis<- end + nbTest
allGenomes[beginBis:endBis, 1:nloci] =  genomesTest[,]

return(allGenomes)

}


###########################################################################

removeSNPNotUnique<-function(indSNPDup, Genome){
##### identical SNPs
##### along the chromosome are filtered out, keeping only the first occurence of that
##### SNP on the chromosome
### input : 
### indSNPDup are indices of duplicated SNPs (first occurence not listed)
### Genome : : matrix nbind x nbmarker

if (length(indSNPDup)>0){
	GenomesWithoutDup<-array(0,dim=dim(Genome[,-indSNPDup])) 
	GenomesWithoutDup<-Genome[,-indSNPDup]
} else {
	GenomesWithoutDup<-array(0,dim=dim(Genome[,])) 
	GenomesWithoutDup<-Genome[,]
}

return(GenomesWithoutDup)


}



######################################################


keepSNPSPolymorphic<-function(markNotPoly, Genome){

### remove from the Genome matrix the non polymporphic SNPs
### input : 
### Genome : : matrix nbind x nbmarker
### markNotPoly are indices of not Polymorphic SNPS

if (length(markNotPoly)>0){
	GenomesWithoutPoly<-array(0,dim=dim(Genome[,-markNotPoly])) 
	GenomesWithoutPoly<-Genome[,-markNotPoly]
	} else {
	GenomesWithoutPoly<-array(0,dim=dim(Genome)) 
	GenomesWithoutPoly<-Genome[,]
}

return(GenomesWithoutPoly)

}



#############################################################################

calcCenteredGenomes<-function(Genomes){
### Center Genomes
### input 
###### Genomes : matrix nbind x nbmarker

nbInd<-dim(Genomes)[1]

GenomeCentered<-array(0,dim(Genomes))
meansByColumn<-colMeans(Genomes, na.rm = TRUE, dims = 1)
mm=matrix(rep(meansByColumn, nbInd, each=T), ncol = dim(Genomes)[2], byrow = T)
GenomeCentered <- Genomes - mm

return(GenomeCentered)

}


##########################################################################



#########################################################

computeTheoreticalQuadraticError<-function(GenomeTestCentered,LDCorrectedWithQTLTRN,QTLeffectsMatrix,genomeAtQTLTestCentered,DesignTerm){
### compute the Theoretical Quadratic Error
### Be careful , we assume that environmental variance is equal to one

#### input
### GenomeTestCentered : Genome of Test individual that have been centered, matrix of size nbTest x NbRemainingMarkers
### genomeAtQTLTestCentered :  Genome of Test individual at the QTL  (centered), matrix of size nbTest x nbQTLs
### QTLeffectsMatrix : matrix of size nbQTLs x 1, containing the QTLs effects
### LDCorrectedWithQTLTRN : matrix of size NbRemainingMarkers X nbQTLs, containing the corrected LD for the TRN individuals
### DesignTerm :  design term, it is equal to NewMe/nbTRN

IndivBias<-GenomeTestCentered %*% LDCorrectedWithQTLTRN %*% QTLeffectsMatrix - genomeAtQTLTestCentered %*% QTLeffectsMatrix
nbTest<-dim(GenomeTestCentered)[1]
SecondTerm<-sum(IndivBias^2)/nbTest
TheoreticalQuadraticError<-DesignTerm + SecondTerm + 1

return(TheoreticalQuadraticError)

}


computeTheoreticalAccuracy<-function(GenomeTestCentered, genomeAtQTLTestCentered, LDCorrectedWithQTLTRN, QTLeffectsMatrix, DesignTerm, VarEnv){
### compute the Theoretical Accuracy

#### input
### GenomeTestCentered : Genome of Test individual that have been centered, matrix of size nbTest x NbRemainingMarkers
### genomeAtQTLTestCentered :  Genome of Test individual at the QTL  (centered), matrix of size nbTest x nbQTLs
### QTLeffectsMatrix : matrix of size nbQTLs x 1, containing the QTLs effects
### LDCorrectedWithQTLTRN : matrix of size NbRemainingMarkers X nbQTLs, containing the corrected LD for the TRN individuals
### DesignTerm :  design term, it is equal to NewMe/nbTRN


T1<-GenomeTestCentered %*% LDCorrectedWithQTLTRN %*% QTLeffectsMatrix 
T2<-genomeAtQTLTestCentered %*% QTLeffectsMatrix
NumTheoreticalAccuracy<-mean(T1*T2)

###Denominator

##Compute Variance of GenomeTestCentered %*% LDCorrectedWithQTLTRN %*% QTLeffectsMatrix
VarTerm<-t(QTLeffectsMatrix)%*%t(LDCorrectedWithQTLTRN) %*% var(GenomeTestCentered) %*% LDCorrectedWithQTLTRN %*% QTLeffectsMatrix
## Compute Genetic Variance in Test
VarGenetTest<-var(genomeAtQTLTestCentered %*% QTLeffectsMatrix)
DenomTheoreticalAccuracy<-sqrt( VarGenetTest + VarEnv) *  sqrt(VarEnv*DesignTerm + VarTerm)
TheoreticalAccuracy<-NumTheoreticalAccuracy/DenomTheoreticalAccuracy
return(TheoreticalAccuracy)


}





#######################################################################################

computeProxy <-function(heritability, Me, nbTRN){
### compute the Proxy
## input : heritability, estimated Me, number of Trainings nbTRN 

numProxy<-heritability/sqrt(1-heritability)
denomProxy <- sqrt(Me/nbTRN + heritability/(1-heritability))
Proxy <- numProxy / denomProxy

return(Proxy)

}



calcMeLiJi<-function(NbCut,NbRemainingMarkers,genome){
## compute Me according to Li and Ji
##input 
## NbCut : number of times the chromosome has to be splitted to able to perform Li and Ji
#nbRemainigMarrkers : number of markers which are unique and polymorphics
### genome : matrix of size nb individual X nb markers 
#### nb individual can be larger than nb markers

MeLiJi <- 0
for (cut in 1:NbCut){
TheBeginning<-floor( (cut-1)* NbRemainingMarkers /NbCut) + 1
TheEnd <- floor(cut*NbRemainingMarkers /NbCut)
MeLiJi <- MeLiJi + calcLiAndJi(genome[, TheBeginning:TheEnd])
}

return(MeLiJi)
}



calcLiAndJi <- function(genome){

## Compute ME corresponding to the matrix genome
##input
## genome : matrix of size nb individual X nb markers
#### nb individual has to be larger than nb markers

## compute 
MYCOR<-cor(genome)
MESVALPROP<-eigen(MYCOR)
VALABSOLUMESVALPROP<-abs(MESVALPROP$values)
Me<-sum( (VALABSOLUMESVALPROP>1) + ( VALABSOLUMESVALPROP- trunc(VALABSOLUMESVALPROP)) )

return(Me)

}




distanceBetweenRemainingSNPS<-function(nbSNPS, NbRemainingMarkers, indSNPDup, markNotPoly, step){
#compute genetic distances between remaining SNPS
## input  
#nbSNPS: total number of markers
#nbRemainigMarrkers : number of markers which are unique and polymorphics
# indSNPDUP : indices of markers that are duplicated
## markNotPoly : indices of markers that are not polymorphic (among non duplicated markers)

### index of all SNPS
allSNPS<-seq(1,nbSNPS)
### Remove indices of duplicated markers (keep first instance)
if (length(indSNPDup)>0){
	SNPWithoutDup<-allSNPS[-indSNPDup]
} else {
	SNPWithoutDup<-allSNPS
}

#### Remove indices of not polymorphic SNPS

if (length(markNotPoly)>0){
	SNPWithoutDupPoly<-SNPWithoutDup[-markNotPoly]
	} else {
	SNPWithoutDupPoly<-SNPWithoutDup
}

##SNPWithoutDupPoly are SNP which are not duplicated and polymorphics

########Compute distance between remaining SNPS
tabDistance<-array(0,dim=c(NbRemainingMarkers,NbRemainingMarkers))
for (i in 1:NbRemainingMarkers){
	for (j in 1:NbRemainingMarkers) { 
		tabDistance[i,j]<-abs(SNPWithoutDupPoly[i]-SNPWithoutDupPoly[j]) * step
	}
}
VectDist<-as.vector(tabDistance)

return(VectDist)
}



###############################################

HillWeir.est<-function(R2,Distbp,TailleEchantillon)
{
#Modèle de Hill&Weir
mod_X <-nls(R2~(((10+4*N*Distbp)/((2+4*N*Distbp)*(11+4*N*Distbp)))*(1+((3+4*N*Distbp)*(12+12*4*N*Distbp+(4*N*Distbp)*(4*N*Distbp)))/(TailleEchantillon*(2+(4*N*Distbp))*(11+(4*N*Distbp))))),start=list(N=0.0001), trace=TRUE)
O<-(summary(mod_X))
#paramètre estimé (Taille de la population)
N.est=O$parameters[1]
N.est
}

#Pour calculer la courbe de Hill&Weir
HillWeir.cal<-function(Distbp,TailleEchantillon,N){
Y<-( (10+4*N*Distbp)/ ( (2+4*N*Distbp)*(11+4*N*Distbp) ) )*( 1+( (3+4*N*Distbp)*(12+12*4*N*Distbp+(4*N*Distbp)*(4*N*Distbp)) )/( TailleEchantillon*(2+(4*N*Distbp))*(11+(4*N*Distbp)) ) )
Y}

########################################################















