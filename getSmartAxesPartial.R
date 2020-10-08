## version of getSmartAxes but for a partial LD

####### function for choosing the best axes   
 
getSmartAxesPartial<-function(RankGenomeTRNCenteredWithoutQTL,ThetaTilde,svdMatTrainingsWithoutQTL,p1){

WeightedL2projections<-rep(0,RankGenomeTRNCenteredWithoutQTL)

	for (s in 1:RankGenomeTRNCenteredWithoutQTL){

		IDAxes=s

		NumCompressFirstVersion<- t(ThetaTilde) %*% t(GenomeTRNCentered[,]) %*% GenomeTRNCenteredWithoutQTL %*% t(GenomeTRNCenteredWithoutQTL[,])  %*% matrixTWithoutQTL  %*% svdMatTrainingsWithoutQTL$u[,IDAxes] %*% t(svdMatTrainingsWithoutQTL$u[,IDAxes]) %*% GenomeTRNCentered[,] %*% ThetaTilde

		NumCompressTerm<- NumCompressFirstVersion/nbTRN

		WeightedL2projections[s] <- NumCompressTerm

	}

IndWeightedL2projectionsOrdered<-order(WeightedL2projections,decreasing=TRUE)
 
 
## we are looking for axes that explain p1 of the weighted projection
ratioWeightedL2=0
i<-0
while(ratioWeightedL2 < p1) {
i<-i+1
ratioWeightedL2<-sum(WeightedL2projections[IndWeightedL2projectionsOrdered[1:i]])/sum(WeightedL2projections)

}

IndWeightedL2LastOrdered<-i

return(IndWeightedL2projectionsOrdered[1:IndWeightedL2LastOrdered])


} ### end of function getSmartAxesPartial


##################################################################
## function for Computing compressed accuracy with smart axes
 
ComputeCompressedAccuracySmartAxes<-function(IDAxes,nbTRN,ThetaTilde,svdMatTrainingsWithoutQTL,GenomeTRNCenteredWithoutQTL,matrixTWithoutQTL,lambda,GenomeTRNCentered){
## It computes a compressed accuracy as a function of the rank of the two matrices

## ID Axes are the smart axes

NbAxes<-length(IDAxes)

### Compute Var Term for Partial LD

LDCorrectedWithQTLTRNSmart<-t(GenomeTRNCenteredWithoutQTL[,]) %*% matrixTWithoutQTL[,]  %*% svdMatTrainingsWithoutQTL$u[,IDAxes] %*% t(svdMatTrainingsWithoutQTL$u[,IDAxes]) %*% GenomeTRNCentered

VarCompressTerm<-t(ThetaTilde)%*% t(LDCorrectedWithQTLTRNSmart) %*% var(GenomeTRNCenteredWithoutQTL) %*% LDCorrectedWithQTLTRNSmart %*% ThetaTilde


### Compute design Term for Partial LD

RRWithoutQTLTRNCompressed<-GenomeTRNCenteredWithoutQTL[,] %*% t(GenomeTRNCenteredWithoutQTL[,]) %*% matrixTWithoutQTL[,] %*% svdMatTrainingsWithoutQTL$u[,IDAxes] %*% t(svdMatTrainingsWithoutQTL$u[,IDAxes])
DesignCompressTerm<-sum(RRWithoutQTLTRNCompressed^2)/nbTRN


######################### Focus on numerator for partial LD

NumCompressFirstVersion<- t(ThetaTilde) %*% t(GenomeTRNCentered[,]) %*% GenomeTRNCenteredWithoutQTL %*% t(GenomeTRNCenteredWithoutQTL[,])  %*% matrixTWithoutQTL  %*% svdMatTrainingsWithoutQTL$u[,IDAxes] %*% t(svdMatTrainingsWithoutQTL$u[,IDAxes]) %*% GenomeTRNCentered[,] %*%ThetaTilde

NumCompressTerm<- NumCompressFirstVersion/nbTRN

############################################

MyInterestingTerms = list(NumCompressTerm, DesignCompressTerm, VarCompressTerm)

return(MyInterestingTerms)

}


