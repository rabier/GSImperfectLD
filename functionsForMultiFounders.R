


WrightFisherPopulationMultiFounders<-function(n,nbgeneration,nloci,genomeDefFixedMapWithQTL,nbfounder){
##  Generate individuals according to Wright Fisher
## Be careful at the early beginning we only have two founders
## input : n number of individuals to generate
### nbgeneration : number of generations
#### nloci: number of loci
### nbfounder: number of founders


MyFounders<-array(0,dim=c(nbfounder,nloci))

for (i in 1:nbfounder) { 
MyFounders[i,]=rbern(nloci,0.5)

}

genomesOvertime<-array(0,dim=c(nbgeneration,n,nloci))
#generate genomes of n individuals that are descendents of the 2 founders
# i.e. it is the first generation of descendents
# genomesOvertime[1,i,] genome of individual i , for generation 1
#recombination according to Haldane


for (i in 1:n) { 
	indFounders<-sample(seq(1,nbfounder), 2, replace = FALSE)

	gamete <- hypredRecombine(genomeDefFixedMapWithQTL,
		genomeA = MyFounders[indFounders[1],],
		genomeB = MyFounders[indFounders[2],],
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
