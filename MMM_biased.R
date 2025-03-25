library(ape)
library(diversitree)
library(phangorn)
library(hisse)
library(geiger)
library(TreeSim)

#Heuer et al 2024 data
mammal_data <- read.csv("/Users/fevi4829/Documents/Bee_popgen/BiSSE/mammal_data_reduced.csv", row.names = 2)
mammal_tree_unpruned <- read.tree("/Users/fevi4829/Documents/Bee_popgen/BiSSE/mammal_phylogeny.phy")

name.check(mammal_tree_unpruned, mammal_data)

force.ultrametric<-function(tree,method=c("nnls","extend")){
  method<-method[1]
  if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
                                     rooted=TRUE,trace=0)
  else if(method=="extend"){
    h<-diag(vcv(tree))
    d<-max(h)-h
    ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
               y=tree$edge[,2])
    tree$edge.length[ii]<-tree$edge.length[ii]+d
  } else 
    cat("method not recognized: returning input tree\n\n")
  tree
}

#force tree to be ultrametric
mammal_tree<-force.ultrametric(mammal_tree_unpruned) ## default method
is.ultrametric(mammal_tree)

#force tree to remove polytonies
is.binary(mammal_tree)
is.binary(multi2di(mammal_tree))
mammal_tree_binary = multi2di(mammal_tree)
is.binary(mammal_tree_binary)

#MMM mammals
behavior <-mammal_data[,"morph"]
names(behavior) <- rownames(mammal_data)

is_champion <- as.numeric(behavior == "1")
names(is_champion) <- names(behavior)
is_champion

sampling.f = c((length(which(mammal_data$morph == "0"))/nrow(mammal_data)),(length(which(mammal_data$morph == "1"))/nrow(mammal_data)))

nbModel <- make.bisse(mammal_tree_binary, is_champion, sampling.f = sampling.f)
p <- starting.point.bisse(mammal_tree_binary)
nbModelMLFit <- find.mle(nbModel, p)

estimated = round(coef(nbModelMLFit))

cnbModel <- constrain(nbModel, lambda1 ~ lambda0)
cnbModel <- constrain(cnbModel, mu1 ~ mu0)

cnbModelMLFit <- find.mle(cnbModel, p[c(-1, -3)], method="optim") 
anova(nbModelMLFit, constrained = cnbModelMLFit) 

#MCMC
prior <- make.prior.exponential(1 / (2 * (p[1] - p[3])))

mcmcRun <- mcmc(nbModel, nbModelMLFit$par, nsteps = 25000, prior = prior, w = 0.1, print.every = 1000)
mean(mcmcRun$lambda1)

col <- c("red", "blue")
profiles.plot(mcmcRun[, c("lambda0", "lambda1")], col.line = col, las = 1, legend = "topright", main = "Speciation rate")

profiles.plot(mcmcRun[, c("mu0", "mu1")], col.line = col, las = 1, legend = "topright", main = "Extinction rate")

profiles.plot(mcmcRun[, c("q01", "q10")], col.line = col, las = 1, legend = "topright", main = "Mammal tree (champion = 1)")

postSamples <- mcmcRun[c("lambda1", "mu1")]

profiles.plot(postSamples, col.line = c("red", "blue"), las = 1, legend = "topright")

postSamples$r <- with(mcmcRun, lambda1 - mu1)
postSamples$eps <- with(mcmcRun, mu1/lambda1)

profiles.plot(postSamples[, c("r", "eps")], col.line = c("red", "blue"), las = 1,  legend = "topright")