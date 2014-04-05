library(igraph)
source("./lib/sgmviaIP.R")

A.graph<- erdos.renyi.game(n=20,p.or.m=0.5)
A<- get.adjacency(A.graph)
randpern<- c(1:3, sample(17))

B<- A[randpern,randpern]
matching.exact.it.list <- sgmViaIP(A,B, m=3)
