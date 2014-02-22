library(igraph)
library(R.matlab)
library(clue)

celeg.mats<- readMat("./data/CElegans_3subgraph_pairs.mat")
mat.it<-1:3
nmc <- 400
mat.size<- c(0,0,0)
m.vals.list<-list()
correct.match <- list()
  random.match <- list()
for (it in mat.it)  {
  n <- nrow(celeg.mats[[it]])
  mat.size[it]<- n
  m.vals <- seq(1,n-4,3)
  m.len  <- length(m.vals)
  m.vals.list[[it]] <- m.vals
  
  correct.match <- c(correct.match,list(array(0,dim=c(m.len,nmc))))
  random.match <- c(random.match,list(array(0,dim=c(m.len,nmc))))
  
  rownames(correct.match[[it]]) <- m.vals
  for (m.it in 1:m.len){
    m <- m.vals[m.it]
    init.mat <-matrix(1/(n-m),n-m,n-m)
    for (mc.it in 1:nmc) {
      seed.perm <- sample(n, replace=FALSE)
      
      A <- celeg.mats[[it]][seed.perm,seed.perm]
      B <- celeg.mats[[it+3]][seed.perm,seed.perm]
      
      matching.it.list <- sgm(A,B, m=m, iteration=20, start=init.mat)
      P.hat   <-  matching.it.list[[2]]
      match.hat <- matching.it.list[[1]]
      correct.match[[it]][m.it,mc.it] <- sum(match.hat[,2]==((m+1):n))/(n-m)
     # rownames(correct.match[[it]])=m.it
      print(paste0("correct.match = ",sum(match.hat[,2]==((m+1):n))/(n-m)))
      random.match[[it]][m.it,mc.it] <- 1/(n-m)
    }
  }
}

#melt(correct.match,id.vars=c(1),measure.vars)

library(reshape2)
library(arrayhelpers)
library(ggplot2)

source("~/Documents/projects/JOFC-GraphMatch/lib/simulation_math_util_fn.R")

for (it in mat.it){
correct.match.lf.it <- array2df(correct.match[[it]],levels=list(num.seed=TRUE,mc.it=NA) ,label.x="true.match.frac")
nrows <- nrow(correct.match.lf.it)
correct.match.lf.it$matrix.size <- as.integer(mat.size[it])
if (it==1) {
  correct.match.lf <- correct.match.lf.it
} else{
  correct.match.lf <- rbind(correct.match.lf,correct.match.lf.it)
}
}

corr.summ<-summarySE(data=correct.match.lf,measurevar="true.match.frac",groupvar=c("num.seed","matrix.size"))

corr.summ$num.seed<-as.numeric(levels(corr.summ$num.seed))[corr.summ$num.seed]
corr.summ$chance <- 1/(corr.summ$matrix.size - corr.summ$num.seed)
corr.summ$chance[corr.summ$chance<0]<- NA
corr.summ$matrix.size <- as.factor (corr.summ$matrix.size)


ggplot(corr.summ, aes(x=num.seed, y=true.match.frac, colour=matrix.size)) + 
  geom_errorbar(aes(ymin=true.match.frac-ci, ymax=true.match.frac+ci)) +
  geom_line(aes(x=num.seed,y=chance),linetype=2)+
  #geom_point(size=2)+geom_line(size=1.2)+
  theme_minimal()+
  #theme(text=element_text(size=22)) +
  labs(title="SGM for Celegans Connectome",x=expression(m),
       y="True Match Fraction") +
       #y=((expression(delta^{(m)})))) +
  scale_x_continuous(breaks=seq(0,120,20)) +
  scale_y_continuous(breaks=seq(0,1,.1)) 
#+
 # guides(colour=guide_legend( title ="graph.size",title.hjust=1,title.vjust=-1,label.hjust=1))
ggsave("SGM.for.Celegans.connectome-igraph.sgm.pdf")

corr.summ.frac.seed<- corr.summ
graph.sizes<- as.numeric(levels(corr.summ$matrix.size))[corr.summ$matrix.size]
corr.summ.frac.seed$num.seed<- corr.summ.frac.seed$num.seed/graph.sizes


fig.2<- ggplot(corr.summ.frac.seed, aes(x=num.seed, y=true.match.frac, colour=matrix.size)) + 
  geom_errorbar(aes(ymin=true.match.frac-ci, ymax=true.match.frac+ci)) +
  geom_line(aes(x=num.seed,y=chance),linetype=2)+
  #geom_point(size=2)+geom_line(size=1.2)+
  theme_minimal()+
  #theme(text=element_text(size=22)) +
  labs(title="SGM for Celegans Connectome",x="Fraction of Seed Vertices",
       y="True Match Fraction") +
  #y=((expression(delta^{(m)})))) +
  scale_x_continuous(breaks=seq(0,1,.1)) +
  scale_y_continuous(breaks=seq(0,1,.1)) 

fig.2
ggsave("SGM.for.Celegans.connectome-igraph.sgm-2.pdf")