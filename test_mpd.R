#set.seed(134)
sq.sz<-100 # square number
sq.ln<-sqrt(sq.sz)
samp<-matrix(rpois(sq.sz,1),nrow=sq.ln)

tmp<-1:sq.ln
rownames(samp)<-paste("p",tmp,sep="")
colnames(samp)<-paste("h",tmp,sep="")
dis<-matrix(abs(rnorm(sq.sz)),nrow=sq.ln)
dis[lower.tri(dis)]<-t(dis)[lower.tri(dis)]
diag(dis)<-0.0
rownames(dis)<-paste("h",tmp,sep="")
colnames(dis)<-paste("h",tmp,sep="")
abundance.weighted=T


N <- dim(samp)[1]
mpd <- numeric(N)
for (i in 1:N) {
  sppInSample <- names(samp[i, samp[i, ] > 0])
  if (length(sppInSample) > 1) {
    sample.dis <- dis[sppInSample, sppInSample]
    if (abundance.weighted) {
      sample.weights <- t(as.matrix(samp[i, sppInSample, 
                                         drop = FALSE])) %*% as.matrix(samp[i, sppInSample, 
                                                                            drop = FALSE])
      mpd[i] <- weighted.mean(sample.dis, sample.weights)
    }
    else {
      mpd[i] <- mean(sample.dis[lower.tri(sample.dis)])
    }
  }
  else {
    mpd[i] <- NA
  }
}

y<-ses.mpd(samp,dis,abundance.weighted = T,runs = 1000)

#CUSTOM VERSION FOR COMPARISON
nreps<-1000
x<-NULL #store mpds
Habund<-colSums(samp)
for (i in 1:dim(samp)[1]){ # go through each parasite
  rep.mpd<-0 # store sum of mpds over reps
  HR<-length(which(samp[i,]>0)) # get host richness for this parasite
  for (rep in 1:nreps){ # replicate many times
    grabbed.hosts<-names(sample(Habund,size=HR,prob=Habund,replace=T))
    rep.mpd<-rep.mpd+mean(dis[grabbed.hosts,grabbed.hosts])
  }
  x<-c(x,rep.mpd/nreps)
}
  
plot(y$mpd.obs,x)




