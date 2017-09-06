# change directory
setwd("~/Investigacion/peptidos")
load("~/Investigacion/peptidos/frequencies.RData")

F<-rep(0,1000) #initialize the output
linemax<-count #total number of lines in the file
lagmax<-200 # maximum number of lines to allow for lags
residues<-84 #total number of residues per frame
#two vectors for the loops translate the lag number to the
#row number.
lag2vec<-seq(from=1,to=(linemax-lagmax),by=residues)
lag1vec<-seq(from=residues,to=linemax,by=residues)
#loop 2 if for going frame by frame. This sets the origin
for (lag2 in 1:length(lag2vec)){
  #loop 1 is for shifting the lag for a given origin
  for (lag1 in 1:length(lag1vec)){
    #this loops for residues within a frame
    for (res in 0:residues){
      if (sum(is.na(wtr[(res),3:7]))==10){
        #there is no water. Add a zero to F. Do nothing
      }else{
        F[lag1]<-F[lag1]+
          sum(wtr[(lag2vec[lag2]+res),3:7]==wtr[(lag1vec[lag1]+lag2vec[lag2]+res),3:7],
              na.rm=T)
      }
#print(c("lag1",lag1,"lag2",lag2,"res",res,"row1",lag2vec[lag2]+res,"row2",lag1vec[lag1]+lag2vec[lag2]+res))
#print(c(wtr[(lag2vec[lag2]+res),2],wtr[(lag1vec[lag1]+lag2vec[lag2]+res),2]))
#print(c(lag2vec[lag2]+res,sum(wtr[(lag2vec[lag2]+res),3:7]==wtr[(lag1vec[lag1]+lag2vec[lag2]+res),3:7],na.rm=T)))
#print(rbind(c(lag1vec[lag1]+lag2vec[lag2]+res),wtr[(lag1vec[lag1]+lag2vec[lag2]+res)),3:7]),c(lag1vec[lag1]+lag2vec[lag2]+res),wtr[(lag1vec[lag1]+lag2vec[lag2]+res)),3:7])))
    }
  }
}


