# change directory
setwd("~/Investigacion/peptidos")
load("~/Investigacion/peptidos/frequencies.RData")

F<-rep(0,1000) #initialize the output
linemax<-2000 #total number of lines in the file
lagmax<-200 # maximum number of lines to allow for lags
residues<-84 #total number of residues per frame
#two vectors for the loops translate the lag number to the
#row number.
lag2vec<-seq(from=1,to=(linemax-lagmax),by=residues)
lag1vec<-seq(from=residues,to=linemax,by=residues)

for (res in 0:residues){
  #loop 2 if for going frame by frame. This sets the origin
  for (lag2 in 1:length(lag2vec)){
    #loop 1 is for shifting the lag for a given origin
    for (lag1 in 1:length(lag1vec)){
      #this loops lag
      if (sum(is.na(wtr[(lag2vec[lag2]+res),3:7]))==5 & lag1==1){
        #there is no water in the origin. Do nothing to F
        cuze<-0 #initialize cumulative zeros
        break
      }else{#current origin has water, compare following lag1
        cmp<-sum(wtr[(lag2vec[lag2]+res),3:7]==
                   wtr[(lag1vec[lag1]+lag2vec[lag2]+res),3:7],
                 na.rm=T)
        if (cmp==0){#there is no water in the next lag
          #test further next lag.
          cmp2<-sum(wtr[(lag2vec[lag2]+res),3:7]==
                    wtr[(lag1vec[lag1+1]+lag2vec[lag2]+res),3:7],
                    na.rm=T)
          if (cmp2==0){#not present in the next frame
            #do nothing to F
            cuze<-c(cuze,0)#accumulate a zero
          }else{#water only skipped one frame
            #consider the water contiguous
            cmp<-1 #this resets the water as present
          }
        }else{#there was water in next lag
          F[lag1]<-F[lag1]+cmp
        }
      }
      #print(c("lag1",lag1,"lag2",lag2,"res",res,"row1",lag2vec[lag2]+res,"row2",lag1vec[lag1]+lag2vec[lag2]+res))
      #print(c(wtr[(lag2vec[lag2]+res),2],wtr[(lag1vec[lag1]+lag2vec[lag2]+res),2]))
      #print(c(lag2vec[lag2]+res,sum(wtr[(lag2vec[lag2]+res),3:7]==wtr[(lag1vec[lag1]+lag2vec[lag2]+res),3:7],na.rm=T)))
      #print(rbind(c(lag1vec[lag1]+lag2vec[lag2]+res),wtr[(lag1vec[lag1]+lag2vec[lag2]+res)),3:7]),c(lag1vec[lag1]+lag2vec[lag2]+res),wtr[(lag1vec[lag1]+lag2vec[lag2]+res)),3:7])))
    }
  }
}


