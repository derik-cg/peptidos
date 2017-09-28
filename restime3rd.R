# change directory
#setwd("~/Investigacion/peptidos")
#load("~/Investigacion/peptidos/frequencies.RData")

F<-rep(0,5000) #initialize the output
linemax<-7471610 #total number of lines in the file
lagmax<-5000 # maximum number of lines to allow for lags
residues<-84 #total number of residues per frame
#two vectors for the loops translate the lag number to the
#row number.
lag2vec<-seq(from=1,to=(linemax-lagmax),by=residues)
lag1vec<-seq(from=residues,to=linemax,by=residues)

for (res in 0:residues){
  #loop 2 if for going frame by frame. This sets the origin
  for (lag2 in 1:length(lag2vec)){
    #loop 1 is for shifting the lag for a given origin
      if (sum(is.na(wtr[(lag2vec[lag2]+res),3:7]))==5 & lag1==1){
        #there is no water in the origin. Do nothing to F
        break
      }else{#current origin has water
        #get the vector of waters
        wavc<-wtr[(lag2vec[lag2]+res),
                c(FALSE,FALSE,!is.na(wtr[(lag2vec[lag2]+res),3:7]))]
        #for each element in the vector run the lags
        nm<-length(wavc)#number of water molecules
        for (m in 1:nm){#do lags
          for (lag1 in 1:length(lag1vec)){
            #this is the loops for the lag
            #this does the next lag
          if (is.element(wtr[(lag2vec[lag2]+res),(2+m)],
                 wtr[(lag1vec[lag1]+lag2vec[lag2]+res),3:7])){
            #accumulate the current lag
            F[lag1]<-F[lag1]+1
          }else{#no water in next lag
            #check one further lag
            cmp2<-is.element(wtr[(lag2vec[lag2]+res),(2+m)],
                 wtr[(lag1vec[lag1+1]+lag2vec[lag2]+res),3:7])
            if (cmp2){#there is water one further lag
              #water is assumed to be there continuously
              F[lag1]<-F[lag1]+1
            }else{break}
          }
        }
      }
      #print(c(wtr[(lag2vec[lag2]+res),2],wtr[(lag1vec[lag1]+lag2vec[lag2]+res),2]))
      #print(c(lag2vec[lag2]+res,sum(wtr[(lag2vec[lag2]+res),3:7]==wtr[(lag1vec[lag1]+lag2vec[lag2]+res),3:7],na.rm=T)))
      #print(rbind(c(lag2vec[lag2]+res,as.numeric(wtr[(lag2vec[lag2]+res),3:7])),c(lag1vec[lag1]+lag2vec[lag2]+res,as.numeric(wtr[(lag1vec[lag1]+lag2vec[lag2]+res),3:7]))))
      #print(c("lag1",lag1,"lag2",lag2,"res",res,"row1",lag2vec[lag2]+res,"row2",lag1vec[lag1]+lag2vec[lag2]+res))
    }
  }
}


