#this script computes the distribution of times in a neighborhood

# change directory
setwd("~/Investigacion/peptidos")
#the nature of the file is the length of the string of ids
#open a connection to the file in mode read-only
wtr<-read.csv(file="trimer-1-15.dat",sep=" ",header=F,
              col.names=c("frame","residue","w1","w2","w3","w4","w5"),
              stringsAsFactors = F)
save(file="frequencies.RData",list=c("wtr"))
# The frequency function adds one if the water molecule
# remains near a specific residue in consecutive frames.
# then the origin is moved.

# for each frame the number of residues is the same, and
# the order of the residues is the same

# lag2 shifts the origin in the basic loop. 
#this lag also has a tail that must be taken into account
#in order to have complete sequences of 5000 time lag1
#In order to maintain consistency in the frequencies, a
#maximum number of time lags was set to 5000

# the sizing of the system should leave the last 5000
# even combining three things: 1. lag1, 2. lag2, 3. res

#this is the lag size
lagmax<-100
#this is the number of residues per frame. Should be consistent
residues<-84
#the output structure is a vector with Fs
#keep track of the entries of this vector
F<-rep(0,1000)
#this is the lag2 loop. This shifts the origin of the
#basic loop.   #count-lagmax

# loops -------------------------------------------------------------------

for (lag2 in (seq(from=0,to=101-lagmax,by=residues))){
  ###############################################################
  #this is the basic loop. This finds the F() and accumulates
  #it. This works without the lag2, that is the shifting of 
  #the origin
  for (res in 1:residues){ #residues cannot exeed 84 frame starts at one
    #this part takes all water molecules at once
    #use index as an anchor and increase frames using a counter
    if (sum(is.na(wtr.df[(res),3:12]))==10){
      #there is no water. Add a zero to F. Do nothing
    }else{
      # look lag frames
      for (lag1 in (seq(from=0,to=lagmax,by=84))){ #lag over future frames
        # This finds the number of equal water ids and accumulates
        #to the current lag1
        F[lag1+1]<-F[lag1+1]+
          sum(wtr.df[(res+lag2),3:12]==wtr.df[(res+lag1+lag2),3:12],
              na.rm=T)
#        print(sum(wtr.df[(res+lag2),3:12]==wtr.df[(res+lag1+lag2),3:12],
#                    na.rm=T))
#this is for checking that the same residue appears in
#in all lag1
#        print(c("lag1",lag1,"lag2",lag2,"res",res,"row",lag1+lag2+res))
#        print(c(wtr.df[(res+lag2),2],wtr.df[(res+lag1+lag2),2]))
        print(c(res+lag2),(res+lag1+lag2))
      }
    }
  }############### here ends the basic loop ####################
}
