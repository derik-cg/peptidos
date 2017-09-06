#this script produces a graph of the contacts betweeen 
#a peptide and a bilipd layer.
#Research with Jardon
rm(list=ls())
setwd("~/Investigacion/peptidos")
#read the file line by line
pep<-readLines("pep94_contacts_residue.dat")
pep<-readLines("pep114_contacts_residue.dat")
pep<-readLines("pep193_contacts_residue.dat")
#break all lines into ordered pairs
len<-length(pep)
matpep<-c(NA,NA)
for (i in 1:len)
{
  tmpvec<-unlist(strsplit(pep[i],split=" "))#splits into residues
  lenren<-length(tmpvec)#gets the length
  frame<-rep(tmpvec[1],(lenren-1))#repeats the frame number
  resnm<-tmpvec[2:lenren]#gets the residues
  matframe<-cbind(frame,resnm)#bind them together
  matpep<-rbind(matpep,matframe)
}
matpep<-matpep[2:(len+1),]

#append column with number frame
matpepgr<-as.numeric(matpep[,1])
require(stringr)
#append column with residue code
pepcode<-str_extract(matpep[,2],"[A-Z]+")
#append column with residue position
matpepgr<-cbind(matpepgr,as.numeric(str_extract(matpep[,2],"[1-9]+")))
plot(matpepgr[,1],matpepgr[,2],xlab="Frame",ylab="Residue",pch=".")
