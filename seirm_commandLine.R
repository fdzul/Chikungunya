
#####################################################################
# seirm.R  	                                                     #
# Author:  Kate Cooper                                               #
# Version: 20140925                                                  #
# Goal:    modified implementation of the ODE system i               #
#          Ruiz-Moreno et al, PLoS Neglected Tropical Diseases 2012  #
######################################################################
args <- commandArgs()
country<-paste(args[6])
S<-as.numeric(args[7])
E<-as.numeric(args[8])
I<-as.numeric(args[9])
R<-as.numeric(args[10])
G<-as.numeric(args[11])
temp<-as.numeric(args[12])
D<-L<-As<-Ae<-Ai<-G

library(deSolve)
source('ruizMoreno.R')
N<-S+E+I+R
e<-exp(1)
parms<- returnTempDepParams(temp)
parms<- c(parms,returnParams())
inits<-c(S=S,E=E,I=I,R=R,D=D,G=G,L=L,As=As,Ae=Ae,Ai=Ai)
#dt =  seq(0, 119, by = 7)
dt = seq(1,4, by = 1/7)
SEIRM<- function(t,x,parms){
        with(as.list(c(parms,x)),{
			dG 	<- f*Am-(mG+zG+xG)*G
			dD 	<- zG*G-vD*D
			dL 	<- vD*D-mL*L-xL*L
			dAs	<- xL*L-mAe*As-((b/N)*Tmh*I*As)
			dAe 	<- (b/N)*Tmh*I*As-mAe*Ae-sAe*Ae
			dAi 	<- sAe*Ae-mAi*Ai
			dS 	<- N-mu*(S*I)-(b/N*Thm*Ai*S)
			dE 	<- (b/N)*Thm*Ai*S-(mu+inf+Thm)*E
			dI 	<- inf*E-((mu+Thm+dmu)*I)
			dR 	<- inf*I-mu*R
        		N<-S+E+I+R
			der<-c(dS,dE,dI,dR,dG,dD,dL,dAs,dAe,dAi)
        	list(der)
        })
}
simulation<-as.data.frame(lsoda(inits,dt,SEIRM,parms=parms))
pdf(paste(country,".pdf",sep=""))
plot(simulation$I,main="SEIRM Model",xlab="Time (Days)", ylab="Individual Count",type="l",col="red",ylim=c(0,N))
lines(simulation$S,type="l",col="green",ylim=c(0,N))
lines(simulation$E,type="l",col="orange",ylim=c(0,N))
lines(simulation$R,type="l",col="blue",ylim=c(0,N))
legend('topleft', legend=c('Susceptible','Exposed', 'Infected', 'Recovered'), col=c('green','orange','red','blue'), lwd=2.5)
dev.off()
write.table(simulation$I, file=paste(country,".csv",sep=""),append=FALSE,sep=",",eol=",",na="N/A",row.names=FALSE,col.names=FALSE)
