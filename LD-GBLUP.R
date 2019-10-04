
ld_data=read.table("hapmap_22_out.txt",header=F)

names(ld_data)=colnames(nmat)

map_new=read.table("plink.map",header=F)

ld=ld_data[,names(ld_data) %in% map_new$V2]



library(genetics)
library(reshape)

L1=makeGenotypes(ld)
ld1=LD(L1)






###############################################

p-link --file hapmap --make-bed

#kinship
./gcta64 --bfile plink --make-grm --out kin



#simulation GCTA
./gcta64  --bfile plink  --simu-qt  --simu-causal-loci snp_list.txt  --simu-hsq 0.2 --simu-rep 50   --out test

#variance estimation
./gcta64 --reml --grm kin --pheno test.phen  --mpheno 1 --out var1

#simulation LDAK

./ldak.4.9 --bfile plink --her 0.2 --num-phenos 50 --num-causals 100 --make-phenos test.phen 

#cut weights

./ldak.4.9 --cut-weights  /home/boby/boby/GBLUP/data/Simulation/hapmap/mono/weight/  --bfile plink

# calcaulting weights

./ldak.4.9 --calc-weights  /home/boby/boby/GBLUP/data/Simulation/hapmap/mono/weight/  --section 1 --bfile plink

#join weights

./ldak.4.9 --join-weights  /home/boby/boby/GBLUP/data/Simulation/hapmap/mono/weight/ --bfile plink


# calcaulte LDAK kinship



 ./ldak.4.9 --calc-kins-direct kin.out --kinship-raw YES  --weights /home/boby/boby/GBLUP/data/Simulation/hapmap/mono/weight/weightsALL   --bfile plink


# calcaulte LDAK varaince

./ldak.5.93 --reml reml.out --pheno test.phen --grm kin.out --weights /home/boby/boby/GBLUP/data/Simulation/hapmap/weight/weightsALL  --mpheno 2  --bfile plink




####huamn data


#################33
#SNP list 

map_new=read.table("plink.map",header=F)
ran=sample(1:dim(map_new)[1],370)
snp=map_new$V2[ran]
write.table(snp,"snp_list.txt",quote=F,row.names=F,col.names=F)



system("./gcta64  --bfile plink  --simu-qt  --simu-causal-loci snp_list.txt  --simu-hsq 0.7 --simu-rep 50   --out test")


#######################################################


#######################################################
#GBLUP

#######################################################

#######################################################




#system("./ldak.4.9 --bfile plink --her 0.7 --num-phenos 50 --num-causals 700 --make-phenos test")


#system("mv test.pheno test.phen")

### change here 
her=0.7

library(synbreed)

library(rrBLUP)

dat=read.table("test.phen",header=F)
dat$V1=NULL

mat=read.table("hapmap_3k.txt",header=F)
mat=as.matrix(mat)

mat[mat==1]=0
M=mat-1




A=A.mat(M)
row.names(A)=dat$V2


var_g_van=vector()
var_err_van=vector()

for (r in 2:51) 
   {	 

ans=kin.blup(data=dat,geno="V2",pheno=names(dat)[r],K=A)
	
var_g_van[r]=ans$Vg
var_err_van[r]=ans$Ve
}

#######################################################


#######################################################
#LDAK

#######################################################

#######################################################
library(synbreed)
library(rrBLUP)
library(MASS)
library(Matrix)


dat=read.table("test.phen",header=F)
dat$V1=NULL

kin=read.table("kin_hapmap_05_ld.txt",header=F)
kin=as.matrix(kin)

row.names(kin)=dat$V2


PD=nearPD(kin)
GC=as.matrix(PD$mat)
row.names(GC)=dat$V2


var_g_ldak=vector()
var_err_ldak=vector()

for (r in 2:51) 
   {
ans_ldak=kin.blup(data=dat,geno="V2",pheno=names(dat)[r],K=GC)


var_g_ldak[r]=ans_ldak$Vg
var_err_ldak[r]=ans_ldak$Ve
}


#######################################################
#######################################################
#LD based approch
######################################################


library(synbreed)
library(rrBLUP)
library(MASS)
library(Matrix)
load("hapmap05_ld.Rdata")

dat=read.table("test.phen",header=F)
dat$V1=NULL

mat=read.table("hapmap_3k.txt",header=F)
mat=as.matrix(mat)



mat[mat==1]=0
M=mat-1



	C1=as.matrix(ld1$"D")
	diag(C1)=1
	ind <- lower.tri(C1)
	C1[ind] <- t(C1)[ind] 


# P value
##########
P1=as.matrix(ld1$"P-value")
	diag(P1)=0
	ind <- lower.tri(P1)
	P1[ind] <- t(P1)[ind] 




#F=which(P1>0.002,arr.ind = TRUE)

#C1[F]=0

##########


Kinv=ginv(C1)

#Kinv=diag(2102)




p1=(apply(M,2,sum)+nrow(M))/(nrow(M)*2)
		p=2*(p1-.5)
		P = matrix(p,byrow=T,nrow=nrow(M),ncol=ncol(M))
		Z = as.matrix(M-P)
		b=1-p1
		ps=p1%*%Kinv
		c=ps%*%b
		d=2*c
		ZK=Z%*%Kinv
		S=ZK%*%t(Z)
		U <- (S/d[1])
	row.names(U)=dat$V2

PD=nearPD(U)
GC=as.matrix(PD$mat)
row.names(GC)=dat$V2

var_g_ld=vector()
var_err_ld=vector()




for (r in 2:51) 
   { 

######
# recombination method with scaling

	ans_ld=kin.blup(data=dat,geno="V2",pheno=names(dat)[r],K=GC)
	
	
var_g_ld[r]=ans_ld$Vg
var_err_ld[r]=ans_ld$Ve

}

system("./loop_50.sh")
system('grep "V(G)/V" all.var >gcta.txt')

system("echo $'V1\tV2\tV3' | cat - gcta.txt > gcta_G.txt")

gc=read.table("gcta_G.txt",header=T)



var_g_van
var_g_ldak
var_g_ld

var_err_van
var_err_ldak
var_err_ld



var_g_ldak[2:51]/(var_g_ldak[2:51]+var_err_ldak[2:51])
var_g_ld[2:51]/(var_g_ld[2:51]+var_err_ld[2:51])
var_g_van[2:51]/(var_g_van[2:51]+var_err_van[2:51])

rep=50

Method=rep(c("VanRaden G","LDAK","Obs. LD-G","GCTA"),each=rep)

vanG=(var_g_van[2:51]/(var_g_van[2:51]+var_err_van[2:51]))-her
LDAK=(var_g_ldak[2:51]/(var_g_ldak[2:51]+var_err_ldak[2:51]))-her
LD=(var_g_ld[2:51]/(var_g_ld[2:51]+var_err_ld[2:51]))-her
GC=gc$V2-her

All=c(vanG,LDAK,LD,GC)


res=data.frame(Method=Method,val=All)

library(ggplot2)
library(gridExtra)
library(gtable)


p=ggplot(res, aes(Method, val,fill = factor(Method, levels=unique(Method)))) +ylab("Deviation from the true values")+xlab("Method")+scale_fill_manual(name="Method",values = c("brown3", "grey69", "forestgreen","slategray3"))+theme(legend.position = "NULL")


postscript("rep.eps",horizontal=T,height = 11, width = 9)

p+ stat_boxplot(geom ='errorbar')+geom_boxplot()

dev.off()




##############

#Plot decay
##################33



load("hapmap05_ld.Rdata")
library(reshape)


C=ld1$"R^2"

sel=seq(1,dim(C)[1],by=5)

C1=C[sel,sel]



map1=read.table("map1.txt",header=T)

map=map1[map1$marker %in% row.names(C1),]


LDdata1m1 <- melt(C1)[melt(upper.tri(C1))$value,]
names(LDdata1m1)[names(LDdata1m1) == 'X1'] <- 'marker'
LDdata1m2=merge(LDdata1m1,map,by="marker")
names(map)[names(map) == 'marker'] <- 'marker2'
names(map)[names(map) == 'pos'] <- 'pos2'
names(LDdata1m2)[names(LDdata1m2) == 'X2'] <- 'marker2'
LDdata1m3=merge(LDdata1m2,map,by="marker2")
LDdata1m4 <- cbind ( LDdata1m3 , "Distance"= as.numeric(LDdata1m3$pos2)-as.numeric(LDdata1m3$pos))


library(ggplot2)
library(gridExtra)
library(gtable)
plot_dat=data.frame(dist=LDdata1m4$Distance, val=LDdata1m4$value)

#setEPS()
postscript("hapmap22_ld.eps",horizontal=T,height = 11, width = 9)


p=ggplot(plot_dat,aes(x=dist, y=val))+ylab(expression(paste("Linkage Disqulibrium"," "(R^{2}))))+xlab("Distance in base pairs")

p+geom_point(shape=1,size=0.01)
dev.off()


