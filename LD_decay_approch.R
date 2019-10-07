
library(genetics)
library(reshape)
library(LDheatmap)
library(synbreed)
library(MASS)
library(rrBLUP)


dat_m=read.table("synbreed.txt",header=T,row.names=1,na.strings="NA")
synnew=t(dat_m)
gp =create.gpData(geno=synnew)
gp.coded =codeGeno(gp,impute=F,label.heter=c("AC","CA","AG","GA","AT","TA","CT","TC","CG","GC","TG","GT"))
mat=gp.coded$geno

phe=read.table("pheno.txt",header=T,na.strings="NA")
nmat=mat[row.names(mat) %in% phe$Gid,]
newphe=phe[order(phe$Gid),]


# normal rrBLUP
A=A.mat(nmat-1)
ans=kin.blup(data=newphe,geno="Gid",pheno="Green",K=A)
print("normal rrBLUP")
cor(ans$g,newphe$Green)
ans$Vg
ans$Ve

# decay approach

geno=read.table("map.txt",header=T,na.string="NA")
M=nmat-1
n =nrow(M)
mark = ncol(M)
        
K=matrix(0,nrow=mark,ncol=mark)
map=geno$Pos
map=map/100
chr=geno$Chr
n=1
 for (i in n:mark)
  {
 for (j in n:mark)
 		{

			if(chr[i]-chr[j]==0)
			{
			r=map[j]-map[i]
			
			K[i,j]=exp(-r*26)
 			K[j,i]=K[i,j]
			}
			else{K[j,i]=0}
			
 				
		}
	
         n=n+1

}

		Kinv=ginv(K)
		Z = as.matrix(M)
		ZK=Z%*%Kinv
		U=ZK%*%t(Z)
		



name=newphe$Gid
row.names(U)=name
ans=kin.blup(data=newphe,geno="Gid",pheno="Green",K=U)
print("decay approach")
cor(ans$g,newphe$Green)
ans$Vg
ans$Ve

# decay approach with the scaling
M=nmat-1
p1=(apply(M,2,sum)+nrow(M))/(nrow(M)*2)
		p=2*(p1-.5)
		P = matrix(p,byrow=T,nrow=nrow(M),ncol=ncol(M))
		Z = as.matrix(M-P)
		
		b=1-p1
		c=p1*b
		d=2*(sum(c))
		ZK=Z%*%Kinv
		S=ZK%*%t(Z)
		U <- (S/d)
name=newphe$Gid
row.names(U)=name
ans=kin.blup(data=newphe,geno="Gid",pheno="Green",K=U)
print("decay approach with the scaling")
cor(ans$g,newphe$Green)
ans$Vg
ans$Ve



