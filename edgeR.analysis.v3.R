library(edgeR)
## TODO: 
# Reportar
# aclimatacion
# <dia2> -  <dia1>
# 27 - 12
# clusters FULL: con raw data

## Funcion para encontrar factores de un numero
## lo uso para layout
get_all_factors <- function(n){
  prime_factor_tables <- lapply(n, function(i)
  {
    if(i == 1) return(table(1L))
    table(as.numeric(gmp::factorize(i)))
  })
  lapply(prime_factor_tables, function(pft)
  {
    powers <- mapply(
      function(name, value) as.numeric(name) ^ seq.int(0L, value),
      names(pft), 
      pft
    )
    power_grid <- expand.grid(powers)
    sort(apply(power_grid, 1, prod)) 
  })
}
# get_all_factors(c(1L, 7L, 60, 2520))


if(FALSE){
  source("rnaSeq_fun.R")  
  myCPM<-function(y,blog2=FALSE,offset=0){
    x<-y$counts
    lsize<-apply(y$samples[,2:3],1,function(xx){xx[1]*xx[2]})
    for(j in seq_along(x[1,])){
     x[,j]<-x[,j]/lsize[j]*1e6
    }
    if(blog2){
     return(log2(x+offset))
    }else{
     return(x+offset)
    }
   }
}

#import counts
if(FALSE){
  temp<-c(12,17,22,27)
  ppheno<-c()
  reads<-c()
  for(i in seq_along(temp)){
    at    <- paste("at_",temp[i],"_",sep="")
    (load(paste("counts/counts",temp[i],".RData",sep="")))
    cn<-colnames(counts$gene.counts)[grep(at,colnames(counts$gene.counts))]
    if(i==1){
      reads<-counts$gene.counts
      ppheno <- data.frame(cn,temp=rep(temp[i],length(cn)),matrix(unlist(strsplit(sub(at,"",cn),"_",fixed=TRUE)),byrow=TRUE,ncol=2))
    }else{
      reads<-cbind(reads,counts$gene.counts[,cn])
      ppheno <- rbind(ppheno,data.frame(cn,temp=rep(temp[i],length(cn)),matrix(unlist(strsplit(sub(at,"",cn),"_",fixed=TRUE)),byrow=TRUE,ncol=2)))
    }
    rm(counts)
    gc()
  }
  colnames(ppheno)<-c("sample","temperature","time","replicate")
  meta<-reads[,1:5]
  reads<-reads[,-c(1:5)]

  pheno<-data.frame(ppheno,hpi=rep(rep(c(1,5,9,13,17,21,25,29,33,37,41,45),each=2),4),
                         light=rep(rep(c(1,0),each=6),8),
                         daytime=rep(rep(c(1,5,9,13,17,21),each=2),8))

  if(FALSE){
    save(pheno,reads,meta,file="reads96.Rdata")
  }
}else{
 (load("reads96.Rdata"))
}

(load("/home/arabinov/doctorado/programacion/rna_seq/viejo/reads96.Rdata"))
##
## edgeR

 #cuantos reads tengo en total?
 (N<-sum(reads))
#  1115399630

#consideramos un gen expresado si recibio al menos minCountPM cuentas por millon de reads  a lo largo de todo el dataset
k <- c()
m <- seq(5, 100, 5)
minCountPM = 5
n     <- minCountPM*N/(1e6*ncol(reads))
keep  <- apply(reads,1,sum) >= n

for(i in m){
  print(i)
  minCountPM = i
  n     <- minCountPM*N/(1e6*ncol(reads))
  keep  <- apply(reads,1,sum) >= n
  k <- c(k, sum(keep))
}
plot(m, k)
table(keep)
# keep
# FALSE  TRUE 
# 10487 23115 

#Parece poco sensible al minCountPM así que voy a usar los mismos criterios de filtrado que uso para las redes.
#Los que pasan el criterio de redes pasan también CMP pero son bastante menos, so...
i1    <- apply(reads,1,function(x){mean(x) > 10})
i22   <- apply(reads/meta[,"effWidth"],1,function(x){mean(x) > 0.05})
keep <- ipass <- i1 & i22
table(ipass)
# ipass
# FALSE  TRUE 
# 19313 14289  

cond   <- factor(apply(pheno,1,function(x){paste(x[2:3],collapse=".")}))
#design <- model.matrix(~ temperature * time, data=pheno[,2:4])

pheno$temperature<-relevel(as.factor(pheno$temperature),"22")
pheno$time<-factor(pheno$time,levels=c(1,2,3,4,5,6,7,8,9,10,11,12))
design <- model.matrix(~ time + time:temperature, data=pheno[,2:4])
cond   <- factor(paste(pheno$temperature,pheno$time,sep="."))
cond   <- relevel(cond,"22.1")
 
 
#calculo modelo 
bcalcul=FALSE
if(bcalcul){ 


y <- DGEList(counts=reads[keep,],group=cond,genes=meta[keep,])

cpm <- apply(y$counts,2, function(x) (x/sum(x))*1000000) 
cpm[1, 1]
countsPerMillion <- cpm(y)
countsPerMillion[1, 1]
summary(countsPerMillion)
countCheck <- countsPerMillion > 1
head(countCheck)
keep <- which(rowSums(countCheck) >= 2)
y <- y[keep,]
summary(cpm(y))
y <- calcNormFactors(y,method=c("TMM","RLE","upperquartile")[1])
y$samples

#data exploration
if(FALSE){
a<-plotMDS(y,main="BCV distance")

if(FALSE)dev.copy(pdf,file="edgeR.MDS.BCVdistance.pdf",width=8,height=8);dev.off()
}

# Detecto genes que se mueven a 22grados
if(FALSE){

i22<-which(pheno$temperature==22) 
#design22 <- model.matrix(~ 0+time, data=pheno[i22,]) ## 2015/02/12: OJO con esto: este modelo no parece servir para ver movimiento
design22 <- model.matrix(~ time, data=pheno[i22,])

#consideramos un gen expresado si recibio al menos minCountPM cuentas por millon de reads  a lo largo de todo el dataset
n22 <- minCountPM*N/(1e6*ncol(reads[,i22]))
minCountPM = 5
keep22 <- apply(reads[,i22],1,sum)>=n22   ## 2015/02/12: OJO con esto! Faltaba el i22
#   keep   keep22
# 3787 19328 0
y22 <- DGEList(counts=reads[keep22,as.character(pheno$sample[i22])],
               group=pheno$hpi[i22],genes=meta[keep22,])

y22 <- calcNormFactors(y22,method=c("TMM","RLE","upperquartile")[1])
y22$samples
y22 <- estimateGLMTrendedDisp(y22,design22)
y22 <- estimateGLMTagwiseDisp(y22,design22)
fit22<-glmFit(y22,design22)
lrt22<-glmLRT(fit22,coef=2)
a<-topTags(lrt22,n=length(lrt22$table[,1])) 
sum(a$table[,"FDR"]<0.01)

}



#y <- estimateGLMCommonDisp(y, design,verbose=TRUE) 
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)
if(FALSE){
plotBCV(y)
dev.copy(png,file="plotBCV.png",width=800,height=800)
dev.off()
} 

# el modelo que fitea es
# log(mu_{gi}) = (x_i)^t beta_g + log N_i
# con 
#  mu_{gi}: valor medio de counts gen g en sample i
#  N_i    : tamanio libreria i normalizada
# log:logaritmo Nepperiano
fit<-glmFit(y,design)

if(FALSE) save(fit,y,file="yfit96.Rdata")
}else{
(load("yfit96.Rdata"))
}

 
 
 if(FALSE){
  lrt<-glmLRT(fit,coef=2:12)
  a<-topTags(lrt,n=length(lrt$table[,1])) 
  lrt.time<-lrt

  lrt<-glmLRT(fit,coef=13:24)
  b<-topTags(lrt,n=length(lrt$table[,1])) 
  lrt.12<-lrt

  lrt<-glmLRT(fit,coef=25:36)
  b<-topTags(lrt,n=length(lrt$table[,1])) 
  lrt.17<-lrt
  
  lrt<-glmLRT(fit,coef=37:48)
  c<-topTags(lrt,n=length(lrt$table[,1])) 
  lrt.27<-lrt
  if(FALSE) write.table(a,file="edgeR.log2FC.csv",sep="\t",quote=FALSE)
 }
 
## Calculo 
##   lrt Fisher para algun tiempo
##   llrt para cada factor time:temp
##   llrt.ac: contrastes de aclimatacion (time25:TX - time1:TX, time29:Tx - time5:Tx,...)
## pvalues de estos ultimos
## notar que los test de cjtos de coef y el de c/cef por separado no son
## exactamente equivalentes...hago esto ultimo para 'identificar' que punto temporal
## 'disparo' la dif observada.
bcalcul<-FALSE
if(bcalcul){
  ttemps<-as.character(c(12,17,27))
  lcoefs<-list(13:24,25:36,37:48) #indices de: colnames(fit$design)
  names(lcoefs)<-ttemps
  lttag<-lttagAcli<-lttagTimePoint<-lmpval<-list()
  for(k in seq_along(ttemps)){   
   lrt<-glmLRT(fit,coef=lcoefs[[ttemps[k]]])
   x<-topTags(lrt,n=length(lrt$table[,1])) 
   llrt <-llrt.ac<-list()
   for(ij in seq_along(x$comparison)){
    cat("Calculating",x$comparison[ij],"\n")
    a<-glmLRT(fit,coef=x$comparison[ij])
    llrt[[x$comparison[ij]]]<-topTags(a,n=length(a$table[,1]))
    if(ij <= length(x$comparison)/2){
     mycont<-rep(0,ncol(design))
     mycont[colnames(design)==x$comparison[ij+6]]<-  1
     mycont[colnames(design)==x$comparison[ij]]  <- -1    
     a<-glmLRT(fit,contrast=mycont)
     llrt.ac[[x$comparison[[ij]]]]<-topTags(a,n=length(a$table[,1])) 
    }
   }
   
   lttag[[ttemps[k]]]         <- x
   lttagTimePoint[[ttemps[k]]]<- llrt
   lttagAcli[[ttemps[k]]]<- llrt.ac
#    assign(paste("ttag",ttemps[k],sep="."),topTags(lrt,n=length(lrt$table[,1])))
#    assign(paste("ttagTimePoint",ttemps[k],sep="."),topTags(llrt,n=length(llrt$table[,1])))
#    assign(paste("mpval",ttemps[k],sep="."),mpval) 
  }
  
  #contrastes t27 -  t12
  mycont<-c()
  for(i in 1:12){
   a<-rep(0,length=length(colnames(fit$coefficients)))
   a[(13:24)[i]]<- -1
   a[(37:48)[i]]<- 1
   mycont<-cbind(mycont,a)
  } 
  colnames(mycont)<-paste("time",1:2,sep="")
  a<-glmLRT(fit,contrast=mycont)
  ttag2712<-topTags(a,n=length(a$table[,1])) 
  
  
  
  if(FALSE) save(lttag,lttagTimePoint,lttagAcli,ttag2712,file="lttag96.Rdata")
 }else{
  (load("lttag96.Rdata"))
}

if(FALSE){
}

##Analizo un poco los genes que se movieron globalmente 
##respecto a T22
pvlim <- 1e-6 #1e-4
lfc   <- log2(2) 
lres<-list()
u<-c()
for(i in seq_along(lttag)){
 tt<-lttag[[i]]$table
 respv <- as.numeric(tt["FDR"]<pvlim)
 clfc<-grep("logFC",colnames(tt),fixed=TRUE)
 if(length(clfc)==1){
  resfc<-abs(tt[,clfc])>lfc
 }else{
  resfc<-abs(apply(tt[,clfc],1,max))>lfc
 }
 lres[[i]]<-which(respv * resfc != 0)
 u<-union(u,names(lres[[i]]))
}
mres<-cbind(u%in%names(lres[[1]]),
            u%in%names(lres[[2]]),
            u%in%names(lres[[3]]))
colnames(mres)<-paste("T",names(lttag),sep="")
rownames(mres)<-u
vennDiagram(mres)


## Visualizacion

ycpm <-cpm(y,log=FALSE,normalized=TRUE)
if(FALSE){
 require(vsn)
 ycpm.vsn<-justvsn(ycpm)    #cpm variancStabilized
 save(ycpm.vsn,file="ycpm.vsn.Rdata")
}else{
 load("ycpm.vsn.Rdata")
}
if(FALSE){
 ycpm2<-2^cpm(y,log=TRUE,normalized=TRUE,prior.count=2)
 plot(ycpm[,1],ycpm2[,1])
} 

ccol<-match(pheno$temperature,unique(pheno$temperature))
ccol<-rainbow(4)[ccol]

#Elijo genes
#genes con curvas diferentes a T22
genes<- rownames(mres)[!mres[,1] & mres[,2] & !mres[,3]]
saux <- "curvas diferentes"
#genes con efecto aclimatacion T12 time1
tt<-lttagAcli[["12"]][["time1:temperature12"]]$table
genes<- rownames(tt)[tt$FDR<0.001 & abs(tt$logFC)>log2(3)]
saux<-"aclimatacion T12 time1"
#genes con efecto aclimatacion T27 time1
tt<-lttagAcli[["27"]][["time1:temperature27"]]$table
genes<- rownames(tt)[tt$FDR<0.001 & abs(tt$logFC)>log2(2)]
saux<-"aclimatacion T27 time1"

#2 genes del reloj
genes<-c("AT1G01060","AT2G46830","AT2G02570")
saux<-"1+2 genes del reloj"

#genes 2712
tt<-ttag2712$table
fc<-apply(tt[,6:17],1,function(x){any(x>log2(3))})
genes<- rownames(tt)[tt$FDR<0.001 & fc]
length(genes)
saux<-"T27-T12"
genes2712<-genes


#
genes<-c("AT5G12010","AT5G49910")
saux<-"foo"

#Por que los dos genes del reloj no estan en el dataset de aclimatacion?
#Rta: esta bien que no este bien. Estaba visualizando sin normalizar por librerias lo que confundia todo....
if(FALSE){
  g <-c("AT1G01060","AT2G46830")[2]
  
  ##empiezo a analizar contrastes de T22
  iT22<-49:72  
  i   <-14
  lrt2<-list()
  for(i in 2){
   lrt2[[i]]<-glmLRT(fit,coef=i)
   names(lrt2)[i]<-lrt2[[i]]$comparison
  } 
  log2(mean(lrt2$fitted.values[g,2*(i-1)+1:2+min(iT22)-1])/mean(lrt2$fitted.values[g,1:2+min(iT22)-1]))# -5.718048
  
  if(FALSE){ #Gordon dixit
# Dear Xinwei,
# 
# 
# No, not quite.
# 
# 
# It would be best to try to think in terms of linear comparisons on the
# log-scale, because that is how glms works. A-B-C+D means
# 
# 
# (average log CPM in A)
# minus
# (average log CPM in B)
# minus
# (average log CPM in C)
# add
# (average log CPM in D)
# 
# 
# You can think of it as twice
# 
# 
# (average log CPM in A and D)
# minus
# (average log CPM in B and C).
# 
# 
# or you could think of it as
# 
# 
# (log fold change from A to B)
# plus
# (log fold change from D to C)
# 
# 
# but none of these are quite the same as what you wrote. Your expression
# interchanges averages and logs in a way that is not correct.
# 
# 
# Note that the logFC you get from the contrast c(1,-1,-1,1) is exactly
# equal to the logFC you get from c(1,-1,0,0) plus what you get from
# (0,0,-1,1). It just adds up in the natural way.
# 
# 
# Best wishes
# Gordon
# 
# 
# 
# 
# On Tue, 2 Apr 2013, Xinwei Han wrote:
# 
# HI Gordon,
# 
# Thanks. Does that mean in the second case, basically logFC = log2 (
# (average CPM in D+ average CPM in A) / (average CPM in B +average CPM in
# C) ), except for other adjustments from negative binomial?
# 
# Thanks,
# Xinwei
# 
# On Tue, Apr 2, 2013 at 7:03 PM, Gordon K Smyth wrote:
# 
# Dear Xinwei,
# 
# Well, in the first case you get C-A (C minus A). In the second case you
# get A-B-C+D.
# 
# Best wishes
# Gordon
# 
# Date: Mon, 1 Apr 2013 23:10:42 -0400
# From: Xinwei Han <xinwei.han@duke.edu>
# To: <bioconductor at stat.math.ethz.**ch <bioconductor@stat.math.ethz.ch>>
# Subject: [BioC] logFC in GLM of egdeR
# 
# Hi,
# 
# I fit a GLM with edgeR like this:
# 
# colnames(fit)
# A B C D
# 
# For contrast = c(-1,0,1,0), I guess logFC = log2 (average CPM in C/
# average
# CPM in A).
# 
# But for contrast = c(1,-1,-1,1), how is logFC calculated?
# 
# Thanks,
# Xinwei
# ______________________________**______________________________**__________
# The information in this email is confidential and intended solely for the
# addressee.
# You must not disclose, forward, print or use it without the permission of
# the sender.
# ______________________________**______________________________**__________  
  
  }
  

  layout(matrix(1:2,ncol=2))
  plot(cpm(y)[g,iT22],ylim=c(0,1000))
  points(cpm(y)[g,iT22-48],pch=18)
  plot(cpm(y)[g,iT22],ylim=c(1,1000),log="y")
  points(cpm(y)[g,iT22-48],pch=18,log="y")
   
  abline(v=12.5,col="gray",lty=2)
  
  if(FALSE){
   layout(matrix(1:2,ncol=2))
   plot(fit$counts[g,iT22],ylim=c(0.1,20000))
   points(fit$fitted.values[g,iT22],pch=20,col=2)
   points(fit$counts[g,iT22-48],pch=18)
   points(fit$fitted.values[g,iT22-48],pch=20,col=4)
   abline(v=12.5,col="gray",lty=2)

   plot(fit$counts[g,iT22],ylim=c(0.1,20000),log="y")
   points(fit$fitted.values[g,iT22],pch=20,col=2)
   points(fit$counts[g,iT22-48],pch=18)
   points(fit$fitted.values[g,iT22-48],pch=20,col=4)
   abline(v=12.5,col="gray",lty=2)
  }  
  
  #analizo que son los logFC reportados para g
  # me fijo en el contraste "time2"
   (r1<-lrt2[["time2"]]$table[g,"logFC"])  #esto es lo que reporta
  
  #log2 de count per million (con librerias normalizadas)
  (a<-cpm(y,log=TRUE,normalized=TRUE)[g,iT22])
  (r2<-mean(a[3:4])-mean(a[1:2]))
  r1-r2<1e-4  ##!ok  entonces logFC equivale a la diferencia lineal entres log2(cpm)'s 
  
  # me fijo en el contraste "time1:temperature12"
  (r1<-lrt2[["time1:temperature12"]]$table[g,"logFC"])
  (a<-c(cpm(y,log=TRUE,normalized=TRUE)[g,1:2],cpm(y,log=TRUE,normalized=TRUE)[g,iT22[1:2]]))
  (r2<-mean(a[1:2])-mean(a[3:4]))
  r1-r2<1e-4  ##!ok  entonces logFC equivale a la diferencia lineal entres log2(cpm)'s 
  
  # me fijo en el contraste "Coefficient:  -1*time1:temperature12 1*time7:temperature12"
  (r1<-lttagAcli[[1]][[1]][g,"logFC"])
  (a<-c(cpm(y,log=TRUE,normalized=TRUE)[g,1:2],cpm(y,log=TRUE,normalized=TRUE)[g,iT22[1:2]],
        cpm(y,log=TRUE,normalized=TRUE)[g,13:14],cpm(y,log=TRUE,normalized=TRUE)[g,iT22[13:14]] ))
  (r2<-mean(a[1:2])-mean(a[3:4])- (mean(a[5:6])-mean(a[7:8])) )
  r1$table$logFC-r2<1e-4  ##!ok  entonces logFC equivale a la diferencia lineal entres log2(cpm)'s 

  blog=FALSE
  bnorm=FALSE
  (a00<-c(cpm(y,log=blog,normalized=bnorm)[g,1:2],cpm(y,log=blog,normalized=bnorm)[g,iT22[1:2]],
        cpm(y,log=blog,normalized=bnorm)[g,13:14],cpm(y,log=blog,normalized=bnorm)[g,iT22[13:14]] ))
  blog=FALSE
  bnorm=TRUE
  (a01<-c(cpm(y,log=blog,normalized=bnorm)[g,1:2],cpm(y,log=blog,normalized=bnorm)[g,iT22[1:2]],
        cpm(y,log=blog,normalized=bnorm)[g,13:14],cpm(y,log=blog,normalized=bnorm)[g,iT22[13:14]] ))

  
}

#yy<-log(reads[genes,]+.5)
#yy<-reads[genes,]
yy<-ycpm[genes,]


#promedio replicas y ploteo
# el modelo que fitea es
# log(mu_{gi}) = (x_i)^t beta_g + log N_i
# con 
#  mu_{gi}: valor medio de counts gen g en sample i
#  N_i    : tamanio libreria i normalizada
# a ver coeficientes: log(mu_{gi}/N_i)
a<-apply(yy,1,function(x){return(apply(matrix(x,nrow=2),2,mean))})
colnames(a)<-genes
layout(matrix(1:2,1,2),width=c(.6,.4))
#readline("press <enter> to start ")
for(i in seq_along(a[1,])){
  plot(a[,i],main=genes[i],typ="b",col="red",pch=20)
  mtext(saux)
  vl<-which(diff(as.numeric(pheno$temperature))!=0)/2+.5
  abline(v=c(vl),col="gray",lty=2)
  
  ma<-matrix(a[,i],nrow=12)
  matplot(1:12,ma[,c(3,1,2,4)],typ="l",col=c("gray","blue","violet","red"),lty=1,lwd=c(15,5,5,5))
#  matplot(1:6,ma[1:6,c(3,1,2,4)],xlim=c(1,12),typ="l",col=c("gray","blue","violet","red"),lty=1,lwd=c(15,5,5,5))
#  matplot(7:12,ma[7:12,c(3,1,2,4)],typ="l",col=c("gray","blue","violet","red"),lty=1,lwd=c(15,5,5,5),add=TRUE)
  mtext(genes[i])
  abline(v=6.5,lty=2,col=1)
  
  readline(paste(i,"/",ncol(a), "press <enter> to continue","\n"))
}

 
## PAra los genes de aclimatacion produzco figura de perfiles
## porque cisualizacion por clusters me parece que no va
if(FALSE){
# a ver coeficientes: log(mu_{gi}/N_i)

for(j in seq_along(lttagAcli)){
   tt<-lttagAcli[[j]][[1]]$table
   genes<-rownames(tt)[tt$FDR<0.001 & abs(tt$logFC)>log2(3)]
   saux<-paste("T",names(lttagAcli)[j],"T22.time1time6",sep="")
   if(length(genes)<1) next
     
   yy<-ycpm[genes,]
   a<-apply(yy,1,function(x){return(apply(matrix(x,nrow=2),2,mean))})
   colnames(a)<-genes
  
  

   #readline("press <enter> to start ")
   pdf(file=paste("aclimat",saux,".FDR1e-3.logFC3.pdf",sep=""),width=8,height=11)
   layout(matrix(1:6,3,2,byrow=TRUE),width=c(.6,.4))
   for(i in seq_along(a[1,])){
     plot(a[,i],main=genes[i],typ="b",col="red",pch=20)
     mtext(saux)
     vl<-which(diff(as.numeric(pheno$temperature))!=0)/2+.5
     abline(v=c(vl),col="gray",lty=2)
  
     ma<-matrix(a[,i],nrow=12)
     matplot(1:12,ma[,c(3,1,2,4)],typ="l",col=c("gray","blue","violet","red"),lty=1,lwd=c(15,5,5,5))
#  matplot(1:6,ma[1:6,c(3,1,2,4)],xlim=c(1,12),typ="l",col=c("gray","blue","violet","red"),lty=1,lwd=c(15,5,5,5))
#  matplot(7:12,ma[7:12,c(3,1,2,4)],typ="l",col=c("gray","blue","violet","red"),lty=1,lwd=c(15,5,5,5),add=TRUE)
     mtext(genes[i])
     abline(v=6.5,lty=2,col=1)
  
     #readline(paste(i,"/",ncol(a), "press <enter> to continue","\n"))
   }
   dev.off()
}



## Analisis Tipo 2A (usando edgeR):
## Para los genes de la Temp T que presentaron diferencias con respecto a la curva control (22)
## (y ademas maxFC >/log2FC), clusterear perfiles (T22,T) (en R^24)
## ATENNCION: 
## uso cpm o cpmvsn
if(FALSE){
 require(dynamicTreeCut)
 
 #voy a seguir con la normalizacion dada por varianceStabilizingTransformation 
 # pero uso la info de expresion diferencial de edgeR!!
 #assay(vsd)
 pvlim    <- 1e-4
 logfclim <- log2(2)

 ds1<-c(0,1,2,3,4)[1]
 ds2<-c(0,1,2,3,4)[3]
 method1<-method2<-"complete"
 bwrite1<-FALSE
 bwrite2<-TRUE
 bplot1 <-FALSE
 bplot2 <-TRUE
 scode<-paste(ds1,ds2,method1,sep="")

 ccol=c("gray","blue","violet","red")
 names(ccol)<-c(22,12,17,27)

 
 for(iTemp in seq_along(lttag)){
  fname<-paste("cpmvsnT22","T",names(lttag)[iTemp],sep="")
 
  scol <- colnames(lttag[[iTemp]]$table)[grep("logFC",colnames(lttag[[iTemp]]$table))]
  maxfc<- apply(lttag[[iTemp]]$table[,scol],1,max)
  ig   <- which(lttag[[iTemp]]$table[,"FDR"]<pvlim & abs(maxfc)>logfclim) #ACA filtro
  genes<- rownames(lttag[[iTemp]]$table)[ig]              

  #columnas de referencia (T22) + T de analisis
  i22<- which(pheno[,"temperature"]==22)
  iT <- which(pheno[,"temperature"]==names(lttag)[iTemp])
    
  geneX <- t(apply(ycpm.vsn[genes,c(i22,iT)],1,function(x){return(apply(matrix(x,nrow=2),2,mean))}))
  rownames(geneX)<-genes
  
  
  #step1
  d<-1-(1+cor(t(geneX)))/2
  h<-hclust(d=as.dist(d),method=method1)

  minClusterSize=20
  deepSplit     =ds1

  ct<-cutreeDynamic(h,minClusterSize=minClusterSize,distM=d,deepSplit=deepSplit)
  ta<-table(ct)
  ua<-sort(unique(as.numeric(ct)))
  z <-t(apply(geneX,1,function(x){return((x-mean(x))/sd(x))}))
  colnames(z)<-c(paste("T22_",c(1:12),sep=""),paste("T",names(lttag)[iTemp],"_",c(1:12),sep=""))
 
  npanel<- signif(length(ua)/10,0)*10
  par(mar=c(2,2,2,2))
  if(bplot1) layout(matrix(1:npanel,ncol=10,byrow=TRUE))
  mgenex<-c()
  for(i in 1:min(length(ua),npanel)){
   ia<-which(ct==ua[i])
   aux<-length(ia)
   xx<-apply(z[ia,],2,mean)
 
   mgenex<-rbind(mgenex,apply(z[ia,],2,mean))

   if(bwrite1==TRUE) {
    write.table(z[ia,],sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE,file=paste(fname,"precluster",i,"csv",sep="."))
   }
   if(bplot1){

    matplot(1:12 ,t(z[ia,1:12]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),xlim=c(1,24),ylim=range(t(z[ia,])))
    matplot(13:24,t(z[ia,13:24]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),add=TRUE)
    
    mtext(paste("preCluster",i," - ",aux),cex=0.7)
    lines(1:12,xx[1:12],lwd=7,col="gray")
    lines(13:24,xx[13:24],lwd=7,col=ccol[names(lttag)[iTemp]])

    if(i%%npanel==0) readline("Sacale una foto y Press <enter>\n")
    if(bwrite2){
     dev.copy(png,file=paste(fname,scode,"preclusters","png",sep="."),width=1200,height=800)
     dev.off()
    }
   }
  }
 
  

  ## step 2
  #dendrograma de expresiones medias
  dd<-0.5*(1-cor(t(mgenex)))
  hh<-hclust(d=as.dist(dd),method=method2)
  if(FALSE) plot(hh)

  #cutree sobre mean profiles
  mminClusterSize=2
  ddeepSplit     =ds2

  cct<-cutreeDynamic(hh,minClusterSize=mminClusterSize,distM=dd,deepSplit=ddeepSplit)
  names(cct)<-seq_along(cct)

  if(FALSE){
   xx<-aux<-c()
   for(i in c(1,2,3,6)){
    aux<-c(aux,i)
    ia<-which(ct%in%aux)
    xx<-c(xx,sqrt(sum(apply(z[ia,],2,function(x){var(x)})^2)))
   }  
  }
 
  tta<-table(cct)
  uua<-sort(unique(as.numeric(cct)))

  rcRatio <- 3/4
  ncol <- signif(sqrt(length(uua)/rcRatio)+0.5,0)
  nraw <- signif(length(uua)/ncol+.5,0)
  npanel<- ncol*nraw
  par(mar=c(2,2,2,2))
  if(bplot2) layout(matrix(1:npanel,ncol=ncol,byrow=TRUE))
  for(i in 1:length(uua)){
   iia<-which(cct==uua[i])
   ia<-which(ct%in%iia)
   aux<-length(ia)
   xx<-apply(z[ia,],2,mean)
   mgenex<-rbind(mgenex,xx)

   if(bwrite2==TRUE) {
    write.table(z[ia,],sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE,file=paste(fname,scode,"cluster",i,"csv",sep="."))
   }
   if(bplot2){
#     matplot(t(z[ia,]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4))
#     mtext(paste("Cluster",i," - ",aux),cex=0.7)
#     lines(1:3,xx[1:3],lwd=7,col="red")
#     lines(4:6,xx[4:6],lwd=7,col="violet")
    
    matplot(1:12 ,t(z[ia,1:12]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),xlim=c(1,24),ylim=range(t(z[ia,])))
    matplot(13:24,t(z[ia,13:24]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),add=TRUE)
    mtext(paste("Cluster",i," - ",aux),cex=0.7)
    
    lines(1:12,xx[1:12],lwd=7,col="gray")
    lines(13:24,xx[13:24],lwd=7,col=ccol[names(lttag)[iTemp]])
    
    
    if(i%%npanel==0) readline("Sacale una foto y Press <enter>\n")
   }
  }
  if(bwrite2){
   dev.copy(png,file=paste(fname,scode,"clusters","png",sep="."),width=1200,height=800)
   dev.off()
  }

  
#  varnorm<-function(x){(x-mean(x))/sd(x)}
#  xx<-t(apply(datX1,1,varnorm))
#  heatmap(xx) 
#   
  
 }

}

## Analisis Tipo 2B (usando DEseq):
## Para los genes de la Temp T que presentaron diferencias con respecto a la curva control (22)
## (y ademas maxFC >/log2FC), clusterear perfiles (T22,T) (en R^24)
## ATENNCION: 
## uso la transformacion de estabilizacion de varianza de DESeq2 pero 
## la info de expresion diferencial de edgeR!!
if(FALSE){
 require(dynamicTreeCut)
 #Empiezo evaluando las maneras de estabilizar varianzas que aparecen en DESeq2
 #analisis DESeq2
 library(DESeq2)
 library("BiocParallel")
 library("vsn")
 if(bcalcul){
  register(MulticoreParam(4))
  dds<-DESeqDataSetFromMatrix(countData = reads[keep,], 
                             colData= pheno,
                             design = ~ time + time:temperature)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)

  rld<-rlog(dds,fast=TRUE)
  vsd <- varianceStabilizingTransformation(dds)
  if(FALSE) save(dds,rld,vsd,file="deseq2.96.Rdata")
 }else{
  require(GenomicRanges)
  (load("deseq2.96.Rdata"))
 }
 #aver...
 if(FALSE){
  par(mfrow=c(1,3))
  notAllZero <- (rowSums(counts(dds))>0)
  meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1))
  mtext("log2")
  meanSdPlot(assay(rld[notAllZero,]))
  mtext("rlog")
  meanSdPlot(assay(vsd[notAllZero,]))
  mtext("vsd")
 }
 
 #voy a seguir con la normalizacion dada por varianceStabilizingTransformation 
 # pero uso la info de expresion diferencial de edgeR!!
 #assay(vsd)
 pvlim    <- 1e-4
 logfclim <- log2(2)

 ds1<-c(0,1,2,3,4)[1]
 ds2<-c(0,1,2,3,4)[3]
 method1<-method2<-"complete"
 bwrite1<-FALSE
 bwrite2<-TRUE
 bplot1 <-FALSE
 bplot2 <-TRUE
 scode<-paste(ds1,ds2,method1,sep="")

 ccol=c("gray","blue","violet","red")
 names(ccol)<-c(22,12,17,27)
 for(iTemp in seq_along(lttag)){
  fname<-paste("T22","T",names(lttag)[iTemp],sep="")
 
  scol <- colnames(lttag[[iTemp]]$table)[grep("logFC",colnames(lttag[[iTemp]]$table))]
  maxfc<- apply(lttag[[iTemp]]$table[,scol],1,max)
  ig   <- which(lttag[[iTemp]]$table[,"FDR"]<pvlim & abs(maxfc)>logfclim) #ACA filtro
  genes<- rownames(lttag[[iTemp]]$table)[ig]              

  #columnas de referencia (T22) + T de analisis
  i22<- which(colData(vsd)[,"temperature"]==22)
  iT <- which(colData(vsd)[,"temperature"]==names(lttag)[iTemp])
    
  geneX <- t(apply(assay(vsd)[genes,c(i22,iT)],1,function(x){return(apply(matrix(x,nrow=2),2,mean))}))
  rownames(geneX)<-genes
  
  
  #step1
  d<-1-(1+cor(t(geneX)))/2
  h<-hclust(d=as.dist(d),method=method1)

  minClusterSize=20
  deepSplit     =ds1

  ct<-cutreeDynamic(h,minClusterSize=minClusterSize,distM=d,deepSplit=deepSplit)
  ta<-table(ct)
  ua<-sort(unique(as.numeric(ct)))
  z <-t(apply(geneX,1,function(x){return((x-mean(x))/sd(x))}))
  colnames(z)<-c(paste("T22_",c(1:12),sep=""),paste("T",names(lttag)[iTemp],"_",c(1:12),sep=""))
 
  npanel<- signif(length(ua)/10,0)*10
  par(mar=c(2,2,2,2))
  if(bplot1) layout(matrix(1:npanel,ncol=10,byrow=TRUE))
  mgenex<-c()
  for(i in 1:min(length(ua),npanel)){
   ia<-which(ct==ua[i])
   aux<-length(ia)
   xx<-apply(z[ia,],2,mean)
 
   mgenex<-rbind(mgenex,apply(z[ia,],2,mean))

   if(bwrite1==TRUE) {
    write.table(z[ia,],sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE,file=paste(fname,"precluster",i,"csv",sep="."))
   }
   if(bplot1){

    matplot(1:12 ,t(z[ia,1:12]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),xlim=c(1,24),ylim=range(t(z[ia,])))
    matplot(13:24,t(z[ia,13:24]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),add=TRUE)
    
    mtext(paste("preCluster",i," - ",aux),cex=0.7)
    lines(1:12,xx[1:12],lwd=7,col="gray")
    lines(13:24,xx[13:24],lwd=7,col=ccol[names(lttag)[iTemp]])

    if(i%%npanel==0) readline("Sacale una foto y Press <enter>\n")
    if(bwrite2){
     dev.copy(png,file=paste(fname,scode,"preclusters","png",sep="."),width=1200,height=800)
     dev.off()
    }
   }
  }
 
  

  ## step 2
  #dendrograma de expresiones medias
  dd<-0.5*(1-cor(t(mgenex)))
  hh<-hclust(d=as.dist(dd),method=method2)
  if(FALSE) plot(hh)

  #cutree sobre mean profiles
  mminClusterSize=2
  ddeepSplit     =ds2

  cct<-cutreeDynamic(hh,minClusterSize=mminClusterSize,distM=dd,deepSplit=ddeepSplit)
  names(cct)<-seq_along(cct)

  if(FALSE){
   xx<-aux<-c()
   for(i in c(1,2,3,6)){
    aux<-c(aux,i)
    ia<-which(ct%in%aux)
    xx<-c(xx,sqrt(sum(apply(z[ia,],2,function(x){var(x)})^2)))
   }  
  }
 
  tta<-table(cct)
  uua<-sort(unique(as.numeric(cct)))

  rcRatio <- 3/4
  ncol <- signif(sqrt(length(uua)/rcRatio)+0.5,0)
  nraw <- signif(length(uua)/ncol+.5,0)
  npanel<- ncol*nraw
  par(mar=c(2,2,2,2))
  if(bplot2) layout(matrix(1:npanel,ncol=ncol,byrow=TRUE))
  for(i in 1:length(uua)){
   iia<-which(cct==uua[i])
   ia<-which(ct%in%iia)
   aux<-length(ia)
   xx<-apply(z[ia,],2,mean)
   mgenex<-rbind(mgenex,xx)

   if(bwrite2==TRUE) {
    write.table(z[ia,],sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE,file=paste(fname,scode,"cluster",i,"csv",sep="."))
   }
   if(bplot2){
#     matplot(t(z[ia,]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4))
#     mtext(paste("Cluster",i," - ",aux),cex=0.7)
#     lines(1:3,xx[1:3],lwd=7,col="red")
#     lines(4:6,xx[4:6],lwd=7,col="violet")
    
    matplot(1:12 ,t(z[ia,1:12]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),xlim=c(1,24),ylim=range(t(z[ia,])))
    matplot(13:24,t(z[ia,13:24]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),add=TRUE)
    mtext(paste("Cluster",i," - ",aux),cex=0.7)
    
    lines(1:12,xx[1:12],lwd=7,col="gray")
    lines(13:24,xx[13:24],lwd=7,col=ccol[names(lttag)[iTemp]])
    
    
    if(i%%npanel==0) readline("Sacale una foto y Press <enter>\n")
   }
  }
  if(bwrite2){
   dev.copy(png,file=paste(fname,scode,"clusters","png",sep="."),width=1200,height=800)
   dev.off()
  }

  
#  varnorm<-function(x){(x-mean(x))/sd(x)}
#  xx<-t(apply(datX1,1,varnorm))
#  heatmap(xx) 
#   
  
 }

}


## Analisis Tipo 3:
## Encuentro cluwsters para cada temperatura
if(FALSE){

## Analisis Tipo 4 (usando edgeR):
## Para los genes que presentaron:
##    curvas diferentes con respecto a la curva control (22)
##   o
##    para genes que se movieron entre T27T12 (ver mas arriba)
##   o
##    genes con efecto aclimatacion T12 time1
## (y ademas maxFC >/log2FC), clusterear PERFILES COMPLETOS (en R^96)
## ATENNCION: 
## uso cpm o cpmvsn
if(FALSE){
 require(dynamicTreeCut)
 
 #voy a seguir con la normalizacion dada por varianceStabilizingTransformation 
 # pero uso la info de expresion diferencial de edgeR!!
 #assay(vsd)
 pvlim    <- 1e-4
 logfclim <- log2(2)

 ds1<-c(0,1,2,3,4)[1]
 ds2<-c(0,1,2,3,4)[3]
 method1<-method2<-"complete"
 bwrite1<-TRUE
 bwrite2<-TRUE
 bplot1 <-TRUE
 bplot2 <-TRUE
 scode<-paste(ds1,ds2,method1,sep="")

 ccol=c("gray","blue","violet","red")
 names(ccol)<-c(22,12,17,27)

 
 ## ELEGIR UNA OPCION
 if(FALSE){  #para genes que se movieron en algun tiempo para alguna Temp respecto a T22
  fname<-paste("cpmvsnFullProfile")
  genes<-c()
  for(iTemp in seq_along(lttag)){
  scol <- colnames(lttag[[iTemp]]$table)[grep("logFC",colnames(lttag[[iTemp]]$table))]
  maxfc<- apply(lttag[[iTemp]]$table[,scol],1,max)
  ig   <- which(lttag[[iTemp]]$table[,"FDR"]<pvlim & abs(maxfc)>logfclim) #ACA filtro
  genes<- union(genes,rownames(lttag[[iTemp]]$table)[ig])              
 }
 }
 if(FALSE){ #para genes que se movieron entre T27T12 (ver mas arriba)
  fname<-paste("T27T12genes.cpmvsnFullProfile")  
  tt<-ttag2712$table
  fc<-apply(tt[,6:17],1,function(x){sum(abs(x)>log2(2))>=2})
  genes<- rownames(tt)[tt$FDR<0.0001 & fc]
 }
 if(FALSE){  #genes con efecto aclimatacion T12 time1 # NO SE VE BIEN...NO USAR
  genes<-c()
  for(i in seq_along(lttagAcli)){
   tt<-lttagAcli[[i]][[1]]$table
   aa<-rownames(tt)[tt$FDR<0.001 & abs(tt$logFC)>log2(3)]
   cat(length(aa),"\n")
   genes<- c(genes,aa)
  } 
  saux<-"aclimatacion T12 time1"
  fname<-paste("Aclimat.cpmvsnFullProfile")
 }
 
  geneX <- t(apply(ycpm.vsn[genes,],1,function(x){return(apply(matrix(x,nrow=2),2,mean))}))
  a<-colnames(ycpm.vsn)  
  a<-a[grep("_A",a)];a<-sub("_A","",sub("at_","T",a));
  colnames(geneX)<-a
  rownames(geneX)<-genes
  
  
  #step1
  d<-1-(1+cor(t(geneX)))/2
  h<-hclust(d=as.dist(d),method=method1)

  minClusterSize=20
  deepSplit     =ds1

  ct<-cutreeDynamic(h,minClusterSize=minClusterSize,distM=d,deepSplit=deepSplit)
  ta<-table(ct)
  ua<-sort(unique(as.numeric(ct)))
  z <-t(apply(geneX,1,function(x){return((x-mean(x))/sd(x))}))
  
  #npanel<- signif(length(ua)/10,0)*10
  npanel<-length(ua)
  par(mar=c(2,2,2,2))
  
  #voy a tratar de ordenar todos estos
  #multidimensional scaling + grilla
  ncolWish<-10
  orderedClusters<-seq_along(unique(ua))
  if(npanel<ncolWish){
   ncol=npanel
  }else{
   aux=get_all_factors(length(ua))[[1]]
   ncol<-aux[which.min(abs(aux-ncolWish))] 
   if(length(aux)==2){
    ncol=ncolWish   #si es primo
    cat("me cacho que es primo!\n")
   } 
   if(length(aux)>2){
   xx<-c()
   for(i in 1:min(length(ua),npanel)){
    ia<-which(ct==ua[i])
    aux<-length(ia)
    xx<-rbind(xx,apply(z[ia,],2,mean)) 
   }
   dxx<-1-(1+cor(t(xx)))/2
   xy<-cmdscale(dxx,2)
   if(FALSE){
    plot(c(-1,1),c(-1,1),typ="n")
    text(xy[,1],xy[,2],1:70)
   } 
   rxy<-cbind(rank(xy[,1]),rank(xy[,2]))
   
  
   orderedClusters<-c()
   for(irow in (nrow(xy)/ncol):1){
    imin<-1+(irow-1)*ncol 
    imax<-ncol+(irow-1)*ncol 
    ii<-which(rxy[,2]>=imin & rxy[,2]<=imax)
    orderedClusters<-c(orderedClusters,ii[order(xy[ii,1],decreasing=FALSE)])
   }
  }
  }
  
  npanel2<- signif(length(ua)/ncol,0)*ncol + ncol* (length(ua)%%ncol!=0)
  if(bplot1) layout(matrix(1:npanel2,ncol=ncol,byrow=TRUE))
  mgenex<-c()
  for(i in 1:min(length(ua),npanel)){
   ia<-which(ct==orderedClusters[ua[i]])
   aux<-length(ia)
   xx<-apply(z[ia,],2,mean)
 
   mgenex<-rbind(mgenex,apply(z[ia,],2,mean))

   if(bwrite1==TRUE) {
    write.table(z[ia,],sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE,file=paste(fname,"precluster",i,"csv",sep="."))
    write.table(y$counts[rownames(z)[ia],],sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE,file=paste(fname,"precluster",i,"raw.csv",sep="."))
   }
   if(bplot1){

#     matplot(1:12 ,t(z[ia,1:12]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),xlim=c(1,24),ylim=range(t(z[ia,])))
#     matplot(13:24,t(z[ia,13:24]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),add=TRUE)
#     
#     mtext(paste("preCluster",i," - ",aux),cex=0.7)
#     lines(1:12,xx[1:12],lwd=7,col="gray")
#     lines(13:24,xx[13:24],lwd=7,col=ccol[names(lttag)[iTemp]])
    matplot(t(z[ia,]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),xlim=c(1,48),ylim=range(t(z[ia,])))
    if(FALSE){
     layout(matrix(1:3,1,3))
     matplot(t(z[ia,]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),xlim=c(1,48),ylim=range(t(z[ia,])))
     matplot(t(y$counts[rownames(z)[ia],]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),ylim=range(y$counts[rownames(z)[ia],]))
     matplot(t(y$counts[rownames(z)[ia],])+0.01,typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),ylim=c(1,max(y$counts[rownames(z)[ia],])),log="y")
    }    

    mtext(paste("preCluster",i," - ",aux),cex=0.7)
    lines(xx,lwd=4,col="green")
    abline(v=12.5*c(1:3),col="black")

    if(i%%npanel==0) readline("Sacale una foto y Press <enter>\n")
    if(bwrite1){
     if(fname=="Aclimat.cpmvsnFullProfile"){
      dev.copy(png,file=paste(fname,scode,"preclusters","png",sep="."),width=1200,height=200)
     }else{
      dev.copy(png,file=paste(fname,scode,"preclusters","png",sep="."),width=1200,height=800)
     }
     dev.off()
    }
   }
  }
 
  

  ## step 2
  #dendrograma de expresiones medias
  dd<-0.5*(1-cor(t(mgenex)))
  hh<-hclust(d=as.dist(dd),method=method2)
  if(FALSE) plot(hh)

  #cutree sobre mean profiles
  mminClusterSize=2
  ddeepSplit     =ds2

  cct<-cutreeDynamic(hh,minClusterSize=mminClusterSize,distM=dd,deepSplit=ddeepSplit)
  names(cct)<-seq_along(cct)

  if(FALSE){
   xx<-aux<-c()
   for(i in c(1,2,3,6)){
    aux<-c(aux,i)
    ia<-which(ct%in%aux)
    xx<-c(xx,sqrt(sum(apply(z[ia,],2,function(x){var(x)})^2)))
   }  
  }
 
  tta<-table(cct)
  uua<-sort(unique(as.numeric(cct)))

  rcRatio <- 3/4
  

  ncol <- signif(sqrt(length(uua)/rcRatio)+0.5,0)
  nraw <- signif(length(uua)/ncol+.5,0)
  npanel<- ncol*nraw
  par(mar=c(2,2,2,2))
  if(bplot2) layout(matrix(1:npanel,ncol=ncol,byrow=TRUE))
  for(i in 1:length(uua)){
   iia<-which(cct==uua[i])
   ia<-which(ct%in%orderedClusters[iia])
   aux<-length(ia)
   xx<-apply(z[ia,],2,mean)
   mgenex<-rbind(mgenex,xx)

   if(bwrite2==TRUE) {
    write.table(z[ia,],sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE,file=paste(fname,scode,"cluster",i,"csv",sep="."))
    write.table(y$counts[rownames(z)[ia],],sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE,file=paste(fname,scode,"cluster",i,"raw.csv",sep="."))   
   }
   if(bplot2){
#     matplot(t(z[ia,]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4))
#     mtext(paste("Cluster",i," - ",aux),cex=0.7)
#     lines(1:3,xx[1:3],lwd=7,col="red")
#     lines(4:6,xx[4:6],lwd=7,col="violet")
    
#     matplot(1:12 ,t(z[ia,1:12]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),xlim=c(1,24),ylim=range(t(z[ia,])))
#     matplot(13:24,t(z[ia,13:24]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),add=TRUE)
#     mtext(paste("Cluster",i," - ",aux),cex=0.7)
#    
#     lines(1:12,xx[1:12],lwd=7,col="gray")
#     lines(13:24,xx[13:24],lwd=7,col=ccol[names(lttag)[iTemp]])

    matplot(t(z[ia,]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),xlim=c(1,48),ylim=range(t(z[ia,])))
    mtext(paste("Cluster",i," - ",aux),cex=0.7)
    
    lines(xx,lwd=4,col="red")
    abline(v=12.5*c(1:3),col="black")
    
    
    if(i%%npanel==0) readline("Sacale una foto y Press <enter>\n")
   }
  }
  if(bwrite2){
   dev.copy(png,file=paste(fname,scode,"clusters","png",sep="."),width=1200,height=800)
   dev.off()
  }

  
#  varnorm<-function(x){(x-mean(x))/sd(x)}
#  xx<-t(apply(datX1,1,varnorm))
#  heatmap(xx) 
#   
  

}

 
 # Detecto genes que se mueven a c/temperatura
 if(FALSE){
  
  bcalcul=FALSE
  if(bcalcul){
   vTemp<-as.character(c(12,17,22,27))
   lttagTemp<-list()
   for(iT in seq_along(vTemp)){
    cat(paste("T",vTemp[iT]),"\n")
    i22<-which(pheno$temperature==vTemp[iT]) 
    design22 <- model.matrix(~ time, data=pheno[i22,])
    n22 <- minCountPM*N/(1e6*ncol(reads[,i22]))
    minCountPM = 5
    keep22 <- apply(reads,1,sum)>=n22
#   keep   keep22
# 1706 21409 0
    y22 <- DGEList(counts=reads[keep22,as.character(pheno$sample[i22])],
                 group=pheno$hpi[i22],genes=meta[keep22,])
 
    y22 <- calcNormFactors(y22,method=c("TMM","RLE","upperquartile")[1])
    y22$samples
    y22 <- estimateGLMTrendedDisp(y22,design22)
    y22 <- estimateGLMTagwiseDisp(y22,design22)
    fit22<-glmFit(y22,design22)
    lrt22<-glmLRT(fit22,coef=2:12)
    a<-topTags(lrt22,n=length(lrt$table[,1])) 
    lttagTemp[[vTemp[iT]]]<-a
   # sum(a$table[,"FDR"]<0.0001)
   }
   if(FALSE) save(lttagTemp,file="lttagTemp.Rdata")
  }else{
   load("lttagTemp.Rdata")
  }
 }


 require(dynamicTreeCut)
 #Empiezo evaluando las maneras de estabilizar varianzas que aparecen en DESeq2
 #analisis DESeq2
 library(DESeq2)
 library("BiocParallel")
 library("vsn")
 bcalcul=FALSE
 if(bcalcul){
  register(MulticoreParam(4))
  dds<-DESeqDataSetFromMatrix(countData = reads[keep,], 
                             colData= pheno,
                             design = ~ time + time:temperature)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)

  rld<-rlog(dds,fast=TRUE)
  vsd <- varianceStabilizingTransformation(dds)
  if(FALSE) save(dds,rld,vsd,file="deseq2.96.Rdata")
 }else{
  (load("deseq2.96.Rdata"))
 }
 #aver...
 if(FALSE){
  par(mfrow=c(1,3))
  notAllZero <- (rowSums(counts(dds))>0)
  meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1))
  mtext("log2")
  meanSdPlot(assay(rld[notAllZero,]))
  mtext("rlog")
  meanSdPlot(assay(vsd[notAllZero,]))
  mtext("vsd")
 }
 
 #voy a seguir con la normalizacion dada por varianceStabilizingTransformation 
 #assay(vsd)
 pvlim    <- 1e-4
 logfclim <- log2(2)

 ds1<-c(0,1,2,3,4)[1]
 ds2<-c(0,1,2,3,4)[3]
 method1<-method2<-"complete"
 bwrite1<-FALSE
 bwrite2<-TRUE
 bplot1 <-FALSE
 bplot2 <-TRUE
 scode<-paste(ds1,ds2,method1,sep="")

 ccol=c("gray","blue","violet","red")
 names(ccol)<-c(22,12,17,27)
 temps<-c("12","17","22","27")
 lclus<-list()
 for(iTemp in seq_along(temps)){
 
  sTemp <- temps[iTemp]
  fname<-paste("T",sTemp,sep="")
 
  x<-lttagTemp[[sTemp]]
  scol <- colnames(x$table)[grep("logFC",colnames(x$table))]
  maxfc<- apply(x$table[,scol],1,max)
  ig   <- which(x$table[,"FDR"]<pvlim & abs(maxfc)>logfclim)
  genes<- rownames(x$table)[ig]  
#   
#   
#   if(sTemp=="22"){
#    (load("yfit96.Rdata"))
#    lrt<- glmLRT(fit,coef=2:12)
#    x  <- topTags(lrt,n=length(lrt$table[,1])) 
#    scol <- colnames(x$table)[grep("logFC",colnames(x$table))]
#    maxfc<- apply(x$table[,scol],1,max)
#    ig   <- which(x$table[,"FDR"]<pvlim & abs(maxfc)>logfclim)
#    genes<- rownames(x$table)[ig]               
#   }else{  
#    scol <- colnames(lttag[[sTemp]]$table)[grep("logFC",colnames(lttag[[sTemp]]$table))]
#    maxfc<- apply(lttag[[sTemp]]$table[,scol],1,max)
#    ig   <- which(lttag[[sTemp]]$table[,"FDR"]<pvlim & abs(maxfc)>logfclim)
#    genes<- rownames(lttag[[sTemp]]$table)[ig]              
#   }
  
  #columnas  T de analisis
  
  iT <- which(data.frame(colData(vsd))[,"temperature"]==sTemp)
    
  geneX <- t(apply(assay(vsd)[genes,c(iT)],1,function(x){return(apply(matrix(x,nrow=2),2,mean))}))
  rownames(geneX)<-genes
  
  
  #step1
  d<-1-(1+cor(t(geneX)))/2
  h<-hclust(d=as.dist(d),method=method1)

  minClusterSize=20
  deepSplit     =ds1

  ct<-cutreeDynamic(h,minClusterSize=minClusterSize,distM=d,deepSplit=deepSplit)
  ta<-table(ct)
  ua<-sort(unique(as.numeric(ct)))
  z <-t(apply(geneX,1,function(x){return((x-mean(x))/sd(x))}))
  colnames(z)<-c(paste("T",sTemp,"_",c(1:12),sep=""))
 
  npanel<- signif(length(ua)/10,0)*10
  par(mar=c(2,2,2,2))
  if(bplot1) layout(matrix(1:npanel,ncol=10,byrow=TRUE))
  mgenex<-c()
  for(i in 1:min(length(ua),npanel)){
   ia<-which(ct==ua[i])
   aux<-length(ia)
   xx<-apply(z[ia,],2,mean)
 
   mgenex<-rbind(mgenex,apply(z[ia,],2,mean))

   if(bwrite1==TRUE) {
    write.table(z[ia,],sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE,file=paste(fname,"precluster",i,"csv",sep="."))
   }
   if(bplot1){

    matplot(1:12 ,t(z[ia,1:12]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),xlim=c(1,24),ylim=range(t(z[ia,])))
    matplot(13:24,t(z[ia,13:24]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),add=TRUE)
    
    mtext(paste("preCluster",i," - ",aux),cex=0.7)
    lines(1:12,xx[1:12],lwd=7,col="gray")
    lines(13:24,xx[13:24],lwd=7,col=ccol[sTemp])

    if(i%%npanel==0) readline("Sacale una foto y Press <enter>\n")
    if(bwrite2){
     dev.copy(png,file=paste(fname,scode,"preclusters","png",sep="."),width=1200,height=800)
     dev.off()
    }
   }
  }
 
  

  ## step 2
  #dendrograma de expresiones medias
  dd<-0.5*(1-cor(t(mgenex)))
  hh<-hclust(d=as.dist(dd),method=method2)
  if(FALSE) plot(hh)

  #cutree sobre mean profiles
  mminClusterSize=2
  ddeepSplit     =ds2

  cct<-cutreeDynamic(hh,minClusterSize=mminClusterSize,distM=dd,deepSplit=ddeepSplit)
  names(cct)<-seq_along(cct)

  if(FALSE){
   xx<-aux<-c()
   for(i in c(1,2,3,6)){
    aux<-c(aux,i)
    ia<-which(ct%in%aux)
    xx<-c(xx,sqrt(sum(apply(z[ia,],2,function(x){var(x)})^2)))
   }  
  }
 
  tta<-table(cct)
  uua<-sort(unique(as.numeric(cct)))

  rcRatio <- 3/4
  ncol <- signif(sqrt(length(uua)/rcRatio)+0.5,0)
  nraw <- signif(length(uua)/ncol+.5,0)
  npanel<- ncol*nraw
  par(mar=c(2,2,2,2))
  if(bplot2) layout(matrix(1:npanel,ncol=ncol,byrow=TRUE))
  for(i in 1:length(uua)){
   iia<-which(cct==uua[i])
   ia<-which(ct%in%iia)
   aux<-length(ia)
   xx<-apply(z[ia,],2,mean)
   mgenex<-rbind(mgenex,xx)

   if(bwrite2==TRUE) {
    write.table(z[ia,],sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE,file=paste(fname,scode,"cluster",i,"csv",sep="."))
   }
   if(bplot2){
#     matplot(t(z[ia,]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4))
#     mtext(paste("Cluster",i," - ",aux),cex=0.7)
#     lines(1:3,xx[1:3],lwd=7,col="red")
#     lines(4:6,xx[4:6],lwd=7,col="violet")
    
    matplot(1:12 ,t(z[ia,1:12]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),xlim=c(1,12),ylim=range(t(z[ia,])))
    #matplot(13:24,t(z[ia,13:24]),typ="l",ylab="",xlab="",col=rgb(0.2,0.2,0.2,0.4),add=TRUE)
    mtext(paste("T",sTemp,"- Cluster",i," - ",aux),cex=0.7)
    
    lines(1:12,xx[1:12],lwd=7,col="gray")
    #lines(13:24,xx[13:24],lwd=7,col=ccol[sTemp])
    
    
    if(i%%npanel==0) readline("Sacale una foto y Press <enter>\n")
   }
   lclus[[paste("T",temps[iTemp],"c",i,sep="")]]<-rownames(z)[ia]

  }
  if(bwrite2){
   dev.copy(png,file=paste(fname,scode,"clusters","png",sep="."),width=1200,height=800)
   dev.off()
  }

  
#  varnorm<-function(x){(x-mean(x))/sd(x)}
#  xx<-t(apply(datX1,1,varnorm))
#  heatmap(xx) 
#   
  
 }
 save(lclus,file="lclus.Rdata")
}



if(FALSE){
 #elijo algunos
 tselection=1
 ttemp<-c("12","17","27")
 temp<-paste("lrt",ttemp[tselection],sep=".")
 imres<-apply(mres,1,function(x){all(x==(ttemp[tselection]==ttemp))})   ##PENSADO PARA  los 3 CASOS PUROS 
 genes<-names(imres)[imres]
 assign("lrt",get(temp))
 assign("mpval",get(paste("mpval",ttemp[tselection],sep=".")))
 
 #yy<-log(reads[rownames(lrt$table)[ires],]+.5)
 yy<-log(reads[genes,]+.5)
 yy<-reads[genes,]

 a<-apply(yy,1,function(x){return(apply(matrix(x,nrow=2),2,mean))})
 layout(matrix(1:2,1,2),width=c(.6,.4))
 for(i in seq_along(a[1,])){
  plot(a[,i],main=genes[i],typ="b",col="red",pch=20)
  vl<-which(diff(as.numeric(pheno$temperature))!=0)/2+.5
  abline(v=c(vl),col="gray",lty=2)
  
  ma<-matrix(a[,i],nrow=12)
  matplot(ma[,c(3,1,2,4)],typ="l",col=c("gray","blue","violet","red"),lty=1,lwd=c(15,5,5,5))
  mtext(genes[i])
  abline(v=6.5,lty=2,col="gray")

  colnames(lrt$table)
  ipv <- mpval[genes[i],]<0.01
  ilfc<- abs(lrt$table[genes[i],1:12]) > log2(1.5)
  ii  <- ipv & ilfc
  ccex<-0.5 + ii *4
  points(rep(1:12),a[(1:12)+(tselection-1)*12,genes[i]],pch=20,cex=ccex)
  
  readline(paste(i,"/",ncol(a),"  press <enter> to continue"))
  
 }

 

  ## comparacion edgeR - cuffdiff
  u<-union(ires,ide)
  vennDiagram(cbind(edgeR=u%in%ires,cuffdidff=u%in%ide))

  plot(M[u],lrt$table[u,"logFC"],xlab="cuffdiff",ylab="edgeR")
  ired<-ires
  ired<-ires[!(ires%in%ide)]
  points(M[ired],lrt$table[ired,"logFC"],pch=20,col=3)
  ired<-ide[!(ide%in%ires)]
  points(M[ired],lrt$table[ired,"logFC"],pch=20,col=2)
  abline(0,1)
  legend("topleft",c("only edgeR","only cuffdiff"),col=c(3,2),pch=20,inset=0.02)


 ## heatmap
 gnames<-rownames(datX1)[ires]

 
 #datX1 <-log2(cpm(y,TRUE)[gnames,]+1e-6)
 datX1<-log2(myCPM(y)[gnames,]+1e-6)
 dd1<-1-(1+cor(t(datX1)))/2
 hc1<-hclust(as.dist(dd1))

 varnorm<-function(x){(x-mean(x))/sd(x)}
 xx<-t(apply(datX1,1,varnorm))
 heatmap(xx) 




 ##
 ## genero FIGURA figura de perfiles de expresion
 fname0<-paste(dir,"/clust.",iclus,".png",sep="")
 cclust<-clus0$labels 
#  >  par("mar")
# [1] 5.1 4.1 4.1 2.1
 par(mar=c(2,2,2,2)) 
 layout(matrix(1:15,byrow=TRUE,5,3))
 for(iclus in 1:15){
  ii<-rownames(datX1)[which(cclust==iclus)]
  nn   <- ncol(datX1)
  nn_2 <- ncol(datX1)/2
  
  aa1<-matrix(apply(matrix(t(datX1[ii,1:nn_2]),byrow=TRUE,ncol=nrep),1,mean),byrow=TRUE,ncol=nn_2/nrep)
  aa2<-matrix(apply(matrix(t(datX1[ii,nn_2+(1:nn_2)]),byrow=TRUE,ncol=nrep),1,mean),byrow=TRUE,ncol=nn_2/nrep)
 
  a1 <- t(apply(cbind(aa1,aa2),1,function(x){x[1:(nn_2/nrep)]/max(abs(x))}))
  a2 <- t(apply(cbind(aa2,aa1),1,function(x){x[1:(nn_2/nrep)]/max(abs(x))}))
  
  
  b1 <-apply(a1,2,mean)
  b2 <-apply(a2,2,mean)
  b  <-rbind(b1,b2)
 
  xxlab<-yylab<-""
  if(iclus>12) xxlab="time"
  if(iclus%in%c(1,4,7,10,13)) yylab="standarized log2CPM"
  matplot(1:(nn/nrep),rbind(t(a1),t(a2)),typ="b",lty=1,pch=20,
       #main=paste("cluster",iclus," size",length(ii)),
       ylim=c(-1,1),xlab=xxlab,ylab=yylab,axes=FALSE)  
     #axis(1,1:(nn/nrep),rep(c(2,6,10,14,18,22),2))
  mtext(paste("cluster",iclus," size",length(ii)),cex=0.8 )
  if(iclus>12){     
     axis(1,1:(nn_2/nrep),c(2,6,10,14,18,22),cex=0.8)#2)
     axis(1,(1:(nn_2/nrep))+(nn_2/nrep),c(2,6,10,14,18,22),cex=0.6)
  }   
  if(iclus%in%c(1,4,7,10,13))   axis(2)         
   points(1:(nn_2/nrep),b[1,],typ="b",lwd=5,col=1,pch=18)
   points((nn_2/nrep)+(1:(nn_2/nrep)),b[2,],typ="b",lwd=5,col=2)  
   abline(v=nn_2/nrep+.5)
   box()
 }
 if(FALSE) dev.copy(pdf,file="clusterProfiles.pdf",width=8,height=11)


}



##
## DEseq
{
library(DESeq)

a<-matrix(unlist(strsplit(as.character(cond),"_",fixed=TRUE)),byrow=TRUE,ncol=2)
a<-data.frame(genomic=a[,1],time=a[,2])
rownames(a)<-colnames(counts.table)
a$time<-factor(as.character(a$time),levels=as.character(c(0,seq(2,22,4))))

cds <- newCountDataSet(counts.table,cond,phenoData=new("AnnotatedDataFrame",data=a))

#normalization
cds <- estimateSizeFactors(cds)
plot(sizeFactors(cds),typ="b",ylim=c(0,2),main="sizeFactors")
abline(h=1,col="gray")

#variance estimation
cds <- estimateDispersions(cds,method="pooled",sharingMode="maximum",fitType="parametric")

#check
if(FALSE){
 plotDispEsts(cds)
 disp.deseq<-head(fData(cds))
} 

fit1 <- fitNbinomGLMs(cds,count ~ genomic * time)
fit0 <- fitNbinomGLMs(cds,count ~ genomic + time)

pvalsGLM <- nbinomGLMTest(fit1,fit0)
padj <- p.adjust(pvalsGLM,method="BH")
names(padj)<-rownames(fit1)

#estos mostraron interaccion genomic:time
ipv <- which(padj < 0.01)  
names(ipv)

##
## variance stabilization for visualizacion & clustering
##
cdsBlind<-estimateDispersions(cds,method="blind")
vsd     <-varianceStabilizingTransformation(cdsBlind)
if(FALSE){
 trsf <- function(x){
  log((x+sqrt(x*x+1))/2)
 }

 library("vsn")
 par(mfrow=c(1,2))
 notAllZero <- (rowSums(counts(cds))>0)
 meanSdPlot(log2(counts(cds)[notAllZero,]+1), ylim=c(0,2.5))
 meanSdPlot(vsd[notAllZero,], ylim=c(0,2.5))
}
}

}