#
###################################
# Useful functions
###################################
#
caminho<-"/Users/Daiane/Documents/Daiane/Doutorado/03.Projeto/04. Mapeamento QTL pedigree/GAW17/"
#
install.package(Rccp)
library(Rcpp)
sourceCpp(file=paste(caminho,"cppKruskall.cpp",sep=""))
################
# Crámeer's coefficient
#
cv.test = function(x,y) {
  CV = sqrt(chisq.test(x, y, correct=FALSE)$statistic /
    (length(x) * (min(length(unique(x)),length(unique(y))) - 1)))
  print.noquote("CramÈr V / Phi:")
  return(as.numeric(CV))}
#
################
### draw a discrete value
#
rDiscreta<-function(p){
 u<-runif(1)
 P<-cumsum(p)
 val<-sum(P<u)+1
 val}
#
################
### Haldane function (in Morgans)
#
Haldane<-function(dist){
 # dist = dist‚ncia em cM para calcular a taxa de recombinaÁ„o
 recomb<-0.5*(1-exp(-2*dist))
 recomb}
#
################
### identifica a origem de cada alelo do genótipo
#
identalel_pai_mae<-function(genpai,genmae,genfilho){
 # Alelo = 2 representa alelo dominante e Alelo = 1 representa alelo recessivo
 if (genfilho==1) {
  afp<-2
  afm<-2} else {
  if (genfilho==-1){
   afp<-1
   afm<-1} else {
   if (genpai==1) {
    afp<-2
    afm<-1} else {
    if (genpai==-1) {
     afp<-1
     afm<-2} else {
     if (genmae==1) {
      afp<-1
      afm<-2} else {
      if (genmae==-1) {
       afp<-2
       afm<-1} else {
       aux<-rbinom(1,1,0.5)
       afp<-(2^aux)*(1^(1-aux))
       afm<-(2^(1-aux))*(1^aux)}}}}}}
 list(afp,afm)}
#
################
### calcula a probabilidade de seleção dos marcadores via Kruskal-wallis e seleciona um marcador
#
prob.selec.marc<-function(dados,res1,qtls,loc.marc){
  # dados = matrix which has the phenotype in first column and markers' genotype in remaining column
  # qtls = vector with QTLs' positions (in CMs) in current model
  if (length(qtls)==0) krusk<-kruskalcpp(dados) else krusk<-kruskalcpp(cbind(res1,dados[,2:ncol(dados)]))  
  if (sum(qtls<=loc.marc[2])>0) krusk[1]<-0
  for (i in 2:(length(loc.marc)-1)) if ((sum(qtls>=loc.marc[i-1] & qtls<=loc.marc[i])>0)&(sum(qtls>=loc.marc[i] & qtls<=loc.marc[i+1])>0)) krusk[i]<-0
  if (sum(qtls>=loc.marc[length(loc.marc)-1])>0) krusk[length(loc.marc)]<-0
  prob<-krusk/sum(krusk)
  marc<-rDiscreta(prob)
  list(marc,log(prob[marc]),krusk)}
#
################
### calcula a probabilidade de exclusão de um QTL via soma dos seus efeitos e seleciona um QTL
#
prob.excl.marc<-function(num.QTLs,vet.coef){
 efeitos<-numeric(num.QTLs)
 for (i in 1:num.QTLs) efeitos[i]<-1/sum(c(abs(vet.coef[(2*i),1]),abs(vet.coef[((2*i)+1),1])))
 proba<-efeitos/sum(efeitos)
 qtl<-rDiscreta(proba)
 list(qtl,log(proba[qtl]))}
#
################
### calcula a posição do QTL aos redores do marcador escolhido
#
pos.qtl<-function(qtls,loc.marc,marc,krusk){
 if (marc==1 | sum(qtls>=loc.marc[marc-1] & qtls<=loc.marc[marc])>0){
  marc.ini<-marc
  marc.fim<-marc+1} else {
  if (marc==length(loc.marc) | sum(qtls>=loc.marc[marc] & qtls<=loc.marc[marc+1])>0){
   marc.ini<-marc-1
   marc.fim<-marc} else {
   if ((sum(qtls>=loc.marc[marc-1] & qtls<=loc.marc[marc])==0)&(sum(qtls>=loc.marc[marc] & qtls<=loc.marc[marc+1])==0)){
    marc.ini<-marc-1
    marc.fim<-marc+1}}}
 marc.pos<-(loc.marc[marc.ini:marc.fim]-loc.marc[marc.ini])/(loc.marc[marc.fim]-loc.marc[marc.ini])
 estas<-krusk[marc.ini:marc.fim]
 mi.pos<-(t(marc.pos)%*%estas)/sum(estas)
 parA<-mi.pos/(1-mi.pos)
 parB<-1
 uger<-rbeta(1,parA,parB)
 loc.qtl<-loc.marc[marc.ini]+(loc.marc[marc.fim]-loc.marc[marc.ini])*uger
 dens<-dbeta(uger,parA,parB,log = TRUE)
 list(loc.qtl,dens)}
#
dens.pos.qtl<-function(QTLmg,qtls,loc.marc,marc,krusk){
 if (marc==1 | sum(qtls>=loc.marc[marc-1] & qtls<=loc.marc[marc])>0){
  marc.ini<-marc
  marc.fim<-marc+1} else {
  if (marc==length(loc.marc) | sum(qtls>=loc.marc[marc] & qtls<=loc.marc[marc+1])>0){
   marc.ini<-marc-1
   marc.fim<-marc} else {
   if ((sum(qtls>=loc.marc[marc-1] & qtls<=loc.marc[marc])==0)&(sum(qtls>=loc.marc[marc] & qtls<=loc.marc[marc+1])==0)){
    marc.ini<-marc-1
    marc.fim<-marc+1}}}
 marc.pos<-(loc.marc[marc.ini:marc.fim]-loc.marc[marc.ini])/(loc.marc[marc.fim]-loc.marc[marc.ini])
 estas<-krusk[marc.ini:marc.fim]
 mi.pos<-(t(marc.pos)%*%estas)/sum(estas)
 parA<-mi.pos/(1-mi.pos)
 parB<-1
 uger<-(QTLmg-loc.marc[marc.ini])/(loc.marc[marc.fim]-loc.marc[marc.ini])
 dens<-dbeta(uger,parA,parB)
 dens}
#
################
### calcula a posição dos QTLs uniformemente
#
priori.loc.qtls<-function(num.qtls,loc.marc){
 loc.qtls<-numeric(num.qtls)
 num.marc<-length(loc.marc)
 prob.loc<-0
 if (num.qtls>0){
 for (i in 1:num.qtls){
   loc.qtls[i]<-runif(1,min=loc.marc[sum(loc.marc<=loc.qtls[i-1])+1],max=loc.marc[num.marc-(num.qtls-i)])
   prob.loc<-prob.loc+dunif(loc.qtls[i],min=loc.marc[sum(loc.marc<=loc.qtls[i-1])+1],max=loc.marc[num.marc-(num.qtls-i)],log=TRUE)}}
 list(loc.qtls,prob.loc)}
#
dens.priori.loc.qtls<-function(loc.qtls,loc.marc){
 num.qtls<-length(loc.qtls)
 num.marc<-length(loc.marc)
 prob.loc<-0
 if (num.qtls>0){
 for (i in 1:num.qtls){
   prob.loc<-prob.loc+dunif(loc.qtls[i],min=loc.marc[sum(loc.marc<=loc.qtls[i-1])+1],max=loc.marc[num.marc-(num.qtls-i)],log=TRUE)}}
 prob.loc}
#
################
### calculam a probabilidade do genótipo de um QTL com base no genótipo dos marcadores flanqueadores
#
# gen = 1 se marcador homozigoto dominante
# gen = 0 se marcador heterozigoto
# gen = -1 se marcador homozigoto recessivo
gen.igual<-function(recomb,gen1) if (gen1==0) {(recomb*recomb)+(1-recomb)**2} else {(1-recomb)**2}
gen.dif1<-function(recomb) 2*(1-recomb)*recomb
gen.dif2<-function(recomb) recomb*recomb
#
calc.prob.gen<-function(gen1,gen2,recomb){
 if (abs(gen1-gen2)==0) {prob<-gen.igual(recomb,gen1)} else {
  if (abs(gen1-gen2)==1) {if (gen1==0) {prob<-gen.dif1(recomb)/2} else {prob<-gen.dif1(recomb)}}
  else {prob<-gen.dif2(recomb)}}
 prob}
#
################
### gera valores da posteriori de sigma 2
#
poster.sigma2<-function(neta.a,neta.b,residuos){
 alpha<-(length(residuos)/2)+neta.a
 beta<-(sum(residuos^2)/2)+neta.b
 sigma2<-1/(rgamma(1,alpha,beta))
 dens<-dgamma((1/sigma2),alpha,beta,log = TRUE)
 list(sigma2,dens)}
#
################
### gera valores da posteriori da média geral mi
#
poster.mi<-function(media.mi,sigma2.mi,sigma2,residuos,mi.anterior){
 res.mi<-residuos+mi.anterior
 quo<-(length(residuos)/sigma2)+(1/sigma2.mi)
 media<-((sum(res.mi)/sigma2)+(media.mi/sigma2.mi))/quo
 variancia<-1/quo
 mi<-rnorm(1,media,sqrt(variancia))
 dens<-dnorm(mi,media,sqrt(variancia),log = TRUE)
 list(mi,dens)}
#
################
### gera valores da posteriori do efeito aditivo alphaj
#
poster.alpha<-function(media.alpha,sigma2.alpha,sigma2,residuos,alpha.anterior,gen.QTL){
 res.alpha<-residuos+(alpha.anterior*gen.QTL)
 quo<-(sum(gen.QTL^2)/sigma2)+(1/sigma2.alpha)
 media<-((sum(gen.QTL*res.alpha)/sigma2)+(media.alpha/sigma2.alpha))/quo
 variancia<-1/quo
 alpha<-rnorm(1,media,sqrt(variancia))
 dens<-dnorm(alpha,media,sqrt(variancia),log = TRUE)
 list(alpha,dens)}
#
################
### gera valores da posteriori do efeito de dominância deltaj
#
poster.delta<-function(media.delta,sigma2.delta,sigma2,residuos,delta.anterior,gen.QTL.dom){
 res.delta<-residuos+(delta.anterior*gen.QTL.dom)
 quo<-(sum(gen.QTL.dom^2)/sigma2)+(1/sigma2.delta)
 media<-((sum(gen.QTL.dom*res.delta)/sigma2)+(media.delta/sigma2.delta))/quo
 variancia<-1/quo
 delta<-rnorm(1,media,sqrt(variancia))
 dens<-dnorm(delta,media,sqrt(variancia),log = TRUE)
 list(delta,dens)}
#
################
### seleciona entre a inclusão ou exclusão de um QTL
#
dec.sp.mg<-function(num.QTLs,num.marc){
 if (num.QTLs==0) {psplit<-1; pmerge<-0} else {if (num.QTLs==(num.marc-1)) {psplit<-0;pmerge<-1} else {psplit<-pmerge<-1/2}}
 prob<-c(psplit,pmerge)
 ind.sp.mg<-rDiscreta(prob)
 list(ind.sp.mg,log(prob))}
#
#################
### gera inclusão de um QTL no modelo e calcula a probabilidade dessa inclusão
#
gera.inclusao.QTL<-function(dados,residuos,mat.delinea,vet.coef,pos.qtls,loc.marc,sigma2.vig,alpha.vig,delta.vig,media.mi,sigma2.mi,media.alpha,sigma2.alpha,media.delta,sigma2.delta,neta.a,neta.b,
mat.ped.pais,founders,nonfounders,alelo.paterno,alelo.materno,Sp,Sm,mat.alelo.qtls){
 #
 marcador<-prob.selec.marc(dados,residuos,pos.qtls,loc.marc)
 qtl<-pos.qtl(pos.qtls,loc.marc,marcador[[1]],marcador[[3]])
 #
 dQTL<-qtl[[1]]
 Marc1<-sum(loc.marc<=dQTL)
 Marc2<-length(loc.marc)-(sum(loc.marc>=dQTL)-1)
 dM1<-loc.marc[Marc1]
 dM2<-loc.marc[Marc2]
 r12<-Haldane(abs(dM2-dM1))
 r1<-Haldane(abs(dQTL-dM1))
 r2<-Haldane(abs(dQTL-dM2))
 matrec1<-matrix(c(1-r1,r1,r1,1-r1),2,2,byrow=TRUE)
 matrec2<-matrix(c(1-r2,r2,r2,1-r2),2,2,byrow=TRUE)
 #
 ## gera o genótipo e define o alelo paterno e materno do QTL dos founders
 #
 matriz.prob.gen.QTL<-NULL
 probQTL<-numeric(3)
 for (j in 1:length(founders)){
  gen<-c(-1,0,1)
  for (i in 1:3) probQTL[i]<-(calc.prob.gen(dados[j,(Marc1+1)],gen[i],r1)*calc.prob.gen(gen[i],dados[j,(Marc2+1)],r2))/calc.prob.gen(dados[j,(Marc1+1)],dados[j,(Marc2+1)],r12)
  matriz.prob.gen.QTL<-rbind(matriz.prob.gen.QTL,probQTL)}
 #
 dados.QTL<-matrix(0,length(founders),1)
 for (i in 1:length(founders)) dados.QTL[i,1]<-gen[rDiscreta(matriz.prob.gen.QTL[i,])]
 monoto<-length(table(dados.QTL[,1]))
 if (monoto==1){
  indiv<-sample(seq(1:length(founders)),1)
  if (dados.QTL[indiv,1]==-1) dados.QTL[indiv,1]<-0
  if (dados.QTL[indiv,1]==0) dados.QTL[indiv,1]<-1
  if (dados.QTL[indiv,1]==1) dados.QTL[indiv,1]<-0}
 #
 ale.pai<-matrix(99,nrow=length(founders),1)
 ale.mae<-matrix(99,nrow=length(founders),1)
 probale<-numeric(2)
 #
 for (k in 1:length(founders)){
  if (dados.QTL[k,1]==-1){
   ale.mae[k,1]<-1
   ale.pai[k,1]<-1}
  if (dados.QTL[k,1]==1){
   ale.mae[k,1]<-2
   ale.pai[k,1]<-2}
  if (dados.QTL[k,1]==0){
   for (ale in 1:2) probale[ale]<-matrec1[alelo.paterno[k,Marc1],ale]*matrec2[ale,alelo.paterno[k,Marc2]]
   probale<-probale/sum(probale)
#   aux2<-which(probale==max(probale))-1
   aux2<-rDiscreta(probale)-1
   ale.mae[k,1]<-1^(aux2)*2^(1-aux2)
   ale.pai[k,1]<-2^(aux2)*1^(1-aux2)}}
 dados.alelo<-cbind(ale.pai,ale.mae)
 #
 ## gera os alelos do QTL dos nonfounders
 #
 matriz.prob.S.pat<-NULL
 matriz.prob.S.mat<-NULL
 probSpat<-numeric(2)
 probSmat<-numeric(2)
 for (j in 1:length(nonfounders)){
  for (i in 1:2){
   probSmat[i]<-matrec1[Sm[j,Marc1],i]*matrec2[i,Sm[j,Marc2]]
   probSpat[i]<-matrec1[Sp[j,Marc1],i]*matrec2[i,Sp[j,Marc2]]}
  matriz.prob.S.mat<-rbind(matriz.prob.S.mat,(probSmat/sum(probSmat)))
  matriz.prob.S.pat<-rbind(matriz.prob.S.pat,(probSpat/sum(probSpat)))}
 Smqtl<-numeric(length(nonfounders))
 Spqtl<-numeric(length(nonfounders))
 for (j in 1:length(nonfounders)){
  Smqtl[j]<-rDiscreta(matriz.prob.S.mat[j,])
  Spqtl[j]<-rDiscreta(matriz.prob.S.pat[j,])}
 #
 # Sp==1 ou Sm==1 indica que o alelo È do avÙ
 # Sp==2 ou Sm==2 indica que o alelo È da avÛ
 #
 for (i in 1:length(nonfounders)) dados.alelo<-rbind(dados.alelo,c(dados.alelo[mat.ped.pais[i,1],Spqtl[i]],dados.alelo[mat.ped.pais[i,2],Smqtl[i]]))
 #
 ## gera o genótipo do QTL de todos individuos
 #
 dados.QTL<-matrix(0,nrow(dados.alelo),1)
 for (i in 1:nrow(dados.alelo)){
  if (dados.alelo[i,1]==1 & dados.alelo[i,2]==1) dados.QTL[i,1]<--1
  if (dados.alelo[i,1]==2 & dados.alelo[i,2]==2) dados.QTL[i,1]<-1
  if (dados.alelo[i,1]!= dados.alelo[i,2]) dados.QTL[i,1]<-0}
 mat.delinea<-cbind(mat.delinea,dados.QTL)
 #
 #
 alpha<-poster.alpha(media.alpha,sigma2.alpha,sigma2.vig,residuos,alpha.vig,dados.QTL[,1])
 vet.coef<-rbind(vet.coef,alpha[[1]])
 predito<-mat.delinea%*%vet.coef
 residuos<-dados[,1]-predito
 #
 dados.QTL<-cbind(dados.QTL,(1-abs(dados.QTL[,1])))
 delta<-poster.delta(media.delta,sigma2.delta,sigma2.vig,residuos,delta.vig,dados.QTL[,2])
 mat.delinea<-cbind(mat.delinea,dados.QTL[,2])
 vet.coef<-rbind(vet.coef,delta[[1]])
 predito<-mat.delinea%*%vet.coef
 residuos<-dados[,1]-predito
 #
 mi<-poster.mi(media.mi,sigma2.mi,sigma2.vig,residuos,vet.coef[1,1])
 vet.coef[1,1]<-mi[[1]]
 predito<-mat.delinea%*%vet.coef
 residuos<-dados[,1]-predito
 sigma2<-poster.sigma2(neta.a,neta.b,residuos)
 #
 list(qtl[[1]],mat.delinea,vet.coef,sigma2[[1]],marcador[[2]],qtl[[2]],dados.alelo,alpha[[2]],delta[[2]],mi[[2]],sigma2[[2]],c(pos.qtls,qtl[[1]]),cbind(mat.alelo.qtls,dados.alelo))}
#
#################
### gera exclusão de um QTL no modelo e calcula a probabilidade dessa exclusão
#
gera.exclusao.QTL<-function(pos.qtls,vet.coef,mat.delinea,sigma2.vig,media.mi,sigma2.mi,neta.a,neta.b,mat.alelo.qtls){
 num.QTLs<-length(pos.qtls)
 qtl<-prob.excl.marc(num.QTLs,vet.coef)
 #
 pos.qtls<-pos.qtls[-qtl[[1]]]
 mat.alelo.qtls<-mat.alelo.qtls[,-c((2*qtl[[1]])-1,(2*qtl[[1]]))]
 mat.delinea<-matrix(mat.delinea[,-c((2*qtl[[1]]),(2*qtl[[1]]+1))],nrow=nrow(dados))
 vet.coef<-matrix(vet.coef[-c((2*qtl[[1]]),(2*qtl[[1]]+1))],ncol=1)
 predito<-mat.delinea%*%vet.coef
 residuos<-dados[,1]-predito
 #
 mi<-poster.mi(media.mi,sigma2.mi,sigma2.vig,residuos,vet.coef[1,1])
 vet.coef[1,1]<-mi[[1]]
 predito<-mat.delinea%*%vet.coef
 residuos<-dados[,1]-predito
 sigma2<-poster.sigma2(neta.a,neta.b,residuos)
 #
 list(qtl[[1]],mat.delinea,vet.coef,sigma2[[1]],qtl[[2]],mi[[2]],sigma2[[2]],sigma2[[2]],sigma2[[2]],sigma2[[2]],sigma2[[2]],pos.qtls,mat.alelo.qtls)}
#
#################
### gera junção de dois QTL no modelo
#
gera.juncao.QTL<-function(dados,num.QTLs,mat.delinea,vet.coef,pos.qtls,media.alpha,sigma2.alpha,sigma2.vig,media.delta,sigma2.delta,media.mi,sigma2.mi,neta.a,neta.b,mat.alelo.qtls){
 cramer<-numeric(num.QTLs-1)
 for (i in 1:(num.QTLs-1)) cramer[i]<-abs(cv.test(mat.delinea[,2*i],mat.delinea[,2*(i+1)])) # cramer entre os QTLs vizinhos
 prob_par<-cramer/sum(cramer)
 junta<-rDiscreta(prob_par)
 par_QTL<-c(junta,junta+1)
 efeitos<-c(1/sum(abs(vet.coef[2*junta,1]),abs(vet.coef[2*junta+1,1])),1/sum(abs(vet.coef[2*(junta+1),1]),abs(vet.coef[2*(junta+1)+1,1])))
 prob<-efeitos/sum(efeitos)
 gera<-rDiscreta(prob) # escolhe o QTL que vai sumir
 qtl<-par_QTL[gera]
 #
 pos.qtls.c<-pos.qtls[-qtl]
 mat.alelo.qtls.c<-mat.alelo.qtls[,-c((2*qtl)-1,(2*qtl))]
 mat.delinea.c<-matrix(mat.delinea[,-c((2*qtl),(2*qtl+1))],nrow=nrow(dados))
 vet.coef.c<-matrix(vet.coef[-c((2*qtl),(2*qtl+1))],ncol=1)
 predito<-mat.delinea.c%*%vet.coef.c
 residuos<-dados[,1]-predito
 #
 alpha<-poster.alpha(media.alpha,sigma2.alpha,sigma2.vig,residuos,vet.coef.c[(2*par_QTL[[1]]),1],mat.delinea.c[,(2*par_QTL[[1]])])
 vet.coef.c[(2*par_QTL[[1]]),1]<-alpha[[1]]
 residuos<-dados[,1]-(mat.delinea.c%*%vet.coef.c)
 delta<-poster.delta(media.delta,sigma2.delta,sigma2.vig,residuos,vet.coef.c[(2*par_QTL[[1]])+1,1],mat.delinea.c[,(2*par_QTL[[1]])+1])
 vet.coef.c[(2*par_QTL[[1]])+1,1]<-delta[[1]]
 residuos<-dados[,1]-(mat.delinea.c%*%vet.coef.c)
 #
 mi<-poster.mi(media.mi,sigma2.mi,sigma2.vig,residuos,vet.coef.c[1,1])
 vet.coef.c[1,1]<-mi[[1]]
 predito<-mat.delinea.c%*%vet.coef.c
 residuos<-dados[,1]-predito
 sigma2<-poster.sigma2(neta.a,neta.b,residuos)
 prob_merge<-log(prob[gera])+log(prob_par[junta])
 #
 list(qtl,mat.delinea.c,vet.coef.c,sigma2[[1]],prob_merge,mi[[2]],sigma2[[2]],alpha[[2]],delta[[2]],par_QTL[which(par_QTL!=qtl)],par_QTL,pos.qtls.c,mat.alelo.qtls.c)}
#
#################
### calcula a probabilidade de aceitação do nascimento
#
prob.aceitacao<-function(residuossp,residuos,sigma2,sigma2sp,misp,mi,media.mi,sigma2.mi,neta.a,neta.b,alphasp,media.alpha,sigma2.alpha,deltasp,media.delta,
sigma2.delta,loc.marc,num.QTLs,num.QTLssp,psplit,pmerge,pmarc,plambdasp,plambdamg,postmi,postsigma2,postalphasp,postdeltasp,postmisp,postsigma2sp,pri.pos.sp,pri.pos){
 vero<-sum(dnorm(residuossp,0,sqrt(sigma2sp),log=TRUE))-sum(dnorm(residuos,0,sqrt(sigma2),log=TRUE))
 priori<-dnorm(misp,media.mi,sqrt(sigma2.mi),log=TRUE)-dnorm(mi,media.mi,sqrt(sigma2.mi),log=TRUE)+
         dgamma((1/sigma2sp),neta.a,neta.b,log=TRUE)-dgamma((1/sigma2),neta.a,neta.b,log=TRUE)+
         dnorm(alphasp,media.alpha,sqrt(sigma2.alpha),log=TRUE)+dnorm(deltasp,media.delta,sqrt(sigma2.delta),log=TRUE)+
         pri.pos.sp-pri.pos
 trans<-pmerge+plambdamg+postmi+postsigma2-psplit-pmarc-plambdasp-postalphasp-postdeltasp-postmisp-postsigma2sp
 prob.ace<-exp(vero+priori+trans)
 prob.ace}
#
#################
### calcula a probabilidade de aceitação do split de um QTL
#
prob.aceitacao.split<-function(residuossp,residuos,sigma2,sigma2sp,misp,mi,media.mi,sigma2.mi,neta.a,neta.b,alphamg,alphasp1,alphasp2,media.alpha,sigma2.alpha,deltamg,deltasp1,deltasp2,media.delta,
sigma2.delta,loc.marc,num.QTLs,num.QTLssp,psplit,pmerge,pmarc,plambdasp,plambdamg,postalpha,postdelta,postmi,postsigma2,postalphasp1,postalphasp2,postdeltasp1,postdeltasp2,postmisp,postsigma2sp,pri.pos.sp,pri.pos){
 vero<-sum(dnorm(residuossp,0,sqrt(sigma2sp),log=TRUE))-sum(dnorm(residuos,0,sqrt(sigma2),log=TRUE))
 priori<-dnorm(misp,media.mi,sqrt(sigma2.mi),log=TRUE)-dnorm(mi,media.mi,sqrt(sigma2.mi),log=TRUE)+
         dgamma((1/sigma2sp),neta.a,neta.b,log=TRUE)-dgamma((1/sigma2),neta.a,neta.b,log=TRUE)+
         dnorm(alphasp1,media.alpha,sqrt(sigma2.alpha),log=TRUE)+dnorm(deltasp1,media.delta,sqrt(sigma2.delta),log=TRUE)+
         dnorm(alphasp2,media.alpha,sqrt(sigma2.alpha),log=TRUE)+dnorm(deltasp2,media.delta,sqrt(sigma2.delta),log=TRUE)-
         (dnorm(alphamg,media.alpha,sqrt(sigma2.alpha),log=TRUE)+dnorm(deltamg,media.delta,sqrt(sigma2.delta),log=TRUE))+
         pri.pos.sp-pri.pos
 trans<-pmerge+plambdamg+postalpha+postdelta+postmi+postsigma2-psplit-pmarc-plambdasp-postalphasp1-postalphasp2-postdeltasp1-postdeltasp2-postmisp-postsigma2sp
 prob.ace<-exp(vero+priori+trans)
 prob.ace}