#
##################
##################
# running the DDRJ
##################
##################
#
neta.a<-1
neta.b<-1
sigma2.mi<-100
media.mi<-0
sigma2.alpha<-100
media.alpha<-0
sigma2.delta<-100
media.delta<-0
#
####### starting the model
#
set.seed(510)
num.marc<-ncol(dados)-1
residuos<-dados[,1]
sigma2.vig<-poster.sigma2(neta.a,neta.b,residuos)[[1]]
#
mi<-poster.mi(media.mi,sigma2.mi,sigma2.vig,residuos,0)[[1]]
#
pos.qtls<-NULL # posição dos QTLs no modelo. Se tivermos QTLs no modelo inicial, precisamos definir seus genÛtipos e efeitos e alterar os resÌduos
num.QTLs<-length(pos.qtls)
mat.alelo.qtls<-NULL # colunas Ìmpares s„o os alelos paternos e colunas pares s„o os alelos maternos
mat.delinea<-matrix(1,nrow(dados),1)
vet.coef<-matrix(mi,1,1)
predito<-mat.delinea%*%vet.coef
residuos<-dados[,1]-predito
sigma2.vig<-poster.sigma2(neta.a,neta.b,residuos)[[1]]
#
amostrasfin<-100
burnin<-100
saltos<-5
AmostrasTotal<-burnin+amostrasfin*saltos
#
indSpMgtotal<-numeric(AmostrasTotal) # 1 - candidato com inclusão de QTL; 2 - candidato com exclusão de QTL
probacetotal<-numeric(AmostrasTotal)
probacemerge<-numeric(AmostrasTotal)
num.QTL.total<-numeric(AmostrasTotal)
indrejtotal<-numeric(AmostrasTotal)
#
install.package(compiler)
library(compiler)
enableJIT(3)
#
for (int in (1:AmostrasTotal)){
  cat('\n', cromos, int)
  cand.sp.mg<-dec.sp.mg(num.QTLs,num.marc)
  indSpMgtotal[int]<-cand.sp.mg[[1]]
  #
  ##############
  ###### QTL birth step
  ##############
  #
  if (indSpMgtotal[int]==1) {candidato<-gera.inclusao.QTL(dados,residuos,mat.delinea,vet.coef,pos.qtls,loc.marc,sigma2.vig,0,0,media.mi,sigma2.mi,media.alpha,sigma2.alpha,media.delta,sigma2.delta,neta.a,neta.b,
mat.ped.pais,founders,nonfounders,alelo.paterno,alelo.materno,Sp,Sm,mat.alelo.qtls)
   num.QTL.total[int]<-num.QTLs+1
   residuossp<-dados[,1]-(candidato[[2]]%*%candidato[[3]])
   sigma2<-sigma2.vig
   sigma2sp<-candidato[[4]]
   mi<-vet.coef[1,1]
   misp<-candidato[[3]][1,1]
   alphasp<-candidato[[3]][(nrow(candidato[[3]])-1),1]
   deltasp<-candidato[[3]][(nrow(candidato[[3]])),1]
   num.QTLssp<-num.QTLs+1
   psplit<-cand.sp.mg[[2]][1]
   pmerge<-dec.sp.mg(num.QTLssp,num.marc)[[2]][2]
   pmarc<-candidato[[5]]
   plambdasp<-candidato[[6]]
  #
   efeito<-NULL
   if (num.QTLs==0) efeito<-0
   if (num.QTLs>0) for (i in 1:num.QTLs) efeito[i]<-1/(abs(vet.coef[2*i,1])+abs(vet.coef[(2*i)+1,1]))
   efeitosp<-1/(abs(alphasp)+abs(deltasp))
   plambdamg<-log(efeitosp/(efeitosp+sum(efeito)))
  #
   res.mi<-residuos+mi
   quo<-(length(residuos)/sigma2sp)+(1/sigma2.mi)
   media<-((sum(res.mi)/sigma2sp)+(media.mi/sigma2.mi))/quo
   variancia<-1/quo
   postmi<-dnorm(mi,media,sqrt(variancia),log = TRUE)
  #
   aa1<-(length(residuos)/2)+neta.a
   bb1<-(sum(residuos^2)/2)+neta.b
   postsigma2<-dgamma((1/sigma2),aa1,bb1,log = TRUE)
  #
   postalphasp<-candidato[[8]]
   postdeltasp<-candidato[[9]]
   postmisp<-candidato[[10]]
   postsigma2sp<-candidato[[11]]
   pos.qtls.sp<-sort(candidato[[12]])
   pri.pos.sp<-dens.priori.loc.qtls(pos.qtls.sp,loc.marc)
   pri.pos<-dens.priori.loc.qtls(pos.qtls,loc.marc)
  #
   probace<-prob.aceitacao(residuossp,residuos,sigma2,sigma2sp,misp,mi,media.mi,sigma2.mi,neta.a,neta.b,alphasp,media.alpha,sigma2.alpha,deltasp,media.delta,
sigma2.delta,loc.marc,num.QTLs,num.QTLssp,psplit,pmerge,pmarc,plambdasp,plambdamg,postmi,postsigma2,postalphasp,postdeltasp,postmisp,postsigma2sp,pri.pos.sp,pri.pos)}
  #
  ##############
  ###### QTL death step
  ##############
  #
  if (indSpMgtotal[int]==2) {candidato<-gera.exclusao.QTL(pos.qtls,vet.coef,mat.delinea,sigma2.vig,media.mi,sigma2.mi,neta.a,neta.b,mat.alelo.qtls)
   num.QTL.total[int]<-num.QTLs-1
   residuosmg<-dados[,1]-(candidato[[2]]%*%candidato[[3]])
   sigma2<-sigma2.vig
   sigma2mg<-candidato[[4]]
   mi<-vet.coef[1,1]
   mimg<-candidato[[3]][1,1]
   num.QTLsmg<-num.QTLs-1
   pmerge<-cand.sp.mg[[2]][2]
   psplit<-dec.sp.mg(num.QTLsmg,num.marc)[[2]][1]
   postmimg<-candidato[[6]]
   postsigma2mg<-candidato[[7]]
   pos.qtls.mg<-sort(candidato[[12]])
   pri.pos<-dens.priori.loc.qtls(pos.qtls,loc.marc)
   pri.pos.mg<-dens.priori.loc.qtls(pos.qtls.mg,loc.marc)
   plambdamg<-candidato[[5]]
   alpha<-vet.coef[(2*candidato[[1]]),1]
   delta<-vet.coef[(2*candidato[[1]]+1),1]
  #
   krusk<-prob.selec.marc(dados,residuosmg,pos.qtls.mg,loc.marc)[[3]]
   prob<-krusk/sum(krusk)
   QTLmg<-pos.qtls[candidato[[1]]]
   Marc1<-sum(loc.marc<=QTLmg)
   Marc2<-num.marc-sum(loc.marc>=QTLmg)+1
   loc.Marc1<-dens.pos.qtl(QTLmg,pos.qtls.mg,loc.marc,Marc1,krusk)
   loc.Marc2<-dens.pos.qtl(QTLmg,pos.qtls.mg,loc.marc,Marc2,krusk)
   plambdasp<-log(prob[Marc1]*loc.Marc1+prob[Marc2]*loc.Marc2)
   pmarc<-0
  #
   quo<-(sum(mat.delinea[,2*candidato[[1]]]^2)/sigma2mg)+(1/sigma2.alpha)
   media<-((sum(mat.delinea[,2*candidato[[1]]]*residuosmg)/sigma2mg)+(media.alpha/sigma2.alpha))/quo
   variancia<-1/quo
   postalpha<-dnorm(alpha,media,sqrt(variancia),log = TRUE)
  #
   quo<-(sum(mat.delinea[,2*candidato[[1]]+1]^2)/sigma2mg)+(1/sigma2.delta)
   media<-((sum(mat.delinea[,2*candidato[[1]]+1]*residuosmg)/sigma2mg)+(media.delta/sigma2.delta))/quo
   variancia<-1/quo
   postdelta<-dnorm(delta,media,sqrt(variancia),log = TRUE)
  #
   res.mi<-residuos+vet.coef[1,1]
   quo<-(length(residuos)/sigma2mg)+(1/sigma2.mi)
   media<-((sum(res.mi)/sigma2mg)+(media.mi/sigma2.mi))/quo
   variancia<-1/quo
   postmi<-dnorm(mi,media,sqrt(variancia),log = TRUE)
  #
   aa1<-(length(residuos)/2)+neta.a
   bb1<-(sum(residuos^2)/2)+neta.b
   postsigma2<-dgamma((1/sigma2),aa1,bb1,log = TRUE)
  #
   probace<-1/prob.aceitacao(residuos,residuosmg,sigma2mg,sigma2,mi,mimg,media.mi,sigma2.mi,neta.a,neta.b,alpha,media.alpha,sigma2.alpha,delta,media.delta,sigma2.delta,
loc.marc,num.QTLsmg,num.QTLs,psplit,pmerge,pmarc,plambdasp,plambdamg,postmimg,postsigma2mg,postalpha,postdelta,postmi,postsigma2,pri.pos,pri.pos.mg)}
  #
  #####################
  #  updating parameters
  #####################
  #
  # update of QTL's number
  #
  probacetotal[int]<-probace
  aux2<-runif(1)
  if (aux2<probace){
   indrejtotal[int]<-0
   pos.qtls<-candidato[[12]]
   num.QTLs<-length(pos.qtls)
   mat.delinea<-candidato[[2]]
   vet.coef<-candidato[[3]]
   sigma2.vig<-candidato[[4]]
   mat.alelo.qtls<-candidato[[13]]
   if (num.QTLs>0){
    posicao<-order(pos.qtls)
    pos.qtls<-pos.qtls[posicao]
    posicao2<-1
    for (i in 1:length(posicao)) posicao2<-c(posicao2,posicao[i]*2,posicao[i]*2+1)
    mat.alelo.qtls<-mat.alelo.qtls[,posicao2[-1]-1]
    mat.delinea<-mat.delinea[,posicao2]
    vet.coef<-matrix(vet.coef[posicao2,],ncol(mat.delinea),1)}
   residuos<-dados[,1]-(mat.delinea%*%vet.coef)}
  #
  if (aux2>=probace) indrejtotal[int]<-1
  #
  #### QTLs merge step
  #
  if (num.QTLs>1){
   candidato<-gera.juncao.QTL(dados,num.QTLs,mat.delinea,vet.coef,pos.qtls,media.alpha,sigma2.alpha,sigma2.vig,media.delta,sigma2.delta,media.mi,sigma2.mi,neta.a,neta.b,mat.alelo.qtls)
   residuosmg<-dados[,1]-(candidato[[2]]%*%candidato[[3]])
   sigma2<-sigma2.vig
   sigma2mg<-candidato[[4]]
   mi<-vet.coef[1,1]
   mimg<-candidato[[3]][1,1]
   num.QTLsmg<-num.QTLs-1
   pmerge<-0
   psplit<-0
   postmimg<-candidato[[6]]
   postsigma2mg<-candidato[[7]]
   postalphamg<-candidato[[8]]
   postdeltamg<-candidato[[9]]
   pos.qtls.mg<-sort(candidato[[12]])
   pri.pos<-dens.priori.loc.qtls(pos.qtls,loc.marc)
   pri.pos.mg<-dens.priori.loc.qtls(pos.qtls.mg,loc.marc)
   plambdamg<-candidato[[5]]
   alphamg<-candidato[[3]][(2*candidato[[11]][1]),1]
   deltamg<-candidato[[3]][(2*candidato[[11]][1]+1),1]
   alphasp1<-vet.coef[(2*candidato[[1]]),1]
   deltasp1<-vet.coef[(2*candidato[[1]]+1),1]
   alphasp2<-vet.coef[(2*candidato[[10]]),1]
   deltasp2<-vet.coef[(2*candidato[[10]]+1),1]
   #
   krusk<-prob.selec.marc(dados,residuosmg,pos.qtls.mg,loc.marc)[[3]]
   prob<-krusk/sum(krusk)
   QTLmg<-pos.qtls[candidato[[1]]]
   Marc1<-sum(loc.marc<=QTLmg)
   Marc2<-num.marc-sum(loc.marc>=QTLmg)+1
   loc.Marc1<-dens.pos.qtl(QTLmg,pos.qtls.mg,loc.marc,Marc1,krusk)
   loc.Marc2<-dens.pos.qtl(QTLmg,pos.qtls.mg,loc.marc,Marc2,krusk)
   plambdasp<-log(prob[Marc1]*loc.Marc1+prob[Marc2]*loc.Marc2)
   pmarc<-0
   #
   quo<-(sum(mat.delinea[,2*candidato[[1]]]^2)/sigma2mg)+(1/sigma2.alpha)
   media<-((sum(mat.delinea[,2*candidato[[1]]]*residuosmg)/sigma2mg)+(media.alpha/sigma2.alpha))/quo
   variancia<-1/quo
   postalpha1<-dnorm(alphasp1,media,sqrt(variancia),log = TRUE)
   #
   quo<-(sum(mat.delinea[,2*candidato[[1]]+1]^2)/sigma2mg)+(1/sigma2.delta)
   media<-((sum(mat.delinea[,2*candidato[[1]]+1]*residuosmg)/sigma2mg)+(media.delta/sigma2.delta))/quo
   variancia<-1/quo
   postdelta1<-dnorm(deltasp1,media,sqrt(variancia),log = TRUE)
   #
   vet_coef_parc<-rbind(candidato[[3]],vet.coef[candidato[[1]]*2,1],vet.coef[candidato[[1]]*2+1,1])
   mat_delinea_parc<-cbind(candidato[[2]],mat.delinea[,candidato[[1]]*2],mat.delinea[,candidato[[1]]*2+1])
   vet_coef_parc<-matrix(vet_coef_parc[-(2*candidato[[11]][1]),],ncol=1)
   mat_delinea_parc<-mat_delinea_parc[,-(2*candidato[[11]][1])]
   residuosparc<-dados[,1]-mat_delinea_parc%*%vet_coef_parc
   quo<-(sum(mat.delinea[,2*candidato[[10]]]^2)/sigma2mg)+(1/sigma2.alpha)
   media<-((sum(mat.delinea[,2*candidato[[10]]]*residuosparc)/sigma2mg)+(media.alpha/sigma2.alpha))/quo
   variancia<-1/quo
   postalpha2<-dnorm(alphasp2,media,sqrt(variancia),log = TRUE)
   #
   vet_coef_parc<-matrix(vet.coef[-(2*candidato[[10]]+1)],ncol=1)
   mat_delinea_parc<-mat.delinea[,-(2*candidato[[10]]+1)]
   residuosparc<-dados[,1]-mat_delinea_parc%*%vet_coef_parc
   quo<-(sum(mat.delinea[,2*candidato[[10]]+1]^2)/sigma2mg)+(1/sigma2.delta)
   media<-((sum(mat.delinea[,2*candidato[[10]]+1]*residuosparc)/sigma2mg)+(media.delta/sigma2.delta))/quo
   variancia<-1/quo
   postdelta2<-dnorm(deltasp2,media,sqrt(variancia),log = TRUE)
   #
   res.mi<-residuos+vet.coef[1,1]
   quo<-(length(residuos)/sigma2mg)+(1/sigma2.mi)
   media<-((sum(res.mi)/sigma2mg)+(media.mi/sigma2.mi))/quo
   variancia<-1/quo
   postmi<-dnorm(mi,media,sqrt(variancia),log = TRUE)
   #
   aa1<-(length(residuos)/2)+neta.a
   bb1<-(sum(residuos^2)/2)+neta.b
   postsigma2<-dgamma((1/sigma2),aa1,bb1,log = TRUE)
   #
   probace<-1/prob.aceitacao.split(residuos,residuosmg,sigma2mg,sigma2,mi,mimg,media.mi,sigma2.mi,neta.a,neta.b,alphamg,alphasp1,alphasp2,media.alpha,sigma2.alpha,deltamg,deltasp1,deltasp2,media.delta,
sigma2.delta,loc.marc,num.QTLsmg,num.QTLs,psplit,pmerge,pmarc,plambdasp,plambdamg,postalphamg,postdeltamg,postmimg,postsigma2mg,postalpha1,postalpha2,postdelta1,postdelta2,postmi,postsigma2,pri.pos,pri.pos.mg)
   #
   probacemerge[int]<-probace
   aux2<-runif(1)
   if (aux2<probace){
    pos.qtls<-candidato[[12]]
    num.QTLs<-length(pos.qtls)
    mat.alelo.qtls<-candidato[[13]]
    mat.delinea<-candidato[[2]]
    vet.coef<-candidato[[3]]
    sigma2.vig<-candidato[[4]]
    residuos<-dados[,1]-(mat.delinea%*%vet.coef)}}  #  #### update QTLs' position - Metropolis Hasting
  #
  if (num.QTLs>0){
   for (i in 1:num.QTLs){
    M.esq<-sum(loc.marc<=pos.qtls[i])
    M.dir<-num.marc-(sum(loc.marc>=pos.qtls[i]))+1
    loc.cand<-runif(1,min=loc.marc[M.esq],max=loc.marc[M.dir])
    #
    vet.delinea<-mat.delinea
    mat.alelo.qtlsc<-mat.alelo.qtls
    r1c<-Haldane(abs(loc.cand-loc.marc[M.esq]))
    r2c<-Haldane(abs(loc.cand-loc.marc[M.dir]))
    r1<-Haldane(abs(pos.qtls[i]-loc.marc[M.esq]))
    r2<-Haldane(abs(pos.qtls[i]-loc.marc[M.dir]))
    r12<-Haldane(abs(loc.marc[M.dir]-loc.marc[M.esq]))
    matrec1<-matrix(c(1-r1,r1,r1,1-r1),2,2,byrow=TRUE)
    matrec2<-matrix(c(1-r2,r2,r2,1-r2),2,2,byrow=TRUE)
    matrec1c<-matrix(c(1-r1c,r1c,r1c,1-r1c),2,2,byrow=TRUE)
    matrec2c<-matrix(c(1-r2c,r2c,r2c,1-r2c),2,2,byrow=TRUE)
    probQTLc<-numeric(3)
    probQTL<-numeric(3)
    logdens<-numeric(3)
    probacum<-0
    probacumc<-0
    densacum<-0
    densacumc<-0
    numer<-denom<-0
    #
    # atualizando o genótipo e os alelos do QTL para os founders e calculando as probabilidades para os alelos e genótipo vigentes
    #
    gen<-c(-1,0,1)
    for (j in 1:length(founders)){
     for (l in 1:3){
      probQTLc[l]<-log((calc.prob.gen(dados[j,(M.esq+1)],gen[l],r1c)*calc.prob.gen(gen[l],dados[j,(M.dir+1)],r2c))/calc.prob.gen(dados[j,(M.esq+1)],dados[j,(M.dir+1)],r12))
      probQTL[l]<-log((calc.prob.gen(dados[j,(M.esq+1)],gen[l],r1)*calc.prob.gen(gen[l],dados[j,(M.dir+1)],r2))/calc.prob.gen(dados[j,(M.esq+1)],dados[j,(M.dir+1)],r12))
      dom<-1-abs(gen[l])
      vet.delinea[j,2*i]<-gen[l]
      vet.delinea[j,(2*i)+1]<-dom
      logdens[l]<-dnorm(dados[j,1],vet.delinea[j,]%*%vet.coef,sqrt(sigma2.vig),log=TRUE)}
     probc<-exp(probQTLc+logdens-max(probQTLc+logdens))/sum(exp(probQTLc+logdens-max(probQTLc+logdens)))
     prob<-exp(probQTL+logdens-max(probQTL+logdens))/sum(exp(probQTL+logdens-max(probQTL+logdens)))
     ger.gen<-rDiscreta(probc)
     vet.delinea[j,2*i]<-gen[ger.gen]
     vet.delinea[j,(2*i)+1]<-1-abs(vet.delinea[j,2*i])
     probacumc<-probacumc+log(probc[ger.gen])
     densacumc<-densacumc+logdens[ger.gen]
     numer<-numer+probQTLc[ger.gen]
     gen.at<-which(gen==mat.delinea[j,2*i])
     probacum<-probacum+log(prob[gen.at])
     densacum<-densacum+logdens[gen.at]
     denom<-denom+probQTL[gen.at]}
    #
    ale.pai<-matrix(99,nrow=length(founders),1)
    ale.mae<-matrix(99,nrow=length(founders),1)
    probale<-numeric(2)
    for (k in 1:length(founders)){
     if (vet.delinea[k,2*i]==-1){
      ale.mae[k,1]<-1
      ale.pai[k,1]<-1}
     if (vet.delinea[k,2*i]==1){
      ale.mae[k,1]<-2
      ale.pai[k,1]<-2}
     if (vet.delinea[k,2*i]==0){
      for (ale in 1:2) probale[ale]<-matrec1[alelo.paterno[k,M.esq],ale]*matrec2[ale,alelo.paterno[k,M.dir]]
      probale<-probale/sum(probale)
      aux2<-rDiscreta(probale)-1
      ale.mae[k,1]<-1^(aux2)*2^(1-aux2)
      ale.pai[k,1]<-2^(aux2)*1^(1-aux2)}
     mat.alelo.qtlsc[k,(2*i-1)]<-ale.pai[k,1]
     mat.alelo.qtlsc[k,(2*i)]<-ale.mae[k,1]}
    #
    matriz.prob.S<-NULL
    matriz.prob.Sc<-NULL
    probSpat<-numeric(2)
    probSmat<-numeric(2)
    probSpatc<-numeric(2)
    probSmatc<-numeric(2)
    #
    # atualizando os alelos e genótipos do QTL para os nonfounders e calculando as probabilidades para os alelos e genótipos vigentes
    #
    for (j in 1:length(nonfounders)){
     for (l in 1:2){
      probSmat[l]<-matrec1[Sm[j,M.esq],l]*matrec2[l,Sm[j,M.dir]]
      probSpat[l]<-matrec1[Sp[j,M.esq],l]*matrec2[l,Sp[j,M.dir]]
      probSmatc[l]<-matrec1c[Sm[j,M.esq],l]*matrec2c[l,Sm[j,M.dir]]
      probSpatc[l]<-matrec1c[Sp[j,M.esq],l]*matrec2c[l,Sp[j,M.dir]]}
     probSpat<-probSpat/sum(probSpat)
     probSmat<-probSmat/sum(probSmat)
     probSpatc<-probSpatc/sum(probSpatc)
     probSmatc<-probSmatc/sum(probSmatc)
     matriz.prob.S<-c(probSpat*probSmat[1],probSpat*probSmat[2])
     matriz.prob.Sc<-c(probSpatc*probSmatc[1],probSpatc*probSmatc[2])
     #
     # calculando as probabilidades para os alelos e genÛtipo vigentes
     #
     mat.possc<-cbind(matrix(c(1,1,2,1,1,2,2,2),4,2,byrow=TRUE),matrix(0,4,2)) # col1 = heranÁa paterna (Sp), col2 = heranÁa materna (Sm)
     for (l in 1:4){
      mat.possc[l,3]<-mat.alelo.qtls[mat.ped.pais[j,1],2*i-(2-mat.possc[l,1])]
      mat.possc[l,4]<-mat.alelo.qtls[mat.ped.pais[j,2],2*i-(2-mat.possc[l,2])]}
     mat.possc<-cbind(mat.possc,matrix(matriz.prob.S,4,1))
     if ((sum(mat.possc[,3])==4 | sum(mat.possc[,3])==8) & (sum(mat.possc[,4])==4 | sum(mat.possc[,4])==8)){
      mat.possc[1,5]<-1
      mat.possc<-matrix(mat.possc[-c(2,3,4),],1,5)} else {
      if (sum(mat.possc[,3])==4 | sum(mat.possc[,3])==8){
       mat.possc[1,5]<-mat.possc[1,5]+mat.possc[2,5]
       mat.possc[3,5]<-mat.possc[3,5]+mat.possc[4,5]
       mat.possc<-mat.possc[-c(2,4),]} else {
       if (sum(mat.possc[,4])==4 | sum(mat.possc[,4])==8){
        mat.possc[1,5]<-mat.possc[1,5]+mat.possc[3,5]
        mat.possc[2,5]<-mat.possc[2,5]+mat.possc[4,5]
        mat.possc<-mat.possc[-c(3,4),]}}}
     mat.possc<-cbind(mat.possc,matrix(0,nrow(mat.possc),2))
     for (l in 1:nrow(mat.possc)){
      if (mat.possc[l,3]==1 & mat.possc[l,4]==1) mat.possc[l,6]<--1
      if (mat.possc[l,3]==2 & mat.possc[l,4]==2) mat.possc[l,6]<-1
      if (mat.possc[l,3]!= mat.possc[l,4]) mat.possc[l,6]<-0}
     mat.possc[,7]<-1-abs(mat.possc[,6])
     #
     logdens<-numeric(nrow(mat.possc))
     for (l in 1:nrow(mat.possc)){
      vet.delinea[length(founders)+j,2*i]<-mat.possc[l,6]
      vet.delinea[length(founders)+j,(2*i)+1]<-mat.possc[l,7]
      logdens[l]<-dnorm(dados[length(founders)+j,1],vet.delinea[length(founders)+j,]%*%vet.coef,sqrt(sigma2.vig),log=TRUE)}
     probc<-exp(log(mat.possc[,5])+logdens-max(log(mat.possc[,5])+logdens))/sum(exp(log(mat.possc[,5])+logdens-max(log(mat.possc[,5])+logdens)))
     gen<-which(mat.possc[,3]==mat.alelo.qtls[length(founders)+j,(2*i)-1] & mat.possc[,4]==mat.alelo.qtls[length(founders)+j,2*i])
     probacum<-probacum+log(probc[gen])
     densacum<-densacum+logdens[gen]
     denom<-denom+log(mat.possc[gen,5])
     #
     # calculando as probabilidades para os alelos e genótipo candidatos
     #
     mat.possc<-cbind(matrix(c(1,1,2,1,1,2,2,2),4,2,byrow=TRUE),matrix(0,4,2)) # col1 = heranÁa paterna (Sp), col2 = heranÁa materna (Sm)
     for (l in 1:4){
      mat.possc[l,3]<-mat.alelo.qtlsc[mat.ped.pais[j,1],2*i-(2-mat.possc[l,1])]
      mat.possc[l,4]<-mat.alelo.qtlsc[mat.ped.pais[j,2],2*i-(2-mat.possc[l,2])]}
     mat.possc<-cbind(mat.possc,matrix(matriz.prob.Sc,4,1))
     if ((sum(mat.possc[,3])==4 | sum(mat.possc[,3])==8) & (sum(mat.possc[,4])==4 | sum(mat.possc[,4])==8)){
      mat.possc[1,5]<-1
      mat.possc<-matrix(mat.possc[-c(2,3,4),],1,5)} else {
      if (sum(mat.possc[,3])==4 | sum(mat.possc[,3])==8){
       mat.possc[1,5]<-mat.possc[1,5]+mat.possc[2,5]
       mat.possc[3,5]<-mat.possc[3,5]+mat.possc[4,5]
       mat.possc<-mat.possc[-c(2,4),]} else {
       if (sum(mat.possc[,4])==4 | sum(mat.possc[,4])==8){
        mat.possc[1,5]<-mat.possc[1,5]+mat.possc[3,5]
        mat.possc[2,5]<-mat.possc[2,5]+mat.possc[4,5]
        mat.possc<-mat.possc[-c(3,4),]}}}
     mat.possc<-cbind(mat.possc,matrix(0,nrow(mat.possc),2))
     for (l in 1:nrow(mat.possc)){
      if (mat.possc[l,3]==1 & mat.possc[l,4]==1) mat.possc[l,6]<--1
      if (mat.possc[l,3]==2 & mat.possc[l,4]==2) mat.possc[l,6]<-1
      if (mat.possc[l,3]!= mat.possc[l,4]) mat.possc[l,6]<-0}
     mat.possc[,7]<-1-abs(mat.possc[,6])
     #
     logdens<-numeric(nrow(mat.possc))
     for (l in 1:nrow(mat.possc)){
      vet.delinea[length(founders)+j,2*i]<-mat.possc[l,6]
      vet.delinea[length(founders)+j,(2*i)+1]<-mat.possc[l,7]
      logdens[l]<-dnorm(dados[length(founders)+j,1],vet.delinea[length(founders)+j,]%*%vet.coef,sqrt(sigma2.vig),log=TRUE)}
     probc<-exp(log(mat.possc[,5])+logdens-max(log(mat.possc[,5])+logdens))/sum(exp(log(mat.possc[,5])+logdens-max(log(mat.possc[,5])+logdens)))
     ger.gen<-rDiscreta(probc)
     vet.delinea[length(founders)+j,2*i]<-mat.possc[ger.gen,6]
     vet.delinea[length(founders)+j,(2*i)+1]<-mat.possc[ger.gen,7]
     mat.alelo.qtlsc[length(founders)+j,(2*i-1)]<-mat.possc[ger.gen,3]
     mat.alelo.qtlsc[length(founders)+j,(2*i)]<-mat.possc[ger.gen,4]
     probacumc<-probacumc+log(probc[ger.gen])
     densacumc<-densacumc+logdens[ger.gen]
     numer<-numer+log(mat.possc[ger.gen,5])}
    #
    #
    paceit<-exp(densacumc+numer+probacum-densacum-denom-probacumc)
    aux3<-runif(1)
    if (aux3<paceit){
     pos.qtls[i]<-loc.cand
     mat.delinea[,2*i]<-vet.delinea[,2*i]
     mat.delinea[,(2*i)+1]<-vet.delinea[,(2*i)+1]
     mat.alelo.qtls[,(2*i-1)]<-mat.alelo.qtlsc[,(2*i-1)]
     mat.alelo.qtls[,(2*i)]<-mat.alelo.qtlsc[,(2*i)]}}}
  #
  #### updating QTLs genotype
  #
  if (num.QTLs>0){
   for (i in 1:num.QTLs){
    M.esq<-sum(loc.marc<=pos.qtls[i])
    M.dir<-num.marc-(sum(loc.marc>=pos.qtls[i]))+1
    r1<-Haldane(abs(pos.qtls[i]-loc.marc[M.esq]))
    r2<-Haldane(abs(pos.qtls[i]-loc.marc[M.dir]))
    r12<-Haldane(abs(loc.marc[M.dir]-loc.marc[M.esq]))
    matrec1<-matrix(c(1-r1,r1,r1,1-r1),2,2,byrow=TRUE)
    matrec2<-matrix(c(1-r2,r2,r2,1-r2),2,2,byrow=TRUE)
    probQTL<-numeric(3)
    logdens<-numeric(3)
    matriz.prob.S<-NULL
    probSpat<-numeric(2)
    probSmat<-numeric(2)
    gen<-c(-1,0,1)
    #
    # atualizando o genótipo e os alelos do QTL para os founders
    #
    for (j in 1:length(founders)){
     for (l in 1:3){
      probQTL[l]<-log((calc.prob.gen(dados[j,(M.esq+1)],gen[l],r1)*calc.prob.gen(gen[l],dados[j,(M.dir+1)],r2))/calc.prob.gen(dados[j,(M.esq+1)],dados[j,(M.dir+1)],r12))
      dom<-1-abs(gen[l])
      vet.delinea<-mat.delinea[j,]
      vet.delinea[2*i]<-gen[l]
      vet.delinea[(2*i)+1]<-dom
      logdens[l]<-dnorm(dados[j,1],vet.delinea%*%vet.coef,sqrt(sigma2.vig),log=TRUE)}
     prob<-exp(probQTL+logdens-max(probQTL+logdens))/sum(exp(probQTL+logdens-max(probQTL+logdens)))
     mat.delinea[j,2*i]<-gen[rDiscreta(prob)]
     mat.delinea[j,(2*i)+1]<-1-abs(mat.delinea[j,2*i])}
    monoto<-length(table(mat.delinea[1:length(founders),2*i]))
    if (monoto==1){
     indiv<-sample(seq(1:length(founders)),1)
     if (mat.delinea[indiv,2*i]==-1) mat.delinea[indiv,2*i]<-0
     if (mat.delinea[indiv,2*i]==0) mat.delinea[indiv,2*i]<-1
     if (mat.delinea[indiv,2*i]==1) mat.delinea[indiv,2*i]<-0
     mat.delinea[indiv,(2*i)+1]<-1-abs(mat.delinea[indiv,2*i])}
    #
    ale.pai<-matrix(99,nrow=length(founders),1)
    ale.mae<-matrix(99,nrow=length(founders),1)
    probale<-numeric(2)
    for (k in 1:length(founders)){
     if (mat.delinea[k,2*i]==-1){
      ale.mae[k,1]<-1
      ale.pai[k,1]<-1}
     if (mat.delinea[k,2*i]==1){
      ale.mae[k,1]<-2
      ale.pai[k,1]<-2}
     if (mat.delinea[k,2*i]==0){
      for (ale in 1:2) probale[ale]<-matrec1[alelo.paterno[k,M.esq],ale]*matrec2[ale,alelo.paterno[k,M.dir]]
      probale<-probale/sum(probale)
      aux2<-rDiscreta(probale)-1
      ale.mae[k,1]<-1^(aux2)*2^(1-aux2)
      ale.pai[k,1]<-2^(aux2)*1^(1-aux2)}
     mat.alelo.qtls[k,(2*i-1)]<-ale.pai[k,1]
     mat.alelo.qtls[k,(2*i)]<-ale.mae[k,1]}
    #
    # atualizando os alelos e genÛtipos do QTL para os nonfounders
    #
    for (j in 1:length(nonfounders)){
     for (l in 1:2){
      probSmat[l]<-matrec1[Sm[j,M.esq],l]*matrec2[l,Sm[j,M.dir]]
      probSpat[l]<-matrec1[Sp[j,M.esq],l]*matrec2[l,Sp[j,M.dir]]}
     probSpat<-probSpat/sum(probSpat)
     probSmat<-probSmat/sum(probSmat)
     matriz.prob.S<-c(probSpat*probSmat[1],probSpat*probSmat[2])
     mat.possc<-cbind(matrix(c(1,1,2,1,1,2,2,2),4,2,byrow=TRUE),matrix(0,4,2)) # col1 = heranÁa paterna (Sp), col2 = heranÁa materna (Sm)
     for (l in 1:4){
      mat.possc[l,3]<-mat.alelo.qtls[mat.ped.pais[j,1],2*i-(2-mat.possc[l,1])]
      mat.possc[l,4]<-mat.alelo.qtls[mat.ped.pais[j,2],2*i-(2-mat.possc[l,2])]}
     mat.possc<-cbind(mat.possc,matrix(matriz.prob.S,4,1))
     if ((sum(mat.possc[,3])==4 | sum(mat.possc[,3])==8) & (sum(mat.possc[,4])==4 | sum(mat.possc[,4])==8)){
      mat.possc[1,5]<-1
      mat.possc<-matrix(mat.possc[-c(2,3,4),],1,5)} else {
      if (sum(mat.possc[,3])==4 | sum(mat.possc[,3])==8){
       mat.possc[1,5]<-mat.possc[1,5]+mat.possc[2,5]
       mat.possc[3,5]<-mat.possc[3,5]+mat.possc[4,5]
       mat.possc<-mat.possc[-c(2,4),]} else {
       if (sum(mat.possc[,4])==4 | sum(mat.possc[,4])==8){
        mat.possc[1,5]<-mat.possc[1,5]+mat.possc[3,5]
        mat.possc[2,5]<-mat.possc[2,5]+mat.possc[4,5]
        mat.possc<-mat.possc[-c(3,4),]}}}
     mat.possc<-cbind(mat.possc,matrix(0,nrow(mat.possc),2))
     for (l in 1:nrow(mat.possc)){
      if (mat.possc[l,3]==1 & mat.possc[l,4]==1) mat.possc[l,6]<--1
      if (mat.possc[l,3]==2 & mat.possc[l,4]==2) mat.possc[l,6]<-1
      if (mat.possc[l,3]!= mat.possc[l,4]) mat.possc[l,6]<-0}
     mat.possc[,7]<-1-abs(mat.possc[,6])
     #
     logdens<-numeric(nrow(mat.possc))
     vet.delinea<-mat.delinea[length(founders)+j,]
     for (l in 1:nrow(mat.possc)){
      vet.delinea[2*i]<-mat.possc[l,6]
      vet.delinea[(2*i)+1]<-mat.possc[l,7]
      logdens[l]<-dnorm(dados[length(founders)+j,1],vet.delinea%*%vet.coef,sqrt(sigma2.vig),log=TRUE)}
     prob<-exp(log(mat.possc[,5])+logdens-max(log(mat.possc[,5])+logdens))/sum(exp(log(mat.possc[,5])+logdens-max(log(mat.possc[,5])+logdens)))
     ger.gen<-rDiscreta(prob)
#     ger.gen<-which(prob==max(prob))
     mat.delinea[length(founders)+j,2*i]<-mat.possc[ger.gen,6]
     mat.delinea[length(founders)+j,(2*i)+1]<-mat.possc[ger.gen,7]
     mat.alelo.qtls[length(founders)+j,(2*i-1)]<-mat.possc[ger.gen,3]
     mat.alelo.qtls[length(founders)+j,(2*i)]<-mat.possc[ger.gen,4]}}}
  residuos<-dados[,1]-(mat.delinea%*%vet.coef)
  #
  #### updating mu - general average
  #
  vet.coef[1,1]<-poster.mi(media.mi,sigma2.mi,sigma2.vig,residuos,vet.coef[1,1])[[1]]
  residuos<-dados[,1]-(mat.delinea%*%vet.coef)
  #
  #### updating alphas and deltas
  #
  if (num.QTLs>0){
   for (i in 1:num.QTLs){
    vet.coef[(2*i),1]<-poster.alpha(media.alpha,sigma2.alpha,sigma2.vig,residuos,vet.coef[(2*i),1],mat.delinea[,(2*i)])[[1]]
    residuos<-dados[,1]-(mat.delinea%*%vet.coef)
    vet.coef[(2*i)+1,1]<-poster.delta(media.delta,sigma2.delta,sigma2.vig,residuos,vet.coef[(2*i)+1,1],mat.delinea[,(2*i)+1])[[1]]
    residuos<-dados[,1]-(mat.delinea%*%vet.coef)}}
  #
  #### updating sigma2
  #
  sigma2.vig<-poster.sigma2(neta.a,neta.b,residuos)[[1]]
  #
  ############################## guarda as informações para inferência
  #
  if (int>burnin & int%%saltos==0){
   cat('',num.QTLs,file=paste(caminho,"numero_QTLs_q1_dd2_crom",cromos,".txt",sep=""),append=T)
   cat('',pos.qtls,file=paste(caminho,"posicao_QTLs_q1_dd2_crom",cromos,".txt",sep=""),append=T)
   cat('',vet.coef,file=paste(caminho,"vetor_coeficientes_q1_dd2_crom",cromos,".txt",sep=""),append=T)
   cat('',sigma2.vig,file=paste(caminho,"sigma2_q1_dd2_crom",cromos,".txt",sep=""),append=T)}
}
cat('',indrejtotal,file=paste(caminho,"indrej_q1_dd2_crom",cromos,".txt",sep=""),append=T)
cat('',round(probacetotal,2),file=paste(caminho,"prob_rej_q1_dd2_crom",cromos,".txt",sep=""),append=T)
cat('',round(probacemerge,2),file=paste(caminho,"prob_rej_mg_q1_dd2_crom",cromos,".txt",sep=""),append=T)