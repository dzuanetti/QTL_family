#
gen_snps<-scan(file=paste(caminho,"dadosgaw17.txt",sep=""))
#
nindiv<-697
nvar<-13791
#
dados.compl<-matrix(gen_snps,nrow=nindiv,ncol=nvar,byrow=TRUE)
dados.compl<-dados.compl[,-c(1,2,6,7)]
# layout base: id, pai, mae, snps
#
fenot<-read.csv(file=paste(caminho,"fam_phen_1.txt",sep=""))
fenot<-fenot[,-c(2,3,4,8)]
# layout base: id, q1, q2, q4
#
num_snps<-c(1136,800,575,541,565,856,594,596,757,883,1447,842,208,394,465,426,654,332,990,313,137,273)
#
#####################################
########## organizing family structure
#####################################
#
organiza<-numeric()
for (i in 1:nrow(dados.compl)){
  if (dados.compl[i,2]==0 & dados.compl[i,3]==0) organiza<-c(organiza,i)}
dados.aux<-dados.compl[organiza,]
dados.aux<-rbind(dados.aux,dados.compl[-organiza,])
#
fenot.aux<-fenot[organiza,]
fenot.aux<-rbind(fenot.aux,fenot[-organiza,])
#
dados.compl<-dados.aux
fenot<-fenot.aux
familia<-dados.compl[,1:3]
#
for (i in 1:nrow(familia)){
  pai<-familia[i,2]
  mae<-familia[i,3]
  if (pai!=0){
    for (j in 1:nrow(familia)) if (pai==familia[j,1]) familia[i,2]<-j}
  if (mae!=0){
    for (j in 1:nrow(familia)) if (mae==familia[j,1]) familia[i,3]<-j}}
#
#####################################
########## specify chromossome 1
#####################################
cromos<-1
num_snps_acum<-c(0,cumsum(num_snps))
snps<-dados.compl[,c(1,2,3,(num_snps_acum[cromos]+4):(num_snps_acum[cromos+1]+3))]
for (j in (4:(num_snps[cromos]+3))){
  for (i in (1:697)){
    if (snps[i,j]==0) snps[i,j]<--1 else {
      if (snps[i,j]==1) snps[i,j]<-0 else {
        if (snps[i,j]==2) snps[i,j]<-1}}}}
#
mapa_marc<-read.table(file=paste(caminho,"mapa_crom",cromos,".txt",sep=""), sep="\t", header=TRUE)
loc.marc<-c(mapa_marc[,2]/100)
loc.marc<-loc.marc-loc.marc[1]
#
dados<-cbind(fenot[,2],snps[,4:ncol(snps)])
gen.marc<-dados[,-1]
#
########## exclude SNPS with low variability
#
bfreq1<-NULL
for (i in 1:ncol(gen.marc)) bfreq1[i]<-round(max(table(gen.marc[,i])/697*100),1)
excl<-which(bfreq1>90)
gen.marc<-gen.marc[,-excl]
nmarc<-ncol(gen.marc)
loc.marc<-loc.marc[-excl]
loc.marc<-loc.marc-loc.marc[1]
dados<-dados[,-(excl+1)]
#
  tira<-which(is.na(loc.marc))
  if (length(tira)>0){
    gen.marc<-gen.marc[,-tira]
    nmarc<-ncol(gen.marc)
    loc.marc<-loc.marc[-tira]
    loc.marc<-loc.marc-loc.marc[1]
    dados<-dados[,-(tira+1)]}
  #
  tira_i_pos<-rep(0,length(loc.marc))
  for (i in 1:(length(loc.marc)-1)){
    if (loc.marc[i]==loc.marc[i+1]) tira_i_pos[i+1]<-1}
  tira_i_pos<-which(tira_i_pos>0)
  if (length(tira_i_pos)>0){
    gen.marc<-gen.marc[,-tira_i_pos]
    nmarc<-ncol(gen.marc)
    loc.marc<-loc.marc[-tira_i_pos]
    loc.marc<-loc.marc-loc.marc[1]
    dados<-dados[,-(tira_i_pos+1)]}
#
nmarc<-ncol(gen.marc)
nfoun<-sum(familia[,2]==0 & familia[,3]==0)
nfil<-nrow(gen.marc)-nfoun
founders<-seq(1:nfoun)
nonfounders<-seq(1:nfil)+nfoun
mat.ped.pais<-familia[(nfoun+1):nrow(dados),2:3]
#
marc.dist<-numeric(length(loc.marc))
marc.dist[1]<-0
for (j in 2:length(loc.marc)) marc.dist[j]<-loc.marc[j]-loc.marc[j-1]
#
tx.rec<-numeric(length(loc.marc))
for (j in 1:length(loc.marc)) tx.rec[j]<-Haldane(marc.dist[j])
tx.rec.c<-1-tx.rec
#
################ determina os alelos maternos e paternos dos founders
#
set.seed(3498)
alelo.materno.fou<-matrix(99,nfoun,nmarc)
alelo.paterno.fou<-matrix(99,nfoun,nmarc)
#
# first marker
for (j in 1:length(founders)){
 if (gen.marc[j,1]==-1){
  alelo.materno.fou[j,1]<-1
  alelo.paterno.fou[j,1]<-1}
 if (gen.marc[j,1]==1){
  alelo.materno.fou[j,1]<-2
  alelo.paterno.fou[j,1]<-2}
 if (gen.marc[j,1]==0){
  aux<-rDiscreta(c(1/2,1/2))-1
  alelo.materno.fou[j,1]<-1^(aux)*2^(1-aux)
  alelo.paterno.fou[j,1]<-2^(aux)*1^(1-aux)}}
#
# remaining markers
for (j in 2:nmarc){
 mat.prob2<- matrix(c(tx.rec.c[j],tx.rec[j],tx.rec[j],tx.rec.c[j]),2,2,byrow=TRUE)# esse vetor tem a Pr(marcador do pai ! marcador anterior do pai), ou seja, usa a dist‚ncia entre os dois marcadores e depois corrige as probabilidades considerando o alelo que vem do marcador do filho
 for (k in 1:length(founders)){
  if (gen.marc[k,(j)]==-1){
   alelo.materno.fou[k,j]<-1
   alelo.paterno.fou[k,j]<-1}
  if (gen.marc[k,(j)]==1){
   alelo.materno.fou[k,j]<-2
   alelo.paterno.fou[k,j]<-2}
  if (gen.marc[k,(j)]==0){
    aux<-rDiscreta(mat.prob2[alelo.paterno.fou[k,j-1],])-1
    alelo.materno.fou[k,j]<-1^(aux)*2^(1-aux)
    alelo.paterno.fou[k,j]<-2^(aux)*1^(1-aux)}}}
#
################ determina os alelos maternos e paternos dos nonfounders
#
alelo.materno<-matrix(99,length(nonfounders),nmarc)
alelo.paterno<-matrix(99,length(nonfounders),nmarc)
for (j in 1:length(nonfounders)){
 for (k in 1:nmarc){
  aux<-identalel_pai_mae(gen.marc[mat.ped.pais[j,1],k],gen.marc[mat.ped.pais[j,2],k],gen.marc[nonfounders[j],k])
  alelo.paterno[j,k]<-aux[[1]]
  alelo.materno[j,k]<-aux[[2]]}}
#
alelo.materno<-rbind(alelo.materno.fou,alelo.materno)
alelo.paterno<-rbind(alelo.paterno.fou,alelo.paterno)
#
#table(alelo.materno,alelo.paterno,gen.marc)
#
############# montando as bases com indicadores de meiose dos individuos nonfounders
#
# Sp==1 ou Sm==1 indica que o alelo e do avooo
# Sp==2 ou Sm==2 indica que o alelo e da avoaa
#
Sp<-matrix(99,length(nonfounders),nmarc)
Sm<-matrix(99,length(nonfounders),nmarc)
for (j in 1:length(nonfounders)){
 for (k in 1:nmarc){
  if (alelo.paterno[nonfounders[j],k]==alelo.paterno[mat.ped.pais[j,1],k]) Sp[j,k]<-1 else Sp[j,k]<-2
  if (alelo.materno[nonfounders[j],k]==alelo.paterno[mat.ped.pais[j,2],k]) Sm[j,k]<-1 else Sm[j,k]<-2}}