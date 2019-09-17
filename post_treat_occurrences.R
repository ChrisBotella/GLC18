# initié 24/10/17
# post-traitement des extractions GBIF en vue de GeoLifeCLef 2018
setwd('C:/Users/Christophe/hubiC/Documents/0_These/Github/R/_base/')
source('functions.R')

data_dir = 'C:/Users/Christophe/hubiC/Documents/0_These/data/geolifecLef/v4.0/'
data_dir = "D:/GeoLifeClef/v4.0/correction csv"
setwd(data_dir)
print(list.files())

#####
# 1) Filtre: occurrences répétées (fusion des extractions)
#####

n_filtres = 2 
for(f in 1:n_filtres){
  #na = paste('DL_gbif_filtre',f,'.csv',sep="")
  na = paste('rm_locality',f,'.csv',sep="")
  
  if(f==1){tab = read.csv(na,sep=";",header=T)
  }else{tab=rbind(tab,read.csv(na,sep=";",header=T))}
}
colnames(tab)

ids = unique(tab$gbifid)
res = NULL
#rem = data.frame(w=1:dim(tab)[1] , id=tab$gbifid)
#sto=NULL
deb=Sys.time()
k=1
for(id in ids){
  #w = which(rem$id==id) 
  #sto = c(sto,w) 
  #res[k] = rem$w[w]
  #if( k/1000==round(k/1000)){ rem=rem[-sto,] ; sto=NULL }
  res[k] = which(tab$gbifid==id)[1]
  k=k+1
  cat('\r  Process ...',100*k/length(ids),' %                               \r')
  flush.console()
}
print(Sys.time()-deb)

tabo = tab[res,]
setwd(data_dir)

write.table(tabo,"tabo.csv",sep=";",col.names=T,row.names=F)

#####
# 2) Virer: carré lat/lon & Lat/Lon inversées
#####

tab = read.csv('tabo.csv',sep=";",header=T)
tabo= tab
w = NULL
for(i in 1:dim(tab)[1]){
  lon = tabo$decimallongitude[i]
  lat = tabo$decimallatitude[i]
  if( (lon>40 & lon<52) & (lat>-6 & lat<11)){
    tabo$decimallatitude[i] = lon
    tabo$decimallongitude[i] = lat
    w=c(w,i)
  }else if( (lat>40 & lat<52) & (lon>-6 & lon<11) ){
    w=c(w,i)
  }
  cat('\r  Process ...',100*i/dim(tabo)[1],' %                               \r')
  flush.console()
}

tabo = tabo[w,]
write.table(tabo,'taboo.csv',sep=";",row.names=F,col.names=T)
dim(tabo)

# Combien de changements de coordonnées?
sum( c(tabo$decimallatitude,
       tabo$decimallongitude)==c(tab$decimallatitude[w],
        tab$decimallongitude[w]) )/2

#####
# 3) Filtre: territoire FR, 
# virer floraine,
# virer matching incertain
#####

setwd(data_dir)
tab = read.csv('taboo.csv',sep=";",header=T)

# premier tri spatial grossier
condi = tab$decimallongitude>-6 & tab$decimallongitude<10 & tab$decimallatitude>41 & tab$decimallatitude<52
tabo = tab[condi,]
colnames(tabo)[colnames(tab)%in%c('decimallongitude','decimallatitude')]=c('Latitude','Longitude')

# On vire Floraine car
# les occ sont aggrégées 
# sur une maille 5x5 
condi = tabo$institutioncode!="Floraine"
tabo =tabo[condi,]

#Eliminer les occurrences dont le matching est fait à un rang supérieur à l'espèce et celles dont le match à l'espèce est fuzzy.
tolerated=c('SPECIES','SUBSPECIES','VARIETY')
tabo = tabo[tabo$taxonrank%in%tolerated,]
tabo=tabo[!(tabo$taxonrank=="SPECIES" & (regexpr('TAXON_MATCH_FUZZY',tabo$issue)>0)),]

# Virer si hors France
tabo = filter_and_correct(tabo)
setwd(data_dir)
write.table(tabo,'tabooo.csv',sep=";",row.names=F,col.names=T)

#####
# 4) Filtre: Trouver le nom retenu d'espèce dans bdtfx et virer sinon
#####

setwd(data_dir)
tab = read.csv('tabooo.csv',sep=";",header=T)

setwd('C:/Users/Christophe/hubiC/Documents/0_These/data/Taxonomy/référentiels/bdtfx4_01du2017_03_15')

ref = read.csv('ref_bdtfx4_01du2017_03_15.csv',sep=";",header=T)

table(ref$Code.rang)


# Code.rang = 
# 290 espèce (dont synomymes)
# 220 genre
# 180 famille

tab$espece_retenue_bdtfx=NA
gbif_names = as.character(unique(tab$scientificname))
signals = data.frame(scientificname=gbif_names,
                     species_bdtfx=NA,
                     false_auth=F,
                     not_retained=F,
                     multi_match=F,
                     no_match=F,
                     subsp=F,
                     higher_rank=F,
                     bug=F)

for(i in 1:dim(signals)[1]){
  na = as.character(signals$scientificname[i])
  print(na)
  # matching dans ref
  id = which(ref$Nom.avec.auteur==na)
  # Y a t'il eu un matching sur Nom+auteur?
  if(length(id)==0){
    print('no match with author')
    # Non, on tente Nom seul
    id = which(ref$Nom.sans.auteur==as.character(tab$species[tab$scientificname==na][1]))
    signals$false_auth[i]=T
  }
  
  # y a t'il au moins un nom retenu parmis Noms qui matchent?
  rete=ref$N.nomenclatural[id] == ref$N.nomenclatural.nom.retenu[id]
  if(sum(rete,na.rm=T)>0){id = id[which(rete)]# oui, on garde ceux là
  }else{signals$not_retained[i]=T    # non
  }
  
  # y-a t'il plusieurs nom.retenus/synonymes matchés?
  if(length(id)>1){
    # oui, on ne garde que celui du niveau espèce 
    # & le plus récent
    print(paste('multi match for',na))
    annees = ref$annee.publi[id]
    codes = ref$Code.rang[id]
    if(290 %in% codes){
      id=id[which(annees==max(annees) & codes==290)]
    }else{
      id=id[which(annees==max(annees))]
    }
    signals$multi_match[i]=T
  }
  
  # n'y a t'il eu aucun nom matché?
  if(length(id)==0){
    # non aucun: le nom est inconnu du référentiel
    print(paste('No match at all for',na))
    valid_name = NA
    signals$no_match[i]=T
  # si 1: le nom matché a un rang inférieur ou égal à l'espèce?
  }else if(ref$Code.rang[id]>=290){
    # oui: on récupère le nom retenu de l'espèce associée
    r.id = reach.higher.rank(id,ref,290)
    if(r.id!=id){print(paste(na,'match rank',ref$Code.rang[id]));signals$subsp[i]=T}
    valid_name = as.character(ref$Nom.retenu.avec.auteur[r.id])
  }else{
    #non, il est au dessus, alors on ne sait pas retrouver l'espèce
    print(paste(ref$Code.rang[id],' rank match for',na))
    valid_name=NA
    signals$higher_rank[i]=T
  }
  if(length(valid_name)>1){signals$bug[i]=T}
  tab$espece_retenue_bdtfx[tab$scientificname==na]=valid_name
  signals$species_bdtfx[i]=valid_name
}

signals$species_bdtfx=sapply(1:dim(signals)[1],function(i) tab$espece_retenue_bdtfx[tab$scientificname==signals$scientificname[i]][1])
# 227 scientificnames ne matchent pas

tabo = tab[!is.na(tab$espece_retenue_bdtfx),]

setwd(data_dir)
write.table(tabo,'taboooo.csv',sep=";",row.names=F,col.names=T)

#####
# 5) Filtre: On vire les doublons espXgeoloc
#####

setwd(data_dir)
tab = read.csv('taboooo.csv',sep=";",header=T)

check = aggregate( list(count=rep(1,dim(tab)[1])), by=list(lon=tab$Longitude,lat=tab$Latitude,esp=tab$espece_retenue_bdtfx),sum)
# cas à problèmes
checko= check[check$count>1,]
vire=NULL
for(i in 1:dim(checko)[1]){
  ids = which(tab$Longitude==checko$lon[i] & tab$Latitude==checko$lat[i] & tab$espece_retenue_bdtfx==checko$esp[i])
  vire = c(vire,sample(ids,length(ids)-1) )
  #print(tab[ids,c('Longitude','Latitude','espece_retenue_bdtfx','datasetkey')])
  cat('\r Process... ',100*i/dim(checko)[1],'%                  \r')
  flush.console()
}
print(length(vire))
tabo=tab[-vire,]
rm(check,checko)
gc(reset=T)

setwd(data_dir)
write.table(tabo,'tabooooo.csv',sep=";",row.names=F,col.names=T)


#####
# 6) Virer les espèces <5 observations
#####

setwd(data_dir)
tab = read.csv('tabooooo.csv',sep=";",header=T)

esp = table(tab$espece_retenue_bdtfx)
esp = as.character(names(esp[esp>=5]))
length(esp)
tabo = tab[tab$espece_retenue_bdtfx%in%esp,]

write.table(tabo,'taboooooo.csv',sep=";",col.names=T,row.names=F)

#####
# 7) Extraction d'un test set
#####

setwd(data_dir)
tab = read.csv('taboooooo.csv',sep=";",header=T)

esp = unique(as.character(tab$espece_retenue_bdtfx))
# calculer les distances par espèce à l'avance
d = rep(0,dim(tab)[1])
deb=Sys.time()
tol_metre = 100
for(e in esp){
  cat('\r Process species ',100*which(esp==e)/length(esp),'          \r')
  flush.console()
  ids = which(tab$espece_retenue_bdtfx==e)
  tabi = tab[ids,]
  i=1
  for(w in ids){
    d_lon = abs(tabi$Longitude[-i]-tabi$Longitude[i])*(pi*6.371e6)/180
    idss = which.min(d_lon)
    if(min(d_lon[idss])>=tol_metre){
      d_lat = abs(tabi$Latitude[-i][idss] -tabi$Latitude[i])*(pi*6.371e6*cos(tabi$Latitude[-i][idss]*pi/180))/180
      d[w] = min( sqrt(d_lat^2 + d_lon^2) )
    }
    i=i+1
  }
}
# nombre d'occurrences potentielles
sum(d>=tol_metre)
pot = which(d>=tol_metre)

# 1/4 des occ
# Garder des obs train pour chaque espèce
# Isoler spatialement les test
tab$test= F
train = 1:dim(tab)[1]
test = NULL
p=0.25*dim(tab)[1]/sum(d>=tol_metre)
while (length(test) <= dim(tab)[1]/4){
  # on tire une occurrence au hasard parmis les elligibles
  i = sample(pot,1)
  
  # on la retient avec proba .25
  if(runif(1,0,1)<p){
    # si retenue, on garde si il reste des obs train pour la même espèce 
    ids = train[train!=i]
    #esp en train
    keep2 = tab$espece_retenue_bdtfx[i]%in%tab$espece_retenue_bdtfx[ids]
    if(keep2){
      test=c(test,i)
      train=setdiff(train,i)
      pot = setdiff(pot,i)
      print(length(test))
    }
  }
}

tab$test[test]=T
tab$patch_id[tab$test]=1:sum(tab$test)
tab$patch_id[!tab$test]=1:sum(!tab$test)



setwd(data_dir)
# occurrences totales
write.table(tab,'occurrences.csv',sep=";",row.names=F,col.names=T)


#####
# 8) extraires variables ponctuelles pour occurrences.csv
#####

setwd(data_dir)
tab=read.csv('occurrences.csv',sep=";",header=T)
variables = as.character(load_w_variables()$variables)

tabo=get_variables(variables,tab,dpath = "C:/Users/Christophe/hubiC/Documents/0_These/data/",d_NA = 4000)

species=as.character(unique(tabo$espece_retenue_bdtfx))
tabo$species_glc_id=sapply(1:dim(tabo)[1],function(i) which(species==tabo$espece_retenue_bdtfx[i]) )

setwd(data_dir)
# full
write.table(tabo,'occurrences.csv',sep=";",row.names=F,col.names=T)
# train 
write.table(tabo[!tabo$test,],'occurrences_train.csv',sep=";",row.names=F,col.names=T)
#test 
te.tabo = tabo[tabo$test,]
te.tabo$glc_id = 1:dim(te.tabo)[1]
given = c("glc_id","datasetkey","publishingorgkey","Latitude","Longitude","coordinateuncertaintyinmeters","coordinateprecision","elevation","elevationaccuracy","depth","depthaccuracy","eventdate","day","month","year","basisofrecord","institutioncode","collectioncode","catalognumber","identifiedby",
          "license","rightsholder","recordedby","typestatus","establishmentmeans","lastinterpreted","mediatype",variables)
write.table(te.tabo[,given],'occurrences_test.csv',sep=";",row.names=F,col.names=T)

# ground truth
gt_file = te.tabo[,c('glc_id','species_glc_id')]
write.table(gt_file,'gt_file.csv',sep=";",row.names=F,col.names=F)

#####
# 9) Création des patchs
#####

setwd(data_dir)
tab=read.csv('occurrences.csv',sep=";",header=T)

variables = load_w_variables()

# train patch images of size 64*64
writeAndSort(tab[(!tab$test) & tab$added,],variables,size=64,dirnames = "patch_dirname",tiff_id="patch_id",opath='D:/GeoLifeClef/v4.0/patchTrain/')

# test patch 
writeAndSort(tab[tab$test & tab$added,],variables,size=64,dirnames = "patch_dirname",tiff_id="patch_id",opath='D:/GeoLifeClef/v4.0/patchTest/')

train = read.csv('occurrences_train.csv',sep=";",header=T)
colnames(train)
sum(train$test)

#####
# Definir les patchs train à créer (Version réduite => non retenue)
#####

setwd(data_dir)
tab=read.csv('occurrences_train.csv',sep=";",header=T)

setwd('C:/Users/Christophe/hubiC/Documents/0_These/data/0_mydata/proxi_eau_fast/')
list.files()
r = raster('proxi_eau_fast.tif')

# on donne un même id a toutes les occurrences qui tombent dans la même case de ce raster
xmn=extent(r)[1]
xmx=extent(r)[2]
ymn=extent(r)[3]
res_x = resolution_x(r)
res_y=resolution_y(r)



row = (tab$Latitude - ymn) %/% res_y
col = (tab$Longitude - xmn) %/% res_x
ncell = ncol(r)*row + col 
patch = data.frame(patch_id=1:length(unique(ncell)),cell=unique(ncell))
print(dim(patch)[1])

tab$patch_id = NA
for(i in 1:dim(patch)[1]){
  tab$patch_id[ncell==patch$cell[i]]=patch$patch_id[i]
}
print(length(unique(tab$patch_id)))
sumup = table(tab$patch_id)
print(sumup[order(sumup,decreasing=T)][1:50])

patch$Longitude = sapply(1:dim(patch)[1],function(i) mean(tab$Longitude[tab$patch_id==patch$patch_id[i]]))
patch$Latitude = sapply(1:dim(patch)[1],function(i) mean(tab$Latitude[tab$patch_id==patch$patch_id[i]]))

setwd(data_dir)
write.table(tab,'occurrences_train.csv',sep=";",row.names=F,col.names=T)
write.table(patch,'patch_train.csv',sep=";",row.names=F,col.names=T)


#####
# Stats espèces
#####

setwd(data_dir)
tab=read.csv('occurrences.csv',sep=";",header=T)

# species
esp = table(tab$espece_retenue_bdtfx)
library('ggplot2')
p=ggplot( as.data.frame(esp),aes(x=1:length(esp),y=as.integer(esp[order(esp,decreasing=T)])))
p=p+geom_line()+ theme(panel.grid.minor = element_line(colour="white", size=0.5))
p=p+scale_y_continuous(minor_breaks=seq(0,max(esp),100),breaks =seq(0,max(esp),200))
p=p+labs(x='Species',y='Number of observations.')
print(p)
plot(1:length(esp),as.integer(esp[order(esp,decreasing=T)]),type="l")

# datasets
bdd = table(tab$datasetkey)
print(bdd[order(bdd,decreasing = T)])

length(unique(tab$patch_id))

# map
library('RgoogleMaps')

?RgoogleMaps
center = 
GetMap()



geoloc = aggregate(list(count=rep(1,dim(occ)[1])),by=list(lat = occ$Latitude,lon = occ$Longitude),sum)
print(geoloc[order(geoloc$count,decreasing=T),][1:100,])

geoEsp = aggregate(list(count=rep(1,dim(occ)[1])),by=list(lat = occ$Latitude,lon = occ$Longitude,esp=occ$espece_retenue_bdtfx),sum)
sum(geoEsp$count>1)
max(geoEsp$count)

print(geoEsp[order(geoEsp$count,decreasing=T),][1:100,],row.names=F)

geoEspX = aggregate(list(count=geoEsp$count),by=list(lon=geoEsp$lon,lat=geoEsp$lat),max)
sum(geoEspX$count>1)

newOcc = aggregate(list(count=rep(1,dim(occ)[1])),by=list(lat = occ$Latitude,lon = occ$Longitude,esp=occ$espece_retenue_bdtfx,dat=occ$datasetkey),sum)


#####

setwd(data_dir)
tab = read.csv('taboooooo.csv',sep=";",header=T)

setwd("C:/Users/Christophe/hubiC/Documents/0_These/data/geolifecLef/v4.0/")
list.files()
tab2 = read.csv('occurrences.csv',sep=";",header=T)
id = which(! (paste(tab$Longitude,tab$Latitude) %in% paste(tab2$Longitude,tab2$Latitude)) )


tab = rbind(tab,  tab2[!tab2$gbifid%in%tab$gbifid ,colnames(tab)] )

tests = which(tab2$test)
idtest = sapply(tests,function(i) which(tab$gbifid==tab2$gbifid[i]))
length(idtest)==sum(tab2$test)
tab$test=F
tab$test[idtest]=T

tab$patch_id=NA
for(i in 1:dim(tab2)[1]){
  tab$patch_id[tab$gbifid==tab2$gbifid[i]] = tab2$patch_id[i]
}

# calculer les distances par espèce à l'avance
d = rep(0,dim(tab)[1])
deb=Sys.time()
tol_metre = 100
for(e in esp){
  cat('\r Process species ',100*which(esp==e)/length(esp),'          \r')
  flush.console()
  ids = which(tab$espece_retenue_bdtfx==e)
  tabi = tab[ids,]
  i=1
  for(w in ids){
    d_lon = abs(tabi$Longitude[-i]-tabi$Longitude[i])*(pi*6.371e6)/180
    idss = which.min(d_lon)
    if(min(d_lon[idss])>=tol_metre){
      d_lat = abs(tabi$Latitude[-i][idss] -tabi$Latitude[i])*(pi*6.371e6*cos(tabi$Latitude[-i][idss]*pi/180))/180
      d[w] = min( sqrt(d_lat^2 + d_lon^2) )
    }
    i=i+1
  }
}
# nombre d'occurrences potentielles
sum(d>=tol_metre)
pot = which(d>=tol_metre & tab$added)

# 1/4 des occ
# Garder des obs train pour chaque espèce
# Isoler spatialement les test
train = which(!tab$test)
test = which(tab$test)
p=0.25*dim(tab)[1]/sum(d>=tol_metre)
while (length(test) <= dim(tab)[1]/4){
  # on tire une occurrence au hasard parmis les elligibles
  i = sample(pot,1)
  
  # on la retient avec proba .25
  if(runif(1,0,1)<p){
    # si retenue, on garde si il reste des obs train pour la même espèce 
    ids = setdiff(train,i)
    #esp en train
    keep2 = tab$espece_retenue_bdtfx[i]%in%tab$espece_retenue_bdtfx[ids]
    if(keep2){
      test=c(test,i)
      train=setdiff(train,i)
      pot = setdiff(pot,i)
      print(length(test))
    }
  }
}

tab$test[test]=T



n=max(tab$patch_id[tab$test],na.rm=T)
tab$patch_id[tab$test & is.na(tab$patch_id)]=(n+1):sum(tab$test)
n=max(tab$patch_id[!tab$test],na.rm=T)
tab$patch_id[!tab$test & is.na(tab$patch_id)]=(n+1):sum(!tab$test)

setwd(data_dir)
# occurrences totales
write.table(tab,'occurrences.csv',sep=";",row.names=F,col.names=T)

variables = as.character(load_w_variables()$variables)

tab=get_variables(variables,tab,dpath = "C:/Users/Christophe/hubiC/Documents/0_These/data/",fix=F)

species=as.character(unique(tab$espece_retenue_bdtfx))
tab$species_glc_id=sapply(1:dim(tab)[1],function(i) which(species==tab$espece_retenue_bdtfx[i]) )

setwd(data_dir)
# full

tab$patch_dirname = as.character(256* -(-tab$patch_id %/% 256) )
tab=tab[,-which(colnames(tab)=='added')]#virer col added 
write.table(tab,'occurrences.csv',sep=";",row.names=F,col.names=T)

# train 
write.table(tab[!tab$test,],'occurrences_train.csv',sep=";",row.names=F,col.names=T)

# test 
given = c("patch_id",'patch_dirname',"datasetkey","publishingorgkey","Latitude","Longitude","coordinateuncertaintyinmeters","coordinateprecision","elevation","elevationaccuracy","depth","depthaccuracy","eventdate","day","month","year","basisofrecord","institutioncode","collectioncode","catalognumber","identifiedby",
          "license","rightsholder","recordedby","typestatus","establishmentmeans","lastinterpreted","mediatype",as.character(variables$variables))
test = tab[tab$test,given]
write.table(test,"occurrences_test.csv",sep=";",row.names=F,col.names=T)


