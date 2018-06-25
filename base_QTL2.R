##############################
## QTL analysis base script ##
##############################

##### Experiment setup ####
source("mpQTL_fun.R") #file containing most of the functions

jobname<-"base_QTL2.R" #name that will be used for the Logfile and for the Result RDS file
link_threshold<-.85 #linkage threshold for the True Positive 
herits<-c(0.2,0.5,0.8) #heritabilities to test
config<-4 #number of genetic configurations
no_cores<-7 #number of cores for parallel computing, 
QTLs<-c(1349,5056) #extreme chromosome, high density, low density. QTLs to be tested

totallele<-"PedigreeSIM/Parents/Total_pop.txt" #adress of the file containing the ancestral alleles of each parental chromosome
mapfile<-"PedigreeSIM/Potato/Potato.map"
#Names of the alledosage and founderallele files
crossnames<-paste0("cross",sprintf("%03.f",0:20))
crossfounder<-paste0("PedigreeSIM/NAM_crosses/1_ancestral/",crossnames,"_founderalleles.dat")
crossdosage<-paste0("PedigreeSIM/NAM_crosses/1_ancestral/",crossnames,"_alleledose.dat")

#Logfile to track the development of the program
logfile<-paste0("Storage/",jobname,".txt")
talk<-function(words,file=logfile){
  cat(paste0(words),date(),"\n",file=file,append = T)
}
write(paste0("Experiment ",jobname," started on: ",date(),
             "\nCrosses: ",paste(crossnames,collapse = " "),
             "\nQTLs: ",paste(QTLs,collapse=" "),
             "\nAllelic configurations: ",config,
             "\nHeritabilities: ",paste(herits,collapse=" ")),file=logfile)



#### Program run ####
start<-Sys.time()

Result<-list()
for(i in 1:length(crossnames)){
  talk(paste("\nStarting",crossnames[i]))
  
  ### Phenotype ###
  talk("Phenotypes sim")
  #Read the founder allele file, count the number of segregating alleles, and then "tidying them" (cocatenate per individual)
  found_alleles<-link_NAM(crossfile=crossfounder[i],totallele=totallele)[QTLs,]
  un_alleles<-apply(found_alleles,1,function(x) length(unique(x)))
  found_alleles<-tidy_allele(found_alleles)
  
  #Design effects, according to the number of configurations. 
  #Each marker has a different number of unique alleles, so we need to calculate the effects per marker
  c<-3.3^(0:(config-1))#the 3.3 has been obtained empirically. Might need to change it
  effects<-lapply(un_alleles,function(un_alleles){
    x<-sapply(c,function(x){
      eff<-((un_alleles-1):0)^x
      max<-max(eff)
      round(eff/max,digits=2)
    })
    colnames(x)<-paste0("eff_",1:config)
    return(x)
  })
  
  seed<-1:(length(herits)*config*length(QTLs))#for each phenotype a different seed.
  phenotypes<-lapply(1:nrow(found_alleles),function(i){
    x<-lapply(1:length(herits),function(j){
      sapply(1:ncol(effects[[i]]),function(x){
        pheno(found_alleles[i,], #alleles
              effects=effects[[i]][,x], #effects
              Evar = 10, #Environmental variance
              h2=herits[j], #heritability
              mu=50,#average phenotype
              seed = seed[x+12*(i-1)+4*(j-1)] #set.seed, one different for each
        )[-1:-un_alleles[i]]#take the first elements of the phenotype vector, as they contain the final values of gen effects
      })
    })
    
    x<-do.call(cbind,x)
    colnames(x)<-sapply(herits,function(h2) paste0("h",h2,"_",1:config)) #set proper colnames
    return(x)
  })
  #Create a single matrix with all phenotypes, in each column
  phenotypes<-do.call(cbind,phenotypes)
  QTLnum<-as.integer(paste0(col(sapply(1:length(QTLs),rep,times=config*length(herits)))))
  colnames(phenotypes)<-paste0("QTL",QTLnum,colnames(phenotypes))
  
  ### QTL analysis ###
  talk("QTL calc start")
  

  #Genotype parametrization
  Res1<-list()
  for(gen.type in c("biallelic","ancestral","parental")){
    if(gen.type=="biallelic"){
      genotypes<-bi.geno(crossdosage[i])
      talk("Biallelic models")
      
    }else if(gen.type=="ancestral"){
      #ancestral genotypes are obtained
      genotypes<-link_NAM(crossfile=crossfounder[i],totallele=totallele)
      talk("Ancestral models")
      
    }else if(gen.type=="parental"){
      #parental chromosome genotypes are obtained
      genotypes<-as.matrix(data.table::fread(crossfounder[i],header=T)[,-1])
      talk("Parental models")
    }
    map<-data.table::fread(mapfile)
    #pvalue calculation with 4 models
    pval<-lapply(1:4,function(x){
      if(x==1){mixed=NULL;pop=NULL}
      if(x==2){mixed=NULL;pop=T}
      if(x==3){mixed=T;pop=NULL}
      if(x==4){mixed=T;pop=T}
      pval<-map.QTL(phenotypes,
                    genotypes,
                    K=mixed,
                    Q=pop,
                    dosage=as.matrix(data.table::fread(crossdosage[1])[,-1]),
                    map=map)
      
      return(pval)})
    
    names(pval)<-c("linear","linearQ","mixed","mixedQ")
    
    Res1[[gen.type]]<-pval
    
  }
  
  Result[[crossnames[i]]]<-Res1

  saveRDS(Result,file=paste(jobname,".RDS"))
}

end<-Sys.time()
elapsed<-end-start
cat(c("\n\nDone! Time elapsed",round(elapsed,digits=2),attr(elapsed,"units")),file=logfile,append=T)
