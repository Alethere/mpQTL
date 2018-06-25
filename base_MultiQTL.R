####################################
## MULTI-ANCESTRAL QTL base script##
####################################

##### Experiment setup ####
source("mpQTL_fun.R") #file containing most of the functions

jobname<-"base_MultiQTL" #name that will be used for the Logfile and for the Result RDS file
link_threshold<-.85 #linkage threshold for the True Positive 
herits<-c(0.2,0.5,0.8) #heritabilities to test
no_cores<-7 #number of cores for parallel computing, 
gen.type<-"biallelic" #can be biallelic, ancestral or parental
mod.type<-"linear" #can be linear or mixed



totallele<-"PedigreeSIM/Parents/Total_pop.txt" #adress of the file containing the ancestral alleles of each parental chromosome
mapfile<-"PedigreeSIM/Potato/Potato.map"
#Names of the alledosage and founderallele files
ancestrals<-c(3,7,10)
pop.number<-25 #number of populations to test for each ancestral group
crossfounder<-c();crossdosage<-c()
for(i in ancestrals){
  cross<-paste0("cross",sprintf("%03.f",0:(pop.number-1)+100*(i-1)))
  crossfounder<-c(crossfounder,
                  paste0("PedigreeSIM/NAM crosses/",i," ancestral/",cross,"_founderalleles.dat"))
  crossdosage<-c(crossdosage,
                 paste0("PedigreeSIM/NAM crosses/",i," ancestral/",cross,"_alleledose.dat"))
}
crossnames<-substr(crossfounder,nchar(crossfounder)-26,nchar(crossfounder)-19)

#Matrix that identifies the ancestral group of each allele
anc_alleles<-sapply(1:10,function(x) 0:39+(x-1)*40)
colnames(anc_alleles)<-paste0("A",0:9)

#Logfile to track the development of the program
logfile<-paste0(jobname,".txt")
talk<-function(words,file=logfile){
  cat(paste0(words),date(),"\n",file=file,append = T)
}
write(paste0("Experiment ",jobname," started on: ",date(),
             "\nCrosses: ",paste(crossnames,collapse = " "),
             "\nQTLs: ",paste(QTLs,collapse=" "),
             "\nAllelic configurations: ",config,
             "\nHeritabilities: ",paste(herits,collapse=" "),
             "\nModel: ",gen.type," ",mod.type,"\n\n"),file=logfile)



#### Program run ####
start<-Sys.time()

Result<-list()
for(i in 1:length(crossnames)){
  talk(paste("\nStarting",crossnames[i]))
  
  ## Detection of ancestral groups ##
  found_alleles<-link_NAM(crossfounder[i])
  n.anc<-sapply(1:nrow(found_alleles),function(k){
    a<-unique(found_alleles[k,])
    b<-sapply(a,function(a) colSums(a==anc_alleles) )
    n.anc<-sum(rowSums(b)!=0)
    return(n.anc)
  })
  QTL_num<-unique(n.anc)#in our data, the number of ancestral groups at each locus is the same. 
  #So here is only one n.anc, and QTL_num wil lbe equal to the number of QTLs to simulate
  
  #sample the QTLs randomly from the genome
  map<-data.table::fread("PedigreeSIM/Potato/Potato.map")
  QTLs<-sample.QTL(QTL_num,map,seed=7)
  found_alleles<-found_alleles[QTLs$index,] #we just keep the QTLs we will use for simulation
  
  #will calculate the percentage of founder allele linkage of each individual
  link<-linkage(QTLs$index,crossfounder[i],totallele = totallele) 
  tp<-apply(link,2,function(x){
    which(x>=link_threshold)
  })
  
  ### Phenotype ###
  talk("Phenotypes sim")
  effects<-multi.effects(found_alleles,anc_alleles,n=1)#obtain effects for the alleles at each locus
  un_alleles<-apply(found_alleles,1,function(x) length(unique(x)))
  found_alleles<-tidy_allele(found_alleles)
  
  seed<-1:length(herits)
  phenotypes<-sapply(1:length(herits),function(j){
    multi.pheno(found_alleles,#all QTLs at once
                effects,herits[j],
                seed=seed[j])
  })
  colnames(phenotypes)<-herits
  
  ### QTL analysis ###
  talk("QTL calc start")
  
  if(gen.type=="biallelic"){
    dosages<-as.matrix(data.table::fread(crossdosage[i]))
    markers<-dosages[,1]; inds<-colnames(dosages[,-1])
    dosages<-dosages[,-1] #take out marker column
    segregants<-sapply(1:nrow(dosages),function(x) segregation<-length(unique(dosages[x,-1])) )
    dosages<-matrix(as.numeric(dosages[which(segregants!=1),]),ncol=ncol(dosages))
    colnames(dosages)<-inds; rownames(dosages)<-markers[which(segregants!=1)]
    dosX<-T #whether to apply the dosage.X function or not
    
  }else if(gen.type=="ancestral"){
    dosages<-link_NAM(crossfile=crossfounder[i],totallele=totallele)
    dosX<-F
    
  }else if(gen.type=="parental"){
    dosages<-as.matrix(data.table::fread(crossfounder[i],header=T)[,-1])
    dosX<-F
  }else{
    stop("Wrong genetic model type")
  }
  
  if(mod.type=="linear"){
    #In case a linear model is selected, the following procedure is done
    cluster<-parallel::makeCluster(no_cores)
    parallel::clusterExport(cl=cluster,c("phenotypes","dosages","dosage.X","dosX"))
    pval<-parallel::parSapply(cl=cluster,1:nrow(dosages),function(i){
      #each lm calculates residuals for all phenotypes at once. We can calculate the pval using the SSR and SSM
      test<-lm(phenotypes~dosage.X(dosages[i,],dosages = dosX))
      SSR<-colSums((test$residuals)^2) #Sum of Squares Residuals
      SSM<-sapply(1:ncol(phenotypes),function(x){ #Sum of Squares of Model
        mean<-mean(phenotypes[,x])
        sum((test$fitted.values[,x]-mean)^2)
      })
      
      Ftest<-SSM/(SSR/test$df.residual)
      pf(Ftest,1,test$df.residual,lower.tail = F)
      
    }) 
    parallel::stopCluster(cluster)

  }else if(mod.type=="mixed"){
    #in a mixed model approach we need to calculate the K and Hinv matrices. 
    #according to Giorgio, using all markers for a K calculation can lead to overcorrection, therefore we must
    #sample a limited number of markers. 
    relat<-as.matrix(data.table::fread(crossdosage[i])[,-1])
    map<-data.table::fread(mapfile,header=T)
    relat<-sample.cM(relat,map)
    
    K<-calc.K(t(relat)) #markers should be in columns
    Z<-diag(nrow(phenotypes)) #nrow(phenotypes) is equal to number of individuals
    Hinv<-calc.Hinv(phenotypes,X=matrix(rep(1,nrow(phenotypes))),Z,K) #X=vecor of 1s, to apply P3D/EMMAX algorithm
    
    cl<-parallel::makeCluster(no_cores)
    parallel::clusterExport(cl,list("phenotypes","Z","K","Hinv",
                                    "dosages","mm.solve","dosage.X","dosX"))
    
    pval<-parallel::parSapply(cl,1:ncol(phenotypes),function(i){#parallely, over each phenotype
      sapply(1:nrow(dosages),FUN=function(k){#calculate the pvalue for each marker
        X<-dosage.X(dosages[k,],dosages = dosX) #obtain the X matrix for that marker
        if(ncol(X)>2){X<-X[,-2]}#the -2 is to take out one of the columns of X, to prevent singularity when solving
        mm.solve(phenotypes[,i],X,Z,K,Hinv[[i]])$pval
      })
    })
    
    parallel::stopCluster(cl)
    
  }else{stop("Wrong statistical model specificaton")}
  
  talk("Saving results")
  tp<-rep(tp,each=config*length(herits)) #this is done just so that Result[[i]][1] (pvals) ==result[[1]][2] (tp)
  Result[[i]]<-list(pval,tp)
  saveRDS(Result,file=paste(jobname,".RDS"))
}

end<-Sys.time()
elapsed<-end-start
cat(c("\n\nDone! Time elapsed",round(elapsed,digits=2),attr(elapsed,"units")),file=logfile,append=T)
