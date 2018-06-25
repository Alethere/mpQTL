##############################
## QTL analysis base script ##
##############################

##### Experiment setup ####
source("mpQTL_fun.R") #file containing most of the functions

jobname<-"base_QTL" #name that will be used for the Logfile and for the Result RDS file
link_threshold<-.85 #linkage threshold for the True Positive 
herits<-c(0.2,0.5,0.8) #heritabilities to test
config<-4 #number of genetic configurations
no_cores<-7 #number of cores for parallel computing, 
QTLs<-c(8,1349,5056) #extreme chromosome, high density, low density. QTLs to be tested
gen.type<-"ancestral" #can be biallelic, ancestral or parental
mod.type<-"linear" #can be linear or mixed

totallele<-"PedigreeSIM/Parents/Total_pop.txt" #adress of the file containing the ancestral alleles of each parental chromosome
mapfile<-"PedigreeSIM/Potato/Potato.map"
#Names of the alledosage and founderallele files
crossnames<-paste0("cross",sprintf("%03.f",0:99))
crossfounder<-paste0("PedigreeSIM/NAM crosses/1 ancestral/",crossnames,"_founderalleles.dat")
crossdosage<-paste0("PedigreeSIM/NAM crosses/1 ancestral/",crossnames,"_alleledose.dat")

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
  
  #will calculate the percentage of founder allele linkage of each individual
  link<-linkage(QTLs,crossfounder[i],totallele = totallele) 
  tp<-apply(link,2,function(x){
    which(x>=link_threshold)
  })
  
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
  Result[[crossnames[i]]]<-list(pval,tp)
  saveRDS(Result,file=paste(jobname,".RDS"))
}

end<-Sys.time()
elapsed<-end-start
cat(c("\n\nDone! Time elapsed",round(elapsed,digits=2),attr(elapsed,"units")),file=logfile,append=T)
