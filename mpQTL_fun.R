#These functions require the following packages:
#data.table, stringr

##### Pedigree SIM ######
makedir<-function(folder){if(!dir.exists(folder)){dir.create(folder)}}

#' Make PedSim Parameter File
#'
#' @param output adress of the output to be created (without extension),
#' @param chrom adress of the chromosome file. If not specified, it will assume the same name as output.
#' @param pedfile adress of the pedigree (.ped) file. If not specified, it will assume the same name as output.
#' @param mapgen adress of the .map file. If not specified, it will assume the same name as output.
#' @param founder adress of the .gen file. If not specified, it will assume the same name as output.
#' @param outname adress of the output files. If not specified, it will assume the same name as output.
#' @param ploidy ploidy level to be used
#' @param mapfun mapping function to be used. Either "HALDANE" or "KOSAMBI"
#' @param miss character indivating missing data
#' @param pairing pairing type to be used. 1 if probabilistic, 0 if determined by .chrom file. 
#'
#' @return creates a parameter for PedigreeSim.
#' @export
#'
#' @examples
makePar <- function(output,chrom=output,pedfile=output,mapgen=output,founder=output, 
                    outname=output,ploidy=4,mapfun="HALDANE", miss="NA",pairing=1){
  #Writes the file
  write(file=paste0(output,".par"),
        c(paste0("PLOIDY = ",ploidy),
          paste0("MAPFUNCTION = ",mapfun),
          paste0("MISSING = ",miss),
          paste0("CHROMFILE = ",chrom,".chrom"),
          paste0("PEDFILE = ",pedfile,".ped"),
          paste0("MAPFILE = ",mapgen,".map"),
          paste0("FOUNDERFILE = ",founder,".gen"),
          paste0("OUTPUT = ",outname),
          paste0("NATURALPAIRING = ",pairing)))
}

#' Make NAM pedigree file
#' @description Generates a NAM pedigree file. Names of the offspring are coded as 
#' O(ffspring)number_crossnumber. 
#'
#' @param parents character vector of parent names
#' @param offspring number of offpsring per NAM 
#' @param output character string for pedigree adress
#' @param prefix prefix to add to offspring names
#'
#' @return A .ped file of NAM structure for PedigreeSim
#' @export
#'
#' @examples
makeNAMPed<- function(
  parents, 
  offspring, 
  output="Pedigree", 
  prefix=NULL
){
  #Create a header
  write(file=paste0(output,".ped"),c("Name Parent1 Parent2"))
  
  #Create the rows regarding parents
  for(i in 1:length(parents)){
    write(file=paste0(output,".ped"),paste0(parents[i]," NA NA"),append=T)}
  
  #Create the rows regarding the offspring. Each offspring will be a cross of P1 with the rest
  for(i in 2:length(parents)){
    for(j in 1:offspring){
      write(file=paste0(output,".ped"),paste0(prefix,"O",sprintf("%02.f",j-1),"_",i-1," ",parents[1]," ",parents[i]),append=T)}
  }
}

#' Make Random Pedigree File
#' @description Generates a pedigree file of a random mixing population. A set of founders 
#' are created labelled as G0# (generation 0), and they are randomly paired, without selfing,
#' to produce the following generations.
#'
#' @param output name of the file to be created
#' @param nfounder number of founder individuals
#' @param popsize population size of each generation
#' @param ngener nubmer of generations to be made
#'
#' @return A PedigreeSim pedigree file following a random mating process.
#' @export
#'
#' @examples
makeRanPed<- function(
  output,
  nfounder=10, 
  popsize=100, 
  ngener=10 
){
  #lists to store the different generations
  Name <- list(); Parent1 <- list(); Parent2 <- list()
  
  #Create names for as many founders as specified (P1, P2...)
  Name[[1]] <- paste0("G0",1:nfounder)
  #create empty lists with slots for as many founders as specified
  Parent1[[1]] <- rep(NA, nfounder)
  Parent2[[1]] <- Parent1[[1]]
  
  #loop that generates offsprings for each generation, in Name the individual names are stored. 
  #Name[1] is a vector with parents, Name[2] is a vector with first generation (G1), etc. 
  #Parent1 and 2 are sampled from previous generation, with replacement
  for (g in 2:(ngener+1)) {
    Name[[g]] <- paste0("G", g-1, "_", sprintf("%03d",1:popsize))
    Parent1[[g]] <- sample(Name[[g-1]], popsize, replace=TRUE)
    Parent2[[g]] <- sample(Name[[g-1]], popsize, replace=TRUE)
    #to prevent selfings, if Parent1 and parent2 are the same, parent 2 is resampled.
    repeat {
      self <- Parent1[[g]] == Parent2[[g]] #self=vector with T or F if they are the same.
      if (sum(self) == 0) break() #if all are F, sum(self)==0
      Parent2[[g]][self] <- sample(Name[[g-1]], sum(self), replace=TRUE) #only those elements with self=true are resampled
    }
  }
  #a dataframe is made that contains the pedigree
  ped <- data.frame(Name=unlist(Name), Parent1=unlist(Parent1), Parent2=unlist(Parent2)) 
  #the df is written into the output file
  write.table(ped,paste0(output,".ped"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
}

#' Make a PedigreeSim chromosome file
#'
#' @param output name of the file to be created
#' @param chms vector of chromosome names
#' @param cM numeric vector of chromosome lengths, in centimorgans
#' @param centros vector of centromere positions
#' @param pref preferential pairing value (all chromosomes will have the same value)
#' @param quad quadrivalent fraction (all chromosomes will have the same value)
#'
#' @return A PedigreeSim chromosome file.
#' @export
#'
#' @examples
makeChrom <- function(
  output, #name or adress of the file to be created
  chms, #vector of chromosome names
  cM, #vector of chromosome length
  centros, #vector of centromere positions
  pref, #preferential pairing (all chr the same value)
  quad #quadrivalent fraction (all chromosomes the same)
){
  
  #check that all vectors have the same length
  if(!length(chms)&&length(cM)&&length(centros)==(length(chms)))
  {stop("Unequal input vector lengths")}
  
  ## Create header:
  write(file=paste0(output,".chrom"),
        c("chromosome length centromere prefPairing quadrivalents"))
  
  #print each row with each vector
  for(i in 1:length(chms)){
    write(file=file.path(folder,paste0(output,".chrom")),
          paste(chms[i],cM[i],centros[i],pref,quad),append=TRUE)
  }
}

#' Make PedigreeSim map file
#'
#' @param output name of the file to be created
#' @param marker vector with marker names
#' @param chr vector with chromosome of each marker
#' @param pos vector with marker positions (in each chromosome)
#' @param head logical value indicating whether the header "marker chromosome position" 
#' should be written or not. Useful for writing maps by chromosome, instead of all at once.
#' @param body logical value indicating whether the body (each marker) should be written or not. 
#' Useful for writing maps by chromosome, instead of all at once.
#'
#' @return A PedigreeSim genetic map file
#' @export
#'
#' @examples
makeMap <- function(
  output, #name of file to be written
  marker, #vector with marker names
  chr, #vector with the chr of each marker
  pos, #vector with marker positions
  head=T, #write the head ("marker chromosome position") or not
  body=T #logical value for whether to write the body or not
){
  
  #Create header. depends on head.
  if(head){write(file=paste0(output,".map"),
                 c("marker chromosome position"))}
  
  #writes rows with marker name, chr to which they belong and position. Depends on body.
  if(body){
    for(i in 1:length(marker)){
      write(file=paste0(output,".map"),
            paste(marker[i],chr[i],pos[i]),append=TRUE)
    }
  }
}

#' Make a PedigreeSim .gen file
#' @description Formats a data.frame into a .gen file as required by PedigreeSim. The data.frame passed
#' needs to have a first column named "markers" which contains marker names.
#' 
#' @param output name of the file to be created
#' @param data data frame where each row is a marker (first column must be named "markers" and contain
#' marker names) and each column is a parental homologue, #' labelled as Parent_1, 
#' Parent_2 ... Parent_ploidy. Each column will contain an allele state (numbers,
#' letters, whatever code).
#'
#' @return A PedigreeSim gen file
#' @export
#'
#' @examples
makeGen<- function(
  output, #output file name
  data #data frame. Each row is a marker. First column is named "markers" and has marker names. Rest of columns contain the alleles of each parental homologue, labelled such as Parent_1, Parent_2 etc..
){
  
  #check right object file cause otherwise write.table will not work
  if(!is.data.frame(data)){print("data file is not dataframe"); break}
  
  #needed for PedigreeSim, checks that first column is called markers (important header in the .gen file)
  names<-colnames(data)
  if(names[1]!="marker"){print("First column of data must be named 'markers', check column names"); break}
  
  write.table(file=paste0(output,".gen"),data,quote=F,row.names=F)
}

#' Random PedigreeSim chromosome file
#' @description For testing purposes. Will generate a random chromosome file with
#' specific chromosome number, random lengths of chromosomes, optional random centromere positions,
#' with specified preferential and quadrivalent probabilities.
#'
#' @param output name of the file to be created
#' @param chms number of chromosomes to be simulated
#' @param CM average chromosome length
#' @param sdCM standard deviation to be applied to chromosome length
#' @param centr average or fixed centromere position (2 is 1/2, 3 is 1/3 and so on)
#' @param rancentr logical value indicating whether centromere position should be 
#' random, or fixed.
#' @param pref probability of preferential pairing.
#' @param quad probability of quadrivalent pairing.
#'
#' @return A PedigreeSim chromosome file with random chromosome structure.
#' @export
#'
#' @examples
ranChrom <- function(
  output, #name of the file
  chms=1, #number of chromosomes
  CM=100, #average chr length
  sdCM=CM/2, #st dev of chr length
  centr=2, #average centromere position (2=half, 3= third and so on) or centromere position (in cM) if rancentr=F
  rancentr=T, #(T/F) whether centr are random
  pref=0, #probability of preferential pairing,
  quad=0 #probability of quadrivalents
){
  
  
  #each chromosome has a letter name. If there are more than 26 chromosomes, double letter names will be used (AA, AB...)
  lets<-c()
  for(i in 1:length(LETTERS)){
    lets<-c(lets,paste0(LETTERS[i],LETTERS))
  }
  chms<-seq(1:chms)
  if(length(chms)>26){chms<-lets[chms]}else{ chms<-LETTERS[chms]}
  
  
  #make vector with chr length. If random, with mean CM and sd=sdCM, otherwise just all same length
  repeat{
    chrlength<-c()
    chrlength<-c(rnorm(length(chms),mean=CM,sd=sdCM))
    if(all(chrlength>0)){break}
  }
  #order CM according to length
  chrlength<-chrlength[order(chrlength,decreasing = T)]
  
  #Make random centromeres vector. Centros are always on the top half of the chromosome length,
  #and are prevented to be negative. If random, mean=chr length/centr sd=chr length/6
  centros<-c()
  if(rancentr){
    for(i in 1:length(chms)){
      repeat{
        centros[i]<-rnorm(1,mean=chrlength[i]/centr,sd=chrlength[i]/6)
        if(centros[i]>0){break}
      }
      if((centros[i]-chrlength[i]/2)>0){centros[i]<-chrlength[i]-centros[i]}
    }
  }else{centros<-rep(centr,length(chms))}
  
  #make the file
  makeChrom(output=output,chms=chms,cM=chrlength,centros=centros,pref=pref,quad=quad)
}

#' Random PedigreeSim map file
#' @description For testing purposes. A random map file is generated, based on a chromosome file.
#' Markers are randomly position but at a specified density along the chromosomes
#'
#' @param chrom .chrom file adress of the file to base the genetic map on.
#' @param output name of the output file. Defaults to the same name as the chromosome file, 
#' with different termination.
#' @param density marker density to simulate, in markers/cM
#'
#' @return A random genetic map file, based on a chromosome file.
#' @export
#'
#' @examples
ranMap<-function(
  chrom, #chrom file name 
  output=chrom, #output name
  density=1 #marker density in markers/cM
){
  #read chrom file to get the chromosome names and the chromosome lengths
  table<-read.delim(chrom,sep=" ")
  
  chromname<-as.vector(table[,1])
  chromlength<-as.vector(table[,2])
  
  #a file is made, but only the header is written
  makeMap(output,body=F)
  
  #for each chromosome
  for(i in 1:length(chromname)){
    #we generate an interval vector containing interval division points.
    interval<-seq(from=0,to=chromlength[i],by=1/density)
    marker<-c();markchrom<-c();markname<-c()
    
    #for each interval, a single marker is positioned randomly (runif), a vector with the chromosome name is made
    #and a marker name is generated such as A001, first marker of chromosome A.
    for(j in 1:(length(interval)-1)){
      marker[j]<-runif(1,interval[j],interval[j+1])
      markchrom[j]<-chromname[i]
      markname[j]<-paste0(chromname[i],sprintf("%03d",j))
      
    }
    
    #for each chromosome, the markers are written in the marker file (without header), and the loop continues.
    makeMap(output,marker=markname,chr=markchrom,pos=marker,head=F)
  }
}

#' Make a random PedigreeSim .gen file
#' @description For testing purposes. Based on a genetic map file, a number of parental
#' genotypes is simulated such that the genotypes originate from a flat distribution. Ploidy
#' of the parents can be specified.
#' 
#' @param map genetic map to base the simulation on (with extension)
#' @param parents number of parents to be simulated
#' @param output name of the output file. Defaults to same name as map, with different extension
#' @param ploidy ploidy number of the parents
#'
#' @return A random PedigreeSim .gen file for a specified number of parents.
#' @export
#'
#' @examples
ranGen<- function(
  map, #name of the map file to be used (with extension)
  parents, #number of parents
  output=map, #name of the output
  ploidy=4
){
  #read map file to obtain marker names. 
  markers<-read.delim(map,sep=" ")
  
  #create data frame where the genotypes will be stored. First column contains marker names.
  gen<-data.frame("marker"=markers[,1])
  #the rest of columns are created, P01_1, P01_2, etc.
  parname<-c()
  for(i in parents){
    parname<-c(parname,paste0(i,"_",seq(1:ploidy)))
  }
  gen[,parname]<-NA
  
  #to obtain a vector with all homologue names.
  homol<-names(gen)
  homol<-homol[2:length(homol)]
  
  #for each marker, in each homologue a sampling procedure decides between A or B
  for(i in 1:length(markers[,1])){
    for(j in 1:length(homol)){
      gen[i,homol[j]]<-sample(c("A","B"),1)
    }
  }
  
  #write it in a file using makeGen function. 
  makeGen(output,gen)
}

#' Run PedigreeSim from R studio
#' @description Run PedigreeSim from R studio by specifying a parameter file and the adress of
#' PedigreeSim. 
#' 
#' @param parfile name of the parameter file.
#' @param path name of the PedigreeSim program.
#'
#' @return
#' @export
#'
#' @examples
run.PedigreeSim <- function(
  parfile, #name of the parfile to be used
  path="PedigreeSIM/PedigreeSim.jar" #specify location of PedigreeSim file
) {
  
  #send the command to the cmd
  system2(command = "java",
          args = c("-jar",
                   path,
                   parfile))
}

#' Parental sampling
#' @description Generates a chromosome file based on sampling from a set of files containing parental
#' chromosomes (.dat files of PedigreeSim). To specify number of ancestrals, one can chose to 
#' specify the number of ancestral groups, from how many different files should the parents be sampled.
#' By default, all files given will be used.
#'
#' @param files table files containing parental chromosomes, labelled as Parent_1, 
#' Parent_2 ... Parent_Ploidy. Each file will constitute an "ancestral group"
#' @param ancgroups number of ancestral groups to use in the sampling. By default, all given ancestral groups
#' will be used. 
#' @param parents number of parents to sample.
#' @param ploidy ploidy of the parents.
#'
#' @return
#' @export
#'
#' @examples
sampleChrom<-function(#Give a vector of file names, containing parental chromosomes. Returns a matrix
  #It generates only the most balanced situations (10 parents from 3 pops will give 3+3+4)
  files,
  ancgroups=length(files), 
  parents=10,
  ploidy=4){
  
  repeat{
    chosen<-sample(files,ancgroups)
    if(length(unique(chosen))==ancgroups){break}
  }
  chosen<-sort(c(rep(chosen,parents%/%ancgroups),head(chosen,parents%%ancgroups)))
  
  #candidate parents of the chosen ancestral pools
  candparents<-sapply(unique(chosen),function(file){
    names<-as.matrix(data.table::fread(file,nrows=1,header=F))[-1]
    names<-unique(substr(names,1,nchar(names)-2))
  })
  
  #choose parents of each ancestral pool (without repeating)
  repeat{
    chosenparents<-apply(candparents[,chosen],2,sample,1)
    if(length(unique(chosenparents))==parents){break}
  }
  
  #Now we use the chosenparents vector to read the corresponding chromosome columns and paste them into one file
  result<-lapply(1:length(chosenparents),function(i){
    chroms<-paste0(chosenparents[i],"_",1:ploidy)
    data.table::fread(chosen[i],select=chroms,header=T)
  })
  result<-do.call(cbind,result)
  markers<-data.table::fread(chosen[1],select="marker",header=T)
  return(cbind(markers,result))
}


##### Mixed models ######

#' Calculation of realized distance matrix (K)
#' @description Using dosage scores, a distance matrix is calculated
#' such that the average distance of an individual with itself is 1, and 
#' the average with an unrelated individual is 0. Based on GWASpoly package.
#' @param matrix Numeric matrix. Individuals on rows and markers in columns
#'
#' @return A numeric matrix nxn where n is the number of rows.
#' @export
#'
#' @examples
calc.K<-function(
  matrix 
){
  #Substract mean dosage for each marker
  M<-apply(matrix,2,function(x) x-mean(x)) 
  K<-M%*%t(M) #Calculate distance matrix
  K<-K/mean(diag(K)) #average the center
  colnames(K)<-rownames(matrix)
  rownames(K)<-rownames(matrix)
  return(K)
}

#' Calculates Hat inverse matrix
#' @description Ridge-regression-based method for calculating Hat inverse matrix.
#' Extracted from GWASpoly & rrBLUP package. In EMMAX/P3D approach to mixed-models 
#' for QTL analysis the variance components are estimated previously. Calculating
#' the Hat inverse matrix is estimating the variance components. The outcome
#' of this function can be recycled on a marker-per-marker analysis,
#' speeding up considerably computation time.
#' 
#' @param y Response vector or matrix. Each column is taken as a different response.
#' @param X Fixed effects design matrix
#' @param Z Random effects design matrix
#' @param K Variance-covariance matrix (Genetic distance matrix). If specified, must
#' be square
#' @param bounds Bounds of the Ridge Regression parameter
#' @param method "REML" or "ML". If "ML", Hinv is calculated using maximum likelihood
#' equations. If REML, Restricted Maximum Likelihood.
#'
#' @return If y is a vector, returns an Hinv matrix. If y is a matrix, returns a list
#' where each element corresponds to the Hinv of each column of y.
#' @export
#'
#' @examples
calc.Hinv<-function(
  y, 
  X, 
  Z, 
  K=NULL, 
  bounds=c(1e-09,1e+09), 
  method="REML" 
){
  if(!dim(K)[1]==dim(K)[2]){stop("K must be square matrix")}
  
  if(is.matrix(y)|is.data.frame(y)){ #to align the dimensions of y and X
    Xmat<-comp.vec(dim(y),dim(X))
  }else if(is.vector(y)){
    Xmat<-comp.vec(length(y),dim(X))
  }
  if(!any(Xmat)){stop("y and X dimensions do not match")}
  if(!Xmat[1]){X<-t(X)}
  
  if(is.matrix(y)|is.data.frame(y)){ #to align the dimensions of y and Z
    Zmat<-comp.vec(dim(y),dim(Z))
  }else if(is.vector(y)){
    Zmat<-comp.vec(length(y),dim(Z))
  }
  if(!any(Zmat)){stop("y and Z dimensions do not match")}
  if(!Zmat[1]){Z<-t(Z)}
  
  if(is.matrix(y)|is.data.frame(y)){ #to align the dimensions of y and K
    Kmat<-comp.vec(dim(y),dim(K))
  }else if(is.vector(y)){
    Kmat<-comp.vec(length(y),dim(K))
  }
  if(!any(Kmat)){stop("y and K dimensions do not match")}
  
  if(!is.null(dim(y))){n<-nrow(y)}else{n<-length(y)}
  p<-ncol(X)
  m<-ncol(Z)
  
  Xinv<-solve(crossprod(X))
  S<-diag(n)-tcrossprod(X%*%Xinv,X)
  
  #In case n<(m+p) the following part is done with eigenvalue decomposition. Otherwise Cholesky decomposition is usable. Because of our dataset, only the first method is implemented. Both require positive semidefinite matrices (all eigenvalues are >=0)
  offset<-sqrt(n)
  if(is.null(K)){
    Hb<-tcrossprod(Z)+offset*diag(n)
  }else{
    Hb<-tcrossprod(Z%*%K,Z)+offset*diag(n)
  }
  
  Hb.system <- eigen(Hb, symmetric = TRUE)
  phi <- Hb.system$values - offset
  if (min(phi) < -1e-06) {
    stop("K not positive semi-definite.")
  }
  U <- Hb.system$vectors
  SHbS <- S %*% Hb %*% S
  SHbS.system <- eigen(SHbS, symmetric = TRUE)
  theta <- SHbS.system$values[1:(n - p)] - offset
  Q <- SHbS.system$vectors[, 1:(n - p)]
  
  Hinv<-lapply(1:(length(y)/n),function(i){#for each number of phenotypes we will recycle Q, allowing us to obtain the Hinv much faster
    suby<-y[((i-1)*n+1):(i*n)]
    
    omega.sq<-crossprod(Q,suby)^2
    
    if (method == "ML") {
      f.ML <- function(lambda, n, theta, omega.sq, phi) {
        n * log(sum(omega.sq/(theta + lambda))) + sum(log(phi+lambda))
      }
      soln <- optimize(f.ML, interval = bounds, n, theta, omega.sq, phi)
      lambda.opt <- soln$minimum
      df <- n
    }
    else {
      f.REML <- function(lambda, n.p, theta, omega.sq) {
        n.p * log(sum(omega.sq/(theta + lambda))) + sum(log(theta + lambda))
      }
      soln <- optimize(f.REML, interval = bounds, n - p, theta, omega.sq)
      lambda.opt <- soln$minimum
      df <- n - p
    }
    Hinv <- U %*% (t(U)/(phi + lambda.opt))
  })
  if(length(Hinv)==1){Hinv<-matrix(unlist(Hinv),ncol=ncol(Hinv[[1]]))}
  
  return(Hinv)
}

#' Remove non-segregants from a dosage matrix
#' @description Markers that are not segregating in the population 
#' are removed from the dosage matrix.
#' @param crossfile ADRESS of a Pedsim doseallele file.
#'
#' @return dosage matrix with segregants removed
#' @export
#'
#' @examples
bi.geno<-function(
  crossfile
){
  result<-as.matrix(data.table::fread(crossfile))
  markers<-result[,1]; inds<-colnames(result[,-1])
  result<-result[,-1] #take out marker column
  segregants<-sapply(1:nrow(result),function(x) segregation<-length(unique(result[x,-1])) )
  result<-matrix(as.numeric(result[which(segregants!=1),]),ncol=ncol(result))
  colnames(result)<-inds; rownames(result)<-markers[which(segregants!=1)]
  return(result)
  
}

#' Calculate a Q design matrix based on a vector
#'
#' @param pop vector where each element is a population identifiyer
#' @param names optionally, a set of names for each element
#'
#' @return a matrix of ncol=unique(pop) identifying each individual belonging to 
#' one population.
#' @export
#'
#' @examples
Q.mat<-function(
  pop, #vector identifying a population
  names=NULL #optionally, a vector defining the names
){
  match<-1*sapply(unique(pop),function(x) x==pop)
  if(!is.null(names)){rownames(match)<-names}
  return(match[,-1]) #we take out the first column, to avoid singularity
}

#' Vector comparison function
#' @description Element per element comparison between two vectors.
#' @param vec1 vector of values
#' @param vec2 vector of values
#'
#' @return logical matrix where cols are vec1 and rows are vec 2, the value 
#' indicating if row and col are equal.
#' @export
#'
#' @examples
comp.vec<-function(vec1,vec2){
  mat<-sapply(vec1,function(x) x==vec2)
  colnames(mat)<-vec1
  rownames(mat)<-vec2
  return(mat)}

#' Mixed Model Solver
#' @description Mixed model solver based on the rrBLUP package mixed.solve function
#'
#' @param y Numeric vector of response variable
#' @param X Fixed effect design matrix. Must include an intercept.
#' @param Z Random effects design matrix.
#' @param K Variance-covariance structure. If not specified, identity matrix is used
#' @param Hinv Hat inverse matrix, will be calculated if not provided
#' @param random wether to return random effects also (slows down significantly)
#' @param no.test Number of parameters in the X matrix that will not be taken into 
#' account in the pvalue calculation. In practice, this means that the first no.test parameters
#' will not be used for the pvalue calculation. So, if 3 cofactors are added to a model, no.test
#' should be 4 (to ignore intercept and the 3 cofactors, place on the 3 first columns after intercept)
#'
#' @return a list with: beta (estimates of parameters), Fstat (estimate of F value),
#' pval (estimate of p-value), se (standard error of the model)
#' @export
#'
#' @examples
mm.solve<-function(y,X,Z,K,Hinv=NULL,random=F,no.test=1){
  if(!is.null(dim(y))){n<-nrow(y)}else{n<-length(y)}
  p<-ncol(X)-no.test #number of genetic parameters. Only works if there is 1 intercept and the rest are genetic parameters
  df<-n-ncol(X)
  
  if(no.test==0){
    X2<-Z%*%X
  }else{ ##### WARNING DOES THIS WORK WITH MORE COFACTORS??
    X2<-cbind(X[,1:no.test],Z%*%X[,-1:-no.test])#do not multiply cofactors
  }
  
  if(is.null(Hinv)){Hinv<-calc.Hinv(y=y,X=X,Z=Z,K=K)}
  W<-crossprod(X2, Hinv %*% X2)
  Winv<-solve(W)
  beta <- Winv %*% crossprod(X2, Hinv %*% y) #estimates
  resid <- y - X2 %*% beta #residual
  s2 <- as.double(crossprod(resid, Hinv %*% resid))/df #variance
  
  Q<-s2*Winv[(1+no.test):ncol(X),(1+no.test):ncol(X)]#We only take the Winv col and row of the parameters we want to test (not intercept, so ncol(X)-1)
  Tt <- solve(Q) #Estimates/se
  Fstat <- crossprod(beta[(1+no.test):ncol(X)], Tt %*% beta[(1+no.test):ncol(X)])/(p) #F-statistic of the genetic model?
  pval<-pf(Fstat,p,df,lower.tail=F) #pvalue of genetic model?
  
  if(random){
    KZt <- tcrossprod(K, Z)
    u <- crossprod(KZt %*% Hinv,(y - X %*% beta))#random effects
    return(list(beta=beta,random=u,se=s2,Fstat=Fstat,pval=pval))
  }
  
  return(list(beta=beta,Fstat=Fstat,pval=pval,se=sqrt(s2)))
}

#' Calculates design matrix of X
#'
#' @param genotypes vector of dosages per individual, or matrix of genotypes per 
#' chromosome
#' @param ploidy 
#'
#' @return if genotypes is vector, a matrix will be returned with first column having an 
#' intercept (all 1's), second column having the dosages. If genotype is a matrix,
#' a design matrix is return with ncol=unique()
#' @export
#'
#' @examples
dosage.X<-function(genotypes,ploidy=4){
  if(all(unique(genotypes)%in%0:ploidy)){
    result<-c(rep(1,length(genotypes)))
    result<-matrix(c(result,genotypes),ncol=2)
    return(result)
  }
  
  #we obtain the different alleles present
  unals<-unique(genotypes)
  #we obtain a design matrix indicating the allele of each chromosome
  match<-sapply(unals,function(x) x==genotypes)
  #we count the number of each allele for each individual
  alcount<-sapply(1:(length(genotypes)/ploidy),
                  function(x) colSums(match[1:ploidy+(x-1)*ploidy,]))
  
  #add an intercept
  result<-t(rbind(rep(1,ncol(alcount)),alcount)) 
  inds<-unique(substr(names(genotypes),1,nchar(names(genotypes))-2))
  rownames(result)<-inds
  colnames(result)<-c("Intercept",unals)

  return(result)
}

#' Sample markers every X centimorgans
#'
#' @param genotypes named vector or matrix where names or rownames are marker names
#' @param map table containing at least a "chromosome" and a "position" column
#' @param cM Size of bin to sample a marker from
#'
#' @return A table as genotypes, with as many rows as can be obtained by sampling
#' markers every X centimorgans
#' @export
#'
#' @examples
sample.cM<-function(
  genotypes,
  map, 
  cM=1
){
  maxlen<-sapply(unique(map$chromosome),function(x) max(map$position[which(map$chromosome==x)]))
  cumlen<-c(0,cumsum(maxlen))
  pos<-map$position+cumlen[map$chromosome]
  
  selpos<-seq(0,max(pos),by=cM) #the cM at which to sample
  
  markers<-sapply(selpos,function(k){
    dif<-abs(k-pos)#what is the closest position to this marker
    a<-which(dif==min(dif))#the index of the closest marker
    if(length(a)>1){#in case there is more than one marker at that position, we select randomly 1
      a<-sample(a,1)
    }
    return(a)
  })
  
  if(is.vector(genotypes)){
    result<-genotypes[markers]
  }else{
    result<-genotypes[markers,]
  }
  
  return(result)
}

#' QTL mapping of a matrix of phenotypes
#'
#' @param phenotypes A matrix of phenotypes, 
#' rows are individuals and columns are different phenotypes. 
#' Not updated for vector phenotypes yet. Must follow same
#' order of individuals as genotypes.
#' @param genotypes A matrix of genotypes, rows are markers and
#' columns may be either (1) individual dosages, with column
#' names coinciding with phenotype individual names or (2) chromosome alleles of each
#' individual. Must follow same order of individuals as phenotypes.
#' @param K NULL, T or distance matrix. If NULL, no relatedness
#' matrix will be used (i.e, a linear model will be applied). If T
#' a K distance matrix will be calculated. A distance matrix may also
#' be directly specified.
#' @param Q NULL, T or vector identifying populations. If NULL, no Q will
#' be included in the model (i.e, a model without Q correction). If T, a pco
#' decomposition will be used to estimate population differentiation. If a
#' vector specifying population of each individual is passed, 
#' it will be used to construct a Q matrix. Vector may contain numerical 
#' or character.
#' @param dosage A dosage matrix, markers on rows, individuals on columns. 
#' If not specified, it will try to read it from "genotypes". If it fails it
#' will not be able to calculate K distance matrix.
#' @param map A table with a genetic map containing at least a "chromosome"
#' and a "position" columns, specifying chromosome number and cM position. 
#' Used to sample marker dosages every 1 cM
#' AN ADRESS of a Pedsim map file: at least with a "chromosome"
#' and "position" columns, specifying chromosome number and cM position 
#' of the marker at that chromosome
#' @param cM Numeric. K distance matrix will be calculated using markers every
#' cM centimorgans. Defaults to 1.
#' @param k Numeric. Number of axis to be used in pco-based Q estimation. Defaults to 2.
#' @param no_cores Numeric. Number of cores to be used in parallel computing
#' of the p-values. Defaults to number of cores -1.
#' 
#' @return a pvalue matrix containing the pvalue of each marker with each phenotype passed.
map.QTL<-function(
  phenotypes,
  genotypes, #genotype matrix
  K=NULL, #distance matrix
  Q=NULL, #population effect matrix
  Z=NULL,
  dosage=NULL, #dosage matrix
  map, #genetic map table
  cM=1, #
  k=2, #number of axis used for pco decomposition
  no_cores=parallel::detectCores()-1
){
  
  if(is.null(K)){
    linear<-T
  }else{
    linear<-F
  }
  
  markers<-rownames(genotypes)
  
  #DEFINITION OF K
  if(all(K==T)){ #if K is T, we calculate K distance matrix
    
    if(is.null(dosage)){
      if(all(rownames(phenotypes)==colnames(genotypes))){
        dosage<-genotypes
      }else{
        stop("If genotypes are not marker dosages, dosages must be specified") 
      }
    }
    
    #first we sample homogeneously markers along the genome
    K<-sample.cM(dosage,map,cM = cM)
    #then we calculate the distance
    K<-calc.K(t(K))
    
  }else if(!all(ncol(K)==dim(X)[1],nrow(K)==dim(X)[1])){
    #in case K is already defined, we check dimensions
    stop("K matrix does not have adequate dimensions")
  }
  
  #DEFINITION OF Q
  if(!is.null(Q)){
    if(Q==T){
      if(is.null(K)){
        if(all(rownames(phenotypes)==colnames(genotypes))){
          dosage<-genotypes
        }else if(is.null(dosage)){
          stop("If genotypes are not marker dosage, dosage must be specified")
        }
        
        #first we sample homogeneously markers along the genome
        K<-sample.cM(dosage,map,cM = cM)
        #then we calculate the distance
        K<-calc.K(t(K))
      }
      Q <- cmdscale(1-K,k=k, eig = F, add = FALSE, x.ret = FALSE)
      
    }else if(is.vector(Q)){#if Q is a vector identifying populations, a Q matrix is designed
      Q<-Q.mat(Q)
    }else{ stop("Wrong Q specification")}
  }
  
  
  if(linear){ 
    #First set up cluster
    
    cluster<-parallel::makeCluster(no_cores)
    
    if(is.null(Q)){
      export<-c("phenotypes","genotypes","dosage.X")
      print("Linear model will be used")
    }else{
      export<-c("phenotypes","genotypes","dosage.X","Q")
      print("Linear model with Q correction will be used")
    }
    
    parallel::clusterExport(cl=cluster,export,
                            envir=environment()) #needed again?
    
    #For each row of genotypes, calculate lm for all phenotypes
    pval<-parallel::parSapply(cl=cluster,1:nrow(genotypes),function(k){
      X<-dosage.X(genotypes[k,])
      nparX<-ncol(X)
      
      
      #if Q is NULL, X will not change
      X<-cbind(X,Q)
      
      test1<-lm(phenotypes~X)
      SSR1<-colSums((test1$residuals)^2) #Sum of Squares Residuals Full
      
      test3<-lm(phenotypes~X[,2:nparX]) #Sum of Squares of Markers
      SSR3<-colSums((test3$residuals)^2)
      
      SST<-apply(phenotypes,2,function(x){
        sum((x-mean(x))^2)
      })
      SSM<-SST-(SSR3-SSR1)-SSR1 #Calculate the amount of variation explained by the markers
      #in a model with a Q structure
      
      df2<-apply(test1$coefficients,2,function(x) nparX-sum(is.na(x)))
      
      Ftest<-(SSM/df2)/(SSR1/test1$df.residual)
      return(pf(Ftest,df2,test1$df.residual,lower.tail = F))
      
    }) 
    parallel::stopCluster(cluster)
    pval<-t(pval)
    
  }else{ #ergo, K must be defined, we must use mixed models
    if(is.null(Z)){
      Z<-diag(nrow(phenotypes)) 
    }
    
    Hinv<-calc.Hinv(phenotypes,
                    #X=vecor of 1s, to apply P3D/EMMAX algorithm
                    X=matrix(rep(1,nrow(phenotypes))),
                    Z,K) 
    
    cl<-parallel::makeCluster(no_cores)
    if(is.null(Q)){
      export<-c("phenotypes","Z","K","Hinv","genotypes","mm.solve","dosage.X","no.test")
      no.test<-1
      print("Mixed model will be used")
    }else{
      export<-c("phenotypes","Z","K","Hinv","genotypes","mm.solve","dosage.X","Q","no.test")
      no.test<-1+ncol(Q)
      print("Mixed model with Q correction will be used")
    }
    
    parallel::clusterExport(cl,export,
                            envir=environment()) #I dunno why, but without this it doesnt work
    
    pval<-parallel::parSapply(cl,1:ncol(phenotypes),function(w){#parallely, over each phenotype
      sapply(1:nrow(genotypes),FUN=function(k){#calculate the pvalue for each marker
        X<-dosage.X(genotypes[k,]) #obtain the X matrix for that marker
        if(ncol(X)>2){X<-X[,-2]}#the -2 is to take out one of the columns of X, to prevent singularity when solving
        X<-cbind(Q,X)
        
        mm.solve(phenotypes[,w],X,Z,K,Hinv[[w]],no.test = no.test)$pval
      })
    })
    
    parallel::stopCluster(cl)
  }
  
  colnames(pval)<-colnames(phenotypes)
  rownames(pval)<-markers

  
  return(pval)
}


##### Phenotyping ######
linkage<-function(#for a series of populations, calculate the linkage of all markers with one marker.
  #Dependant on link_NAM
  marker, #marker to analyze the linkage with
  crossfile, #crossfile
  totallele="PedigreeSIM/Parents/Total_pop.txt",
  no_cores=7 #careful with changing this
){
  cross_anc<-link_NAM(crossfile=crossfile,parental=T,totallele=totallele) #ancestral alleles of the parents
  cross_pa<-data.table::fread(crossfile,header=T) #cross expressed in parental alleles
  cross_anc<-as.matrix(cross_anc)
  markers<-unlist(cross_pa[,1,with=F])
  cross_pa<-as.matrix(cross_pa[,-1,with=F]) #take out the marker column
  
  #This is a little geniality that makes this program way faster. It simply treats the cross_pa matrix as an index vector, but each row is added as many columns as cross_anc, 
  #so that items on the first row will translate with items from first row of cross_anc, second with second and so on. Very fast.
  indexvector<-cross_pa+(1:nrow(cross_pa)-1)*40+1
  result<-matrix(t(cross_anc)[indexvector],ncol=ncol(cross_pa))
  colnames(result)<-colnames(cross_pa);rownames(result)<-rownames(cross_pa)
  
  linked<-sapply(marker,function(mark){
    match<-t(t(result)==result[mark,])
    link<-rowSums(match)/ncol(match)
  })
  
  colnames(linked)<-markers[marker] #the ones given by the users
  rownames(linked)<-markers #the complete list, read from cross_pa
  
  return(linked)
}

link_NAM<-function(#function to connect ancestral founder alleles and NAM alleles
  crossfile, #founderallele file of the cross
  marker=NULL, #index of marker or markers to get
  totallele="PedigreeSIM/Parents/Total_pop.txt", #pedigreeSim-type table with chromosomes and individuals, specifying alleles for each parental chromosome
  parental=F, #F:return the whole cross NAM, T:returns only the parental alleles
  parents=10, #number of parents
  ploidy=4
){
  #Obtain chromosome names of parents
  headcross<-unlist(data.table::fread(crossfile,nrow=1,header=F)) #first column is "marker"
  par_names<-headcross[(1+(1:(parents*ploidy)))]#first columns must contain parent names
  markers<-data.table::fread(crossfile,select="marker",header=T)
  
  #find genotypes in datafile with true alleles
  par_founder<-data.table::fread(totallele,select=par_names,header=T)
  par_founder<-as.matrix(par_founder)
  
  if(!parental){
    #read the cross genotypes
    cross_geno<-data.table::fread(crossfile,header=T)
    
    cross_geno<-as.matrix(cross_geno[,-1,with=F]) #turn into matrix for speed
    #use chromosome numbers in cross_geno as indeces to on par_founder, which contains the true alleles ordered so that they coincide with parental chromosomes
    result<-t(sapply(1:nrow(markers),function(x) par_founder[x,cross_geno[x,]+1]))
    rownames(result)<-unlist(markers)
    colnames(result)<-headcross[-1] #change names of matrix, exclude the "marker column"
  }else if(parental){
    rownames(par_founder)<-marker
    return(par_founder)
  }else{stop("Parentals must be either T or F")}
  
  if(!is.null(marker)) result<-result[marker,] #if certain markers are selected, only return those
  
  return(result)
}

multi.effects<-function( #custom function. Requires same number of rows in gen, as ancestral groups present. Returns effects for the alleles at each position asked.
  gen, #genotypes (per chromosome of individual)
  anc_alleles, #matrix defining which alleles belong to which ancestral group
  size=NULL, #vector of effects for each locus (each row in gen)
  n=1 #number of alleles that are not 0
){
  
  if(is.null(size)){size<-(nrow(gen):1)}
  
  anc_mat<-lapply(1:nrow(gen),function(k){
    a<-unique(gen[k,]) #unique alleles present in locus k
    b<-sapply(a,function(a) colSums(a==anc_alleles)) #to which ancestral does each allele correspond
    c<-rowSums(b) #number of alleles presents of each ancestral group
    d<-which(c!=0)[k] #we take the kth ancestral group
    alleles<-a[as.logical(b[d,])] #we obtain the alleles of the kth ancestral group
    effects<-rep(0,length(a))#create a vector of empty effects for all alleles
    if(n>length(alleles)){n<-length(alleles)}
    effects[which(a%in%alleles)]<-c(rep(size[k],n),#n alleles with an effect
                                    rep(0,times=(length(alleles)-n)))#the rest of alleles with no effect
    names(effects)<-a
    return(effects)
  })
  
  return(anc_mat)
}

multi.pheno<-function(#Obtain a phenotype based on the QTLs defined in the matrix gen. Each column should correspond to the concatenated allels of a specific locus
  gen, #genotypes in a tidy manner
  effects, #list of effects. Each element of the list should have as many alleles as there are in each locus
  herit, #total heritability of the character
  partial.herit=NULL, #vector of partial heritabilities for each QTL. APPROXIMATE RESCALING
  mu=50,
  Evar=1,
  nchar=3, #number of characters of each allele
  seed=7){
  
  pattern<-paste0("(.{",nchar,"})")
  als<-apply(gen,1,function(x){
    x<-gsub(pattern, "\\1 ", x)
  })
  #alleles found in the data. Check that there are as many effects as alleles
  unals<-apply(als,2,function(x){
    unique(unlist(strsplit(as.character(x),split=" ")))
  })
  check1<-lapply(1:length(unals),function(i) length(unals[[i]])==length(effects[[i]]))
  if(!all(unlist(check1))){
    stop("Effects and number of alleles do not match in at least 1 locus")
  }
  
  #We obtain a dosage list, where each element is one locus
  dosages<-lapply(1:length(unals),function(x) {
    a<-sapply(unals[[x]],function(y) stringr::str_count(als[,x],y))
    rownames(a)<-rownames(als)
    colnames(a)<-paste(x,colnames(a))
    return(a)
  })
  
  pheno<-lapply(1:length(dosages),function(i){dosages[[i]]%*%effects[[i]]})
  pheno<-do.call(cbind,pheno)#this has a column for the effect of each locus
  
  if(!is.null(partial.herit)){#if partial heritabilities are defined, the genetic effects are rescaled
    varmat<-var(pheno)
    beta<-optim(rep(200,nrow(gen)),lower=0,method="L-BFGS-B",function(beta){
      sum(sapply(1:nrow(gen),function(i){
        (beta[i]^2*varmat[i,i]+sum(beta[-i]*beta[i]*varmat[i,-i])-partial.herit[i])^2
      }))
    })
    
    pheno<-lapply(1:length(dosages),function(i){dosages[[i]]%*%effects[[i]]*beta$par[i]}) #update phenotypes
    pheno<-do.call(cbind,pheno)
    
    effects<-apply(1:length(effects),function(x) effects[[x]]*beta$par[i]) #update effects
  }

  set.seed(seed)
  env<-rnorm(nrow(pheno),mean=mu,sd=Evar)
  Sg<-sum(var(pheno))
  Se<-var(env)
  
  #Now we will rescale the genetic effects
  alpha2<-Se*herit/(Sg*(1-herit))#calculation of alpha squared
  effects2<-lapply(effects,function(x) x*sqrt(alpha2))#rescaling of genetic effects
  pheno2<-lapply(1:length(dosages),function(i){dosages[[i]]%*%effects2[[i]]})
  pheno2<-do.call(cbind,pheno2)#the trait has heritability of 0.7 now
  
  result<-rowSums(pheno2)+env#add up genetic and environmental effects
  names(result)<-colnames(gen)#names of the individuals
  return(result)
}

pheno<-function(#obtain a phenotype table with first rows, allele effects, and phenotype for each individual.
  genotypes, #named vector of genotypes
  effects, #should be equal to cj*alpha
  length=3, #all alleles should be equal size
  Evar=1, #environmental variance to be used
  mu=50, #the average phenotype
  h2=0.3, #the heritability to be simulated
  seed=NULL
){
  
  
  pattern<-paste0("(.{",length,"})")
  als<-sapply(genotypes,function(x){
    x<-gsub(pattern, "\\1 ", x)
  })
  
  #genotypes found in the data. Check that there are as many cjs as genotypes
  unals<-unique(unlist(strsplit(as.character(als),split = " ")))
  if(length(unals)!=length(effects)){stop("Number of alleles and effects do not match")}
  
  #add as many columns as alleles, count the number of occurences
  dosages<-sapply(unals,function(x) stringr::str_count(als,x) )
  rownames(dosages)<-names(als)
  
  #calculate the Dosage Sum of Sqares (SSD), and the total
  dosages<-cbind(dosages,sapply(1:length(unals), function(x){
    avg<-mean(dosages[,unals[x]])
    (effects[x]*(dosages[,unals[x]])-avg)^2
  }))
  colnames(dosages)<-c(unals,paste("SSD",unals))
  
  SST<-rowSums(dosages[,paste("SSD",unals)])^2
  SST<-sum(SST)
  alpha<-sqrt(h2*Evar*(nrow(dosages)-1)/((1-h2)*SST)) #the equation
  
  phenotype<-sapply(1:length(unals), function(x) dosages[,unals[x]]*effects[x]*alpha)
  set.seed(seed)
  phenotype<-rowSums(phenotype)+rnorm(nrow(phenotype),mean=0,sd=Evar)+mu
  
  result<-c(round(effects*alpha,digits=4),phenotype)
  return(result)
}


sample.QTL<-function( #returns a selection of unlinked QTL positions based on a map table
  QTL_number, #the number of QTLs to simulate
  map, #a table with a "marker" a "chromosome" and a "position" (in cM) columns.
  seed=NULL #an option to add a seed
){
  set.seed(seed)
  QTL_chr<-sample(unique(map$chromosome),QTL_number)
  
  set.seed(seed)
  QTLs<-sapply(QTL_chr,function(chr){ #which markers to sample
    a<-sample(map$marker[which(map$chromosome==chr)],1)
    index<-which(a==map$marker)
    return(cbind(map[index,],index=index))
  })
  
  QTLs<-t(QTLs)
  QTLs<-data.frame(marker=unlist(QTLs[,1]),chromosome=unlist(QTLs[,2]),
                   position=unlist(QTLs[,3]),index=unlist(QTLs[,4]))
  return(QTLs)
}

tidy_allele<-function(#gets a matrix with chromosomal alleles and returns a concatenation of the alleles
  #of the same individual. Expects chromosomes to be labelled as name_1, name_2; and will group all alleles
  #from the same "name".
  alleles, #matrix with a chromosome in each column
  nchar=3, #characters of the alleles
  ploidy=4
){
  if(is.vector(alleles)){
    names<-names(alleles)
    alleles<-matrix(alleles,nrow=1)
    colnames(alleles)<-names
    }
  alleles<-t(alleles) #turn into vertical
  format<-paste0("%0",nchar,".f")
  
  #get parent names
  parents<-sapply(rownames(alleles),function(x) substr(x,1,nchar(x)-2))  
  #paste together all alleles belonging to same individual

  compalleles<-sapply(unique(parents),function(x){
    y<-parents==x
    z<-split(sprintf(format,alleles[y,]),rep(1:ncol(alleles),each=ploidy))
    sapply(z,paste0,collapse="")
  })
  
  return(compalleles)
}

#### Result analysis ####
#Some plotting functions, and Liji Correction

LiJi <- function(m,alpha) {
  ## m, matrix of dosages with markers in rows and individuals in columns
  ## alpha, significance level
  
  # step 1: correlation matrix between tests (markers)
  cm <- cor(t(m))   # correlation matrix between markers
  
  # step 2: estimate the effective number of independent tests (Meff)
  eval <- eigen(cm, symmetric=T, only.values=T)   # eigenvalues
  epos <- eval$values[eval$values>=0]  # only positive eigenvalues
  epos2 <- as.numeric(epos>=1) + (epos - floor(epos))  # Li Ji method
  Meff <- sum(epos2)  # number of independent tests
  
  # step 3: adjust alpha
  alpha.adj <- 1-(1-alpha)^(1/Meff)  # similar to (alpha/Meff)
  
  return(alpha.adj)
}

ggd.qqplot<-function( #computes the QQplot for a set of pvalues
  pvector, #vector of pvalues
  main=NULL, #title
  ... #other plot parameters
  ){
  o = -log10(sort(pvector,decreasing=F))#observed
  e = -log10( 1:length(o)/length(o) )#expected
  plot(e,o,pch=19,cex=0.7, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  lines(e,e,col="red")
}

colplotQTL<-function(
  pvals, #vector or matrix of pvalues
  pos, #cM postiion of the markers
  alpha=0.05,
  main=NULL,#title of the plot
  QTLs=NULL,
  mean=T, #whether or not to plot the rowmean
  ...
){
  testnum<-dim(pvals)[2]
  if(is.null(testnum)){testnum<-1}
  if(testnum==1){mean=F}
  #QTLplot
  plot(rep(pos,testnum),-log10(pvals),
       ylab="Significance",xlab="cM",
       cex=0.5,col="darkblue",
       main=main,...)
  #average pvalue line
  if(mean){
    avg<-rowMeans(-log10(pvals))
    lines(pos,avg,col="skyblue")
  }
  
  #the QTL lines
  if(!is.null(QTLs)){
    abline(v=QTLs,lty=3)
  }
  #threshold line
  abline(h=-log10(alpha),col="red")
}

colQQplot<-function(
  pvals,
  main=NULL,
  ...
){
  testnum<-dim(pvals)[1]
  if(is.null(testnum)){testnum<-length(pvals)}
  
  #first plot
  e<- -log10(1:testnum/testnum)
  plot(e,e,type="l",col="red",main=main,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)),ylim=c(0,max(-log10(pvals))),...)
  #points superposition
  apply(pvals,2,function(x){
    o<- -log10(sort(x,decreasing=F))
    e<- -log10(1:length(o)/length(o))
    points(e,o,pch=19,cex=0.5,col=rgb(0.7,0.7,0.7,0.5))
  })
  #make expected line on top
  lines(e,e,col="red",cex=1.2)
}

comp.skyplot<-function(#computes the row maximum of each dataframe in a list and plots it.
  list, #list with elements to be plotted. Each element is expected to be a table with columns having QTL analysis
  x,#values to use for the X axis
  main=NULL,
  xlab="Genomic cM position",
  ylab=expression(Maximum~~-log[10](italic(p))),
  QTLs=NULL, #optional QTL positions to be drawn (as vertical lines)
  legend=names(list), #names to be used in the legend
  index=T, #value indicating column indexes for the tables. If not specified, all.
  h=c(120,240), #hue for the colours
  coltype="divergent", #type of colour palette
  colors,
  legtitle="Models",
  threshold=7e-5, #threshold line
  c=100,
  min=T,
  xindex=T
){
  n<-length(list)
  
  #First, we calculate the ylimit of the main plot
  ylim<-sapply(list,function(x){min(x[xindex,index],na.rm =T)})
  ylim<-min(ylim)
  #Then, we create the colour palette we are going to use
  if(coltype=="divergent"){
    #creates a divergent palette from a central neutral color
    cols<-colorspace::diverge_hcl(n+1,h=h,c=c,power=0.7)
    half<-(n+1)%/%2+(n+1)%%2
    cols<-c(rev(cols[1:(half-1)]),cols[-1:-(half-1)])
    cols<-cols[-c(half)]#to take out the central neutral color
  }else if(coltype=="sequential"){
    #creates a sequential palette based on one hue
    cols<-colorspace::sequential_hcl(n,h=h,c=c)
  }else if(coltype=="rainbow"){
    cols<-colorspace::rainbow_hcl(n,c=c,l=60)
  }else if(coltype=="spec"){
    cols<-colors
  }
  
  
  #We will generate two plots, the first one containing the skyplot, the second the legend.
  #layout splits the device in two
  layout(matrix(c(1,2),ncol=2),
         widths = c(12,2),heights = c(4,1))
  #we change the margins to fit the two plots well together
  par(mar=c(4,4,4,0)+0.1)
  plot(x,x,type="n",ylim=c(0,-log10(ylim)),
       xlab=xlab,ylab=ylab,main=main)
  #we draw the maximum skyplot lines
  if(!min){x<-rep(x,ncol(list[[1]][,index]))}
  for(i in 1:n){
    result<-list[[i]][xindex,index]
    if(min){y<-apply(result,1,min,na.rm=T)}
    if(!min){
      y<-result
      }
    points(x,-log10(y),col=cols[i],cex=0.5,pch=19)
  }
  abline(v=QTLs,lty=3)
  abline(h=-log10(threshold),col="red",lwd=2)
  
  #Then the legend
  x0<-rep(0.5,n)
  x1<-rep(1,n)
  y0<-n:1/n
  y1<-n:1/n
  
  par(mar=c(7,0,7,0)) #set up adequate margins. Why does it look so terrible
  plot(1,1,type="n",axes=F,xlab="",ylab="",
       xlim=c(0,3),ylim=c(0,1))
  segments(x0,y0,x1,y1,col=cols,lwd=5) #make the lines
  text(rep(2.2,6),y0,legend) #add the names
  text(1.5,1.2,legtitle,xpd=T) #add the title "Models"
  
}

comp.QQplot<-function( #compares QQplots of columns of a matrix
  pvals,
  main=NULL,
  coltype="rainbow",
  legtitle="Models",
  h=c(120,240),
  c=100,
  l=70,
  legend=colnames(pvals),
  lim=NULL
){
  n<-dim(pvals)[1]
  e<--log10(1:n/n)
  
  if(is.null(lim)){
    xmax<-max(e)
    ymax<-max(-log10(pvals))
  }else{
    xmax<-lim
    ymax<-lim
  }
  
  layout(matrix(c(1,2),ncol=2),
         widths = c(12,2),heights = c(4,1))
  par(mar=c(4,4,4,0)+0.1)
  plot(e,e,type="l",col="red",main=main,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,xmax),ylim=c(0,ymax))
  
  m<-ncol(pvals)
  
  if(coltype=="divergent"){
    #creates a divergent palette from a central neutral color
    cols<-colorspace::diverge_hcl(m+1,h=h,c=c,power=0.7)
    half<-(m+1)%/%2+(m+1)%%2
    cols<-c(rev(cols[1:(half-1)]),cols[-1:-(half-1)])
    cols<-cols[-c(half)]#to take out the central neutral color
  }else if(coltype=="sequential"){
    #creates a sequential palette based on one hue
    cols<-colorspace::sequential_hcl(m,h=h,c=c)
  }else if(coltype=="rainbow"){
    cols<-colorspace::rainbow_hcl(m,c=c,l=60)
  }else if(coltype=="spec"){
    cols<-colors
  }
  
  sapply(1:m,function(x){
    o<- -log10(sort(pvals[,x],decreasing=F))
    e<- -log10(1:length(o)/length(o))
    points(e,o,pch=19,cex=0.5,col=cols[x])
  })
  
  x0<-rep(0.5,m)
  x1<-rep(1,m)
  y0<-m:1/m
  y1<-m:1/m
  
  par(mar=c(7,0,7,0)) #set up adequate margins. Why does it look so terrible
  plot(1,1,type="n",axes=F,xlab="",ylab="",
       xlim=c(0,3),ylim=c(0,1))
  segments(x0,y0,x1,y1,col=cols,lwd=5) #make the lines
  text(rep(2.2,6),y0,legend) #add the names
  text(1.5,1.2,legtitle,xpd=T) #add the title "Models"
  
}

mosaic.list<-function(
  list,
  x=NULL, #values of x
  h=240,
  ...
){
  
  if(is.null(x)){
    x<-1:length(list)
  }
  
  plot(1,1,type="n",
       xlim=c(0,max(x)),ylim=c(0,1),...)
  
  for(i in 1:length(list)){
    lin<-list[[i]]
    x0<-rep(x[i],length(lin))
    x1<-x0
    y1<-cumsum(lin/sum(lin))
    y0<-c(0,y1[-length(y1)])
    data<-cbind(x0,y0,x1,y1)
    cols<-colorspace::sequential_hcl(length(lin),c=100,h=h)
    segments(data[,1],data[,2],data[,3],data[,4],col=cols)
  }
}

