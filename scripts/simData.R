
#################################################################################
#### Step 1, run netwSim_Datanew.R #############################################
#################################################################################
## This is a bit different and I think more efficient from making random data from netwSim_datanew.R
#> Differences:
start_time <- Sys.time()
test <- TRUE ## this is a test run
saveit <- TRUE ## save the sample?

# Time difference of 1.000327 mins
# setwd("/home/ag/projects/connectivity/")


## packages installment
## used in sum_fx.r for "graphAM" object
if (!require("BiocManager", quietly = TRUE))  install.packages("BiocManager")
BiocManager::install("Rgraphviz", dependencies=TRUE)
BiocManager::install("RBGL")

pkgs <- c("dplyr","tidyverse","ggplot2", "purrr","broom", "reader",
          "deSolve","parallel",
          "lhs","Cairo","ggrastr","tikzDevice","latex2exp","GGally",
          "igraph","igraph","data.tree",
          "Rgraphviz","RBGL",
          "here",
          # "plyr",
          "future","furrr","purrr") ## multicore computation in Windows:
## "McMasterPandemic" is only used for unpack

# Then we select only the packages that aren't currently installed.
install.lib <- pkgs[!pkgs %in% installed.packages()]
# And finally we install the missing packages, including their dependency.
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
vapply(pkgs, library, logical(1L),character.only = TRUE, logical.return = TRUE)

ncores <- detectCores()-2 ## to be used in parallel computing function such as  mclapply (careful not to use up all cores!)

## multicore computation in Windows and Linux systems:
# future::plan(multicore, workers = ncores)
library(magrittr)
plan(multisession, workers = ncores)

## test igraph
mat <- matrix(0, 3, 3)
diag(mat) <- 0  # no self-loops

# Create a graphAM object
g <- new("graphAM", adjMat = mat, edgemode = "directed")
print(g)

g2 <- graph_from_adjacency_matrix(mat, mode = "directed")

## My functions: -----------------------------------------------------------  
here::here()
source(here::here("scripts", "fun_con.R"))

## Use Fipex functions: --------------------------------------------------------
# fipexpath <- "~/projects/DCI-R-Code-2020/2021 Debug/" ## AG's Linux machine

## clone the forked DCI machinery: git clone git@github.com:aligharouni/DCI-R-Code-2020.git
## go to the directory that you cloned DCI-R-Code-2020 package
fipexpath <- "C:/Ali/projects/DCI-R-Code-2020/2021 Debug/" ## AG's Windows Machine
source(paste0(fipexpath,"get_adj_matrix_from_gis_AG.r"))
source(paste0(fipexpath,"sum_fx.r"))
source(paste0(fipexpath,"dci_calc_fx_AG.r")) ## for getting dci_p and dci_s


# ##############################################################################
# Simulation Inputs
source(here::here("scripts","netwSimInputs.R"))

## Simulation control parameters:
ssize <- 60 ## sample size integer old was 60
nodeRange <-  c(3:10) #c(21,22,23,24,26,27,28,29) #c(25,30,35,40) #c(3:10) #c(11:20),  #c(45,50) ## 
## 3,7,15 are compelete binary graphs c(3:17)
treetypeRange <- as.list(c("linear","binary"))
## Is the network directed? If FALSE, it is assumed the passability up=down which gives a symmetric c_ij for the metapopulation cum pass matrix
directed <- c("Sym","Asym") #c(FALSE , TRUE) ## parameter which determines the symmetric of c_ij passability matrix 
## This is to set D parameter further down.
Dispersal <- c("Short","Long") ## readline(prompt="Enter the dispersal type Short/Long: ")
AllLengths <- "random" ## "unfixed":  for random lengths of reaches

## filename to save samp.l in ############################
## Note that these sample files can be really large and depends on your RAM and cpu could take a long time. The way around it, is to split sampling based on the network size.

fname0 <- "samp3_10" #"samp20_30"#"samp3"#"samp45_50" #"samp25_40" #"samp11_20" ## 

# Define the desired folder path
data_path <- here("data")

# Check if it exists, and create it if it doesn't
if (!dir.exists(data_path)) {
  dir.create(data_path, recursive = FALSE)
}

makefname<-function(fname0, test=TRUE){
  fnameini<-here::here("data",paste0(if(test)"test_",fname0,"_ini.Rdata")) ## initial sample
  fname<-here::here("data",paste0(if(test)"test_",fname0,".Rdata"))
  return(c(fnameini,fname))
  }

fnamevec<-makefname(fname0 = fname0, test = test)
fnameini<-fnamevec[1] ## unprocessed initial sample
fname<- fnamevec[2]## processed sample


rm(samp.lini)
if(exists("samp.lini")){
  cat("1. samp.l is already loaded, 2. make sure this is the right sample, thus check fname")
} else {samp.lini<-list() ## Network samples (lengths and passabilities list)
}

if(AllLengths=="fixed"){
  cat("All reaches have same length")  
  fixedlength <- 1 ## in (0,1] since these are relative lengths
} else{
  fixedlength<-NA
}


## This needs to be run once. if needed mode nodeRange, add to the sample_ini
set.seed(1234) # set the seed for reproducibility
for(nodeNum in nodeRange){
  temp.list <- list()
  ## Number of params to be sampled:
  ## nparams <- 2*nodeNum-1 # number of input parameters: n reaches, n-1 passibilities
  networksize <- nodeNum
  # passabilitysize <- 2*(nodeNum-1) ## if directed network assumed, ie alpha_up \neq alpha_down
  ## if directed==TRUE, for now we assume alpha_d=1 but alpha_u=randomized
  ## if directed==FALSE, alpha_u=alpha_d
  passabilitysize <- ifelse(directed=="Asym", nodeNum-1,nodeNum-1) ## for symmetric passability (ie when alpha_u=alpha_d, we sample only the pass up and equate up and down)
  nparams.all <- networksize+passabilitysize ## total number of parameters # number of input parameters: n reaches, n-1 passibilities
  nparams.rand <-ifelse(AllLengths=="fixed", nodeNum-1,nparams.all) ## num of randomized parameters
  nparams.fixed<- nparams.all-nparams.rand
  ## Note that the sampled passibilities are the product a_u*a*d, thus the sqrt of it should be considered in c_ij matrix (see passibility_cuml2())
  reachNames <-paste(1:nodeNum,"s",sep = "_") ## reaches name (eg: "1_s","2_s",...)
  barNames <- c(paste0(2:nodeNum,"u"),paste0(2:nodeNum,"d")) ## barrier names (Note there is a natural order on the names barID between vertex i-1 and i is tagged as i)
  ## Sampling part:--------------
  # set.seed(1234) # leaving set seed in the loop makes the sample i+1[element] be exactly the same as i[element], with the random 1 reach and passability different
  # sample of reaches' lengths and passabilities
  samp.df<- (list( df_fixed=data.frame(matrix(fixedlength,nrow=ssize,ncol=nparams.fixed)),
                  df_rand=data.frame(randomLHS(ssize, nparams.rand)),
                  df_passdown=data.frame(matrix(1,ssize, passabilitysize)) ## downstream passability
                  )
  %>% bind_cols() ## bind the columns of df_fixed and df_rand
  %>% purrr::set_names(c(reachNames, barNames)) ## segments first and barriers second
  %>% mutate(sample=c(1:ssize),
             n_seg=nodeNum
             # Dispersal=rep(Dispersal, each = ssize), ## AGnew this works for a vectorized Dispersal and directed
             # directed=rep(directed, each = ssize) ## AGnew
             ) 
  # %>% replace(is.na(.), 0) ## not sure if this is needed?
  ## FIXME: scaling the l_i's better to be done out of this loop 
  # %>% mutate(across(.cols = all_of(reachNames), ~ bumpupLength(.x))) ## bumpup lengths to (0.5,1) interval
  %>% mutate(across(.cols = all_of(reachNames), 
                    ~ case_when(AllLengths!="fixed"~ bumpupLength(.x), ## bumpup lengths to (0.5,1) interval
                                AllLengths=="fixed"~ (.x) ) ## If fixed mode
                    )) 
  %>% mutate(L = rowSums(across(.cols = all_of(reachNames)))) ## add scalar column
  )
  
  # Create a vector for Dispersal (Short, Long) and directed (FALSE,TRUE) meaning (Sym, Asym)
  Disp <- rep(Dispersal, each = nrow(samp.df))
  # Create a bonded data frame with copies of df and a column of factors
  samp.df <- cbind(samp.df[rep(seq_len(nrow(samp.df)), length.out = length(Disp)), ], Dispersal = Disp)
  ## directional info adding 
  Dire <- rep(directed, each = nrow(samp.df))
  samp.df <- cbind(samp.df[rep(seq_len(nrow(samp.df)), length.out = length(Dire)), ], directed = Dire)
  
  
  ## consider the same length/pass samples for all treetypeRange with the same # of nodes  
  paramsamp <- list()
  paramsamp <- map(.x=treetypeRange,.f=function(treetype)samp.df) ## make copies of samp.df for the treetypeRange
  names(paramsamp)<-treetypeRange ## make the same sample for both binary and linear topologies
  temp.df<- bind_rows(paramsamp, .id = "Topology")
  
  x<-interaction(unlist(treetypeRange),Dispersal,directed)
  ## split the df samples into list elements::
  
  temp.list<-split(temp.df, list(temp.df$Topology,temp.df$Dispersal, temp.df$directed,temp.df$n_seg,temp.df$sample),"_",drop=T)
  # temp.list <- (split(temp.df, seq(nrow(temp.df)),"par",drop=T)
  #               %>% set_names(unlist(map(.x=treetypeRange,
  #                                        .f=function(treetype)paste(treetype,nodeNum,seq(ssize),sep = "_"))   ))
  #               )
  samp.lini<-append(samp.lini,temp.list)
}

# Remove temp stuff!
rm(samp.df,temp.df,temp.list,paramsamp)
length(names(samp.lini))
length(unique(names(samp.lini)))
##save initial sample; 
if(saveit) saveRDS(samp.lini, file = fnameini)

## Process the initial sample ############################################
if(!saveit){samp.lini<-readRDS(fnameini)}

length(samp.lini)
samp.l<-vector("list", length=length(samp.lini))
samp.l<-samp.lini
## Added on May 22, 2023 by AG to incorporate the pass up=pass down upfront before making the sumtable 
## make passability down=up=sqrt(up) if directed =FALSE (this is symmetric network)
samp.l<- lapply(samp.l, function(samp){
  unpack(as.list(samp))
  ## barrier names (Note there is a natural order on the names barID between vertex i-1 and i is tagged as i)
  # barNames <- c(paste0(2:n_seg,"u"),paste0(2:n_seg,"d")) 
  
  if(directed=="Sym")  samp[paste0(2:n_seg,"d")]<-samp[paste0(2:n_seg,"u")]<-sqrt(samp[paste0(2:n_seg,"u")])
  return(samp)
})

if(saveit) saveRDS(samp.l, file = fname)

# make length vector scaled by L
samp.l<-lapply(samp.l, function(params.df) c(params.df,list(length.s=makelengthvec.s(params.df)) ))

## Make segment_matrix given the Topology and n_seg of network
samp.l<- lapply(samp.l, function(samp){
  unpack(as.list(samp))
  c(samp,list(segment_matrix=make_segment_matrix_igraph(Topology,n_seg)))
})

## Make pass df 
#> note that the pass data frame is created by makepassdf(), which AG made it directional
#> this directional passability matrix will be fed to sum_fx() from FiPex and resulted in the directional sum_table 
samp.l<- lapply(samp.l, function(samp){
  unpack(as.list(samp))
  barNames <- c(paste0(2:n_seg,"u"),paste0(2:n_seg,"d")) ## barrier names
  passdf<- (samp[barNames] ## extract passability vals
            %>% unlist
            %>% t()
            %>% as.data.frame()
            %>% set_names(barNames) ## Barrier names
  )
  c(samp,list(pass=makepassdf(passdf = passdf,n_seg,segment_matrix)))
} )

## Make adj matrix
samp.l<-lapply(samp.l, function(samp){
  unpack(as.list(samp))
  c(samp,list(adj_matrix=makeAdjMat(segment_matrix)))
})


## Make sum_table and sum_tableNull and add it to the list
add2listsumtable <- function(samp) {
  unpack(as.list(samp))
  passNull<-(pass %>% transform(Pass=1))
  ## sorted sum_tables
  sum_table<-segsort(makeSum_table(adj_matrix,pass,length.s))
  sum_tableNull<-segsort(makeSum_table(adj_matrix,passNull,length.s))
  return(c(samp,list( sum_table=sum_table, sum_tableNull=sum_tableNull)))
}

## Make sum_table and sum_tableNull
samp.l<-future_map(samp.l,.f=add2listsumtable)
if(saveit) saveRDS(samp.l, file = fname)

# ################################################################
# Calculate the DCIp and DCIs index 
# ################################################################
system.time(
  samp.ltmp<-lapply(samp.l, 
                      # mc.cores=ncores, ## mine is 8 core
                      FUN=function(samp){
                        unpack(as.list(samp))
                        dci<- dci_calc_fx_AG(sum_table,length.s,all_sections=T) ##calcDCIp(sum_table,length.s,all_sections=T)
                        list(DCIp=dci[[1]][["DCIp"]], DCIs=dci[[2]])
                      })
)
## append the calculated samp.ltmp to the main samp.l
samp.l <- map2(samp.l,samp.ltmp,c)
rm(samp.ltmp)
if(saveit) saveRDS(samp.l, file = fname)
# samp.l<-readRDS(fname)
#####################################
## Adding metapop params to samp.list from ./netwSimInputs.R
## params2 is a global param
samp.l<-lapply(samp.l, function(samp) c(samp,as.list(params2) ))

## Adding scaled metapop params as list ti samp.l
add2listpram.scaled <-function(samp){
  unpack(as.list(samp))
  ## read params2 included in the samp.l
  params.df.init<-data.frame(c=c,e=e,D=D,gamma=gamma,omega=omega,epsilon=epsilon)
  ## sacle params (these are tilda params)
  D.s=ifelse(Dispersal=="Long",L/L,L/(L*n_seg))
  ## scaled params df 
  params.s <-(params.df.init
              %>% mutate(e=e/(L)^omega, ##etilda
                         c=c*(L)^(gamma+epsilon), ##ctilda
                         e=e/c, ## This is delta tildda now, Note this is L dependent so we have delta_i for each sample i
                         c=1, ## ie c/c
                         ## D=runif(ssize, min = 0, max = 1) ## random dispersal distance
                         # D=L/L, ## Long dist dispersal; max cumulative relative length
                         # D=L/(L*n_seg), ## Short Distance dispersal
                         D=D.s,
                         gamma=gamma,
                         omega=omega,
                         epsilon=epsilon) 
              )
  return(c(samp,list(params.s=params.s)))
}

# samp.l <-llply(samp.l, .fun=add2listpram.scaled,.progress= progress_text(char = '+'))
samp.l <-lapply(samp.l, FUN=add2listpram.scaled)

## Adding vector l and f to the sample:
system.time(
  samp.ltmp<- lapply(samp.l,
                       # mc.cores=ncores, ## mine is 8 core
                       FUN=function(samp){
                         unpack(as.list(samp))
                         ## lengths vector unnamed
                         l.s <- length.s %>% select(Shape_Length) %>% unlist() %>% c() %>% unname()
                         f <-rep(1,n_seg) # Habitat quality of each segment, I am assuming f_i =1 for all segments
                         E <-E_ij(params.s,l.s,f) ## Extinction matrix n-by-n
                         list(l.s=l.s,f=f,E=E )
                       })
)
samp.l <- map2(samp.l,samp.ltmp,c)
rm(samp.ltmp)
if(saveit) saveRDS(samp.l, file = fname)

## sum_table dependent indices:
add2listJac <- function(samp){
  unpack(as.list(samp))
  d <- dij_cuml3(sum_table) ## mid of seg j to mid of seg i
  c_ijmetpop <- passibility_cuml3(sum_table)
  ##A Check::  if(!directed) isSymmetric(c_ijmetpop) should be TRUR
  M_land <- M_landscapematrix(params=params.s,d,l.s,f,c_ijmetpop)
  M_mod <- M_model(params=params.s,d,l.s,f,c_ijmetpop)
  Jac<-Jac_fipex2(params.s = params.s,E=E,M_mod=M_mod)
  ## Null model: ie fully passable network model
  c_ijmetpopNull <- matrix(1,nrow = n_seg,ncol = n_seg) ## A fully passable network, identical to passibility_cuml3(sum_tableNul, directed)
  M_landNull <- M_landscapematrix(params=params.s,d,l.s,f,c_ijmetpopNull)
  M_modNull <- M_model(params=params.s,d,l.s,f,c_ijmetpopNull)
  JacNull<-Jac_fipex2(params.s = params.s,E=E,M_mod=M_modNull)
  list(d=d, c_ijmetpop=c_ijmetpop, M_land=M_land, M_mod=M_mod,
       M_landNull=M_landNull,M_modNull=M_modNull,
       Jac=Jac,JacNull=JacNull )
}

samp.ltmp<-future_map(samp.l,.f=add2listJac)
samp.l <- map2(samp.l,samp.ltmp,c)
rm(samp.ltmp)
if(saveit) saveRDS(samp.l, file = fname)

## Lambda
#####################################
# Calculate the lambda, lambdaNull, lambda.s index to the list
system.time(
  samp.ltmp<- lapply(samp.l,
                       # mc.cores=ncores, ## mine is 8 core
                       FUN=function(samp){
                         unpack(as.list(samp))
                         # Jac<-Jac_fipex2(params.s = params.s,E=E,M=M)
                         lambda<-lambda_fipex2(Jac)
                         ## fully passable network model
                         # JacNull<-Jac_fipex2(params.s = params.s,E=E,M=MNull)
                         lambdaNull<-lambda_fipex2(JacNull)
                         ## outputs:
                         list(lambda=lambda,lambdaNull=lambdaNull) #lambda.s=lambda/abs(lambdaNull)
                       })
)
samp.l <- map2(samp.l,samp.ltmp,c)
rm(samp.ltmp)
if(saveit) saveRDS(samp.l, file = fname)

#####################################
## R0 Based on van den Driesshe & Watmough 2002 work:
#####################################

add2listR0<-function(samp){
  unpack(as.list(samp))
  FVi <- FVinv(params.s,E,M_mod)
  FViNull <- FVinv(params.s,E,M_modNull)
  ## a list contains R, left and right eigenvectors as well (we use the R for now)
  R.l<-Rnum2(params.s = params.s,FVinv = FVi) 
  RNull.l<-Rnum2(params.s = params.s,FVinv = FViNull)
  list(R=R.l[["R0"]],RNull=RNull.l[["R0"]] ) ##  we define later R.s=R.l[["R0"]]/ na_if(RNull.l[["R0"]],0)
}

samp.ltmp<-future_map(samp.l,.f=add2listR0)
samp.l <- map2(samp.l,samp.ltmp,c)
rm(samp.ltmp)
if(saveit) saveRDS(samp.l, file = fname)


## Calculate lambda_M from Ovaskainen Hanski 2001 paper
add2listlambda_M<-function(samp){
  unpack(as.list(samp))
  lambda_M<-lambdaM(M_land)
  lambda_MNull<-lambdaM(M_landNull)
  list(lambdaM=lambda_M,lambdaMNull=lambda_MNull) ##,lambdaM.s=lambda_M/lambda_MNull
}

samp.ltmp<-future_map(samp.l,.f=add2listlambda_M)
samp.l <- map2(samp.l,samp.ltmp,c)
rm(samp.ltmp)
if(saveit) saveRDS(samp.l, file = fname)


#####################################
## SSOI steady state occupancy index
#####################################
## Adding the initial state vec and simulation time to samp.l
# Test to remove an element of samp.l (Note thinking of rm ssoi and ssoi Null and recal with higher sim time)

# setting up the simulation sequence::
simtimes <- seq(from=0,to=500,by=2) ## set simulation time

## Using 2 variables passing to function (see https://stackoverflow.com/questions/14427253/passing-several-arguments-to-fun-of-lapply-and-others-apply)

add2listtarjs <- function(samp){
  unpack(as.list(samp))
  ## simulation parameters
  state_init <- c(rep(1,n_seg)) ## all 1
  if(lambda<=0) ## return initial state with vec 0
    trajs <- rbind(c(time=1,rep(0,n_seg) %>% setNames(paste(seq(1:n_seg))) )) %>% data.frame()
  else 
    trajs<-steppingODE2(model=metapop_mod2,params=params.s,M_mod=M_mod,E=E,state=state_init,times = simtimes)
  ## M_mod=M_modNull
  trajsNull<-steppingODE2(model=metapop_mod2,params=params.s,M_mod=M_modNull,E=E,state=state_init,times = simtimes)
  return(list(trajs=trajs, trajsNull=trajsNull, simtimes=simtimes))
  }


samp.ltmp<-future_map(samp.l,.f=add2listtarjs)
samp.l <- map2(samp.l,samp.ltmp,c)
rm(samp.ltmp)
if(saveit) saveRDS(samp.l, file = fname)

## indices from Trajectories
add2listsso<-function(samp){
  unpack(as.list(samp))
  sso<-ssi_readtrajs(trajs)
  ssoNull<-ssi_readtrajs(trajsNull)
  # sum_table_sorted <- segsort(sum_table) ## used only in the naming of the sections
  ## Note sum_table is already sorted
  lsso<-(sso ## local steady state
         %>% select(-ssoi)
         %>% set_names(unique(sum_table$start))
         %>% pivot_longer(cols=everything(),names_to = "sections",values_to = "lss") ## "lss": local steady state
         %>% data.frame()
  )
  lssoNull<-(ssoNull ## local steady state
             %>% select(-ssoi)
             %>% set_names(unique(sum_table$start))
             %>% pivot_longer(cols=everything(),names_to = "sections",values_to = "lss0") ## "lss": local steady state
             %>% data.frame()
  )
  dfmerg<-merge(lsso,lssoNull,by="sections")
  lsso.all<-(dfmerg %>% mutate(lss.s=lss/ na_if(lss0,0))) ##Division by 0 results in NA rather than Inf (see example below)
  return(list(ssoi=sso[["ssoi"]],ssoiNull=ssoNull[["ssoi"]],lsso.all=lsso.all)) ##,lsso=lsso,lssoNull=lssoNull
  }


samp.ltmp<-future_map(samp.l,.f=add2listsso)
samp.l <- map2(samp.l,samp.ltmp,c)
rm(samp.ltmp)
if(saveit) saveRDS(samp.l, file = fname)


#####################################
# ## Local indices: reproductive value of each reach
#####################################
## L2 norm scaled Left eigen vector (max_vec) corresponding to lambda max
## max_vec <- abs(max_vec)/sqrt(sum(max_vec^2)) ## scaled and absolute value of eigenvec
add2listrepro<-function(samp){
  ## Output: df of,sections,repro,reproNull,repro.s(ie= repro/reproNull)
  unpack(as.list(samp))
  repro<-repro_fipex_all2(sum_table=sum_table,Jac=Jac,lambda=lambda)
  ## Null model
  reproNull<-(repro_fipex_all2(sum_table=sum_tableNull,Jac=JacNull,lambda=lambdaNull) 
              %>% rename("reproNull"="repro"))
  ## combine columns and calc repro.scaled (repro.s)  
  repro.df <-merge(repro,reproNull,by="sections")
  return(list(repro.df=repro.df))
}

samp.ltmp<-future_map(samp.l,.f=add2listrepro)
samp.l <- map2(samp.l,samp.ltmp,c)
rm(samp.ltmp)
if(saveit) saveRDS(samp.l, file = fname)

#################################################################################
## Discrete Sensitivity Analysis of Global indicators wrt barriers
#################################################################################
## Specify the barriers to be removed 1-by-1 from the pass matrix and save it as a vector to work on. 
system.time(
  ## Trajectories ode
  samp.ltmp<-lapply(samp.l,
                      # mc.cores=ncores, ## mine is 8 core
                      FUN=function(samp,directed){
                        unpack(as.list(samp))
                        ## In each sample, we specify the barrier IDs to be removed as a vector 
                        ## Create a list of Barriers, and specified u/d? This way one can flexibly change the up/down given a barrier number
                        if(directed=="Asym") ##AG28Jan
                          ## for Asymmetric case, only remove the upward passability:
                          BarIDs_toberem <- map(c(2:n_seg),~paste0(.,c("u"))) %>% purrr::set_names((2:n_seg))
                        else
                          ## for symmetric case, remove both up and down passabilities since they are identical
                          BarIDs_toberem <- map(c(2:n_seg),~paste0(.,c("u","d"))) %>% purrr::set_names((2:n_seg))
                          # BarIDs_toberem <- grep("u",unique(pass[["Bar_ID"]]), value=TRUE)
                          ## I am thinking to keep a mechanism here for future if directed, do one thing otherwise do something else...
                          # if(directed)
                          #   BarIDs_toberem <- grep("u",unique(pass[["Bar_ID"]]), value=TRUE)
                          # else
                          #   BarIDs_toberem <- grep("u",unique(pass[["Bar_ID"]]), value=TRUE)
                          ## directional network: the upstream Bar_IDs will be removed 1-by-1 
                        return(list(BarIDs_toberem=BarIDs_toberem))
                        },directed=directed )
  )

samp.l <- map2(samp.l,samp.ltmp,c)
rm(samp.ltmp)
if(saveit) saveRDS(samp.l, file = fname)

## Based on the barriers specified in the vector BarIDs_toberem, update the passability dependent matrices (sum_table.upd,c_ijmetpop.upd,M.upd)
add2listsensMat<-function(samp){
  unpack(as.list(samp))
  ## In each sample do:
  BarIDs_toberem
  ## local function RETURN 
  ## sum_table.upd , c_ijmetpop.upd, M.upd which has 1 barrier removed based on what in Bar_ID_subset list is.
  update_sumtableWRTbarremoval <- function(BarId){
    ## sum_table.upd , c_ijmetpop.upd, M.upd which has 1 barrier removed based on what in Bar_ID_subset list is.
    c(Bar_ID=BarId, remove_Bar2(samp = samp, Bar_IDrem=BarId)) 
    }
  
  sensMat.l <- list() ## sensativity matrix list wrt to 1 barrier removal specified in BarIDs_toberem
  sensMat.l<- (future_map(BarIDs_toberem, .f=update_sumtableWRTbarremoval)
               %>% purrr::set_names(names(BarIDs_toberem)) 
               )
  return(sensMat.l)
}

## OLD
# add2listsensMat<-function(samp){
#   unpack(as.list(samp))
#   ## In each sample do:
#   BarIDs_toberem
#   sensMat.l <- list() ## sensativity matrix list wrt to 1 barrier removal specified in BarIDs_toberem
#   
#   sensMat.l<- (map(BarIDs_toberem,
#                    .f=function(BarId){
#                      c(Bar_ID=BarId,
#                        ## sum_table.upd , c_ijmetpop.upd, M.upd which has 1 barrier removed based on what in Bar_ID_subset list is.
#                        remove_Bar2(samp = samp, Bar_IDrem=BarId)) })
#                %>% purrr::set_names(names(BarIDs_toberem)) )
#   return(sensMat.l)
# }

## This takes a long time!
samp.ltmp<-future_map(samp.l,.f=add2listsensMat)
samp.l <- map2(samp.ltmp,samp.l,c) 
rm(samp.ltmp)
if(saveit) saveRDS(samp.l, file = fname)

## So far the passability-dependent matrices are updated based on removal of 1 barrier from BarIDs_toberem##

barrierSpecificglob.ind.r <- function(barid){
  ## assumeption:: smap is unpacked
  unpack(barid)
  Jac.r<-Jac_fipex2(params.s,E=E,M_mod=M_mod.upd)
  FVi.r<-FVinv(params.s,E=E,M_mod=M_mod.upd)
  lambda.r= lambda_fipex2(Jac=Jac.r)## lambda with removed 1 barrier
  R.r=Rnum2(params.s,FVinv=FVi.r)[["R0"]] 
  dci<- dci_calc_fx_AG(sum_table=sum_table.upd,length.s,all_sections=FALSE)
  DCIp.r<-dci[[1]][["DCIp"]]
  return(list(lambda.r=lambda.r,
              R.r=R.r,
              DCIp.r=DCIp.r/100,
              Bar_ID=  (unlist(barid[grep("Bar_ID",names(barid))]) %>% unname()) ##barid[["Bar_ID"]]
  )) ## could also include Jac.r=Jac.r,FVi.r=FVi.r
}


add2listglob.ind.r <-function(samp){
  unpack(as.list(samp))
  ## list of barrier-specific updated matrices when that barrier is removed.
  temp.l<-samp[names(BarIDs_toberem)]  ## eg. `2` includes "2u"
  glob.ind.r<-(future_map(temp.l,.f=barrierSpecificglob.ind.r)
               %>% bind_rows()
               %>% data.frame
               %>% mutate(lambda=lambda, ## Keep initial lambda and DCIp from samp.l for comparison
                          DCIp=DCIp/100,
                          # lambda.rr=lambda.r/lambda,DCIp.rr=DCIp.r/DCIp,  ## relative removeal indices
                          R=R
                          # R.rr=R.r/R
                          ))
  return(list(glob.ind.r=glob.ind.r))
  }


samp.ltmp<-future_map(samp.l,.f=add2listglob.ind.r)
samp.l <- map2(samp.l,samp.ltmp,c)
rm(samp.ltmp)
if(saveit) saveRDS(samp.l, file = fname)


end_time <- Sys.time()
end_time - start_time



################################################################################
## Step 2; some visual checking
################################################################################
# # str(samp.l)
# elemL.l <- list()
# 
# # elemL.l <- map(samp.l,~ unlist(.[c("L")]))
# elemL.glob <- map(samp.l,~ unlist(.[c("L", #"params.s",
#                                       "DCIp","lambda","lambdaNull",
#                                    "R","RNull","lambdaM","lambdaMNull","ssoi","ssoiNull")]))
# 
# 
# elemLglob.df<-dplyr::bind_rows(elemL.glob, .id="NetworkId")
# 
# elemLglob.df<-(elemLglob.df 
#            %>% separate(col = NetworkId,
#                          into = c("Topology","Dispersal","directed","Node","Sample"),sep="_", remove = FALSE)
#            %>% mutate(Sample=as.numeric(Sample))
#            %>% unite(col = "TDd", c(Topology,Dispersal,directed),sep="_",remove = FALSE )
#            %>% unite(col = "TD", c(Topology,Dispersal),sep="_",remove = FALSE )
#               )
# 
# ##
# pl.tmp1<-(elemLglob.df
#           %>% ggplot(aes(x=DCIp,y=(ssoi/ssoiNull), group=TDd))
#           + geom_point(aes(color=directed))
#           + facet_grid(TD~Node)
#           + theme_bw()
#           )
# print(pl.tmp1)
# 
# 
# ## df Long 
# elemLglob.dfLong<-(elemLglob.df
#                    %>% pivot_longer(cols=c("DCIp","lambda","lambdaNull",
#                                            "R","RNull","lambdaM","lambdaMNull",
#                                            "ssoi","ssoiNull"), names_to = "Indicator",
#                                     values_to = "value"))
# 
# str(elemLglob.dfLong)
# pl.tmp2<-(elemLglob.dfLong
#           %>% filter(Indicator %in% c("lambdaM","ssoi")) ## "lambda","R",
#           %>% ggplot(aes(x=Sample,y=value, group=TDd))
#           + geom_point(aes(color=Indicator, shape=directed, size=1) )
#           + scale_y_log10()
#           + facet_grid(TD~Node)
#           + theme_bw()
# )
# print(pl.tmp2)
# 





