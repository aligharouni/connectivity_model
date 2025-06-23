# Functions for analysis of the interplay metapopulation and connectivity component of the model 
# By Ali Gharouni (May 13, 2023)

# while (!require(remotes)) {
#   install.packages("remotes")
# }
# ## install development version of bbmle:
# if (!require("bbmle") || packageVersion("bbmle") < "1.0.23.5") {
#   remotes::install_github("bbolker/bbmle")
# }
# ## install the target package and all its dependencies:
# while (!require(McMasterPandemic)) {
#   remotes::install_github("bbolker/McMasterPandemic",
#                           dependencies = TRUE #,
#                           # build_vignettes = TRUE
#   )
# }



# unpack <- McMasterPandemic::unpack
# the unpack func was just used from McMasterPandemic package, no need to install the whole pack!
unpack<- function (x) {
  if (any(names(x) == "")) 
    stop("unnamed elements in x")
  invisible(list2env(as.list(x), envir = parent.frame()))
}

# ####################################################
# Length vector scaling and op
# ####################################################

## scale up the reaches' length samples to (0.5,1)  
bumpupLength <- function(lvec) {
  # l is a vector (row of df l_1,...,l_n)
  lvec <- 0.5+lvec*0.5 ## project reach length to (0.5,1) interval
  return(lvec)
}

# transform the reach abs length to relative length, ie l_i -> l_i/L
scaleLength <- function(lvec) {
  # l is a vector (row of df l_1,...,l_n,L) where L is the scaler
  # class(lvec) <- "params_pansim"
  unpack(as.list(lvec))
  # L <- sum(lvec, na.rm = TRUE)
  lTilda <- (lvec[names(lvec)!=c("L")])/L
  return(lTilda)
}

## vectorize and scale the lengths by L
makelengthvec.s<-function(params.df){
  unpack(as.list(params.df))
  dfout<-(params.df
          %>% select(c(paste(1:n_seg,"s",sep = "_"),L))
          %>% apply(.,1,FUN = scaleLength) ## sacle by L
          %>% t() 
          %>% as.data.frame()
          %>% rename("sink"="1_s")
          %>% pivot_longer(cols = everything() ,names_to = "Seg_ID", values_to = "Shape_Length")
          %>% data.frame()
  )
  return(dfout)
}


# ####################################################
## Functions handling summary table (see Fipex), and translate table to Matrix metapopulation model inputs ----------------
# ####################################################
## Make similar to segments_and_barriers in Fipex but for an igraph This will be just used in making adj matrix
# segments_and_barriers


## Function to remove the loops, create a column of barrier Id
## eg: extract_barrierID(c("2_s","3_s")) is 3 which is the max #, "sink" 
## "Seg_ID"="from","Seg"="to" so it is from 2_s to 3_s, low reach to higher reach, so the passability is upward, gives 3u 
## try this:
# vecin <-c("sink","2_s")
# vecin <- gsub("sink", "1_s", vecin) ## replace "sink" by "1_s"
# tmp <- as.numeric(gsub("\\D", "", vecin))
# tmp
# ifelse(tmp[1]>tmp[2], paste0(max(tmp),"d"),paste0(max(tmp),"u"))

extract_barrierID <- function(vecin){
  ## from <- vecin[1], to <- vecin[2]
  ## input: vector of c(from,to)
  ## output: barrier_ID as integer
  vecin <- gsub("sink", "1_s", vecin) ## replace "sink" by "1_s"
  tmp <- as.numeric(gsub("\\D", "", vecin)) ##remove underscore "_"
  Barrier_ID <- ifelse(tmp[1]>tmp[2], paste0(max(tmp),"d"),paste0(max(tmp),"u"))
  return((Barrier_ID))
}



make_segment_matrix_igraph<-function(Topology,n_seg){
  ## Inputs: Topology =c("binary","linear"), n_seg: is int representing the # of reaches
  childNum <- ifelse(Topology=="binary",2,1)
  g <- graph.tree(n=n_seg, children=childNum, mode = "in")  ##(with no loop, we need loop for Fipex setting later)
  ## Rename vertices compatible with Fipex setup, starting with sink at the mouth of watershed and expand up
  V(g)$name <- c("sink", paste(2:vcount(g),"s",sep = "_"))
  EdgeList <- get.data.frame(g, what="edges") 
  ## TODO: find a better/smarter tidyverse version of combining df1 and df2: 
  # reverse the to-from paths and bind them into df_seg:
  df1 <- get.data.frame(g, what="edges")
  df2 <- transform((df1 %>% filter(to!=from)), from=to,to=from)
  segment_matrix_igraph <- (rbind(df1,df2) %>% rename("Seg_ID"="from","Seg"="to"))
  return(segment_matrix_igraph)
}


## make passibility df required in making pass df, and sum_table 
## (eg of the input samp.l[["binary_2_1"]]) 
makepassdf<-function(passdf, n_seg,segment_matrix){
  # #extract df of barriers' pass
  ## Input passdf: dataframe wide format of passabilities, n_seg: is int representing the # of reaches
  # if( ncol(passdf)-n_seg+1) {print("Error::mismatch# of barriers"); return()}
  dfpass<-(passdf
           %>% pivot_longer(cols = everything() ,names_to = "Bar_ID", values_to = "Pass")
           %>% data.frame()
  )
  df_seg <- (segment_matrix
             %>% filter(Seg_ID!=Seg)
             %>% mutate(Bar_ID= apply(.,MARGIN=1, extract_barrierID),nat_barrier = FALSE) ## MARGIN=1 is row-wise
             %>% unite("section1_2",c(Seg_ID,Seg),remove = F, sep = ",")
             %>% arrange(Bar_ID)
             %>% rename("Seg_1"="Seg_ID","Seg_2"="Seg")
  )
  pass<-merge(df_seg,dfpass,by="Bar_ID")
  return(pass)
}

##  Building tree/graph/network Structure: ################ 
# (eg of the input samp.l[["binary_2_1"]]) which is a sample
makeAdjMat <- function(segment_matrix){
  ## Input: segment_matrix is a data frame in long format col: Seg_ID, Seg
  ## get adj_matrix from segment_mat using Fipex functions 
  adj_matrix <- get_adj_matrix_from_gis_AG(inputname=segment_matrix)
  return(adj_matrix)
}

# This fun returns the cummulative sum of the lengths of a given path by looking up the segments from lengths df
# path<-"2_s,sink,3_s"
cumlength <- function(length.s, path){
  # length is a df with columns: Seg_ID, Shape_Length
  # path is a string of "seg1,seg2,...,segn"
  pathvec<-unlist(strsplit(path,","))
  o <-as.numeric(length.s  
                 %>% filter(Seg_ID %in% pathvec) 
                 #%>% summarise(sum(Shape_Length))
                 ## the lengths are from center to center of patches (specifically the dist between 2 adjacent node is half of sum(l_i,l_j)) 
                 %>% summarise(sum(Shape_Length)-0.5*(Shape_Length[Seg_ID==pathvec[1]]+Shape_Length[Seg_ID==pathvec[length(pathvec)]]))
  )
  
  return(o)
}


## Sum_table
makeSum_table <- function(adj_matrix,pass,length.s){
  # FIPEX function
  sum_table <- sum_fx(adj_matrix=adj_matrix,passability=pass,lengths=length.s) ## all these 3 inputs are already in samp.l
  ## MY addings:: Add the length-along the path to the sum_table 
  cum_length <- apply(sum_table %>% select(path2), 1, FUN = function(par) cumlength(length.s,par))
  sum_table <- (sum_table  %>% mutate(cum_length=cum_length)) ## update sum_table with cummalitive length along reaches of the network
  return(sum_table)
}

## distance along all paths in a network defined by the sum_table from FiPex package
dij_cuml3 <- function(sum_table_sorted){
  # sum_table_sorted <- segsort(sum_table)
  ## output: a matrix (dim=n_seg by n_seg) of distances along network, symmetric
  n_seg <- length(unique(sum_table_sorted$start))
  dij <- matrix(sum_table_sorted$cum_length, ncol = n_seg)
  diag(dij) <- 0 ## Changes Jul 11
  return(dij)
}

# Sort out the sum_table from sink to other segments
## sort is required to bring the sink  at first fir Jacobian, ... 
segsort <- function(sum_table){
  # Input: sum_table, dataframe, created in Ali_DCI0.R ("~/projects/DCI-R-Code-2020/2021 Debug/")
  # sum_table is a dataframe with a column start (refering to start segments)
  o <- (sum_table 
        %>% mutate(across( ## subs sink by 1_s in start and end columns (goal to order start&end 1_s, 2_s,...)
          .cols = c(start,end),
          .fns = ~ dplyr::if_else(stringr::str_detect(.x, "sink"), "0001", .x)
        ))
        %>% mutate(startNum=as.numeric(gsub('_s', '', start)),endNum=as.numeric(gsub('_s', '', end)))
        %>% group_by(startNum, endNum)
        %>% arrange(startNum, endNum, .by_group = TRUE) ## sort
        %>% ungroup()
        %>% mutate(across( ## back substitute 1_s by sink
          .cols = c(start,end),
          .fns = ~ dplyr::if_else(stringr::str_detect(.x, fixed("0001")), "sink", .x)
        ))
        %>% select(-c(startNum, endNum))
        %>% data.frame()
  )
  return(o)
}




#########
## Metapopulation indicators:

## Extinction matrix
E_ij <- function(params,l,f){
  unpack(as.list(params))
  n_seg <- length(l)
  E <- diag(e/(l*f)^omega,nrow = n_seg,ncol=n_seg)
  return(E)
}
## The c_ij from the sum_table is not exactly the c_ij used in the metapopulation model, namely, c_ijtilda in the text. This function derive c_ijtilda from the sum_table
# passibility_cuml3 <- function(sum_table_sorted, directed=c(FALSE,TRUE)){
#   # This is c_ij in Cote 2009, with c_11 is the sink
#   ## output: Matrix of cumulative passibility from patch j to patch i
#   n_seg <- length(unique(sum_table_sorted$start))
#   c_ij <- matrix(sum_table_sorted$pathway_pass, ncol = n_seg) ## c_ij is the cum passability from j to i matching the matrix model notation (This works as long as working with sorted sum table  ## see matrix(c(1:9),ncol=3))
#   if(directed) c_ijmetpop<-c_ij else c_ijmetpop<-sqrt(c_ij) ## Since in sum_table the pathway_path is the a_ui*a_di
#   return(c_ijmetpop)
# }
passibility_cuml3 <- function(sum_table_sorted){
  # This is c_ij in Cote 2009, with c_11 is the sink
  ## output: Matrix of cumulative passibility from patch j to patch i
  n_seg <- length(unique(sum_table_sorted$start))
  c_ij <- matrix(sum_table_sorted$pathway_pass, ncol = n_seg) ## c_ij is the cum passability from j to i matching the matrix model notation (This works as long as working with sorted sum table  ## see matrix(c(1:9),ncol=3))
  c_ijmetpop<-c_ij ## Since in sum_table the pathway_path is the a_ui*a_di
  return(c_ijmetpop)
}

## c_ijmetpop<-passibility_cuml3(sum_table_sorted, directed=TRUE)

## distance along all paths in a network defined by the sum_table from FiPex package
dij_cuml3 <- function(sum_table_sorted){
  # sum_table_sorted <- segsort(sum_table)
  ## output: a matrix (dim=n_seg by n_seg) of distances along network, gives symmetric matrix
  n_seg <- length(unique(sum_table_sorted$start))
  dij <- matrix(sum_table_sorted$cum_length, ncol = n_seg)
  diag(dij) <- 0 ## Changes Jul 11
  return(dij)
}

# d <- dij_cuml3(sum_table_sorted)

## Matrix M, landscape matrix, form Ovaskainan & Hanski 2001 paper
M_landscapematrix <- function(params.s,d,l,f,c_ijmetpop){
  ## Inputs:
  ## params.s: scaled parameters
  ## d <- dij_cuml2(sum_table_sorted) ## Matrix of cumulative distances from seg j to seg i
  ## l: vector of length[i] for each segment
  ## f: f <-rep(1,n_seg) # I am assuming f_i =1 for all segments
  unpack(as.list(params.s))
  # l<-l.s
  n_seg <- length(l)
  M <- matrix(NA,ncol = n_seg, nrow = n_seg)
  dispersal <- exp(-d/D)-diag(n_seg) ##dispersal=  exponential-Identity
  source_sec <- diag((l*f)^epsilon,nrow = n_seg,ncol=n_seg) ## source section contribution
  sink_sec <- diag((l*f)^(gamma+omega),nrow = n_seg,ncol=n_seg) ## sink (or destination) contribution
  
  M <- sink_sec%*%(dispersal*c_ijmetpop)%*%source_sec ## Notice dispersal*c_ij is Hadamard product
  return(M)
}

## model M
M_model <- function(params.s,d,l,f,c_ijmetpop){
  ## Inputs:
  ## params.s: scaled parameters
  ## d <- dij_cuml2(sum_table_sorted) ## Matrix of cumulative distances from seg j to seg i
  ## l: vector of length[i] for each segment
  ## f: f <-rep(1,n_seg) # I am assuming f_i =1 for all segments
  unpack(as.list(params.s))
  # l<-l.s
  n_seg <- length(l)
  M <- matrix(NA,ncol = n_seg, nrow = n_seg)
  dispersal <- exp(-d/D)-diag(n_seg) ##dispersal=  exponential-Identity
  source_sec <- diag((l*f)^epsilon,nrow = n_seg,ncol=n_seg) ## source section contribution
  sink_sec <- diag((l*f)^(gamma),nrow = n_seg,ncol=n_seg) ## sink (or destination) contribution
  M <- sink_sec%*%(dispersal*c_ijmetpop)%*%source_sec ## Notice dispersal*c_ij is Hadamard product
  return(M)
}


## M<-M_landscapematrix(params=params,d,l,f,c_ijmetpop)
lambdaM <- function(M_land){
  # M is M_landscapematrix
  lambdaout <- max(eigen(M_land)$values) ## principle eigenvalue 
  return(lambdaout)
}

Jac_fipex2 <- function(params.s,E,M_mod){
  unpack(as.list(params.s))
  # Note that this M is different from the M in Hanski & Ovaskainan 2001 paper (see FL handwritten notes)
  J <- c*M_mod-E ## Jacobian at p=0
  return(J) # Jacobian Matrix
}

FVinv <- function(params.s,E,M_mod){
  unpack(as.list(params.s))
  out <- c*M_mod %*% solve(E) ## following Watmough 2002 R0 calculation, Note solve() gives inv of a matrix
  return(out)
}
 

# Jac<-Jac_fipex2(params=params,E=E,M=M)

lambda_fipex2 <- function(Jac){
  ## Input: Jacobean matrix
  ## Output the principal eigenvalue
  lambdaout <- max(Re(eigen(Jac)$values)) ## principle eigenvalue of Jaccobian Matrix
  return(lambdaout)
}

# # A function to return R0 and the corresponding eigenvectors -- left and right
Rnum2 <- function(params.s,FVinv){
  unpack(as.list(params.s))
  eigt <- eigen(t(FVinv)) ## transpose of J, to calculate the left eigenvector (Note eig=eigt)
  # return max eigenval and the corresponding left and right eigen vectors:
  lambdamax <- max(Re(eigt$values)) ## principle eigenvalue
  ## Right eigenvector
  j<-which.max(Re(eigt$values)) ## index of the max real part eigen val
  eigvecR <- abs(eigt$vectors[,j])
  ## Left eigenvector (representing reproductive values)
  eigvecL <- abs(eigt$vectors[,j])
  return(list(R0=lambdamax,eigvecR=eigvecR,eigvecL=eigvecL))
}

################################################################################
## Trajectories::

metapop_mod2 <- function(time,state,params,M_mod,E){
  unpack(as.list(params))
  ## first extract the state variables
  p <- state
  ## Colonization rate for patches
  C <- c*M_mod%*%p ## this is a n by 1 vector
  ## Equations ODE, where pi(t) is the probability that patch i is occupied at time t by the focus species
  dpdt <- C*(1-p)-E%*%p ## column vector
  return(list(dpdt))
}

steppingODE2 <- function(model,params,M_mod,E,state,times,tol=1e-6,checkpoint=100){
  unpack(as.list(params))
  ## Inputs:: model: ode specifically metapop_mod(time,state,params,sum_table), 
  ## times is a vector seq(from=0,to=,by=), state is init sate vec (1,1,1,...)
  ## params eg samp.l[["linear_2_1"]][["params.s]] which includes c("c","e","D","gamma","omega","epsilon") 
  ## checkpoint: int 
  # rtol <- 1.0e-4 ## defaults: rtol = 1e-6, atol = 1e-6)
  # atol <- 1.0e-6
  stepsz <- times[2]-times[1] ## stepsize of showing the simulation's result
  # n_seg <- length(unique(sum_table$start)) ## numb of segments
  out <- c() ##output
  timecyc <-floor(length(times)/checkpoint) ## the num of time cycles for forloop
  ## run ode till checkpoint simulationtime
  for(i in 1:timecyc) { ## i is cycle through simtimes vec
    o <- c()
    tstart<-((i-1)*checkpoint) #+1
    tend<-(i*checkpoint)
    if(i==timecyc) tend<-length(times) ## To avoid NSs ate the last cycle
    o <- ode(func=model,
             y=state,
             times= times[tstart:tend],
             parms=params,
             M=M_mod,
             E=E
             # atol = atol, ##rtol, and atol determine the error control performed by the solver. 
             # rtol = rtol
    )
    # print(o)
    if(i==1) out <- o[1,]
    d1<-dist(rbind(o[nrow(o)-1, -1] , o[nrow(o), -1])) ## check the distance of the last 2 sol
    # print(d1)
    out <- rbind(out, o[nrow(o),] )
    if (d1/stepsz < tol) break ##  AG added after talking to FL (stepsz)
    state <- o[nrow(o),-1]
  }
  out<-data.frame(rbind(out))
  rownames(out) <- NULL
  return(out)
}
## trajs<-steppingODE2(model=metapop_mod.test,params=params.s,M=M,E=E,state=state_init,times = simtimes)

# ## The old ssi function input list of trajectories stored intreelistsamp and ave the last row
ssi_readtrajs <- function(trajdf){
  # # input: trajname=c("trajs","trajsNull")
  nodeNum<-ncol(trajdf)-1 ## the first col is Time
  suffix <- paste(seq(1:nodeNum))
  # # ssoi alongside the equilibrium
  o <- data.frame(
    # ssi=mean(as.numeric(trajdf[nrow(trajdf),paste(seq(1:nodeNum))]))*100 ## ave SS occupancy
    trajdf[nrow(trajdf),colnames(trajdf)!="time"], ## steady state variables= the last row of trajdf
    ssoi=mean(as.numeric(trajdf[nrow(trajdf),colnames(trajdf)!="time"])) ## *100 ave SS occupancy
  )
  return(o)
}

# ###########################################################
### Sensativity functions in discrete sense Remove a barrier (make it 1) 
# ###########################################################
## remove_Bar2 is similar to remove_Bar() is updated for directional passability removal (Date of update: May 18, 2023 by AG)

## Note that makeSum_table <- function(adj_matrix,pass,length)
remove_Bar2 <- function(samp, Bar_IDrem){
  ## Inputs:: samp from samp.l
  ## Here just the passability and cumpass are modified not the length 
  ## Bar_IDrem is a vector of "int""u/d" (eg "2u", "2d", or Bar_IDrem <-c("2u","3u") ) in [["pass"]][["Bar_ID"]]. 
  ##Thus you can remove a set of Barriers and return a single sum_table with all Bar_IDrem removed.
  ## Out:: sum_table, c_ijmetpop, M shown by .upd with Bar_removed (=1)
  unpack(as.list(samp))
  if(any(!Bar_IDrem%in%unique(pass[["Bar_ID"]]))) {print("Error::Bar_removed is not in right range");return()}
  sum_table_sorted <- sum_table
  dimsum<-dim(sum_table_sorted)
  ## sum_table.upd: updated
  sum_table.upd <- (matrix(NA,dimsum[1],dimsum[2]) %>% data.frame())  ## space allocation 
  # update passability and sum_table
  pass.upd <-(pass %>% mutate(across(.cols = Pass, ~ ifelse(Bar_ID%in%Bar_IDrem,1,Pass) )))
  sum_table.upd.unsorted <- sum_fx(adj_matrix=adj_matrix,pass=pass.upd,lengths=length.s)
  ## updated elements:
  ## FIXME: try not to use the segsort, match the input and output sum_tables row somehow!
  sum_table.upd <-cbind(segsort(sum_table.upd.unsorted), cum_length=sum_table_sorted[["cum_length"]])
  c_ijmetpop.upd <- passibility_cuml3(sum_table_sorted = sum_table.upd)
  M_land.upd <- M_landscapematrix(params=params.s,d=d,l=l.s,f=f,c_ijmetpop=c_ijmetpop.upd)
  M_mod.upd <- M_model(params=params.s,d=d,l=l.s,f=f,c_ijmetpop=c_ijmetpop.upd)
  return(list( 
    ## updated passability-related matrices
    sum_table.upd=sum_table.upd, 
    c_ijmetpop.upd=c_ijmetpop.upd,
    M_land.upd=M_land.upd,
    M_mod.upd=M_mod.upd
    ## unchanged (so not related to pass but needed for other indicators) params.s=params.s,E=E
    )  
    ) 
}

### calculate a local index: reproduction value
# ###########################################################
# left eigenvec reproduction ----------------------------------------------
# see Caswell book Applied Mathematical Demography page 209 
# Thus, if we take “the contribution of stage i to long-term population size”
# as a reasonable measure of the “value of stage i,” the left eigenvector v1
# gives the relative reproductive values of the stages (Goodman 1968, Keyfitz
#                                                       1968). 
# \cite{goodman1968elementary}

## calculate all sectional reproductive values as a vector,  
repro_fipex_all2 <- function(sum_table,Jac,lambda){
  # input :reach_index integer (this will be the component of the left eigenvector corr to lambda0)
  # Output: eigenvec corresponding to lambda, as df with 2 columns: sections (segment names from sum_table) and repro values
  Jt <- t(Jac)## transpose Jacobian at p=0
  max_indx <- which(Re(eigen(Jac)$values)==lambda)
  max_vec <- eigen(Jt)$vectors[,max_indx] ##left eigenvector
  ## FIXME: abs() FL NOV24 fixed (it was abs(max_vec)/sum(max_vec) )
  max_vec <- abs(max_vec)/sqrt(sum(max_vec^2)) ## scaled and absolute value of eigenvec
  ## Q: why L2 norm and why scaled by Sum of Squares 
  reproduction <- data.frame(sections=unique(sum_table$start), 
                             repro=max_vec)
  return(reproduction)
}

################################################################################
## FUnctions for visualizations and further analysis

## pull any element either a scaler, vector, df, or list
pullelem <-function(samp,elem){
  ## Inputs: samp is a list, elem is a vector of chars included in the samp
  unpack(as.list(samp))
  samp[elem]
}
## test pullelem:
# l1 <- list(a=1, b=data.frame(x=1:4,y=5:8), c=list(q=1,p=2)
#            )
# pullelem(samp = l1, elem = c("a","b"))

# pull indices 1 dim data (such as lambda, DCIp,...) from treelistcomb into a dataframe 
pullInd <- function(samp,elem){
  unpack(as.list(samp))
  # pull elements from this lnested list treelistcomb[[typenode]][[sample]][[elements]]
  return(samp[elem] %>% data.frame %>% set_names(elem) )
}

## binding the lists with id into df of indices
bindtodf <- function(listin){
  dfout<- (bind_rows(lapply(listin,"["), .id = "id")
           %>% data.frame()
           # %>% mutate(sim="old")
  )
  return(dfout)
}

## This function read a samp.l RDS and extarct indexvec from it and returndf
savewhatnecassary<-function(samp.lname,indexvec=globindex){
  # eg: globindex <- c("DCIp","lambda","Dispersal","directed",...) # 
  samp_l <-readRDS(file = samp.lname) 
  ## pull the indices from samp.l 
  # samp_lpart <- mclapply(samp_l, mc.cores=ncores, FUN=function(samp)pullInd(samp,elem=indexvec))
  samp_lpart <- lapply(samp_l, FUN=function(samp)pullInd(samp,elem=indexvec))
  samp.df<-bindtodf(samp_lpart)
  return(samp.df)
}

## This fun pull any form of data from a file and return it as is (sublist a list)
savewhatnecassary2<-function(samp.lname,indexvec){
  ## eg: globindex <- c("DCIp","lambda","Dispersal","directed",...) # 
  samp_l <-readRDS(file = samp.lname)
  unpack(as.list(samp_l))
  ## the vectorized form
  ## samp_l <- (samp.lname %>% map(readRDS) %>% unlist) ## read multiple files at the same time
  ## pull the indices from samp.l 
  # samp_lpart <- mclapply(samp_l, mc.cores=ncores, FUN=function(samp) samp[indexvec]   )
  samp_lpart <- lapply(samp_l, FUN=function(samp) samp[indexvec]   )
  return(samp_lpart)
}


## Read a fname which has a list, and read the elements sublist
load_sublist<-function(fname,element){
  samp_l <-readRDS(file = fname)
  samp_lpart <- lapply(samp_l, FUN=function(samp) samp[element]   )
  return(samp_lpart)
}









