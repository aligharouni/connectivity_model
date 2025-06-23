## By: Ali
## Created in Feb 2024 
#> Origin :C:/Ali/projects/connectivity/scripts/netwSim_local_all.R
#> 
# ##############################################################################
## Analysis: PART 1 of local indicators: barrier-removed effects on global indicators
# ##############################################################################

## A list of all dframes of global indicators when a barrier removed at a time (set =1)
# packages
pkgs <- c("dplyr","tidyverse","ggplot2", "purrr","broom",
          "Cairo","ggrastr","tikzDevice","latex2exp","GGally",
          "parallel",
          "directlabels","plotly",
          "here",
          ## multicore computation
          "future","furrr","purrr","magrittr")
vapply(pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)

## setwd("/home/ag/projects/connectivity/scripts")
ncores <- detectCores()-1 ## to be used in parallel computing function such as  mclapply (careful not to use up all cores!)
plan(multisession, workers = ncores)
## Plotting options::
theme_set(theme_bw())

## Functions required:#########################################
here::here() ## first restart R from session tab! You may need to do this twice.
source(here::here("scripts", "fun_con.R"))
## simulation parameters
source(here::here("scripts", "netwSimInputs.R"))

## DATA upload: ###############################################
list.files(here::here("data"))
## Note that the local1.Rdata is created in local1.R
fnamelocal2 <-here::here("data","local1.Rdata")
samp.comb<-readRDS(file = fnamelocal2)

## extract glob.ind.r from each sample
globInBarRem <- vector("list", length = length(samp.comb))
globInBarRem <- future_map(samp.comb,.f=function(samp){
  unpack(as.list(samp))
  df <- data.frame(glob.ind.r
                   ## adding the dispersal and directionality info as a column
                   %>% mutate(Dispersal=Dispersal,Direction=directed,Topology=Topology, Node=n_seg) 
  )
  return(df)} 
  ) 

## list to df
globInBarRem <- (globInBarRem %>% bind_rows(.id = "id"))

## Stats part:
globInBarRem.stats <- (globInBarRem
                       # Wrangling part:
                       %>% select(-"id")
                       %>% unite(col = "id", c(Topology,Node,Dispersal,Direction),sep="_",remove = TRUE ) ## create new id
                       %>% select(c("id","Bar_ID","DCIp.r","R.r")) ## favorite variables
                       %>% rename("Barrier"="Bar_ID")
                       # %>% mutate(Barrier=as.factor(Barrier))
                       # # Correlation part:
                       %>% pivot_longer(cols = c("DCIp.r"),names_to="x_var", values_to="x_val")
                       %>% pivot_longer(cols = c("R.r"), names_to="y_var", values_to="y_val")
                       %>% nest(data=c(x_val, y_val)) ##you get one row for each group defined by the non-nested columns
                       %>% mutate(cor_test = map(data, ~cor.test(.x$x_val, .x$y_val, method = "spearman",exact = FALSE)),
                                  tidied = map(cor_test, tidy))
                       %>% unnest(tidied)
                       %>% select(c(id,Barrier,x_var,y_var,estimate,p.value)) ## deselect to see all details of correlations
                       %>% data.frame()
                       %>% unite('Comparison',x_var,y_var,sep = "_",remove = FALSE)
                       %>% separate(col=id, into=c("Topology","Node","Dispersal","Direction"), sep="_") ##Node is number of nodes in the system
)
str(globInBarRem.stats,2)


## plot preping:
## Labels for the plots

RsTag<- paste0("(", "$DCI_p$", ", Basic reproduction number ", "$\\mathcal{R}_s$)")
# labvec <- c(DCIp.r_R.r = "corr$({DCI_p},\\mathcal{R}_s)$" )## OLD before May 2025 R: Basic reproduction number
labvec <- c(DCIp.r_R.r = RsTag )## R: Basic reproduction number

labvecDir<-c(Sym="Symmetric", Asym="Asymmetric")
labvecTop <-c(linear ="Linear",binary="Binary")
labvecDis<-c(Long="Global",Short="Local")




## Change the labels
# globInBarRem.stats <- mutate(globInBarRem.stats, across(Comparison, ~ labvec[.]))

BarRange.tmp0 <- (globInBarRem.stats 
                  # %>% filter(Direction=="Sym" & str_detect(Barrier, "u"))
                  %>% select(Barrier) 
                  %>% unlist 
                  %>% unique
) 


BarRange.tmp0
# Sort the vector based on integer first and then alphabet
BarRange.tmp <- str_sort(BarRange.tmp0, numeric = TRUE)
# Extract only the numeric part (if Sym)
BarRange.tmp.num <- str_extract(BarRange.tmp, "\\d+") %>% as.numeric()
# levBar<-paste0("Barrier id: ",sort(BarRange.tmp)) ## levels of id to sorting the plots (Barrier levels)
levBar<-paste0("Barrier id: ",(BarRange.tmp)) ## levels of id to sorting the plots (Barrier levels)

nodeRange.tmp <- unique(globInBarRem.stats[["Node"]])
levNodes<-paste0("Network size=",sort(as.numeric(nodeRange.tmp))) ## levels of id to sorting the plots

## The palette with black:
#> ref: https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
## check
pie(rep(1, length(colorBlindBlack8)), col = colorBlindBlack8)
## reorded to match the old plots in the manuscript
custom_colorBlindBlack8  <- c( "#009E73", "#D55E00", "#000000","#0072B2","#E69F00", "#F0E442",  "#CC79A7","#56B4E9")
#> colorblindness: read https://ggplot2.tidyverse.org/reference/scale_manual.html
#>  see also this paper: https://arxiv.org/abs/1903.06490 refed in the link above. Note scale_<aesthetic>_<datatype>_<colorscale>()


selectedNodes<-c(10,20,30,40,50)


## plot
# pl3_localcor_sym <-( globInBarRem.stats
#                      %>% mutate(Node=as.integer(Node))
#                      %>% filter(Direction=="Sym" & str_detect(Barrier, "u") & Node %in% selectedNodes) ## & Node<35
#                      %>% unite(col="TCD",c("Topology","Comparison","Dispersal","Direction"),sep = ": ",remove = FALSE)
#                      %>% unite(col="TC",c("Topology","Comparison"),sep = ": ",remove = FALSE)
#                      %>% mutate(TCD=as.factor(TCD), TC=as.factor(TC))
#                      %>% mutate(Barrier= str_extract(Barrier, "\\d+") %>% as.integer()) ## com this if you want 2u, 2d
#                       # create facet labels and change the labels
#                      %>% mutate(label=factor(paste0("Network size=",Node),levels = levNodes) )
#                      %>% mutate(across(Comparison, ~ labvec[.]))
#                      %>% mutate(across(Topology, ~labvecTop[.] ))
#                       ## plot part:
#                      %>% ggplot(aes(x=as.integer(Barrier),y=estimate,group=TCD )) ## old
#                       # %>% ggplot(aes(x=factor(Barrier, levels= BarRange.tmp.num),y=estimate,group=TCD ))
#                       + geom_point(aes(colour=TC, shape=Dispersal),size=2,alpha=1)
#                       + geom_line(aes(color=TC,linetype=Dispersal),linewidth=1,alpha=1)
#                       + facet_grid(Topology ~label)
#                       + theme(legend.position="bottom")
#                       + theme(legend.title=element_blank())
#                       + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#                      + scale_colour_manual(values = colorBlindBlack8) ## or select custom_colorBlindBlack8
#                       # + scale_colour_brewer(palette="Set1")
#                       + xlab("Barrier number")
#                       + ylab( "Correlation (Spearman)") #
#                       # + scale_x_continuous(breaks = nodeRange.tmp)
#                       + theme(panel.spacing=grid::unit(0,"lines"))
# )

test <-(globInBarRem.stats
        %>% mutate(Node=as.integer(Node))
        %>% filter(Direction%in%c("Sym","Asym") & str_detect(Barrier, "u") & Node %in% selectedNodes) ## & Node<35
        ## create facet labels and change the labels
        %>% mutate(label=factor(paste0("Network size=",Node),levels = levNodes) )
        %>% mutate(across(Comparison, ~ labvec[.]))
        %>% mutate(across(Direction, ~labvecDir[.] ))
        %>% mutate(across(Topology, ~labvecTop[.] ))
        %>% mutate(across(Dispersal, ~labvecDis[.] ))
        ## group the cols for faceting 
        %>% unite(col="TD",c("Topology","Direction"),sep = ":",remove = FALSE)
        %>% mutate(TD=factor(TD, levels=c("Binary:Symmetric","Linear:Symmetric","Binary:Asymmetric","Linear:Asymmetric")))
        %>% mutate(Barrier= str_extract(Barrier, "\\d+") %>% as.integer()) ## com this if you want 2u, 2d
        ## plot part:
        %>% ggplot(aes(x=as.integer(Barrier),y=estimate )) ## old
        # %>% ggplot(aes(x=factor(Barrier, levels= BarRange.tmp.num),y=estimate,group=TCD ))
        + geom_point(aes(shape=Dispersal),size=2,alpha=1)
        + geom_line(aes(linetype=Dispersal),linewidth=1,alpha=1)
        + facet_grid(TD ~label)
        + theme(legend.position="bottom")
        + theme(legend.title=element_blank())
        + theme(axis.text.x = element_text(angle = 90, hjust = 1))
        + scale_colour_manual(values = colorBlindBlack8) ## or select custom_colorBlindBlack8
        # + scale_colour_brewer(palette="Set1")
        + xlab("Barrier number")
        + ylab( "Correlation (Spearman)") #
        # + scale_x_continuous(breaks = nodeRange.tmp)
        + theme(panel.spacing=grid::unit(0,"lines"))
        )

print(test)


pl3_localcor_sym <-( globInBarRem.stats
                     %>% mutate(Node=as.integer(Node))
                     %>% filter(Direction%in%c("Sym","Asym") & str_detect(Barrier, "u") & Node %in% selectedNodes) ## & Node<35
                     # %>% unite(col="TCD",c("Topology","Comparison","Dispersal","Direction"),sep = ": ",remove = FALSE)
                     # %>% unite(col="TC",c("Topology","Comparison"),sep = ": ",remove = FALSE)
                     %>% unite(col="TD",c("Topology","Direction"),sep = ":",remove = FALSE)
                     %>% mutate(TD=factor(TD, levels=c("Binary:Symmetric","Linear:Symmetric","Binary:Asymmetric","Linear:Asymmetric")))
                     # %>% mutate(TCD=as.factor(TCD), TC=as.factor(TC))
                     %>% mutate(Barrier= str_extract(Barrier, "\\d+") %>% as.integer()) ## com this if you want 2u, 2d
                     # create facet labels and change the labels
                     %>% mutate(label=factor(paste0("Network size=",Node),levels = levNodes) )
                     %>% mutate(across(Comparison, ~ labvec[.]))
                     %>% mutate(across(Direction, ~labvecDir[.] ))
                     %>% mutate(across(Topology, ~labvecTop[.] ))
                     %>% mutate(across(Dispersal, ~labvecDis[.] ))
                     ## plot part:
                     %>% ggplot(aes(x=as.integer(Barrier),y=estimate )) ## old
                     # %>% ggplot(aes(x=factor(Barrier, levels= BarRange.tmp.num),y=estimate,group=TCD ))
                     + geom_point(aes(shape=Dispersal),size=2,alpha=1)
                     + geom_line(aes(linetype=Dispersal),linewidth=1,alpha=1)
                     + facet_grid(TD ~label)
                     + theme(legend.position="bottom")
                     + theme(legend.title=element_blank())
                     + theme(axis.text.x = element_text(angle = 90, hjust = 1))
                     + scale_colour_manual(values = colorBlindBlack8) ## or select custom_colorBlindBlack8
                     # + scale_colour_brewer(palette="Set1")
                     + xlab("Barrier number")
                     + ylab( "Correlation (Spearman)") #
                     # + scale_x_continuous(breaks = nodeRange.tmp)
                     + theme(panel.spacing=grid::unit(0,"lines"))
)

print(pl3_localcor_sym)

pl3_localcor_asym <-( globInBarRem.stats
                       %>% mutate(Node=as.integer(Node))
                       # %>% filter(Direction=="Sym" & str_detect(Barrier, "u") ) ##& Node<35
                      %>% filter(Direction=="Asym" & str_detect(Barrier, "u") & Node %in% selectedNodes)
                       %>% unite(col="TCD",c("Topology","Comparison","Dispersal","Direction"),sep = ": ",remove = FALSE)
                       %>% unite(col="TC",c("Topology","Comparison"),sep = ": ",remove = FALSE)
                       %>% mutate(TCD=as.factor(TCD), TC=as.factor(TC))
                      %>% mutate(Barrier= str_extract(Barrier, "\\d+") %>% as.integer())
                       # create facet labels
                       %>% mutate(label=factor(paste0("Network size=",Node),levels = levNodes) )
                       ## plot part:
                       %>% ggplot(aes(x=as.integer(Barrier),y=estimate,group=TCD )) ## old
                       # %>% ggplot(aes(x=factor(Barrier, levels= BarRange.tmp),y=estimate,group=TCD ))
                       + geom_point(aes(colour=TC, shape=Dispersal), size=2,alpha=1)
                       + geom_line(aes(color=TC,linetype=Dispersal), linewidth=1,alpha=1)
                       + facet_wrap(~label)
                       + theme(legend.position="bottom")
                       + theme(legend.title=element_blank())
                       + theme(axis.text.x = element_text(angle = 90, hjust = 1))
                      + scale_colour_manual(values = colorBlindBlack8) ## or select custom_colorBlindBlack8
                       # + scale_colour_brewer(palette="Set1")
                       + xlab("Barrier")
                       + ylab( "Correlation (Spearman)") #
                       # + scale_x_continuous(breaks = nodeRange.tmp)
                       + theme(panel.spacing=grid::unit(0,"lines"))
)
print(pl3_localcor_asym)


## Linux
# options(tikzLatex = 'C:/Users/Ali.Gharouni/AppData/Local/Programs/MiKTeX/miktex/bin/x64/pdflatex.exe')
## Windows:
options(tikzLatex = 'C:/Users/Ali.Gharouni/AppData/Roaming/TinyTeX/bin/windows/pdflatex.exe')
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}","\\usepackage[T1]{fontenc}",
                               "\\usetikzlibrary{calc}", "\\usepackage{amssymb}" ))

  
pl1name_sym <- here::here("plots","draftplots","localInd2_sym.tex")
pl1name_asym <- here::here("plots","draftplots","localInd2_asym.tex")

ggsave(pl3_localcor_sym
       # + ggtitle("Long vs Short D")
       ,
       device = tikz,
       filename = pl1name_sym,
       standAlone = TRUE,
       width = 8, height = 7, units = "in")

ggsave(pl3_localcor_asym
       # + ggtitle("Long vs Short D")
       ,
       device = tikz,
       filename = pl1name_asym,
       standAlone = TRUE,
       width = 8, height = 7, units = "in")



## IN windows:: (getwd() first;)
getwd()
tools::texi2pdf(pl1name_sym)
tools::texi2pdf(pl1name_asym)













## OLD and TESTS ######################################################################### 


# pl3_localcor_asym <-( globInBarRem.stats
#                       %>% mutate(Node=as.integer(Node))
#                       # %>% filter(Direction=="Sym" & str_detect(Barrier, "u") ) ##& Node<35
#                       %>% filter(Direction=="Asym" & str_detect(Barrier, "u") & Node<55)
#                       %>% unite(col="TCD",c("Topology","Comparison","Dispersal","Direction"),sep = ": ",remove = FALSE)
#                       %>% unite(col="TC",c("Topology","Comparison"),sep = ": ",remove = FALSE)
#                       %>% mutate(TCD=as.factor(TCD), TC=as.factor(TC))
#                       # create facet labels
#                       %>% mutate(label=factor(paste0("Network size=",Node),levels = levNodes) )
#                       ## plot part:
#                       # %>% ggplot(aes(x=as.integer(Barrier),y=estimate,group=TCD )) ## old
#                       %>% ggplot(aes(x=factor(Barrier, levels= BarRange.tmp),y=estimate,group=TCD ))
#                       + geom_point(aes(colour=TC, shape=Dispersal),size=1)
#                       + geom_line(aes(color=TC,linetype=Dispersal))
#                       + facet_wrap(~label)
#                       + theme(legend.position="bottom")
#                       + theme(legend.title=element_blank())
#                       + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#                       + scale_colour_brewer(palette="Set1")
#                       + xlab("Barrier")
#                       + ylab( "Correlation (Spearman)") #
#                       # + scale_x_continuous(breaks = nodeRange.tmp)
#                       + theme(panel.spacing=grid::unit(0,"lines"))
# )
# print(pl3_localcor_asym)






