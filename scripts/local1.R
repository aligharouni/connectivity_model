

## Origin: C:/Ali/projects/connectivity/scripts/netwSim_PlotLocalLSDisp.R

# ##############################################################################
## PART 2 of local indicators: reach-based effects
# ##############################################################################


## Indicators include: DCI_s (structural), and reproductive value of each reach (principal eigenvector) and local steady state occupancy (lss)  
locindex <- c("DCI_s","repro","lss") ## column names in the merged df

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

## Indicators to pull
localindex <- c("BarIDs_toberem","lsso.all","Dispersal","DCIs","repro.df","glob.ind.r")
## adding some other params to be read
localindex <- c(localindex,c("Topology","Dispersal","n_seg","directed"))

## these files include all info from Long/Short dispersal, Sym/Asym network, linear/binary networks
## these files are created ~/projects/connectivity/tests/test_pipeline20240216.R

## Global params
FirstRun <- FALSE
saveplots<-TRUE ## ggsave?
## Target folder for the plots
targetFolder <- "submissionplot1"


fname1<- here::here("data","test_samp3_10.Rdata")
fname2<- here::here("data","test_samp11_20.Rdata")
fname3<- here::here("data","test_samp25_40.Rdata")
fname4<- here::here("data","test_samp45_50.Rdata")
# fnamevec<-c(fname1,fname2,fname3,fname4)
fnamevec<-c(fname1,fname2,fname3,fname4)

## subset elements from list::
if(FirstRun) {
  samp.ltmp<-future_map(fnamevec,.f=function(fname){load_sublist(fname,element = localindex)})
  samp.comb<-as.list(flatten(samp.ltmp))
  rm(samp.ltmp)
  saveRDS(samp.comb,file = fnamelocal2)
}

fnamelocal2 <-here::here("data","local1.Rdata")

samp.comb<-readRDS(file = fnamelocal2)

samp.l.local<-(lapply(samp.comb,
                       FUN = function(samp){
                         unpack(as.list(samp))
                         df<-(list(repro.df,DCIs,lsso.all) ## lsso.all has lss, lss0, lss.s
                              %>% reduce(full_join,by="sections") ## merging 3 dfs with a common col
                              %>% mutate(repro.s=repro/na_if(reproNull,0)) ## include the scaled repro
                         )
                         # print(df)
                         return(df)
                       })
                %>% bind_rows(.id = "id") ## bind the list into a df
                ## %>% select(-c(reproNull)) ## take out unscaled indices (Ali changed NOV24)
                %>% mutate(DCI_s=DCI_as/100,repro=(repro)) ## scale the DCI_s by 100
                # %>% separate(col = id,into = c("Topology","Node","Sample"),sep="_")
                %>% separate(col = id,into = c("Topology","Dispersal","Direction","Node","Sample"),sep="_")
                %>% select(-"Sample")
                # %>% unite(col = "id", c(Topology,Node),sep="_",remove = TRUE )
                %>% unite(col = "id", c(Topology,Node,Direction),sep="_",remove = TRUE )
                # %>% mutate(Dispersal=factor("Long"))
)

anyNA(samp.l.local)

## correlation data calculated for plotting:: 
cordat_loc <- (samp.l.local
               %>% select(-c(reproNull,lss0,DCI_as)) ## take out unscaled indices
               %>% pivot_longer(cols = c("DCI_s"),names_to="x_var", values_to="x_val")
               %>% pivot_longer(cols = c("repro","lss","repro.s","lss.s"), names_to="y_var", values_to="y_val") ## "lss" was added 28 Nov 2022.  
               %>% nest(data=c(x_val, y_val)) ##you get one row for each group defined by the non-nested columns
               %>% mutate(cor_test = map(data, ~cor.test(.x$x_val, .x$y_val, method = "spearman",exact = FALSE)),
                          tidied = map(cor_test, tidy))
               %>% unnest(tidied)
               %>% select(c(id,sections,Dispersal,x_var,y_var,estimate,p.value)) ## deselect to see all details of correlations
               %>% data.frame()
               %>% unite('Comparison',x_var,y_var,sep = "_",remove = FALSE)
               # %>% separate(col=id, into=c("Topology","Node"), sep="_") ##Node is number of nodes in the system
               %>% separate(col=id, into=c("Topology","Node","Direction"), sep="_") ##Node is number of nodes in the system
)

## Labels for the plots

# little p and r for the local, capitals are for global indices
rsTag<- paste0("(", "$DCI_s$", ", Reach reproductive value ", "$\\mathcal{V}_s(i)$)")
psTag<- paste0("(", "$DCI_s$", ", Reach occupancy ", "$\\mathcal{O}_s(i)$)")

labvec <- c(DCI_s_lss = "$DCI_s$,$lsso$",
            DCI_s_repro = "$DCI_s$,$r$" ## r:reproduction value
            , DCI_s_repro.s= rsTag
            # , DCI_s_lss.s= "$DCI_s$,$lsso_s$" ## used before 23 Apr, 2025
            , DCI_s_lss.s= psTag
)

labvecDir<-c(Sym="Symmetric", Asym="Asymmetric")
labvecTop <-c(linear ="Linear",binary="Binary")
labvecDis<-c(Long="Global",Short="Local")


## Change the factors' name
# cordat_loc <- mutate(cordat_loc, across(Comparison, ~ labvec[.]))


nodeRange.tmp <- unique(cordat_loc[["Node"]])
## old Network size= n for all facets
levNodes<-paste0("Network size=",sort(as.numeric(nodeRange.tmp))) ## levels of id to sorting the plots


## The palette with black:
#> ref: https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
## check

## reorded to match the old plots in the manuscript
custom2_colorBlindBlack8  <- c( "#D55E00", "#000000", "#009E73", "#0072B2","#E69F00", "#F0E442",  "#CC79A7","#56B4E9") ## used before May 2025

## If lighter colors is needed for the local indicators compared to the global indicators, do this:
# custom2_colorBlindBlack8  <- c( "#F7944D", # (from ColorBrewer’s light orange range) old was "#D55E00"
#                                 "#7F7F7F", #(a soft, medium-light gray) old was"#000000", 
#                                 "#009E73", "#0072B2","#E69F00", "#F0E442",  "#CC79A7","#56B4E9")



#> colorblindness: read https://ggplot2.tidyverse.org/reference/scale_manual.html
#>  see also this paper: https://arxiv.org/abs/1903.06490 refed in the link above. Note scale_<aesthetic>_<datatype>_<colorscale>()
pie(rep(1, length(custom2_colorBlindBlack8)), col = custom2_colorBlindBlack8)



## New Network size for only the first facet
selectedNodes<-c(10,30,50)
# sortednoderange<-sort(as.numeric(selectedNodes))
# levNodes<-c(paste0("Network size=",sortednoderange[1]),sortednoderange[-1]) ## levels of id to sorting the plots

## plot all for the Appendix 
pl3_localcor1 <- (cordat_loc
                  %>% filter(Node %in% selectedNodes) ## Directed=="Sym"
                  %>% mutate(Segment=as.numeric(gsub('_s|sink', '',sections))) ## for ordering legend create Segment column where sink is going to be NA at first and then replaced by 01
                  %>% mutate(Segment=as.factor(replace_na(Segment,01)))
                  ## change the labels:
                  %>% mutate(across(Direction, ~labvecDir[.] ))
                  %>% mutate(across(Topology, ~labvecTop[.] ))
                  %>% mutate(across(Dispersal, ~labvecDis[.] ))
                  %>% mutate(across(Comparison, ~ labvec[.]))
                  ## sorting out the order of visualization
                  %>% unite(col="TCDD",c("Topology","Comparison","Dispersal","Direction"),sep = ": ",remove = FALSE)
                  %>% unite(col="TC",c("Topology","Comparison"),sep = ":",remove = FALSE)
                  %>% unite(col="TD",c("Topology","Direction"),sep = ":",remove = FALSE)
                  %>% mutate(TCDD=as.factor(TCDD),TC=as.factor(TC))
                  ## You can change the order of levels here
                  %>% mutate(TD=factor(TD, levels=c("Binary:Symmetric","Linear:Symmetric","Binary:Asymmetric","Linear:Asymmetric")))
                  %>% filter(y_var%in%c("repro.s","lss.s")) ## scaled local indices Topology=="linear",
                  ## create facet labels
                  %>% mutate(Nodelabel=factor(paste0("Network size=",Node),levels = levNodes) )
                  %>% ggplot(aes(x=as.integer(Segment),y=estimate))
                  + theme_bw()
                  + geom_point(aes(colour=Comparison, shape=Dispersal),size=2,alpha=1)
                  + geom_line(aes(color=Comparison,linetype=Dispersal),linewidth=1,alpha=1) # linetype=TopCom
                  # + facet_wrap(~ Nodelabel) ## make sure this matches the label for node
                  # + facet_grid(Direction ~ Nodelabel) ## make sure this matches the label for node
                  + facet_grid(TD ~ Nodelabel) ## make sure this matches the label for node
                  + scale_colour_manual(values = custom2_colorBlindBlack8) ## or select custom_colorBlindBlack8
                  # + scale_colour_brewer(palette="Set1")
                  + xlab("Reach number")
                  + ylab( "Correlation (Spearman)") #
                  # + scale_x_continuous(breaks = nodeRange.tmp)
                  # + coord_cartesian(xlim = c(-2, max(selectedNodes))) ## add extra space to the x-axis cut off for the annotations
                  + theme(panel.spacing=grid::unit(0,"lines"), 
                          axis.text.x = element_text(angle = 90, size=10, vjust = 0.5, hjust=1),
                          axis.text.y = element_text(size=10)
                          )
                  ## Legend:
                  # + theme(legend.position="bottom")
                  # + theme(legend.title=element_blank())
                  + theme(legend.position = "bottom",
                          legend.box = "horizontal",
                          legend.box.just = "center",
                          axis.title.x = element_text(size = 12),
                          axis.title.y = element_text(size = 12),
                          legend.text = element_text(size = 12),     # <- adjust text size
                          legend.title = element_text(size = 12)  # <- adjust title size) 
                          ) +
                    guides(colour = guide_legend(ncol = 1, title = "Correlation:", order = 1),
                           shape = guide_legend(ncol = 1, title = "Dispersal:", order = 2, override.aes = list(colour = "#7F7F7F")),  # gray for points
                           linetype = guide_legend(ncol = 1, title = "Dispersal:", order = 2, override.aes = list(colour = "#7F7F7F"))  # gray for lines
                           )
                  ) 

print(pl3_localcor1)

## Linux
# options(tikzLatex = 'C:/Users/Ali.Gharouni/AppData/Local/Programs/MiKTeX/miktex/bin/x64/pdflatex.exe')

## Windows:
# options(tikzLatex = 'C:/Users/Ali.Gharouni/AppData/Roaming/TinyTeX/bin/windows/pdflatex.exe')
options(tikzLatex = 'C:/Users/Ali.Gharouni/AppData/Local/Programs/MiKTeX/miktex/bin/x64/pdflatex.exe')
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}","\\usepackage[T1]{fontenc}",
                               "\\usetikzlibrary{calc}", "\\usepackage{amssymb}" ))

# localInd1 <- here::here("plots","draftplots","localInd1.tex")

pll1name <- here::here("plots",targetFolder,"localInd1.tex")
pll1namePDF <- here::here("plots",targetFolder,"localInd1.pdf")

## Adding benchmark hlines to the correlation plots.
# pl3_localcor1<-(pl3_localcor1 + 
#                 # Add horizontal 
#                 geom_hline(yintercept = 0.4, linetype = "solid", color = "gold",show.legend = FALSE) +
#                 geom_hline(yintercept = 0.6, linetype = "solid", color = "orange",show.legend = FALSE) +
#                 geom_hline(yintercept = 0.8, linetype = "solid", color = "red",show.legend = FALSE)
# )

## Final plot
pl3_localcor1F <- (pl3_localcor1 
  #                  + ## Adding annotations to the hlines
  # geom_text(data = data.frame(
  #   x = -4,
  #   y = 0.4,
  #   label = "Moderate (0.4–0.59)",
  #   Nodelabel = "Network size=10"
  #   # Direction = "Symmetric"
  # ),
  # aes(x = x, y = y, label = label),
  # inherit.aes = FALSE,
  # color = "gold", size = 3, hjust = 0, show.legend = FALSE) +
  # 
  # geom_text(data = data.frame(
  #   x = -4,
  #   y = 0.6,
  #   label = "Strong (0.6–0.79)",
  #   Nodelabel = "Network size=10"
  #   # Direction = "Symmetric"
  # ),
  # aes(x = x, y = y, label = label),
  # inherit.aes = FALSE,
  # color = "orange", size = 3, hjust = 0, show.legend = FALSE) +
  # 
  # geom_text(data = data.frame(
  #   x = -4,
  #   y = 0.8,
  #   label = "Very Strong (0.8–1.0)",
  #   Nodelabel = "Network size=10"
  #   # Direction = "Symmetric"
  # ),
  # aes(x = x, y = y, label = label),
  # inherit.aes = FALSE,
  # color = "red", size = 3, hjust = 0, show.legend = FALSE)
)

print(pl3_localcor1F)


# system("pdflatex pl1_globcor.tex")
# system(paste("pdflatex -output-directory=./test",localInd1))

## IN windows:: (getwd() first;)
# getwd()
# tools::texi2pdf(localInd1)



if(saveplots){
  ggsave(pl3_localcor1F
         # + annotate("text", x=0, y=1, label="text\\textsubscript{subscript}")
         # + ggtitle("Long vs Short D")
         ,
         device = tikz,
         filename = pll1name,
         standAlone = TRUE,
         width = 10, height = 7, units = "in")
  
  # system("pdflatex pl1_globcor.tex")
  # system(paste("pdflatex -output-directory=./test",pl1name))
  ## IN windows:: (getwd() first;)
  getwd()
  setwd(here::here("plots",targetFolder)) 
  tools::texi2pdf(pll1name, clean=TRUE)
  setwd(here::here("scripts"))
}






## PLOTS for APENDIX ## #######################################################
## New Network size for only the first facet
selectedNodes<-c(5,10,15,20,25,30,35,40,45,50)
## plot all for the Appendix 
pl3_localcorAll <- (cordat_loc
                    %>% filter(Node %in% selectedNodes) ## Directed=="Sym"
                    %>% mutate(Segment=as.numeric(gsub('_s|sink', '',sections))) ## for ordering legend create Segment column where sink is going to be NA at first and then replaced by 01
                    %>% mutate(Segment=as.factor(replace_na(Segment,01)))
                    ## change the labels:
                    %>% mutate(across(Direction, ~labvecDir[.] ))
                    %>% mutate(across(Topology, ~labvecTop[.] ))
                    %>% unite(col="TCDD",c("Topology","Comparison","Dispersal","Direction"),sep = ": ",remove = FALSE)
                    %>% unite(col="TC",c("Topology","Comparison"),sep = ":",remove = FALSE)
                    %>% unite(col="TD",c("Topology","Direction"),sep = ":",remove = FALSE)
                    %>% mutate(TCDD=as.factor(TCDD),TC=as.factor(TC))
                    ## You can change the order of levels here
                    # %>% mutate(TD=factor(TD, levels=c("Binary:Sym","Linear:Sym","Binary:Asym","Linear:Asym")))
                    %>% mutate(TD=factor(TD, levels=c("Binary:Symmetric","Linear:Symmetric","Binary:Asymmetric","Linear:Asymmetric")))
                    %>% filter(y_var%in%c("repro.s","lss.s")) ## scaled local indices Topology=="linear",
                    ## create facet labels
                    %>% mutate(label=factor(paste0("Network size=",Node),levels = levNodes) )
                    %>% ggplot(aes(x=as.integer(Segment),y=estimate))
                    + theme_bw()
                    + geom_point(aes(colour=Comparison, shape=Dispersal),size=1,alpha=1)
                    + geom_line(aes(color=Comparison,linetype=Dispersal),linewidth=1,alpha=1) # linetype=TopCom
                    ## + facet_wrap(~ label) ## make sure this matches the label for node
                    ## + facet_grid(Direction ~ label) ## make sure this matches the label for node
                    + facet_grid(TD ~ label) ## make sure this matches the label for node
                    + theme(legend.position="bottom")
                    + theme(legend.title=element_blank())
                    + scale_colour_manual(values = colorBlindBlack8) ## or select custom_colorBlindBlack8
                    ## + scale_colour_brewer(palette="Set1")
                    + xlab("Reach number")
                    + ylab( "Correlation (Spearman)") #
                    ## + scale_x_continuous(breaks = nodeRange.tmp)
                    + theme(panel.spacing=grid::unit(0,"lines"))
                    + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
                    )


print(pl3_localcorAll)

localInd1_all <- here::here("plots","draftplots","localInd1_all.tex")


ggsave(pl3_localcorAll
       # + ggtitle("Long vs Short D")
       ,
       device = tikz,
       filename = localInd1_all,
       standAlone = TRUE,
       width = 10, height = 7, units = "in")


## IN windows:: (getwd() first;)
getwd()
tools::texi2pdf(localInd1_all)









