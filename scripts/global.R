## By: Ali
## Created in Feb 2024 (from: ~/projects/connectivity/scripts/old/netwSim_glob_all.R) to combine all correlation plots and comparing symmetric-Asymmetric network, Long-short dispersal

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

## GLOBAL INDICATORS needed for this analysis: 
globindex <- c("DCIp","lambda","lambdaNull","lambdaM","lambdaMNull","R","RNull","ssoi","ssoiNull", "Dispersal","directed") # 
globindex_unscaled <- c("lambda","lambdaM","R","ssoi")

## these files include all info from Long/Short dispersal, Sym/Asym network, linear/binary networks
## these files are created ~/projects/connectivity/tests/test_pipeline20240216.R

FirstRun <- FALSE #TRUE (this prevent to upload all the files accidentally)
saveplots<-TRUE ## ggsave?

fname1<- here::here("data","test_samp3_10.Rdata")
fname2<- here::here("data","test_samp11_20.Rdata")
fname3<- here::here("data","test_samp25_40.Rdata")
fname4<- here::here("data","test_samp45_50.Rdata")
fname5<- here::here("data","test_samp20_30.Rdata")
fnamevec<-c(fname1,fname2,fname3,fname4,fname5)

if(FirstRun) {
  samp.df <- data.frame()
  # samp.df<-future_map(fnamevec,.f=function(fname)savewhatnecassary(samp.lname = fname,indexvec=globindex) ) 
  for(fname in fnamevec){
    print(fname)
    samp.df <-rbind(samp.df,savewhatnecassary(samp.lname = fname,indexvec=globindex ))
  }
  saveRDS(samp.df,file = here::here("data","global.Rdata"))
}

samp.df <-readRDS(here::here("data","global.Rdata"))

## mutate a transformed column if needed Here:
IndGlobComb <- (samp.df
                    %>% separate(col = id,into =  c("Topology","Dispersal","Direction","Node","Sample"),sep="_")
                    %>% select(-c("Sample", "directed"))
                    %>% transform(DCIp=DCIp/100)
                    # %>% unite(col = "id", c(Topology,Node),sep="_",remove = FALSE )
                    ## mutate scaled indices:
                    %>% mutate(lambda.s=lambda/na_if(lambdaNull,0),
                               lambdaM.s=lambdaM/na_if(lambdaMNull,0),
                               ssoi.s=ssoi/na_if(ssoiNull,0),
                               R.s=R/na_if(RNull,0),
                               # directed=factor(ifelse(directed, "Asym","Sym")),
                               growth.s=exp(lambda-lambdaNull) ## This is annual growth scaled by the Null network
                    )
                    ## mutate ordered factors:
                    %>% mutate(Node=factor(Node, stringr::str_sort(unique(Node), numeric = TRUE)))
                    # %>% mutate(across(c(Dispersal,Topology), as.factor))
)
## IF .tmp is ok, use it:
anyNA(IndGlobComb)

## Check the levels of network size to ensure the order for later on use in the multi var model
levels(IndGlobComb[["Node"]])

test <-(IndGlobComb %>% filter(Node==23))

## mutual correlation between global indices ##################################


#> Some Notes about Spearman cor: A rank-correlation may be used to pick up monotonic association between variates
#> see https://stats.stackexchange.com/questions/132830/is-it-okay-to-plot-a-regression-line-for-ranked-data-spearman-correlation
#> terminology: low/high Spearman correlation show low/high monotonic association between x and y.

head(IndGlobComb,2)
IndGlobComb.stats <- (IndGlobComb 
                      ## Picking the favorite variables
                      # %>% select(-"Sample")
                      %>% select(c("Topology","Node","Dispersal","Direction",
                                   "DCIp","lambda.s","ssoi.s","R.s",
                                   "lambdaM","lambdaM.s","growth.s")) 
                      %>% drop_na() ## remove all NA observations with drop_na() ## FIXEME:: not sure if this is the best way
                      %>% pivot_longer(cols = c("DCIp"),names_to="x_var", values_to="x_val")
                      %>% pivot_longer(cols = c("lambda.s","ssoi.s","R.s","lambdaM","lambdaM.s","growth.s"),
                                       names_to="y_var", values_to="y_val")
                      %>% nest(data=c(x_val, y_val)) ##you get one row for each group defined by the non-nested columns
                      %>% mutate(cor_test.Sp = map(data, ~cor.test(.x$x_val, .x$y_val, method = "spearman",exact = FALSE)), ## for ties (see https://stackoverflow.com/questions/10711395/spearman-correlation-and-ties)
                                 Spearman = map(cor_test.Sp, tidy))
                      %>% mutate(cor_test.Pe = map(data, ~cor.test(.x$x_val, .x$y_val, method = "pearson",exact = FALSE)), ## for ties (see https://stackoverflow.com/questions/10711395/spearman-correlation-and-ties)
                                 Pearson = map(cor_test.Pe, tidy))
                      %>% unnest(c(Spearman,Pearson), names_sep = "_")
                      # %>% select(c(id,x_var,y_var,estimate,p.value)) ## deselect to see all details of correlations
                      %>% select(c(Topology,Node,Dispersal,Direction,x_var,y_var,Spearman_estimate,Pearson_estimate)) ## deselect to see all details of correlations
                      # %>% data.frame()
                      %>% rename("Spearman"="Spearman_estimate","Pearson"="Pearson_estimate")
                      %>% unite('Correlation',x_var,y_var,sep = "_",remove = FALSE)
                      # %>% separate(col=id, into=c("Topology","Node"), sep="_")
                      ## %>% rename("Correlation"="estimate")
                      %>% pivot_longer(cols = c("Spearman","Pearson"), names_to = "Method", values_to = "estimate")
                      %>% mutate(Method=factor(Method)) #Node=as.integer(Node)
                      %>% filter(Method=="Spearman")
                      # %>% arrange(Topology,Node,Dispersal,Direction)
)
head(IndGlobComb.stats,2)

levels(IndGlobComb.stats[["Node"]])
sort(unique(as.numeric(as.character(IndGlobComb.stats[["Node"]]))))

unique(IndGlobComb.stats[["Correlation"]])

# desired_order <- c("DCIp_growth.s", "DCIp_R.s", "DCIp_ssoi.s", 
#                    "DCIp_lambdaM", "DCIp_lambdaM.s", "DCIp_lambda.s")
# 
# # Step 2: Relevel Correlation (before applying labvec)
# IndGlobComb.stats <- IndGlobComb.stats %>%  mutate(Correlation = factor(Correlation, levels = desired_order))

# ##############################################################################
# Plotting Part
# ##############################################################################
## 1. mutual correlation ###########
## Labels for the plots

gsTag <- paste0("(", "$DCI_p$", ", Metapopulation growth rate ", "$\\mathcal{G}_s$)")
RsTag<- paste0("(", "$DCI_p$", ", Basic reproduction number ", "$\\mathcal{R}_s$)")
PsTag<- paste0("(", "$DCI_p$", ", Steady state occupancy ", "$\\mathcal{P}_s$)")


labvec <- c(DCIp_lambda.s = "$(DCI_p,\\lambda_s)$",
            DCIp_ssoi.s = PsTag, # DCIp_ssoi.s = "$(DCI_p,\\mathcal{P}_s)$",
            DCIp_R.s = RsTag, # DCIp_R.s ="$(DCI_p,\\mathcal{R}_s)$",
            DCIp_lambdaM="$(DCI_p,\\lambda_M)$",
            DCIp_lambdaM.s="$(DCI_p,\\lambda_{Ms})$",
            DCIp_growth.s= gsTag ##"$(DCI_p,\\mathcal{G}_s)$"
)

desired_order <- c(gsTag, RsTag, PsTag, 
                   "$(DCI_p,\\lambda_s)$", 
                   "$(DCI_p,\\lambda_M)$", 
                   "$(DCI_p,\\lambda_{Ms})$"
                   )


labvecDir<-c(Sym="Symmetric", Asym="Asymmetric")
labvecTop <-c(linear ="Linear",binary="Binary")
labvecDis<-c(Long="Global",Short="Local")

## Change the factors' name
# IndGlobComb.stats <- mutate(IndGlobComb.stats, across(Correlation, ~ labvec[.]))
str(head(IndGlobComb.stats,2))

## The palette with black:
#> ref: https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
## check
# pie(rep(1, length(colorBlindBlack8)), col = colorBlindBlack8)

## reorded to match the old plots in the manuscript
custom_colorBlindBlack8  <- c( "#009E73", "#000000","#D55E00","#0072B2","#E69F00", "#F0E442",  "#CC79A7","#56B4E9")
#> colorblindness: read https://ggplot2.tidyverse.org/reference/scale_manual.html
#>  see also this paper: https://arxiv.org/abs/1903.06490 refed in the link above. Note scale_<aesthetic>_<datatype>_<colorscale>()

## plot all nodes:
selectedNodes<- sort(unique(as.numeric(as.character(IndGlobComb.stats[["Node"]]))))
# selectedNodes<- c(5,10,15,20,25,30,35,40,45,50) ## all: sort(unique(as.numeric(as.character(IndGlobComb.stats[["Node"]]))))
pl1_globcor <- (IndGlobComb.stats
                # %>% filter(!y_var%in%c("lambdaM","lambdaM.s","lambda.s") 
                #            # & directed=="Sym" 
                #            & Node %in% selectedNodes)
                %>% filter(y_var%in%c("ssoi.s","R.s","growth.s") & Node %in% selectedNodes)
                ## change the labels:
                %>% mutate(across(Correlation, ~ labvec[.]))
                %>% mutate(Correlation = factor(Correlation, levels = desired_order)) ## Reordering the labels
                %>% mutate(across(Direction, ~labvecDir[.] ))
                %>% mutate(across(Topology, ~labvecTop[.] ))
                %>% mutate(across(Dispersal, ~labvecDis[.] ))
                ## sorting out the order of visualization
                %>% mutate(Direction=factor(Direction, levels=c("Symmetric","Asymmetric")) )
                %>% ggplot(aes(x=as.numeric(as.character(Node)),y=estimate,group=interaction(Topology,Correlation,Dispersal)))
                + theme_bw()
                # + facet_grid(directed~Topology) ## In case one wants to compare sym with Asym
                + facet_grid(Direction~Topology) ## TODO:: put sym on top row
                + geom_point(aes(colour=Correlation, shape=Dispersal),size = 2,alpha = 1)
                + scale_shape_manual(values = c(16, 17)) ## Dispersal=Long, Short shape poiunts see https://www.datanovia.com/en/blog/ggplot-point-shapes-best-tips/
                + geom_line(aes(x=as.numeric(as.character(Node)),y=estimate,colour=Correlation,linetype=Dispersal),
                            linewidth=1.,alpha=1)
                + labs(x=("Network size"),
                       y=("Correlation (Spearman)"))
                # + scale_x_continuous(breaks=sort(unique(as.numeric(as.character(IndGlobComb.stats[["Node"]]))))) ## instead of nodeRange from "./netwSimInputs.R"
                # + scale_x_continuous(breaks=sort(selectedNodes) ) ## instead of nodeRange from "./netwSimInputs.R"
                + scale_x_continuous(breaks = c(3,5,10,15,20,25,30,35,40,45,50)) ## Only labels these on x-axis
                # + coord_cartesian(xlim = c(-5, max(selectedNodes))) ## add extra space to the x-axis cut off for the annotations
                + ylim(-0.1,1)
                + scale_colour_manual(values = custom_colorBlindBlack8)
                # + scale_colour_brewer(palette="Dark2")
                + theme(panel.spacing=grid::unit(0,"lines"), 
                        axis.text.x = element_text(angle = 90, size=10, vjust = 0.5, hjust=1),
                        axis.text.y = element_text(size=10)
                )
                ## Legend :
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
                  guides(
                    colour = guide_legend(ncol = 1, title = "Correlation:", order = 1),
                    shape = guide_legend(ncol = 1, title = "Dispersal:", order = 2, 
                                         override.aes = list(colour = "#7F7F7F")),  # gray for points
                    linetype = guide_legend(ncol = 1, title = "Dispersal:", order = 2, 
                                            override.aes = list(colour = "#7F7F7F"))  # gray for lines
                  )
                )

print(pl1_globcor)

## Get tikz working ###########################################################
# tinytex::install_tinytex()
## install Miktex packages, set the global options as follows:
#> Rstudio settings: 
#>                  sweave: sweave, pdflatex, use Tinytex
## add a package to the defaults see Ben's comment here: https://stackoverflow.com/questions/12514612/how-to-annotate-ggplot-with-latex
## in cmd type: H:\>where pdflatex and put the path to the options down.

# options("tikzLatex"='/home/ag/bin/pdflatex')
# options("C:/Users/Ali.Gharouni/AppData/Local/Programs/MiKTeX")
# options(tikzDefaultEngine = 'pdftex') 

## Linux
# options(tikzLatex = 'C:/Users/Ali.Gharouni/AppData/Local/Programs/MiKTeX/miktex/bin/x64/pdflatex.exe')
## Windows:
# options(tikzLatex = 'C:/Users/Ali.Gharouni/AppData/Roaming/TinyTeX/bin/windows/pdflatex.exe')
options(tikzLatex = 'C:/Users/Ali.Gharouni/AppData/Local/Programs/MiKTeX/miktex/bin/x64/pdflatex.exe')

options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}","\\usepackage[T1]{fontenc}",
                               "\\usetikzlibrary{calc}", "\\usepackage{amssymb}" ))

# tikzDevice::tikzTest()
# tikzDevice::tikzTest("\\char77")

pl1name <- here::here("plots","submissionplot1","globalInd1_all.tex")
pl1namePDF <- here::here("plots","submissionplot1","globalInd1_all.pdf")


# pl1_globcor<-(pl1_globcor + 
#                 # Add horizontal 
#                 geom_hline(yintercept = 0.4, linetype = "solid", color = "gold",show.legend = FALSE) +
#                 geom_hline(yintercept = 0.6, linetype = "solid", color = "orange",show.legend = FALSE) +
#                 geom_hline(yintercept = 0.8, linetype = "solid", color = "red",show.legend = FALSE)
#                 )

## Final plot
pl1_globcorF <- (pl1_globcor 
                 # + geom_text(data = data.frame(
                 #   x = -4,
                 #    y = 0.4,
                 #    label = "Moderate (0.4–0.59)",
                 #    Topology = "Binary"
                 #    # Direction = "Symmetric"
                 #  ),
                 #  aes(x = x, y = y, label = label),
                 #  inherit.aes = FALSE,
                 #  color = "gold", size = 3, hjust = 0, show.legend = FALSE) +
                 #  
                 #  geom_text(data = data.frame(
                 #    x = -4,
                 #    y = 0.6,
                 #    label = "Strong (0.6–0.79)",
                 #    Topology = "Binary"
                 #    # Direction = "Symmetric"
                 #  ),
                 #  aes(x = x, y = y, label = label),
                 #  inherit.aes = FALSE,
                 #  color = "orange", size = 3, hjust = 0, show.legend = FALSE) +
                 #  
                 #  geom_text(data = data.frame(
                 #    x = -4,
                 #    y = 0.8,
                 #    label = "Very Strong (0.8–1.0)",
                 #    Topology = "Binary"
                 #    # Direction = "Symmetric"
                 #  ),
                 #  aes(x = x, y = y, label = label),
                 #  inherit.aes = FALSE,
                 #  color = "red", size = 3, hjust = 0, show.legend = FALSE)
                 )
# pl1_globcorF <- (pl1_globcor +
#   geom_text(data = data.frame(
#     x = -4,  #-4
#     y = c(0.4, 0.6, 0.8),
#     label = c("Moderate (0.4–0.59)", "Strong (0.6–0.79)", "Very Strong (0.8–1.0)"),
#     # Direction = "Symmetric",       # Or "Asymmetric", depending on which row
#     Topology = "Binary",            # The column you want it in
#     color = c("yellow","orange", "red")
#   ),
#   aes(x = x, y = y, label = label, color=color),
#   inherit.aes = FALSE,
#   show.legend = FALSE,
#   size = 3, 
#   hjust = 0))

# library(grid)
# 
# # Define the text grob
# outside_label <- grid::textGrob(
#   "Very Strong (0.8,1)", 
#   x = unit(0, "npc") + unit(0.1, "lines"),  # Horizontal: move it to the left margin
#   y = unit(0.8, "native"),                 # Vertical: in data coordinates
#   just = "left",
#   gp = gpar(col = "red", fontsize = 9)
# )
# 
# # Add it to the plot using annotation_custom
# pl1_globcorF <- pl1_globcor + 
#   annotation_custom(grob = outside_label,
#                     xmin = -Inf, xmax = -Inf,  # ensures it's treated as "global" on x-axis
#                     ymin = -Inf, ymax = Inf)   # still constrained vertically
# 
print(pl1_globcorF)
# dev.off()

if(saveplots){
  ggsave(pl1_globcorF
         # + annotate("text", x=0, y=1, label="text\\textsubscript{subscript}")
         # + ggtitle("Long vs Short D")
         ,
         device = tikz,
         filename = pl1name,
         standAlone = TRUE,
         width = 10, height = 7, units = "in")

  # system("pdflatex pl1_globcor.tex")
  # system(paste("pdflatex -output-directory=./test",pl1name))
  ## IN windows:: (getwd() first;)
  getwd()
  setwd(here::here("plots","submissionplot1")) 
  tools::texi2pdf(pl1name, clean=TRUE)
  setwd(here::here("scripts"))
  }


# ## 2nd option
# library(knitr)
# knitr::knit2pdf(input=pl1name ,output=pl1namePDF)

## 3rd option (didn't work)
# system("pandoc -f latex C:/Ali/projects/connectivity/plots/draftplots/globalInd1_all.tex -s -o C:/Ali/projects/connectivity/plots/draftplots/globalInd1_all.pdf")


## TESTS and Others stuff ######################################################
## TESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTT
## TESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTTESTT


## Check: (20240607) ####
## If the R_s and lambda_Ms are the same? Note that the spearman correlation of these indicators and the DCI_P seems identical.


(IndGlobComb
 %>% select(c(Topology, Dispersal, Direction,Node,lambdaM.s,R.s,lambdaM,R))
 %>% mutate(Delta=lambdaM-R,Delta.s=lambdaM.s-R.s)
 )->tmp

any(tmp$Delta<0.0)

(tmp %>% ggplot(aes(x=Node,y=Delta))
  + geom_boxplot(aes(color=Dispersal))
  +facet_grid(Topology~Direction)
  )

### TEST (20240315): ################
## Is DCIp higher or lowr in linear vs binary networks across groups?
## Ans: average DCIp of binary is higher than linear networks, specifically for network size > 4

pltest_dcip<- (IndGlobComb
               # %>% filter(Node%in%c("5"))
               %>% select(Topology,Dispersal,Direction,Node,DCIp)
               %>% filter(Direction%in%c("Sym"))
               ## DCIp is the same for long/short dispersal
               %>% ggplot(aes(x=as.numeric(as.character(Node)),y=DCIp,group=interaction(Node,Topology,Direction) ))
               + geom_boxplot(aes(colour=Topology))
               # + facet_grid(Direction~.)
               + theme_bw()
               + labs(x=("Network size"))
               )
print(pltest_dcip)

if(saveplots){
  pltname_dcip <- here::here("plots","suplots","boxp_dcip.tex")
  ggsave(pltest_dcip,
       device = tikz,
       filename = pltname_dcip,
       standAlone = TRUE,
       width = 8, height = 8, units = "in")
  ## IN windows:: (getwd() first;)
  getwd()
  setwd(here::here("plots","suplots"))
  tools::texi2pdf(pltname_dcip, clean=TRUE)
  }

## Check the differences between the medians as a function of Network size: ####

summedInd <-(IndGlobComb  
             %>% select(c("Topology","Node","Dispersal","Direction",
               "DCIp","ssoi.s","R.s","growth.s")) ## ,"lambda.s","ssoi.s","R.s","lambdaM","lambdaM.s","growth.s"
            %>% drop_na() ## remove all NA observations with drop_na() ## FIXEME:: not sure if this is the best way
            %>% pivot_longer(cols = c("DCIp","ssoi.s","R.s","growth.s"),names_to="x_var", values_to="x_val")
            %>% group_by(Topology,Node,Dispersal,Direction,x_var)
            %>% summarise(
              count = n(),                             # Count observations in each group
              mean_value = mean(x_val),                # Calculate mean value
              median_value = median(x_val),            # Calculate median value
              sd_value = sd(x_val)
              )
            %>% ungroup()
            )

## To check whether the grouping and summarising above is working::
# temp.test<-(IndGlobComb 
#             %>% filter(Topology=="binary", Node==3, Dispersal=="Long", Direction=="Asym")
#             %>% select(DCIp)
#             )
# median(c(temp.test$DCIp))

test<-(summedInd 
  %>% filter(x_var%in%c("DCIp"))
  %>% group_by(Topology,Dispersal,Direction,x_var)
  %>% summarise(df = diff(median_value)
  )
  %>% ungroup()
  # %>% ggplot(aes(x=as.numeric(as.character(Node)),y=median_value, group=x_var) ## c(abs(diff(median_value)),NA)
  #            # ,group=c(Topology,Dispersal,Direction,x_var)
  #            )
  # + geom_point(aes(color=x_var))
  # + geom_line(aes(y=c(0,diff(median_value))))
  # + facet_grid(Direction~interaction(Topology,Dispersal))
  )


2*2*2*length(selectedNodes)


diff(c(0,2,1))










## comparison of ave metrics across groups::
metrics<-c("DCIp","ssoi.s","R.s","growth.s")

# pltest_dcip2<- (IndGlobComb
#                 %>% filter(Node%in%c("3"))
#                 %>% select(all_of(c("Topology","Dispersal","Direction","Node",metrics)))
#                 %>% pivot_longer(cols=metrics, names_to = "Metric", values_to = "value")
#                 # %>% filter(Direction=="Sym")
#                 %>% ggplot(aes(x=as.numeric(as.character(Node)),y=value,group=interaction(Node,Topology,Dispersal,Direction,Metric) ))
#                 + geom_boxplot_jitter(aes(colour=interaction(Dispersal,Metric)))
#                 + facet_grid(Direction~Topology)
#                 + theme_bw()
#                 + labs(x=("Network size"))
# )
# print(pltest_dcip2)


#> Q: What can be inferred from the plot of the distribution (boxplot) of the indicators -- DCIp and metapops -- 
#> grouped by "Topology","Dispersal","Direction","Node"?
#> Ans: DCIp is decreasing ~ 1/L^2 (to be verified) and networks with size > 20 are relatively more similar than
#> the smaller networks. This claim can be visualized by plotting the differences between the node-specific medians 
#> of the consecutive distributions (boxplots). 

#> Self note: A more statistical approach would require a glm multivariate model factored by 
#>  ("Topology","Dispersal","Direction","Node"). If after certain node number (which we visually picked 20), 
#>  the differences were insignificant then we could pick that node number as the cutoff of plots. 
#>  [See below the modelling]
#>  
#>  This similarity (insignificant variability) in indicators for large networks can be explained as follows:
#> larger networks are more alike in terms of passabilities and reach lengths than small-sized networks.  


## replace labels by::
labvec2 <- c(lambda.s = "$\\lambda_s$",  ssoi.s = "$\\mathcal{P}_s$",  R.s ="$\\mathcal{R}_s$",
             growth.s= "$\\mathcal{G}_s$" ,DCIp="$DCI_p$")

pltest_metrics<- (IndGlobComb
                  # %>% filter(Node%in%c("5","10"))
                  %>% select(all_of(c("Topology","Dispersal","Direction","Node",metrics)))
                  %>% pivot_longer(cols=metrics, names_to = "Metric", values_to = "value")
                  %>% mutate(across(Metric, ~labvec2[.]))
                  %>% ggplot2::ggplot(aes(x=as.numeric(as.character(Node)),y=value,
                                 group=interaction(Node,Topology,Dispersal,Direction,Metric) ))
                  + geom_boxplot_jitter(aes(colour=Dispersal))
                  + facet_grid(interaction(Direction,Topology)~Metric)
                  # + facet_wrap(interaction(Direction,Topology)~Metric, scale = "free_y")
                  # + scale_y_continuous(trans='log10')
                  + theme_bw()
                  + labs(x=("Network size"), y="Metrics")
)
print(pltest_metrics)

pltname_metrics <- here::here("plots","suplots","boxp_glob.tex")
# pltname_metrics <- here::here("plots","suplots","boxp_globfreey.tex")
ggsave(pltest_metrics,
       device = tikz,
       filename = pltname_metrics,
       standAlone = TRUE,
       width = 10, height = 10, units = "in")

## IN windows:: (getwd() first;)
getwd()
setwd(here::here("plots","suplots"))
tools::texi2pdf(pltname_metrics, clean=TRUE)



### TEST 2 ############
yfacets2 <- c("R","lambda","ssoi")
yfacets2.s <- c("R.s","lambda.s","ssoi.s","growth.s") ## scaled ind

colnames(IndGlobComb)

df3.tmp<-(IndGlobComb
          # %>% filter(id=="linear_3") ## just checking:
          %>% unite(col="id",c("Topology","Dispersal","Direction","Node" ),sep = ":",remove = FALSE)
          # %>% select("id",all_of(yfacets2),all_of(yfacets2.s),"DCIp",Topology,Node)
          %>% select((c("id",all_of(yfacets2.s),"DCIp","Topology","Node")))
          # %>% transform(ssoi=ssoi/100)
          # %>% separate(col=id, into=c("Topology","Node"), sep="_")
          # %>% mutate(Topology=factor(Topology), Node=factor(Node))
          %>% pivot_longer(cols = c(yfacets2.s), names_to = "y", values_to = "value")
)
## replace labels by::
labvec3 <- c(
  lambda = "$\\lambda$",  ssoi = "$SSO$",  R ="$\\mathcal{R}$",
  lambda.s = "$\\lambda_s$",  ssoi.s = "$\\mathcal{P}_s$",  R.s ="$\\mathcal{R}_s$",
  growth.s= "$\\mathcal{G}_s$" )

## Linux
# options(tikzLatex = 'C:/Users/Ali.Gharouni/AppData/Local/Programs/MiKTeX/miktex/bin/x64/pdflatex.exe')
## Windows:
options(tikzLatex = 'C:/Users/Ali.Gharouni/AppData/Roaming/TinyTeX/bin/windows/pdflatex.exe')
options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}","\\usepackage[T1]{fontenc}",
                               "\\usetikzlibrary{calc}", "\\usepackage{amssymb}" ))

## Change the factors' name
IndGlobComb.stats <- mutate(IndGlobComb.stats, across(Correlation, ~ labvec[.]))

# ## Change the factors' name
# df3.tmp <- mutate(df3.tmp, across(y, ~ labvec3[.]))

## selected Network size thta you would like the scatter plots of:
selected_netsize <- 10 #c(5,10,20,25,30,35,40,45,50)

#> plotted the regression line for rank(x) and rank(y) to visualize the monotonic association between DCIp and the meta population indicatores.
#> see explanation here: https://stats.stackexchange.com/questions/132830/is-it-okay-to-plot-a-regression-line-for-ranked-data-spearman-correlation

pl3.temp <- (df3.tmp
             %>% filter(Node %in% selected_netsize)
             %>% mutate(across(y, ~labvec3[.]))
             %>% ggplot(aes(x = (DCIp), y = (value), group = y))
             ## Why it doesn't make sense to draw 1:1 line in Spearman cor case? see https://stats.stackexchange.com/questions/132830/is-it-okay-to-plot-a-regression-line-for-ranked-data-spearman-correlation
             # + geom_abline(slope = 1, intercept = 0, alpha = 0.7, lty = 2)
             + geom_point(aes(color=id),size=1)
             + stat_smooth(aes(x = (DCIp), y = (value)),method = lm, formula = y ~ poly(x, 2),
                           se = T, alpha=0.5, linewidth=.5)
             + facet_wrap(y~id, scale = "free_y")
             # + facet_grid(y~factor(Node, stringr::str_sort(unique(Node), numeric = TRUE))) ## Sort factors out
             + labs(x=("DCIp"),
                    y=("Metapopulation Indices"))
             + theme(legend.position="bottom")
             + scale_colour_brewer(palette="Dark2")
             + theme(panel.spacing=grid::unit(0,"lines"))
)
print(pl3.temp)

plotname <- here::here("plots","suplots",paste0("scat_glob",selected_netsize,".tex"))
ggsave(pl3.temp
       # + ggtitle("Long vs Short D")
       ,
       device = tikz,
       filename = plotname,
       standAlone = TRUE,
       width = 12, height = 10, units = "in")

## IN windows:: (getwd() first;)
getwd()
setwd(here::here("plots","suplots"))
tools::texi2pdf(plotname, clean=TRUE)
## back to scripts
setwd(here::here("scripts"))







df3test<-(df3.tmp 
          %>% filter(Node %in% selected_netsize, 
                     str_detect(df3.tmp$id, pattern = paste("linear","Long","Asym", sep = ": ")),
                     y %in% c("ssoi.s")
          )
)

# saveRDS(df3test, file="C:/Ali/projects/connectivity/data/linearLongAsym40.Rdata")

## Are they al diff?
length(unique(df3test$value))

str(df3test)

(df3test 
  %>% ggplot(aes(x=DCIp, y=value))
  + geom_point(aes(color=y))
)

cor.test(df3test$DCIp, df3test$value, method = "spearman",exact = FALSE)

q90x<-quantile(df3test$DCIp,0.9)
q90y<-quantile(df3test$value,0.9)

# df3testclean <- df3test[df3test$DCIp<=q90x & df3test$value<=q90y , ]
df3testclean <- df3test[ df3test$value<=q90y , ]

cor.test(df3testclean$DCIp, df3testclean$value,method = "spearman",exact = FALSE)


## TEST:: Statistical Modelling ##################################################  

#> Self note: A more statistical approach would require a glm multivariate model factored by 
#>  ("Topology","Dispersal","Direction","Node"). If after certain node number (which we visually picked 20), 
#>  the differences were insignificant then we could pick that node number as the cutoff of plots. 
#>  [See below the modelling]

str(head(IndGlobComb,5))

sort(unique(as.numeric(as.character(IndGlobComb[["Node"]]))))

library(nlme)
library(mgcv)

## prep data for modelling
dm <-(IndGlobComb
      %>% mutate(across(c("Topology","Dispersal","Direction"), ~ as.factor(.) ),
                 Node.n=as.numeric(as.character(Node))
      )
)

str(head(dm,5))
(levels(dm$Node))

glm.m0 <- glm(DCIp ~ Node, 
              data = dm,
              weights = as.numeric(as.character(Node)),
              family= binomial
)


glm.m1 <- glm(DCIp ~ Topology*Dispersal*Direction*Node, 
              data = dm,
              # weights = as.numeric(as.character(Node)),
              family= binomial
)

glmmodel<-glm.m0
summary(glmmodel)
plot.lme(glmmodel)





m1 <- gam(DCIp ~ s(Node.n, by=interaction(Topology,Dispersal,Direction)), 
          data = dm,
          method="REML",
          # family= binomial
)


summary(m1)





