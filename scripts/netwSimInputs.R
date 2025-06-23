## 0. Metapopulation params
# Metapopulation parameters
params2 <-c( c = 10, ##10
             e = 1,
             D = 1, # was 1## FIXME: D is independent for now, can be a randomly selected # in [l_critical, sum of li's]? 
             # habitat power parameters
             gamma = 1,
             omega = 1,
             epsilon=1
             # n_seg=3 ## added May 7, 2023
)

# gamma = 0, Ovaskainan model?
# #############################################################################
# ## 1. make the samples and list them
# 
# ssize <- 4 ## sample size integer old was 60
# nodeRange <- 3 #c(30,40) ## 3,7,15 are compelete binary graphs c(3,5,7,11,15)
# # nodeRange <- c(5,10,15,20) ## 3,7,15 are compelete binary graphs c(3:17)
# treetypeRange <- as.list(c("linear","binary"))