## SET UP ----

## load libraries 
library(tidyverse)
library(cowplot)
library(jpeg)
library(ggpubr)

## set working to script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("dataFiles/jpegs_45/")

## plot list
plot_list<-list()
plot_names<-NULL
for(i in list.files()) {
  name<-ifelse(grepl("24", i), "LT",
               ifelse(grepl("34", i), "HT", 
                      ifelse(grepl("FA", i,), "FT", 
                             "Control")))
  print(name)
  print(i)
  img<-readJPEG(i)
  p<-ggdraw()+draw_image(img, clip="on")
  plot_list[[i]] <- p
  plot_names<-c(plot_names, name)
}

plot_names<-paste0(plot_names, rep(1:5, 4))

plot_names_t<-c(plot_names[c(1,6,11,16)], 
              plot_names[c(2,7,12,17)],
              plot_names[c(3,8,13,18)],
              plot_names[c(4,9,14,19)],
              plot_names[c(5,10,15,20)])

plot_list_t<-c(plot_list[c(1,6,11,16)], 
               plot_list[c(2,7,12,17)],
               plot_list[c(3,8,13,18)],
               plot_list[c(4,9,14,19)],
               plot_list[c(5,10,15,20)])

## combine
supp <- plot_grid(plotlist=plot_list_t,
                labels=plot_names_t,
                label_size=10,
                ncol=4)

## save
ggsave("../../figuresTables/Figure S[histo].pdf", supp, width=7, height=9)
ggsave("../../figuresTables/Figure S[histo].png", supp, width=7, height=9)

## end script