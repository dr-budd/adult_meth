## SET UP ----

## set working to script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(Gviz)
library(cowplot)

## for plots
ats <- 0.5 ## smaller than default
atf <- 1   ## plain text (not bold)
fst <- 12

## Methylation Cyp19a1a ----

## import methylation data (output from seqmonk)
data_cyp19a1_meth <- read.csv("dataFiles/seqMonk/cyp19a1a_meth.txt",sep= "\t")

## make NaN function for data frames (because it only works for lists)
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

## make NaN's zeros
data_cyp19a1_meth[is.nan(data_cyp19a1_meth)] = 0

## make a data frame containing the chromosome number, start position, end position and strand
df_cyp19a1_meth <- data.frame("chr"=rep(2,nrow(data_cyp19a1_meth)),"start"=data_cyp19a1_meth$Start,
                              "end"=data_cyp19a1_meth$End,"strand"=rep("*",nrow(data_cyp19a1_meth)))

## take the data frame and automatically find the columns that describe genomic ranges using makeGRangesFromDataFrame. 
## (returns them as a GRanges object)
## HT
gr_cyp19a1_meth_34D5 <- makeGRangesFromDataFrame(df_cyp19a1_meth)
gr_cyp19a1_meth_34D1 <- makeGRangesFromDataFrame(df_cyp19a1_meth)
gr_cyp19a1_meth_34D2 <- makeGRangesFromDataFrame(df_cyp19a1_meth)
## CT 
gr_cyp19a1_meth_24D1 <- makeGRangesFromDataFrame(df_cyp19a1_meth)
gr_cyp19a1_meth_24D2 <- makeGRangesFromDataFrame(df_cyp19a1_meth)
gr_cyp19a1_meth_24D3 <- makeGRangesFromDataFrame(df_cyp19a1_meth)
## CT 
gr_cyp19a1_meth_28D1 <- makeGRangesFromDataFrame(df_cyp19a1_meth)
gr_cyp19a1_meth_28D2 <- makeGRangesFromDataFrame(df_cyp19a1_meth)
gr_cyp19a1_meth_28D3 <- makeGRangesFromDataFrame(df_cyp19a1_meth)

## add the methylation values for females to the GRanges object (gr_cyp19a1_meth_F)
## HT 
values(gr_cyp19a1_meth_34D5) <- DataFrame(HT5 = data_cyp19a1_meth$X34D5)
values(gr_cyp19a1_meth_34D1) <- DataFrame(HT1 = data_cyp19a1_meth$X34D1)
values(gr_cyp19a1_meth_34D2) <- DataFrame(HT2 = data_cyp19a1_meth$X34D2)
## CT 
values(gr_cyp19a1_meth_24D1) <- DataFrame(CT1 = data_cyp19a1_meth$X24D1)
values(gr_cyp19a1_meth_24D2) <- DataFrame(CT2 = data_cyp19a1_meth$X24D2)
values(gr_cyp19a1_meth_24D3) <- DataFrame(CT3 = data_cyp19a1_meth$X24D3)
## CT 
values(gr_cyp19a1_meth_28D1) <- DataFrame(CT1 = data_cyp19a1_meth$X28D1)
values(gr_cyp19a1_meth_28D2) <- DataFrame(CT2 = data_cyp19a1_meth$X28D2)
values(gr_cyp19a1_meth_28D3) <- DataFrame(CT3 = data_cyp19a1_meth$X28D3)

## create Annotation Track object from the GRanges objects you made 
## HT
cyp19a1_meth_track_34D5 <- AnnotationTrack(gr_cyp19a1_meth_34D5, name = "HT5")
cyp19a1_meth_track_34D1 <- AnnotationTrack(gr_cyp19a1_meth_34D1, name = "HT1")
cyp19a1_meth_track_34D2 <- AnnotationTrack(gr_cyp19a1_meth_34D2, name = "HT2")
## CT
cyp19a1_meth_track_24D1 <- AnnotationTrack(gr_cyp19a1_meth_24D1, name = "CT1")
cyp19a1_meth_track_24D2 <- AnnotationTrack(gr_cyp19a1_meth_24D2, name = "CT2")
cyp19a1_meth_track_24D2 <- AnnotationTrack(gr_cyp19a1_meth_24D3, name = "CT3")
## Control
cyp19a1_meth_track_28D1 <- AnnotationTrack(gr_cyp19a1_meth_28D1, name = "CT1")
cyp19a1_meth_track_28D2 <- AnnotationTrack(gr_cyp19a1_meth_28D2, name = "CT2")
cyp19a1_meth_track_28D2 <- AnnotationTrack(gr_cyp19a1_meth_28D3, name = "CT3")

## Expression cyp19a1a ----

## read in RNAseq expression data
data_cyp19a1_rna <- read.csv("dataFiles/seqMonk/cyp19a1a_expr.txt",sep= "\t")


df_cyp19a1_rna<- data.frame("chr"=rep(2,nrow(data_cyp19a1_rna)),"start"=data_cyp19a1_rna$Start,
                             "end"=data_cyp19a1_rna$End,"strand"=rep("*",nrow(data_cyp19a1_rna)))


## HT
gr_cyp19a1_rna_34R5 <- makeGRangesFromDataFrame(df_cyp19a1_rna)  
gr_cyp19a1_rna_34R1 <- makeGRangesFromDataFrame(df_cyp19a1_rna)  
gr_cyp19a1_rna_34R2 <- makeGRangesFromDataFrame(df_cyp19a1_rna)  
## CT
gr_cyp19a1_rna_24R1 <- makeGRangesFromDataFrame(df_cyp19a1_rna)  
gr_cyp19a1_rna_24R2 <- makeGRangesFromDataFrame(df_cyp19a1_rna) 
gr_cyp19a1_rna_24R3 <- makeGRangesFromDataFrame(df_cyp19a1_rna) 
## Control
gr_cyp19a1_rna_28R1 <- makeGRangesFromDataFrame(df_cyp19a1_rna)  
gr_cyp19a1_rna_28R2 <- makeGRangesFromDataFrame(df_cyp19a1_rna) 
gr_cyp19a1_rna_28R3 <- makeGRangesFromDataFrame(df_cyp19a1_rna) 

## HT
values(gr_cyp19a1_rna_34R5) <- DataFrame(HT5 = data_cyp19a1_rna$X34R5_CB6VFANX.bam)
values(gr_cyp19a1_rna_34R1) <- DataFrame(HT1 = data_cyp19a1_rna$X34R1_CB6VFANX.bam) 
values(gr_cyp19a1_rna_34R2) <- DataFrame(HT2 = data_cyp19a1_rna$X34R2_CB6VFANX.bam) 
## CT
values(gr_cyp19a1_rna_24R1) <- DataFrame(CT1 = data_cyp19a1_rna$X24R1_CB6VFANX.bam) 
values(gr_cyp19a1_rna_24R2) <- DataFrame(CT2 = data_cyp19a1_rna$X24R2_CB6VFANX.bam) 
values(gr_cyp19a1_rna_24R3) <- DataFrame(CT3 = data_cyp19a1_rna$X24R3_CB6VFANX.bam)
## Control
values(gr_cyp19a1_rna_28R1) <- DataFrame(Control1 = data_cyp19a1_rna$X28R1_CB6VFANX.bam) 
values(gr_cyp19a1_rna_28R2) <- DataFrame(Control2 = data_cyp19a1_rna$X28R2_CB6VFANX.bam) 
values(gr_cyp19a1_rna_28R3) <- DataFrame(Control3 = data_cyp19a1_rna$X28R3_CB6VFANX.bam)

## create Annotation Track object from the GRanges objects you made 

cyp19a1_rna_track_34R5 <- AnnotationTrack(gr_cyp19a1_rna_34R5, name = "HT5")
cyp19a1_rna_track_34R1 <- AnnotationTrack(gr_cyp19a1_rna_34R1, name = "HT1")
cyp19a1_rna_track_34R2 <- AnnotationTrack(gr_cyp19a1_rna_34R2, name = "HT2")

cyp19a1_rna_track_24R1 <- AnnotationTrack(gr_cyp19a1_rna_24R1, name = "CT1")
cyp19a1_rna_track_24R2 <- AnnotationTrack(gr_cyp19a1_rna_24R2, name = "CT2")
cyp19a1_rna_track_24R3 <- AnnotationTrack(gr_cyp19a1_rna_24R3, name = "CT3")

cyp19a1_rna_track_28R1 <- AnnotationTrack(gr_cyp19a1_rna_28R1, name = "Control1")
cyp19a1_rna_track_28R2 <- AnnotationTrack(gr_cyp19a1_rna_28R2, name = "Control2")
cyp19a1_rna_track_28R3 <- AnnotationTrack(gr_cyp19a1_rna_28R3, name = "Control3")

## file contains some geneModel info for cyp19a1a 
gene_model_cyp <- read.csv("dataFiles/cyp19a1_gene_A_barra.csv")

## create genome axis track
gtrack_cyp <- GenomeAxisTrack(showTitle=TRUE, 
                          name=paste("Chr", gene_model_cyp[1,1]), 
                          fontcolor.title="#808080",
                          fontsize.title=10,
                          rotation.title=0)

## create grtrack from file
grtrack_cyp <- GeneRegionTrack(gene_model_cyp, 
                           name="Cyp",
                           background.title=NA,
                           rotation.title=0,
                           fontcolor.title="#808080")

## Plot HT-Control cyp19a1a ----

alltracks_cyp <- list(DataTrack(gr_cyp19a1_meth_28D1, name="Control (M1)\nMethylation (%)", 
                                background.title="#3B528BFF",
                                fontsize.title=fst, cex.axis=ats, font.axis=atf, type=c("p", "h"), 
                                ylim=c(0,100), col="black"),
                      DataTrack(gr_cyp19a1_rna_28R1, name="Control (M1)\nRNA expression", 
                                background.title="#3B528BFF",
                                fontsize.title=fst, cex.axis=ats, font.axis=atf, type="h", 
                                ylim=c(0, 4000), col="black"),
                      DataTrack(gr_cyp19a1_meth_28D2, name="Control (M2)\nMethylation (%)", 
                                background.title="#3B528BFF",
                                fontsize.title=fst, cex.axis=ats, font.axis=atf, type=c("p", "h"), 
                                ylim=c(0,100), col="black"),
                      DataTrack(gr_cyp19a1_rna_28R2, name="Control (M2)\nRNA expression", 
                                background.title="#3B528BFF",
                                fontsize.title=fst, cex.axis=ats, font.axis=atf, type="h", 
                                ylim=c(0, 4000), col="black"),
                      DataTrack(gr_cyp19a1_meth_28D3, name="Control (M1)\nMethylation (%)", 
                                background.title="#3B528BFF",
                                fontsize.title=fst, cex.axis=ats, font.axis=atf, type=c("p", "h"), 
                                ylim=c(0,100), col="black"),
                      DataTrack(gr_cyp19a1_rna_28R3, name="Control (M1)\nRNA expression", 
                                background.title="#3B528BFF",
                                fontsize.title=fst, cex.axis=ats, font.axis=atf, type="h", 
                                ylim=c(0, 4000), col="black"),
                      DataTrack(gr_cyp19a1_meth_34D1, name="HT (T1)\nMethylation (%)", 
                                background.title="#ED7953FF",
                                fontsize.title=fst, cex.axis=ats, font.axis=atf, type=c("p", "h"), 
                                ylim=c(0,100), col="black"),
                      DataTrack(gr_cyp19a1_rna_34R1, name="HT (T1)\nRNA expression", 
                                background.title="#ED7953FF",
                                fontsize.title=fst, cex.axis=ats, font.axis=atf, type="h", 
                                ylim=c(0, 4000), col="black"),
                      DataTrack(gr_cyp19a1_meth_34D2, name="HT (T2)\nMethylation (%)", 
                                background.title="#ED7953FF",
                                fontsize.title=fst, cex.axis=ats, font.axis=atf, type=c("p", "h"), 
                                ylim=c(0,100), col="black"),
                      DataTrack(gr_cyp19a1_rna_34R2, name="HT (T2)\nRNA expression", 
                                background.title="#ED7953FF",
                                fontsize.title=fst, cex.axis=ats, font.axis=atf, type="h", 
                                ylim=c(0, 4000), col="black"),
                      DataTrack(gr_cyp19a1_meth_34D5, name="HT (F1)\nMethylation (%)", 
                                background.title="#ED7953FF",
                                fontsize.title=fst, cex.axis=ats, font.axis=atf, type=c("p", "h"), 
                                ylim=c(0,100), col="black"),
                      DataTrack(gr_cyp19a1_rna_34R5, name="HT (F1)\nRNA expression", 
                                background.title="#ED7953FF",
                                fontsize.title=fst, cex.axis=ats, font.axis=atf, type="h", 
                                ylim=c(0, 4000), col="black"),
                      grtrack_cyp,
                      gtrack_cyp)

pdf(file="figuresTables/Figure S[singlebp_cyp].pdf", width=10, height = 13)
plotTracks(alltracks_cyp, from = 7180303, to = 7184800)
pushViewport(viewport(x=0.5, y=0.5, width=0.5, height=0.5, just=c("center","center")))
grid.rect(gp=gpar(fill="#92D895", alpha=0.3),
          x=unit(-0.06, "npc"), y=unit(0.56, "npc"), 
          width=unit(0.54, "npc"), height=unit(1.85, "npc"), 
          just="center")
popViewport()
dev.off()

png("figuresTables/Figure S[singlebp_cyp].png", width=10*300, height=13*300, res=300)
plotTracks(alltracks_cyp, from = 7180303, to = 7184800)
pushViewport(viewport(x=0.5, y=0.5, width=0.5, height=0.5, just=c("center","center")))
grid.rect(gp=gpar(fill="#92D895", alpha=0.3),
          x=unit(-0.06, "npc"), y=unit(0.56, "npc"), 
          width=unit(0.54, "npc"), height=unit(1.85, "npc"), 
          just="center")
popViewport()
dev.off()

## Methylation dmrt1 ----

# import methylation data (output from seqmonk)
data_dmrt1_meth <- read.csv("dataFiles/seqMonk/dmrt1_meth.txt",sep= "\t")

# make NaN function for data frames (because it only works for lists)
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

# make NaN's zeros
data_dmrt1_meth[is.nan(data_dmrt1_meth)] = 0

# make a data frame containing the chromosome number, start position, end position and strand
df_dmrt1_meth <- data.frame("chr"=rep(13,nrow(data_dmrt1_meth)),"start"=data_dmrt1_meth$Start,
                              "end"=data_dmrt1_meth$End,"strand"=rep("*",nrow(data_dmrt1_meth)))

## HT
gr_dmrt1_meth_34D5 <- makeGRangesFromDataFrame(df_dmrt1_meth)
gr_dmrt1_meth_34D1 <- makeGRangesFromDataFrame(df_dmrt1_meth)
gr_dmrt1_meth_34D2 <- makeGRangesFromDataFrame(df_dmrt1_meth)
## CT 
gr_dmrt1_meth_24D1 <- makeGRangesFromDataFrame(df_dmrt1_meth)
gr_dmrt1_meth_24D2 <- makeGRangesFromDataFrame(df_dmrt1_meth)
gr_dmrt1_meth_24D3 <- makeGRangesFromDataFrame(df_dmrt1_meth)
## CT 
gr_dmrt1_meth_28D1 <- makeGRangesFromDataFrame(df_dmrt1_meth)
gr_dmrt1_meth_28D2 <- makeGRangesFromDataFrame(df_dmrt1_meth)
gr_dmrt1_meth_28D3 <- makeGRangesFromDataFrame(df_dmrt1_meth)

## add the methylation values

## HT 
values(gr_dmrt1_meth_34D5) <- DataFrame(HT5 = data_dmrt1_meth$X34D5)
values(gr_dmrt1_meth_34D1) <- DataFrame(HT1 = data_dmrt1_meth$X34D1)
values(gr_dmrt1_meth_34D2) <- DataFrame(HT2 = data_dmrt1_meth$X34D2)
## CT 
values(gr_dmrt1_meth_24D1) <- DataFrame(CT1 = data_dmrt1_meth$X24D1)
values(gr_dmrt1_meth_24D2) <- DataFrame(CT2 = data_dmrt1_meth$X24D2)
values(gr_dmrt1_meth_24D3) <- DataFrame(CT3 = data_dmrt1_meth$X24D3)
## CT 
values(gr_dmrt1_meth_28D1) <- DataFrame(CT1 = data_dmrt1_meth$X28D1)
values(gr_dmrt1_meth_28D2) <- DataFrame(CT2 = data_dmrt1_meth$X28D2)
values(gr_dmrt1_meth_28D3) <- DataFrame(CT3 = data_dmrt1_meth$X28D3)

## create Annotation Track object from the GRanges objects

##HT
dmrt1_meth_track_34D5 <- AnnotationTrack(gr_dmrt1_meth_34D5, name = "HT5")
dmrt1_meth_track_34D1 <- AnnotationTrack(gr_dmrt1_meth_34D1, name = "HT1")
dmrt1_meth_track_34D2 <- AnnotationTrack(gr_dmrt1_meth_34D2, name = "HT2")
##CT
dmrt1_meth_track_24D1 <- AnnotationTrack(gr_dmrt1_meth_24D1, name = "CT1")
dmrt1_meth_track_24D2 <- AnnotationTrack(gr_dmrt1_meth_24D2, name = "CT2")
dmrt1_meth_track_24D2 <- AnnotationTrack(gr_dmrt1_meth_24D3, name = "CT3")
##Control
dmrt1_meth_track_28D1 <- AnnotationTrack(gr_dmrt1_meth_28D1, name = "CT1")
dmrt1_meth_track_28D2 <- AnnotationTrack(gr_dmrt1_meth_28D2, name = "CT2")
dmrt1_meth_track_28D2 <- AnnotationTrack(gr_dmrt1_meth_28D3, name = "CT3")

## Expression dmrt1 ----

## read in RNAseq expression data
data_dmrt1_rna <- read.csv("dataFiles/seqMonk/dmrt1_expr.txt",sep= "\t") #mydata

df_dmrt1_rna<- data.frame("chr"=rep(13,nrow(data_dmrt1_rna)),"start"=data_dmrt1_rna$Start, #mydata
                            "end"=data_dmrt1_rna$End,"strand"=rep("*",nrow(data_dmrt1_rna)))

## HT
gr_dmrt1_rna_34R5 <- makeGRangesFromDataFrame(df_dmrt1_rna)  
gr_dmrt1_rna_34R1 <- makeGRangesFromDataFrame(df_dmrt1_rna)  
gr_dmrt1_rna_34R2 <- makeGRangesFromDataFrame(df_dmrt1_rna)  
## CT
gr_dmrt1_rna_24R1 <- makeGRangesFromDataFrame(df_dmrt1_rna)  
gr_dmrt1_rna_24R2 <- makeGRangesFromDataFrame(df_dmrt1_rna) 
gr_dmrt1_rna_24R3 <- makeGRangesFromDataFrame(df_dmrt1_rna) 
## Control
gr_dmrt1_rna_28R1 <- makeGRangesFromDataFrame(df_dmrt1_rna)  
gr_dmrt1_rna_28R2 <- makeGRangesFromDataFrame(df_dmrt1_rna) 
gr_dmrt1_rna_28R3 <- makeGRangesFromDataFrame(df_dmrt1_rna) 

## HT
values(gr_dmrt1_rna_34R5) <- DataFrame(HT5 = data_dmrt1_rna$X34R5_CB6VFANX.bam)
values(gr_dmrt1_rna_34R1) <- DataFrame(HT1 = data_dmrt1_rna$X34R1_CB6VFANX.bam) 
values(gr_dmrt1_rna_34R2) <- DataFrame(HT2 = data_dmrt1_rna$X34R2_CB6VFANX.bam) 
## CT
values(gr_dmrt1_rna_24R1) <- DataFrame(CT1 = data_dmrt1_rna$X24R1_CB6VFANX.bam) 
values(gr_dmrt1_rna_24R2) <- DataFrame(CT2 = data_dmrt1_rna$X24R2_CB6VFANX.bam) 
values(gr_dmrt1_rna_24R3) <- DataFrame(CT3 = data_dmrt1_rna$X24R3_CB6VFANX.bam)
## Control
values(gr_dmrt1_rna_28R1) <- DataFrame(Control1 = data_dmrt1_rna$X28R1_CB6VFANX.bam) 
values(gr_dmrt1_rna_28R2) <- DataFrame(Control2 = data_dmrt1_rna$X28R2_CB6VFANX.bam) 
values(gr_dmrt1_rna_28R3) <- DataFrame(Control3 = data_dmrt1_rna$X28R3_CB6VFANX.bam)

## create Annotation Track object from the GRanges objects you made 
dmrt1_rna_track_34R5 <- AnnotationTrack(gr_dmrt1_rna_34R5, name = "HT5")
dmrt1_rna_track_34R1 <- AnnotationTrack(gr_dmrt1_rna_34R1, name = "HT1")
dmrt1_rna_track_34R2 <- AnnotationTrack(gr_dmrt1_rna_34R2, name = "HT2")

dmrt1_rna_track_24R1 <- AnnotationTrack(gr_dmrt1_rna_24R1, name = "CT1")
dmrt1_rna_track_24R2 <- AnnotationTrack(gr_dmrt1_rna_24R2, name = "CT2")
dmrt1_rna_track_24R3 <- AnnotationTrack(gr_dmrt1_rna_24R3, name = "CT3")

dmrt1_rna_track_28R1 <- AnnotationTrack(gr_dmrt1_rna_28R1, name = "Control1")
dmrt1_rna_track_28R2 <- AnnotationTrack(gr_dmrt1_rna_28R2, name = "Control2")
dmrt1_rna_track_28R3 <- AnnotationTrack(gr_dmrt1_rna_28R3, name = "Control3")

## file contains some geneModel info for dmrt1 
gene_model_dmrt <- read.csv("dataFiles/dmrt1_gene_b_barra.csv")

## create genome axis track
gtrack_dmrt <- GenomeAxisTrack(showTitle=TRUE, 
                          name=paste("Chr", gene_model_dmrt[1,1]), 
                          fontcolor.title="#808080",
                          rotation.title=0)

## create grtrack from file
grtrack_dmrt <- GeneRegionTrack(gene_model_dmrt, 
                               name="dmrt1",
                               background.title=NA,
                               rotation.title=0,
                               fontcolor.title="#808080")

## Plot HT-Control dmrt1 ----

alltracks_dmrt1 <- list(DataTrack(gr_dmrt1_meth_28D1, name="Control (M1)\nMethylation (%)", 
                                  fontsize.title=fst, cex.axis=ats, font.axis=atf, type=c("p", "h"), 
                                  ylim=c(0,100), col="black", background.title="#3B528BFF"),
                        DataTrack(gr_dmrt1_rna_28R1, name ="Control (M1)\nRNA expression",
                                  fontsize.title=fst, cex.axis=ats, font.axis=atf, type="h", 
                                  ylim=c(0,1000), col="black", background.title="#3B528BFF"),
                        DataTrack(gr_dmrt1_meth_28D2, name="Control (M2)\nMethylation (%)", 
                                  fontsize.title=fst, cex.axis=ats, font.axis=atf, type=c("p", "h"), 
                                  ylim=c(0,100), col="black", background.title="#3B528BFF"),
                        DataTrack(gr_dmrt1_rna_28R2, name ="Control (M2)\nRNA expression",
                                  fontsize.title=fst, cex.axis=ats, font.axis=atf, type="h", 
                                  ylim=c(0,1000), col="black", background.title="#3B528BFF"),
                        DataTrack(gr_dmrt1_meth_28D3, name="Control (M1)\nMethylation (%)",
                                  fontsize.title=fst, cex.axis=ats, font.axis=atf, type=c("p", "h"), 
                                  ylim=c(0,100), col="black", background.title="#3B528BFF"),
                        DataTrack(gr_dmrt1_rna_28R3, name ="Control (M1)\nRNA expression",
                                  fontsize.title=fst, cex.axis=ats, font.axis=atf, type="h", 
                                  ylim=c(0,1000), col="black", background.title="#3B528BFF"),
                        DataTrack(gr_dmrt1_meth_34D1, name="HT (T1)\nMethylation (%)",
                                  fontsize.title=fst, cex.axis=ats, font.axis=atf, type=c("p", "h"), 
                                  ylim=c(0,100), col="black", background.title="#ED7953FF"),
                        DataTrack(gr_dmrt1_rna_34R1, name ="HT (T1)\nRNA expression",
                                  fontsize.title=fst, cex.axis=ats, font.axis=atf, type="h", 
                                  ylim=c(0,1000), col="black", background.title="#ED7953FF"),
                        DataTrack(gr_dmrt1_meth_34D2, name="HT (T2)\nMethylation (%)",
                                  fontsize.title=fst, cex.axis=ats, font.axis=atf, type=c("p", "h"), 
                                  ylim=c(0,100), col="black", background.title="#ED7953FF"),
                        DataTrack(gr_dmrt1_rna_34R2, name ="HT (T2)\nRNA expression",
                                  fontsize.title=fst, cex.axis=ats, font.axis=atf, type="h", 
                                  ylim=c(0,1000), col="black", background.title="#ED7953FF"),
                        DataTrack(gr_dmrt1_meth_34D5, name="HT (F2)\nMethylation (%)",
                                  fontsize.title=fst, cex.axis=ats, font.axis=atf, type=c("p", "h"), 
                                  ylim=c(0,100), col="black", background.title="#ED7953FF"),
                        DataTrack(gr_dmrt1_rna_34R5, name ="HT (F2)\nRNA expression",
                                  fontsize.title=fst, cex.axis=ats, font.axis=atf, type="h", 
                                  ylim=c(0,1000), col="black", background.title="#ED7953FF"),
                        grtrack_dmrt,
                        gtrack_dmrt)


pdf("figuresTables/Figure S[singlebp_dmrt1].pdf", width=10, height = 13)
plotTracks(alltracks_dmrt1, from = 23722498, to = 23752272)
pushViewport(viewport(x=0.5, y=0.5, width=0.5, height=0.5, just=c("center","center")))
grid.rect(gp=gpar(fill="#92D895", alpha=0.3),
          x=unit(-0.21, "npc"), y=unit(0.56, "npc"), 
          width=unit(0.044, "npc"), height=unit(1.85, "npc"), 
          just="center")
popViewport()
dev.off()

png("figuresTables/Figure S[singlebp_dmrt1].png", width=10*300, height=13*300, res=300)
plotTracks(alltracks_dmrt1, from = 23722498, to = 23752272)
pushViewport(viewport(x=0.5, y=0.5, width=0.5, height=0.5, just=c("center","center")))
grid.rect(gp=gpar(fill="#92D895", alpha=0.3),
          x=unit(-0.21, "npc"), y=unit(0.56, "npc"), 
          width=unit(0.044, "npc"), height=unit(1.85, "npc"), 
          just="center")
popViewport()
dev.off()

## Plot promoter region (+/- 200 bp) for supp ----

## plot
pdf(file="figuresTables/Figure S[singlebp_dmrt1_zoom].pdf", width=10, height = 13)
plotTracks(alltracks_dmrt1, from=23724298, to=23724500) ## (-200 to + 200)
## add highlight using grid functions
pushViewport(viewport(x=0.5, y=0.5, width=0.5, height=0.5, just=c("center","center")))
grid.rect(gp=gpar(fill="#92D895", alpha=0.3),
          x=unit(0.25, "npc"), y=unit(0.56, "npc"), 
          width=unit(0.13, "npc"), height=unit(1.85, "npc"), 
          just="center")
popViewport()
dev.off()

png("figuresTables/Figure S[singlebp_dmrt1_zoom].png", width=10*300, height=13*300, res=300)
plotTracks(alltracks_dmrt1, from=23724298, to=23724500) ## (-200 to + 200)
## add highlight using grid functions
pushViewport(viewport(x=0.5, y=0.5, width=0.5, height=0.5, just=c("center","center")))
grid.rect(gp=gpar(fill="#92D895", alpha=0.3),
          x=unit(0.25, "npc"), y=unit(0.56, "npc"), 
          width=unit(0.13, "npc"), height=unit(1.85, "npc"), 
          just="center")
popViewport()
dev.off()

## Plot less (for main text) ----

## cyp19a1a

tracks_cyp19a1a <- list(DataTrack(gr_cyp19a1_meth_28D1, 
                          name="Methylation (%)", fontsize.title=fst,
                          cex.axis=ats, font.axis=atf, #col.axis="black",
                          type=c("p", "h"), background.title="#3B528BFF",
                          ylim=c(0,100), col="black"),
                DataTrack(gr_cyp19a1_rna_28R1, 
                          name="RNA expression", fontsize.title=fst,
                          cex.axis=ats, font.axis=atf,
                          type="h", background.title="#3B528BFF",
                          ylim=c(0,4000), col="black"),
                DataTrack(gr_cyp19a1_meth_34D5, 
                          name="Methylation (%)",fontsize.title=fst,
                          cex.axis=ats, font.axis=atf,
                          type=c("p", "h"), background.title="#ED7953FF",
                          ylim=c(0,100), col="black"),
                DataTrack(gr_cyp19a1_rna_34R5, 
                          name="RNA expression", fontsize.title=fst,
                          cex.axis=ats, font.axis=atf,
                          type="h", background.title="#ED7953FF",
                          ylim=c(0,4000), col="black"),
                ## NB the below wont show up if your plot is too small
                grtrack_cyp, 
                gtrack_cyp)

## dmrt1
tracks_dmrt1 <- list(DataTrack(gr_dmrt1_meth_28D1, 
                          name="Methylation (%)",fontsize.title=fst,
                          cex.axis=ats, font.axis=atf,
                          type=c("p", "h"), background.title="#3B528BFF",
                          ylim=c(0,100), col="black"),
                DataTrack(gr_dmrt1_rna_28R1, 
                          name="RNA expression", fontsize.title=fst,
                          cex.axis=ats, font.axis=atf,
                          type="h", background.title="#3B528BFF",
                          ylim=c(0,1000), col="black"),
                DataTrack(gr_dmrt1_meth_34D5, 
                          name="Methylation (%)",fontsize.title=fst,
                          cex.axis=ats, font.axis=atf,
                          type=c("p", "h"), background.title="#ED7953FF",
                          ylim=c(0,100), col="black"),
                DataTrack(gr_dmrt1_rna_34R5, 
                          name="RNA expression", fontsize.title=fst,
                          cex.axis=ats, font.axis=atf,
                          type="h", background.title="#ED7953FF",
                          ylim=c(0,1000), col="black"),
                ## NB the below wont show up if your plot is too small
                grtrack_dmrt, 
                gtrack_dmrt)

jpeg(file="figuresTables/Figure [singlebp_cyp_less].jpeg", 
     width=8.25*300, height=11.75*300/3, res=300)
plotTracks(tracks_cyp19a1a, from=7180303, to=7184800)
dev.off()

jpeg(file="figuresTables/Figure [singlebp_dmrt1_less].jpeg", 
     width=8.25*300, height=11.75*300/3, res=300)
plotTracks(tracks_dmrt1, from = 23722498, to = 23752272)
dev.off()

## join

dmrt1j<-jpeg::readJPEG("figuresTables/Figure [singlebp_dmrt1_less].jpeg")
cyp19a1aj<-jpeg::readJPEG("figuresTables/Figure [singlebp_cyp_less].jpeg")

rect_grob_dmrt <- grid::rectGrob(
  x=0.146, y=0.2, width=0.013, height=0.79, just=c("left", "bottom"),
  gp=grid::gpar(fill="#92D895", col=NA, alpha=0.25))

rect_grob_cyp <- grid::rectGrob(
  x=0.092, y=0.2, width=0.26, height=0.79, just=c("left", "bottom"),
  gp=grid::gpar(fill="#92D895", col=NA, alpha=0.25))

red_rect_grob_cyp <- grid::rectGrob(
  x=0.5238, y=0.2, width=0.0163, height=0.79, just=c("left", "bottom"),
  gp=grid::gpar(fill="darkred", col=NA, alpha=0.25))

dmrt1<-ggdraw()+draw_image(dmrt1j, clip="on")+
  draw_text("   Dmrt1", x=0.05, y=0.18, size=7.5, colour="#808080")+
  draw_text("HT female", x=0.019, y=0.4, size=9, angle=90)+
  draw_text("Control male", x=0.018, y=0.8, size=9, angle=90)+
  draw_grob(rect_grob_dmrt)

cyp19a1a<-ggdraw()+draw_image(cyp19a1aj, clip="on")+
  draw_text("   Cyp19a1a", x=0.05, y=0.18, size=7.5, colour="#808080")+
  draw_text("HT female", x=0.019, y=0.4, size=9, angle=90)+
  draw_text("Control male", x=0.018, y=0.8, size=9, angle=90)+
  draw_grob(rect_grob_cyp)+
  draw_grob(red_rect_grob_cyp)

hs <- plot_grid(cyp19a1a, dmrt1,
                labels=c("A", "B"),
                ncol=1)

ggplot2::ggsave(plot=hs, filename="figuresTables/Figure [singlebp].pdf",
       device="pdf", width=8.75, height=8, units=c("in"))

ggplot2::ggsave(plot=hs, filename="figuresTables/Figure [singlebp].png",
                width=8.75, height=8, units=c("in"))

## end script
