###### R script for comparing cell culture media ######

### 1. Comparison of DMEM/DMEM F-12/CDM v1.0 ----------------------------

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(pheatmap)

########## 1.1 Heatmap ##########
Growth.media <- read.table("Growth.media.txt", header = T, sep = "\t")

rownames(Growth.media) <- Growth.media[, 1]
Growth.media <- Growth.media[, -1] # Romove the first col beacuse it's already there
# THMAC.ARG <- THMAC.ARG[which(rowSums(THMAC.ARG) > 0), ]
# new.col <- paste("THMAC_PS", 1: 30, sep = "")            # Column series for PS samples
# THMAC.ARG <- THMAC.ARG[, new.col]
Growth.media[is.na(Growth.media)] <- 0
   # Establish annotation info.

annotation_row = data.frame(
  Components = factor(rep(c("Amino Acid", "Vitamins", "Growth Factors", "Supplements"), c (26,15,6,27))))
rownames(annotation_row) = rownames(Growth.media)
head(annotation_row)


pdf("Growth.media-2.pdf", width=9.7, height=16.2)
a <- pheatmap(Growth.media, scale = 'row', cluster_cols = T, gaps_row = c(26, 41, 47),
                        color = colorRampPalette(c("navy", "white", "firebrick3"))(60),
                        cluster_rows = F, angle_col = c("45"), annotation_row = annotation_row,
                        fontsize = 12, border_color = "black", display_numbers = Growth.media,
              number_color = "blue")# display_numbers = matrix(ifelse(THMAC.ARG > 0.005, "*", ""), nrow(THMAC.ARG)))
dev.off()

pdf("Growth.media-3.pdf", width=9.7, height=16.2)
a <- pheatmap(Growth.media, scale = 'row', cluster_cols = T, cutree_rows = 5,
              color = colorRampPalette(c("navy", "white", "firebrick3"))(60),
              cluster_rows = T, angle_col = c("45"), annotation_row = annotation_row,
              fontsize = 12, border_color = "black", display_numbers = Growth.media,
              number_color = "blue")# display_numbers = matrix(ifelse(THMAC.ARG > 0.005, "*", ""), nrow(THMAC.ARG)))
dev.off()


########## 2. Plackett–Burman design (Design of experiments) ----------------------------
 
 # pb stands for Plackett-Burman. Plackett-Burman designs (Plackett and Burman 1946) are generally used 
 # for screening many variables in relatively few runs, when interest is in main effects only, at least 
 # initially. Different from the regular fractional factorial designs created by function FrF2, they do not 
 # perfectly confound interaction terms with main effects but distribute interaction effects over several 
 # main effects. The designs with number of runs a power of 2 are an exception to this rule: they are just 
 # the resolution III regular fractional factorial designs and are as such not very suitable for screening 
 # because of a high risk of very biased estimates for the main effects of the factors. Where possible, 
 # these are therefore replaced by different designs (cf. below).

 # The resolution III 12 run array or 20 run array have been recommended for screening purposes, 
 # because they avoid complete aliasing of main effects with 2fis (because of GR > 3).

library(FrF2)

 ## Step 2.1 use pb() function
  # For proteins
plan.annotated <- pb(12, seed = 15143, factor.names = list(ITS = c(0, 0.005), FN = c(0, 2), ReHA = c(0, 0.8), FGF2 = c(0, 10),
                                  IGF1 = c(0, 5), EGF = c(0, 10), TGFb1 = c(0, 0.5), VEGF = c(0, 5)),
                     n12.taguchi = TRUE)
  # For non-protein components
plan.annotated2 <- pb(12, seed = 15143, factor.names = list(SS = c(0, 0.01), FC = c(0, 2), ZnSO4 = c(0, 1), Hydro.corti = c(0, 25),
                                                           Ethan. = c(0, 50), ME = c(0, 2.2), HEPES = c(0, 3000), P_LAA = c(0, 20)),
                     n12.taguchi = TRUE)
plan.annotated3 <- pb(12, seed = 15143, factor.names = list(RHA = c(0.5, 0.75), IGF1 = c(10, 20), LA = c(0, 2), CP2 = c(1, 1.5)),
                      n12.taguchi = TRUE)
  # 12.12.2022 updated on 01.13.2023
plan.annotated <- pb(12, seed = 15143, factor.names = list(SS = c(0, 4),
                                                            FC = c(0, 2), 
                                                            ZnSO4 = c(0, 1), 
                                                            HC21 = c(0, 37),
                                                            L_AG = c(0, 0.003),
                                                            Ethan. = c(0, 10), 
                                                            ME = c(0, 2.2), 
                                                            P_LAA = c(0, 25), 
                                                            T3 = c(0, 1.5),
                                                            T4 = c(0, 10),
                                                            Lipmix = c(0, 2)),
                      n12.taguchi = TRUE)
plan.annotated

summary(plan.annotated)

 ## Step 2.2 output design

export.design(plan.annotated, filename = "PB for OT2", type = "all", OutDec = ",")

 ## Step 2.3 Run the designated experiments
 ## Step 2.4 Data interpretation
   
PB_output <- read.csv("Results PB for 0.7.6.csv", header = 1)
asl <- PB_output$results

 # “asl” is average responses from the x repeated measurements
 # asl <- c(5.5, 3.655, 2.8075, 3.62, 4.395, 4.9775, 2.8075, 1.995, 1.515, 2.9175, 5.825, 1.81)

pot.annotated.resp <- add.response(plan.annotated, asl)
summary(lm(pot.annotated.resp))

 ## Step 2.5 plot
GF.MEplot <- MEPlot(pot.annotated.resp, abbrev = 5, cex.xax = 1.6, cex.main = 2,
                    main = paste("Main effects plot for ch-ADC 12-22-2022"))

GF.danielplot <- DanielPlot(pot.annotated.resp, code = TRUE, half = TRUE, alpha = 1, cex.main = 1.8, cex.pch = 1.2, 
           cex.lab = 1.4, cex.fac = 1.4, cex.axis = 1.2)

########## 3. Orthogonal design (Design of experiments) ----------------------------
 #当选出重要因素时，下一步常见是正交试验或响应面分析，用来优选参数。其实正交表的设计原理就是用尽量少的步骤
 #遍历掉因子空间，这样进行一定次数试验就可以发现最优组合。R中的实现基本都在DoE.base 包里，这个包
 #内置了一堆可以直接调用的正交表，可以根据需求进行查询。

 #例如我有6个因素，
 #水平数分别是2，3，3，2，2，6，然后我只打算做不超过54次试验，这时可以直接调用show.oas函数进行查询，
 #给出的正交表随意选一个就可以继续。

show.oas(nruns = c(0, 24), nlevels = c(3, 3, 3, 3, 3), showmetrics = TRUE)
show.oas(nruns = c(0, 12), nlevels = c(2, 2, 2, 2, 2), showmetrics = TRUE)
show.oas(nruns = c(0, 36), nlevels = c(2, 2, 3, 2, 3), showmetrics = TRUE)

## no suitable  resolution IV or more  array found
## 2  orthogonal  arrays found
## name              nruns                  lineage   GR GRind regular SCones   A3  A4  A5  A6   A7 A8
## 31 L24.2.16.3.1    24 2~13;3~1;4~1;:(4~1!2~3;) 3.00  3.00   FALSE  3 37.8 140 326 667 1149 1531
## 33 L24.2.13.3.1.4.1  24                          3.18  3.18   FALSE  0 36.8 124 262 496  765 863

#这里我们选分辨率略高的L24.2.13.3.1.4.1，从名字上看，这是一个24次试验表，可以包含10个两水平，
#8个三水平与1个六水平因子，这也是唯一一个分辨率高于3，可以排除交互作用的设计方法。这里我反复提到分
#辨率，实际上就是一种考察试验设计合理性的指标，分辨率3一般指只能区分没有交互作用的各因子贡献差异，
#高于3就可以区分一定的因子交互作用，可以用GR去计算一个广义分辨率并用oa.design来进一步优化这个设计，
#因为其实符合正交表只是众多选择的一个子集，不过根据优化方法的不同，优化时间也不太一样。
oa.GF <- oa.design(L12.2.11, nlevels = c(2, 2, 2, 2, 2), columns = "min34", factor.names = list(ITS = c(0.5, 1), ReHA = c(0.8, 1.6), FGF2 = c(5, 10),
                                                                                              IGF1 = c(5, 10), TGFb1 = c(0.25, 0.5)))
oa.GF
export.design(oa.GF, filename = "OA for GF", type = "all", OutDec = ",", replace = T)

# oa.Protein <- oa.design(L18.2.1.3.7, nlevels = c(3, 2, 3), columns = "min34", factor.names = list(RHA = c(0.25, 0.5, 0.75),
#                                                                                                  LA = c(0, 2), 
#                                                                                                  CP2 = c(1, 1.5, 2)))
# export.design(oa.GF, filename = "OA for GF", type = "all", OutDec = ",", replace = T)


# asl3 <- c(9.31, 7.98, 10.9, 8.42, 8.28, 7.31, 4.73, 6.43, 10.67, 7.02, 5.47, 7.46)
# pot.annotated.resp <- add.response(oa.GF, asl3)
# summary(lm(pot.annotated.resp))

OA_output <- read.csv("Results OA for 0.7.6.csv", header = 1)
asl <- PB_output$results

pot.annotated.resp <- add.response(oa.GF, asl3)
summary(lm(pot.annotated.resp))

########## 4. Response-surface analysis ---------------------------
library("rsm")

ccd.design <- ccd(3, n0= c(3,3), alpha = "rotatable", randomize = F, inscribed = F, 
                  oneblock = T)
ccd.design.2 <- ccd(3, n0= c(3,3), alpha = "rotatable", randomize = F, inscribed = F, 
                  oneblock = T)

ccd.design$RHA <- 250*ccd.design$x1+750
ccd.design$LA <- 1*ccd.design$x2+ 2.5
ccd.design$CP2 <- 300*ccd.design$x3+ 1500

ccd.design.2$Insulin <- 3*ccd.design.2$x1+5
ccd.design.2$Transf <- 1*ccd.design.2$x2+ 3
ccd.design.2$Se <- 2*ccd.design.2$x3+ 5

options(digits = 3)
View(ccd.design)
View(ccd.design.2)

write.table(ccd.design, file = "rsm design for ITS.csv", sep = ",", row.names = T)

 ### 将实验结果输入.csv文件最后一列（results）后，可以进行分析
 # Data interpretation
rsm_output <- read.csv("rsm design for ITS.csv", header = 1) # results输入csv格式的实验设计

ccd.design$y <- rsm_output$results
SO.ccd <- rsm(y ~ SO(x1,x2,x3), data = ccd.design)
summary(SO.ccd)
 # The second run
 # y.09162022 <- c(6.795,7.76,5.25,6.35,5.84,8.72,4.66,6.43,6.65,6.87,5.69,3.85,
 #                 5.40,5.47,6.8,9.6,6.43,5.84,7.09,5.84)
 # ccd.design.2$y <- y.09162022
 # SO.ccd.2 <- rsm(y ~ SO(x1,x2,x3), data = ccd.design.2)
 # summary(SO.ccd.2)


  # Optimized value calculation 09022022
RHA <- 250*0.2123387+750; RHA  # 0.2xxx is given as optimized x1
LA <- 1*1.0312313+ 2.5; LA
CP2 <- 300*-0.3710273+ 1500; CP2

 # Optimized value 09162022
 # Insulin <- 3*1.585+5; Insulin
 # Transf <- 1*(-10.418)+ 3; Transf
 # Se <- 2*-0.588+ 5; Se

  # plotting
jet.colors <- colorRampPalette( c("red", "dark grey") )
persp(SO.ccd.2, x2~x3, contours = list(z="bottom"), col = jet.colors(150), ticktype = "detailed")
persp(SO.ccd.2, x1~x2, contours = list(z="bottom"), col = jet.colors(150), ticktype = "detailed")
persp(SO.ccd.2, x1~x3, contours = list(z="bottom"), col = jet.colors(150), ticktype = "detailed")

########## 5. CDM data exploration ---------------------------
library(ggplot2)
library(vegan)

CDM <- read.csv("CDM.csv", header = 1)
cdm.GF <- as.data.frame(CDM)[ , 2:11]

# Correlation panel
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y, method = "spearman"), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Customize upper panel
lower.panel<-function(x, y){
  points(x,y, pch=19, col=rgb(0.5, 0, 1, alpha = 0.7))
  r <- round(cor(x, y, method = "spearman"), digits=2)
  txt <- paste0("R = ", r)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  text(0.7, 0.85, txt)
}

cdm.GF[6:29 ,2:10] <- decostand(cdm.GF, 'standardize', 1)     # z-score transform, mean = 0, sd = 1

pairs(cdm.GF[6:29 ,1:10], upper.panel = NULL,
      lower.panel = lower.panel, main = "Correlations between factors, CDM 0.5-0.7")







