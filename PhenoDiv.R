library("Phenoflow")
library("dplyr")
library("gridExtra")
library("ggplot2")
library("RColorBrewer")

path = "data"

flowData <- read.flowSet(path = path, transformation = FALSE, pattern=".fcs")

flowData_transformed <- transform(flowData,`488nm 530/30BP-A`=asinh(`488nm 530/30BP-A`), 
                                  `SSC-A`=asinh(`SSC-A`), 
                                  `488nm 780/60BP-A`=asinh(`488nm 780/60BP-A`), 
                                  `FSC-A`=asinh(`FSC-A`))
param=c("488nm 530/30BP-A", "488nm 780/60BP-A","SSC-A","FSC-A")
flowData_transformed = flowData_transformed[,param]
remove(flowData)

### Create a PolygonGate for denoising the dataset
### Define coordinates for gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)
sqrcut1 <- matrix(c(5.5,5.5,2,2,15,15,
                    3,7,7,15,15,3),ncol=2, nrow=6)
colnames(sqrcut1) <- c("488nm 530/30BP-A","488nm 780/60BP-A")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

###  Gating quality check
xyplot(`488nm 780/60BP-A` ~ `488nm 530/30BP-A`, data=flowData_transformed[5:6], filter=polyGate1,
       scales=list(y=list(limits=c(0,14)),
                   x=list(limits=c(0,16))),
       axis = axis.default, nbin=125, 
       par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)

### Isolate only the cellular information based on the polyGate1
flowData_transformed <- Subset(flowData_transformed, polyGate1)

summary <- fsApply(x=flowData_transformed,FUN=function(x) apply(x,2,max),use.exprs=TRUE)
max = max(summary[,1])
mytrans <- function(x) x/max
flowData_transformed <- transform(flowData_transformed,`488nm 530/30BP-A`=mytrans(`488nm 530/30BP-A`),
                                  `488nm 780/60BP-A`=mytrans(`488nm 780/60BP-A`), 
                                  `SSC-A`=mytrans(`SSC-A`),
                                  `FSC-A`=mytrans(`FSC-A`))

Diversity.ww <- Diversity_rf(flowData_transformed, d=3, param = param, R = 100)

### Extract metadata

meta <-strsplit(gsub(as.character(Diversity.ww$Sample_names), pattern=".fcs", replacement=""),"_")
meta <- data.frame(do.call(rbind,meta))
colnames(meta) <- c("Timepoint", "Pond", "Position","Staining")

### Merge data
results <- data.frame(Diversity.ww, meta)
results2 <- results[results$Staining=="S",]
results2 <- droplevels(results2)

### Perform beta-diversity analysis but only for pond/splitter samples
flowData_beta <- FCS_resample(flowData_transformed[results$Staining=="S" & results$Timepoint!="T0"], replace=TRUE)
fbasis <- flowBasis(flowData_beta, nbin=128, param=param, bw=0.01)
beta <- beta_div_fcm(fbasis, ord.type="PCoA")
var <- round(vegan::eigenvals(beta)/sum(vegan::eigenvals(beta))*100,1)
beta_results <- data.frame(beta$points, meta[meta$Staining=="S" & results$Timepoint!="T0",])

### primary conclusions from the alpha-diversity analysis
# Wastewater has the highest diversity, the pure strains have the lowest diversity,
# as expected (I've plotted these at timepoint 0. There is a shift in alpha-diversity
# visible Pond2 and the splitter start at relatively low diversity but they increase over time,
# with the splitter community having the highest diversity at T3.
# 
a1 <- ggplot(data=results2, aes(x=Timepoint,y=D2,fill=Pond))+
  geom_point(shape=21, size=6)+
  scale_color_brewer(palette="Accent")+
  geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), width=0.05)+
  theme_bw()+
  labs(y="Phenotypic diversity (D2)", x="Timepoint")


### Primary conclusions from the beta-diversity analysis
# The difference between inlet and outlet becomes larger over time. The communities are
# really similar at T1 but differ more at T2/T3. There is a temporal shift noticeable
# as well for both the splitter as the ponds. Pond 2 starts out really similar to the 
# splitter community but becomes much more similar to pond1 at T2/T3. Diversity of pond 1 is quite
# high at T1 but decrease to comparable values with P2/SP at T3.
b1 <- ggplot(data=beta_results, aes(x=X1, y=X2, color=Pond, shape=Position))+
  geom_point(size=8)+
  theme_bw()+
  scale_color_brewer(palette="Accent")+
  labs(x = paste0("Axis1 (",var[1], "%)"), y = paste0("Axis2 (",var[2], "%)"),
       fill="")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
        title=element_text(size=20), legend.text=element_text(size=14))

b2 <- ggplot(data=beta_results[beta_results$Timepoint!="T0",], aes(x=X1, y=X2, color=Pond, shape=Timepoint))+
  geom_point(size=8)+
  theme_bw()+
  scale_color_brewer(palette="Accent")+
  labs(x = paste0("Axis1 (",var[1], "%)"), y = paste0("Axis2 (",var[2], "%)"),
       fill="")+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
        title=element_text(size=20), legend.text=element_text(size=14))

grid.arrange(b1, b2, ncol=2)

### Lets make some contrasts between the ponds
# P1 vs. P2 at timepoint 0
ct1 <- fp_contrasts(fbasis, comp1=beta_results$Pond=="P1" & beta_results$Timepoint=="T1", comp2=beta_results$Pond=="P2" & beta_results$Timepoint=="T1", param=param[1:2],
                    thresh = 0.1)
# P1 vs. P2 at timepoint 3
ct2 <- fp_contrasts(fbasis, comp1=beta_results$Pond=="P1" & beta_results$Timepoint=="T3", comp2=beta_results$Pond=="P2" & beta_results$Timepoint=="T3", param=param[1:2],
                    thresh = 0.1)

### Plot contrasts
### Red/positive values indicate higher density for pond1.
### blue/negative values indicate lower density for pond1.
v1 <- ggplot2::ggplot(ct1, ggplot2::aes(`488nm 530/30BP-A`, `488nm 780/60BP-A`, z = Density))+
  ggplot2::geom_tile(ggplot2::aes(fill=Density)) + 
  ggplot2::geom_point(colour="gray", alpha=0.4)+
  ggplot2::scale_fill_distiller(palette="RdBu", na.value="white", limits=c(-0.6,0.4)) + 
  ggplot2::stat_contour(ggplot2::aes(fill=..level..), geom="polygon", binwidth=0.1)+
  ggplot2::theme_bw()+
  ggplot2::geom_contour(color = "white", alpha = 1)+
  ggplot2::labs(title="Pond1 - Pond2 at T1")+
  ylim(50,85)+
  xlim(50,90)

v2 <- ggplot2::ggplot(ct2, ggplot2::aes(`488nm 530/30BP-A`, `488nm 780/60BP-A`, z = Density))+
  ggplot2::geom_tile(ggplot2::aes(fill=Density)) + 
  ggplot2::geom_point(colour="gray", alpha=0.4)+
  ggplot2::scale_fill_distiller(palette="RdBu", na.value="white", limits=c(-0.6,0.4)) + 
  ggplot2::stat_contour(ggplot2::aes(fill=..level..), geom="polygon", binwidth=0.1)+
  ggplot2::theme_bw()+
  ggplot2::geom_contour(color = "white", alpha = 1)+
  ggplot2::labs(title="Pond1 - Pond2 at T3")+
  ylim(50,85)+
  xlim(50,90)

grid.arrange(v1, v2, ncol=2)
