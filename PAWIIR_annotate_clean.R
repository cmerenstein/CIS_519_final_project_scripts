#pawiir annotate
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(matrixStats)

data.dir = "~/Box/GCB_Courses/519_ML/precip_data/PAWIIR/filtered_feature_bc_matrix"
data = Read10X(data.dir) #readinData

pawiir = CreateSeuratObject(counts = data$`Gene Expression`)
pawiir[['ADT']] = CreateAssayObject(counts = data$`Antibody Capture`)
#extract features
pawiir@assays$ADT %>% rownames() #21 antibodies

adt.counts.pawiir = pawiir@assays$ADT@counts %>% as.matrix() %>% t() %>% as_data_frame()

adt.features = colnames(adt.counts.pawiir) %>% print()


#plot all markers

plot.list = list()
for (i in 1:length(adt.features)){
  print(adt.features[i])
  print(i)
  
  plot.list[[i]] = ggplot(adt.counts.pawiir, aes_string(x=as.name(adt.features[i]))) + 
    geom_histogram() + 
    scale_x_continuous(limits = c(0,100))
}

plot.list[9]

pdf("~/Box/GCB_Courses/519_ML/precip_data/pawiir.expression_scaled.pdf", height = 20, width = 20)
cowplot::plot_grid(plotlist = plot.list)
dev.off()

#set thresholds
adt.features 
thresholds = c(25,25,25,20,100,25, 20, 20, 50, 25, 25, 25, 30, 25, 15, 30, 30, 30, 35, 25, 100, 30)
expression.thresholds = tibble(marker = adt.features,
                               expression.threshold = thresholds) %>% print()

#plot thresholds
plot.list = list()
for (i in 1:length(adt.features)){
  print(adt.features[i])
  print(i)
  
  plot.list[[i]] = ggplot(adt.counts.pawiir, aes_string(x=as.name(adt.features[i]))) + 
    geom_histogram() + geom_vline(xintercept = thresholds[i], color = "red")
    scale_x_continuous(limits = c(0,100))
}

plot.list[1]

pdf("~/Box/GCB_Courses/519_ML/precip_data/PAWIIR_Annotation/pawiir.expression_withthreshold_scale100.pdf", height = 20, width = 20)
cowplot::plot_grid(plotlist = plot.list)
dev.off()

#table of cell x expression
adt.counts.pawiir

expression = adt.counts.pawiir #keep same dimensions

for (i in 1:length(colnames(expression))){
  
  if(colnames(expression[,i]) != expression.thresholds$marker[i]){
    stop()
    print('name not in same order')
  }else{
    
    expression[,i] = expression[,i] > expression.thresholds$expression.threshold[i] #evaluate wheter cell passes expression threshold
  }
  
}

expression$cells = pawiir@assays$ADT %>% colnames() #append cell names
expression


#mutually exclusive pairs
adt.features

t.exclusive = c("CD4.1", "CD8A.1", "CD3", "CD7.1", "CD2.1", "CD5.1", "CD1A.1", "CD94")
b.exclusive = c("CD19.1", "CD22.1")
myeloid.exclusive = c("CD10", "CD123", "CD14.1", "CD33.1", "CD66B")

mutually.exclusive.features = tibble(pair1 = c(rep(t.exclusive, 2),
                                               rep(myeloid.exclusive, 2),
                                               rep(t.exclusive, length(myeloid.exclusive))), 
                                    pair2 =c(rep(b.exclusive[1], length(t.exclusive)),
                                             rep(b.exclusive[2], length(t.exclusive)),
                                             rep(b.exclusive[1], length(myeloid.exclusive)),
                                             rep(b.exclusive[2], length(myeloid.exclusive)),
                                             rep(myeloid.exclusive[1], length(t.exclusive)),
                                             rep(myeloid.exclusive[2], length(t.exclusive)),
                                             rep(myeloid.exclusive[3], length(t.exclusive)),
                                             rep(myeloid.exclusive[4], length(t.exclusive)),
                                             rep(myeloid.exclusive[5], length(t.exclusive))),
                                    pair.name = paste(pair1,pair2,"double",sep="."))

precipitated = tibble(cells= pawiir@assays$ADT %>% colnames())

for (j in 1:nrow(mutually.exclusive.features)){
  
  precipitated[,mutually.exclusive.features$pair.name[j]] = expression[,mutually.exclusive.features$pair1[j]] & expression[,mutually.exclusive.features$pair2[j]]
  
}

precipitated$precip = apply(precipitated[,2:ncol(precipitated)],1,max)
precipitated$precip
#precipitated$concordance = apply(precipitated[,2:(ncol(precipitated)-1)],1,sum)/(ncol(precipitated)-2)

#precipitated$concordance[precipitated$precipitated == T] %>% hist(breaks = 30)
#precipitated$concordance[precipitated$precipitated == T]%>% table()


final.table = adt.counts.pawiir %>% cbind(precipitated) %>% as_tibble()
final.table$precip %>% table()
#make plots
plot.list2 =  plot.list3 = list()
for (j in 1:nrow(mutually.exclusive.features)){
  
  x.axis = mutually.exclusive.features$pair1[j]
  y.axis = mutually.exclusive.features$pair2[j]
  highlight = mutually.exclusive.features$pair.name[j]
  highlight.2 = "precip"
  
  plot.list2[[j]] = ggplot(final.table, aes_string(x=as.name(x.axis), y = as.name(y.axis), color = as.name(highlight))) + geom_point() + 
  scale_x_continuous(limits = c(0,300)) + 
  scale_y_continuous(limits = c(0,300)) + theme(legend.position = "none")
  
  plot.list3[[j]] = ggplot(final.table, aes_string(x=as.name(x.axis), y = as.name(y.axis), color = as.name(highlight.2))) + geom_point() + 
    scale_x_continuous(limits = c(0,300)) + 
    scale_y_continuous(limits = c(0,300)) + theme(legend.position = "none")
  
}


ggplot(data = final.table, aes(x=CD127, y = CD38, color = precip)) + geom_point() + 
  scale_x_continuous(limits = c(0,300)) + 
  scale_y_continuous(limits = c(0,300))

ggplot(data = final.table, aes(x=!! "CD127", y = !! "CD38")) + geom_point() + 
  scale_x_continuous(limits = c(0,300)) + 
  scale_y_continuous(limits = c(0,300))

pdf("~/Box/GCB_Courses/519_ML/precip_data/PAWIIR_Annotation/pawiir.precips.pdf", height = 20, width = 20)
cowplot::plot_grid(plotlist = plot.list2)
cowplot::plot_grid(plotlist = plot.list3)
dev.off()




#save files
setwd("~/Box/GCB_Courses/519_ML/precip_data/PAWIIR_Annotation/")
write.csv(precipitated %>% select(cells, precip), "precipitated_cells.csv")
write.csv(expression, "binarized_expression.csv")
write.csv(adt.counts.pawiir, "raw_adt_counts.csv")



