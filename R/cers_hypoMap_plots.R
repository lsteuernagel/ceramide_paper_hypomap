
##########
### Load single cell data
##########

# hypoMap
# requires to download HypoMap from: https://doi.org/10.17863/CAM.87955
hypoMap = readRDS("hypoMap.rds")

library(Seurat)
library(dplyr)
library(ggplot2)

#colors
short_palette = as.character(palette.colors(palette = "Okabe-Ito"))
short_palette = short_palette[!short_palette %in% c("#999999","#000000")]
getOkabeItoPalette = colorRampPalette(short_palette)

sourcedata_path = "philipp_cers_paper/source_data/"
dir.create(sourcedata_path)

##########
###Overview plot
##########

#Idents(hypoMap) = "Label"
hypoMap_cluster_plot = DimPlot(hypoMap,group.by = "C7_named",order = TRUE,cols = getOkabeItoPalette(7),pt.size = 1,raster = TRUE,raster.dpi = c(2048,2048),repel = TRUE,label = TRUE,label.size = 7)+
  xlab("UMAP")+ylab("UMAP")+theme(text = element_text(size=20),axis.text = element_text(size=20))+ggtitle("HypoMap overview")  +NoLegend() #+NoAxes()#
hypoMap_cluster_plot
ggsave(hypoMap_cluster_plot,filename = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/analysis_projects/philipp_cers_paper/plots/hypoMap_celltype_overview.pdf",width = 300,height = 300,units = "mm",device = "pdf")

# source files
figure1_sourcedata = hypoMap_cluster_plot$data
figure1_sourcedata$Cell_ID = rownames(figure1_sourcedata)

##########
### UMAPs pro cers
##########
library(scales)

cers_genes = c("Cers1","Cers2","Cers3","Cers4","Cers5","Cers6")


for(g in cers_genes){
  print(g)
  hypoMap_cers_plot = FeaturePlot(hypoMap,features = g,order = TRUE,cols = c("grey90","#0b3ebd"),
                                  pt.size = 1.3,raster = TRUE,raster.dpi = c(2048,2048),repel = TRUE,label = F,label.size = 6)+
    xlab("UMAP")+ylab("UMAP")+theme(text = element_text(size=20),axis.text = element_text(size=20))+
    scale_color_gradient(low = "grey90",high = "#0b3ebd",limits=c(0,4),oob=squish) 
   # scale_fill_continuous(limits=c(0,4),oob=squish) #+NoAxes()# +NoLegend()
  #hypoMap_cers_plot
  ggsave(hypoMap_cers_plot,filename = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/projects/analysis_projects/philipp_cers_paper/plots/",g,"_hypoMap_ordered.pdf"),width = 300,height = 300,units = "mm",device = "pdf")
}

## add cers expression
cers_expression = Seurat::FetchData(hypoMap,cers_genes)
figure1_sourcedata = cbind(figure1_sourcedata,cers_expression)
data.table::fwrite(figure1_sourcedata,file=paste0(sourcedata_path,"source_data_1_umap.txt"),sep="\t")

##########
### dotplot 6 x cersX C66 (damit Pomc,Agrp ersichtlich sind)  NUR NEURONEN
##########

target_col = "C66_named"
unsorted = as.numeric(stringr::str_remove(stringr::str_extract(unique(hypoMap@meta.data$C66_named ),"C66-[0-9]+"),"C66-"))
names(unsorted) = unique(hypoMap@meta.data$C66_named )
sorted = sort(unsorted)
hypoMap@meta.data$C66_named = factor(hypoMap@meta.data$C66_named ,levels = names(sorted))

idents_include = as.character(unique(hypoMap@meta.data$C66_named))
idents_include = idents_include[grepl("GABA|GLU",idents_include)]

# plot
Idents(hypoMap) = target_col
dotplot_cers = Seurat::DotPlot(object = hypoMap,
                                  features = cers_genes,
                                  scale = FALSE,
                                  cluster.idents = F,dot.scale = 7,idents = idents_include)
dotplot_cers = dotplot_cers + guides(color=guide_colourbar('Avg. Exp.'),size = guide_legend("Pct. Exp.")) +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size = 15,angle = 90, vjust = 0.35, hjust=0.75),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_color_gradient(low = "grey90",high = "#0b3ebd",limits = c(0,4), oob = scales::squish)+
  coord_flip()
dotplot_cers

ggsave(dotplot_cers,filename = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/projects/analysis_projects/philipp_cers_paper/plots/dotplot_cers_C66neurons.pdf"),width = 450,height = 250,units = "mm",device = "pdf")


##
figure2_sourcedata = dotplot_cers$data %>% dplyr::select(gene = features.plot, cluster = id, avg.exp, pct.exp, -avg.exp.scaled)
data.table::fwrite(figure2_sourcedata,file=paste0(sourcedata_path,"source_data_2_dotplot_c66.txt"),sep="\t")

##########
### dotplot 6 x cersX über Author Classs -- reduzieren
##########

target_col = "Author_Class_Curated"

idents_include = as.character(unique(hypoMap@meta.data$Author_Class_Curated))
idents_include = idents_include[!idents_include %in% c("ParsTuber","Dividing","Erythroid-like","Hypendymal")]

# plot
Idents(hypoMap) = target_col
dotplot_cers = Seurat::DotPlot(object = hypoMap,
                               features = cers_genes,
                               scale = FALSE,
                               cluster.idents = F,dot.scale = 20,idents = idents_include)
dotplot_cers = dotplot_cers + guides(color=guide_colourbar('Avg. Exp.'),size = guide_legend("Pct. Exp.")) +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size = 15,angle = 90, vjust = 0.35, hjust=0.75),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_color_gradient(low = "grey90",high = "#0b3ebd",limits = c(0,4), oob = scales::squish)+
  coord_flip()
dotplot_cers

ggsave(dotplot_cers,filename = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/projects/analysis_projects/philipp_cers_paper/plots/dotplot_cers_classes.pdf"),width = 450,height = 250,units = "mm",device = "pdf")

#
figure3_sourcedata = dotplot_cers$data %>% dplyr::select(gene = features.plot, cluster = id, avg.exp, pct.exp, -avg.exp.scaled)
data.table::fwrite(figure3_sourcedata,file=paste0(sourcedata_path,"source_data_3_dotplot_class.txt"),sep="\t")

##########
### dotplot 6 x cersX über Agrp, Pomc, sf1 positive 
##########

# Group Gpr149, Vcan and Tcf4 from GLU-3 into one Sf1 cluster 

# + POMC + AGRP

hypoMap@meta.data$new_neurons = "ignore"
hypoMap@meta.data$new_neurons[hypoMap@meta.data$C66_named == "C66-19: Pomc.GLU-5" ] = "POMC neurons"
hypoMap@meta.data$new_neurons[hypoMap@meta.data$C66_named == "C66-46: Agrp.GABA-4" ] = "AgRP neurons"
hypoMap@meta.data$new_neurons[hypoMap@meta.data$C66_named %in% c("C66-12: Gpr149.GLU-3","C66-13: Vcan.GLU-3","C66-14: Tcf4.GLU-3") ] = "Sf1 neurons"

target_col = "new_neurons"

idents_include = c("POMC neurons","AgRP neurons","Sf1 neurons")

# plot
Idents(hypoMap) = target_col
dotplot_cers = Seurat::DotPlot(object = hypoMap,
                               features = cers_genes,
                               scale = FALSE,
                               cluster.idents = F,dot.scale = 25,idents = idents_include)
dotplot_cers = dotplot_cers + guides(color=guide_colourbar('Avg. Exp.'),size = guide_legend("Pct. Exp.")) +
  theme(text = element_text(size = 25),
        axis.text.x = element_text(size = 15,angle = 90, vjust = 0.35, hjust=0.75),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_color_gradient(low = "grey90",high = "#0b3ebd",limits = c(0,4), oob = scales::squish)+
  coord_flip()
dotplot_cers

ggsave(dotplot_cers,filename = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/projects/analysis_projects/philipp_cers_paper/plots/dotplot_cers_specific_neurons.pdf"),width = 250,height = 250,units = "mm",device = "pdf")

figure4_sourcedata = dotplot_cers$data %>% dplyr::select(gene = features.plot, cluster = id, avg.exp, pct.exp, -avg.exp.scaled)
data.table::fwrite(figure4_sourcedata,file=paste0(sourcedata_path,"source_data_4_dotplot_neurons.txt"),sep="\t")

##########
### barplot
##########

# out of all cersX positive cells how many are neurons, oligos, astros etc
# 1 bar per cers gene with different colors
# to show that cers6 is mostly in neurons

target_col = "Author_Class_Curated"

idents_include = as.character(unique(hypoMap@meta.data$Author_Class_Curated))
idents_include = idents_include[!idents_include %in% c("ParsTuber","Dividing","Erythroid-like","Hypendymal")]

# count cells per class that I am interested in
gene_expr = Seurat::FetchData(hypoMap,vars = cers_genes)
group_factor = hypoMap@meta.data[,target_col]
# calc pct
gene_expr[gene_expr > 0] <- 1 # set to 1 for occ
cluster_length = table(group_factor)
per_Cluster_occ=apply(gene_expr,2,function(x,group_factor){tapply(x,INDEX=group_factor,FUN=sum)},group_factor = group_factor)

# summarize
per_Cluster_occ_filt = per_Cluster_occ[rownames(per_Cluster_occ) %in% idents_include,]
per_Cluster_occ_filt_long = per_Cluster_occ_filt %>% as.data.frame() %>% dplyr::mutate(cluster = rownames(per_Cluster_occ_filt)) %>%
  tidyr::gather(key="gene",value="count",-cluster) %>% dplyr::group_by(gene) %>% dplyr::mutate(gene_total = sum(count)) %>%
  dplyr::mutate(pct = count / gene_total *100) %>% dplyr::arrange((pct))
# chnage order
per_Cluster_occ_filt_long$cluster = factor(per_Cluster_occ_filt_long$cluster,levels = per_Cluster_occ_filt_long$cluster[per_Cluster_occ_filt_long$gene == "Cers5"])
per_Cluster_occ_filt_long$Celltype = per_Cluster_occ_filt_long$cluster
#plot
set.seed(35)
set.seed(37)
barplot_cers = ggplot(per_Cluster_occ_filt_long,aes(x=gene,y=pct,fill=Celltype))+
  geom_bar(stat="identity")+
 # geom_text(aes(x = gene, y = 20,label=gene_total), hjust =1,size=6)+ # add cell numbers
  scale_fill_manual(values=sample(getOkabeItoPalette(10),10))+
  ylab("Percentage of all expressing cells")+xlab("")+
  coord_flip()+cowplot::theme_cowplot( font_size = 20)

barplot_cers

## with 
ggsave(barplot_cers,filename = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/projects/analysis_projects/philipp_cers_paper/plots/barplot_cers.pdf"),
       width = 300,height = 250,units = "mm",device = "pdf")

ggsave(barplot_cers,filename = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/projects/analysis_projects/philipp_cers_paper/plots/barplot_cers_cellnumbers.pdf"),
       width = 300,height = 250,units = "mm",device = "pdf")

figure5_sourcedata = barplot_cers$data %>% dplyr::select(-Celltype,-.group)
data.table::fwrite(figure5_sourcedata,file=paste0(sourcedata_path,"source_data_5_barplot.txt"),sep="\t")

##########
### All source data as one list
##########

source_list = list(
  umap_data_1 = figure1_sourcedata,
  dotplot_c66_2 = figure2_sourcedata,
  dotplot_class_3 = figure3_sourcedata,
  dotplot_neurons_4 = figure4_sourcedata,
  barplot_5 = figure5_sourcedata
)


WriteXLS::WriteXLS(x = source_list,ExcelFileName = paste0(sourcedata_path,"source_data_hypomap.xlsx"),SheetNames = names(source_list),col.names=TRUE)




