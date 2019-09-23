package sci_commands::gas_seurat;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("gas_seurat");

sub gas_seurat {

@ARGV = @_;
# Defaults

getopts("O:R:X", \%opt);

$die2 = "
scitools gas_seurat [options] [directory containing gene gene_activity_matrixes]




Options:
   -O    [STR]   Output Directory (default is [current working directory]/cicero_output)
   -R    [STR]   Rscript call (def = $Rscript)
   -X           Retain intermediate files (Default = delete)
                  
";


if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = "$ARGV[0]/seurat_output"};

system("mkdir $opt{'O'}");

open R, ">$opt{'O'}/seurat.r";

print R "
#this is the seurat3 version
library(\"Seurat\")
library(dplyr)

#read in unnorm matrix
UNNORM.datam <-read.table(\"$ARGV[0]/cicero.gene_activity_matrix.txt\",sep=\"\\t\",header=T)


#create seurat obj
scATACunnorm.data<-CreateSeuratObject(UNNORM.datam)

#we dont need the mitocondrial plot so just look at feature
png(\"Feturesvscounts.png\",width=12,height=12,units=\"in\",res=600)
VlnPlot(object = scATACunnorm.data, features = c(\"nFeature_RNA\", \"nCount_RNA\"), ncol = 2)
dev.off()

pdf(\"Feturesvscounts.pdf\",width=12,height=12)
VlnPlot(object = scATACunnorm.data, features = c(\"nFeature_RNA\", \"nCount_RNA\"), ncol = 2)
dev.off()



#look at feature scatter vs ncount 
png(\"Feturesvscounts2.png\",width=12,height=12,units=\"in\",res=600)
FeatureScatter(object = scATACunnorm.data, feature1 = \"nCount_RNA\", feature2 = \"nFeature_RNA\")
dev.off()

pdf(\"Feturesvscounts2.pdf\",width=12,height=12)
FeatureScatter(object = scATACunnorm.data, feature1 = \"nCount_RNA\", feature2 = \"nFeature_RNA\")
dev.off()


#filter unnorm data
scATACunnorm.data <- subset(x = scATACunnorm.data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)


#do global scale norm
scATACunnorm.data <- NormalizeData(object = scATACunnorm.data, normalization.method = \"LogNormalize\", scale.factor = 1e4)

#select for most variable genes
scATACunnorm.data <- FindVariableFeatures(object = scATACunnorm.data, selection.method = \'mean.var.plot\', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
length(x = VariableFeatures(object = scATACunnorm.data))

#normalize ncount RNA, much like in monocle
scATACunnorm.data <- ScaleData(object = scATACunnorm.data, features = rownames(x = scATACunnorm.data), vars.to.regress = c(\"nCount_RNA\"))

#run PCA
scATACunnorm.data <- RunPCA(object = scATACunnorm.data, pc.genes = scATACunnorm.data@var.genes,npcs = 100, do.print = FALSE)

#visualize it
png(\"First2PCA.png\",width=12,height=12,units=\"in\",res=600)
VizDimLoadings(object = scATACunnorm.data, dims = 1:2)
dev.off()

#visualize it
pdf(\"First2PCA.pdf\",width=12,height=12)
VizDimLoadings(object = scATACunnorm.data, dims = 1:2)
dev.off()


#DimPlot(object = scATACunnorm.data)


png(\"First12PCAHM.png\",width=12,height=12,units=\"in\",res=600)
DimHeatmap(object = scATACunnorm.data, dims = 1:20, cells = 500, balanced = TRUE)
dev.off()

pdf(\"First12PCAHM.pdf\",width=12,height=12)
DimHeatmap(object = scATACunnorm.data, dims = 1:20, cells = 500, balanced = TRUE)
dev.off()




# NOTE: This process can take a long time for big datasets, comment out for expediency.
# More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
scATACunnorm.data <- JackStraw(object = scATACunnorm.data, num.replicate = 100,dims = 70)
scATACunnorm.data <- ScoreJackStraw(object = scATACunnorm.data, dims = 1:70)
png(\"JackStrawPlo.png\",width=12,height=12,units=\"in\",res=600)
JackStrawPlot(object = scATACunnorm.data, dims = 1:70)
dev.off()
pdf(\"JackStrawPlo.pdf\",width=12,height=12)
JackStrawPlot(object = scATACunnorm.data, dims = 1:70)

dev.off()
png(\"ElbowPlot.png\",width=12,height=12,units=\"in\",res=600)
ElbowPlot(object = scATACunnorm.data,ndims = 100)
dev.off()


pdf(\"ElbowPlot.pdf\",width=12,height=12)
ElbowPlot(object = scATACunnorm.data,ndims = 100)
dev.off()
#chose 100


#you have to choose number of dims here
scATACunnorm.data <- FindNeighbors(object = scATACunnorm.data, dims = 1:100)
scATACunnorm.data <- FindClusters(object = scATACunnorm.data, resolution = 0.4)


scATACunnorm.data <- RunTSNE(object = scATACunnorm.data, dims = 1:100)

# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
png(\"TSNEof100PCA.png\",width=12,height=12,units=\"in\",res=600)
DimPlot(object = scATACunnorm.data, reduction = 'tsne')
dev.off()

png(\"TSNEof100PCA.pdf\",width=12,height=12)
DimPlot(object = scATACunnorm.data, reduction = \'tsne\')
dev.off()



#write out PCA
write.table(t(scATACunnorm.data@reductions$pca@cell.embeddings), quote=F,sep=\"\\t\", file = \"PCACoordinates_PC100.matrix\")



write.table(scATACunnorm.data@active.ident, quote=F,sep=\"\\t\",col.names = F, file = \"groups.annot\")


#find all markers seperating cell types
scATACnorm.markers <- FindAllMarkers(object = scATACunnorm.data, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
scATACnorm.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top10 <- scATACnorm.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)





png(\"DE_accross_clusters_heatmap.png\",width=12,height=12,units=\"in\",res=600)
DoHeatmap(object = scATACunnorm.data, features = top10$gene) + NoLegend()
dev.off()

png(\"DE_accross_clusters_heatmap.pdf\",width=12,height=12)
DoHeatmap(object = scATACunnorm.data, features = top10$gene) + NoLegend()
dev.off()



";

close R;
system("$Rscript $opt{'O'}/seurat.r");


if (!defined $opt{'X'}) {
    system("rm -f $opt{'O'}/seurat.r");
    system("rm -f $opt{'O'}/*.tmp");
}
}
1;
