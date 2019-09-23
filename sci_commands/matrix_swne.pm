package sci_commands::matrix_swne;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_swne");

sub matrix_swne {

@ARGV = @_;
# Defaults
$range_default = "1-15";

getopts("O:XRAC:D:", \%opt);

$die2 = "
scitools matrix-swne [options] [input tf-idf matrix] [k value based on the reconstruction error]
   or    swne
   or    piglet
 Applies SWNE based dim reduction and clustering on given TF-IDF matrix 

Options:
   -O   [STR]   Output prefix (default is [input].k.)
   -X           Retain intermediate files (def = delete)
   -D           Dims to use for pca def: 1-15
   -R   [STR]   Rscript call (def = $Rscript)
   -A           Annotation file, when provided [input].cells.factors.annot is output where factors are added 
   -C			Color Annotation file, when provided [input].colors.annot is output where factor with color back is added 

Note: Requires Seurat and swne R packages so dependencies need to look for those

";


if (!defined $ARGV[0]) {die $die2};
if (!defined $ARGV[1]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (!defined $opt{'D'}) {$opt{'D'} = $range_default};
read_ranges($opt{'D'});
$k=$ARGV[1];


open R, ">$opt{'O'}.SWNE.r";
print R "


#packages
library(Seurat)
library(swne)

norm.matrix<-read.table(file=\"$ARGV[0]\",row.names=1)

#use seurat package for PCA
S_matrix  <- CreateSeuratObject(raw.data = norm.matrix, project = \"SWNE\")

#Create scaled data, which is TF-IDF matrix in our case
S_matrix\@scale.data<-S_matrix\@raw.data
#run SVD much like in normal analysis, this is redundant in future
max_PC=max($range_R_set)+5


#first calculate PCA, calculate first 20 PC usuallym, this is done by the irlba package which computes fast truncated svd (very similar to what we do normally
S_matrix <- RunPCA(object = S_matrix,pcs.compute=max_PC,pc.genes = rownames(S_matrix\@data),reduction.name=\"svd\", do.print = F)

#png(\"$opt{'O'}.PCA.bow.plot.png\")
#PCElbowPlot(S_matrix, num.pc = max_PC)
#dev.off()

#calculate SNN based on PCA/svd and find internal clusters
S_matrix <- FindClusters(object = S_matrix, reduction.type = \"svd\", dims.use = $range_R_set, resolution = 0.6, print.output = 0, save.SNN = TRUE)
#or use this if you already know cell type, add case when cell types read in
#S_matrix <- BuildSNN(S_matrix, dims.use = $range_R_set, k.param = 30, k.scale = 10, prune.SNN = 1/20, force.recalc = T)

#if find clusters run to build SNN
clusters <- S_matrix\@ident; names(clusters) <- S_matrix\@cell.names;

#write out SNN clusters
write.table(as.matrix(clusters),file=\"$opt{'O'}.SNN.Clusters.$k.SWNE.annot\",quote=FALSE,sep=\"\\t\",row.names=TRUE,col.names=FALSE)

norm.counts<-as.matrix(S_matrix\@scale.data,\"dgCMatrix\")

## Unguided NMF
loss <- \"mse\" ## Loss function
n.cores <- 30 ## Number of cores to use
seed <- 32566 ## Set seed for 

# do NMF with k provided
nmf.res <- RunNMF(norm.counts, k = $ARGV[1], alpha = 0, init = \"nnsvd\", n.cores = n.cores, loss = loss)

#calc weight matrix
nmf.res\$W <- ProjectFeatures(norm.counts, nmf.res\$H, loss = \"mse\", n.cores = n.cores)
nmf.scores <- nmf.res\$H

## Run SWNE embedding
snn.matrix <- S_matrix\@snn[colnames(nmf.scores), colnames(nmf.scores)]
swne.embedding <- EmbedSWNE(nmf.scores, snn.matrix, alpha.exp = 1.0, snn.exp = 1, n_pull = 4, dist.use = \"IC\")



#write out dims and factors then combine and write combined dims
write.table(as.matrix(swne.embedding\$sample.coords),file=\"$opt{'O'}.Embedded.cells.$k.SWNE.dims\",quote=FALSE,sep=\"\\t\",row.names=TRUE,col.names=FALSE)

reform<-as.data.frame(cbind(swne.embedding\$H.coords\$x,swne.embedding\$H.coords\$y))
annot_form<-as.data.frame(cbind(swne.embedding\$H.coords\$name,rep(\"factor\",times=length(swne.embedding\$H.coords\$name))))

names(reform)<-c(\"x\",\"y\")
row.names(reform)<-swne.embedding\$H.coords\$name
out_comb<-rbind(as.matrix(swne.embedding\$sample.coords),as.matrix(reform))
write.table(as.matrix(reform),file=\"$opt{'O'}.Embedded.factors.$k.SWNE.dims\",quote=FALSE,sep=\"\\t\",row.names=TRUE,col.names=FALSE)
write.table(out_comb,file=\"$opt{'O'}.Embedded.cell.factors.$k.SWNE.dims\",quote=FALSE,sep=\"\\t\",row.names=TRUE,col.names=FALSE)
write.table(as.matrix(annot_form),file=\"$opt{'O'}.Embedded.factors.$k.SWNE.annot\",quote=FALSE,sep=\"\\t\",row.names=FALSE,col.names=FALSE)

## Associate factors with genes using the gene loadings (W) matrix
gene.loadings <- nmf.res\$W
gene.loadings <- t(apply(nmf.res\$W, 1, function(x) (x - min(x))/(max(x) - min(x))))



write.table(as.matrix(gene.loadings),file=\"$opt{'O'}.Embedded.factors.$k.loadings.SWNE.txt\",quote=FALSE,sep=\"\\t\",row.names=TRUE,col.names=TRUE)
top.factor.genes.df <- SummarizeAssocFeatures(gene.loadings, features.return = 1000)
write.table(as.matrix(top.factor.genes.df),file=\"$opt{'O'}.Embedded.factors.loadings.$k.top1000.loadings.SWNE.txt\",quote=FALSE,sep=\"\\t\",row.names=TRUE,col.names=TRUE)

## Make gene factor association heatmaps
gene.factor.heat <- gene.loadings[unique(top.factor.genes.df\$feature),]

pdf(\"$opt{'O'}.gene_factor_heatmap.$k.SWNE.pdf\", width = 7.5, height = 7.0)
ggHeat(gene.factor.heat, clustering = \"both\", x.lab.size = 14, y.lab.size = 5)
dev.off()


## Associate factors with cell clusters
clusters.list <- UnflattenGroups(clusters)
clusters.matrix <- t(swne:::.genesets_indicator(clusters.list, inv = F, return.numeric = T))
cluster.nmf.assoc <- FactorAssociation(clusters.matrix, nmf.scores, n.cores = n.cores, metric = \"IC\")

pdf(\"$opt{'O'}cluster_factor_heatmap.$k.SWNE.pdf\", width = 7.5, height = 4.5)
ggHeat(cluster.nmf.assoc, clustering = \"both\", x.lab.size = 14, y.lab.size = 14)
dev.off()

";
close R;

system("$Rscript $opt{'O'}.SWNE.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.SWNE.r");
}

#create new annot files with factors
if (defined $opt{'A'}) 
{

system("cat $opt{'A'} $opt{'O'}.Embedded.factors.$k.SWNE.annot > $opt{'O'}.cells.factors.$k.SWNE.annot");

};

#create new color annot files with factors
if (defined $opt{'C'}) 
{

system("cat $opt{'C'} \"factor\t\"black\"\" > $opt{'O'}.Embedded.factors.$k.colors.SWNE.annot");

};

}
1;
