package sci_commands::cds_monocle;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("cds_monocle");

sub cds_monocle {

@ARGV = @_;
# Defaults

getopts("O:R:XD:L:Ti:d:P:G", \%opt);

$die2 = "
scitools cds_monocle [options] [directory containing cds files]
   or    monocle_cds

Prior to running, ensure you have ran scitools matrix-makecds. That function will convert matrix files into the CDS format that Monocle3 requires.
Runs the current version of Monocle3 within a directory containing files for CDS format input. 


Options:
   -O   [STR]   Output Directory (default is [current working directory]/monocle_output)
   -R   [STR]   Rscript call (def = $Rscript)
   -T   [FLAG]   Use internal clustering within Monocle call. If not specified, 
                  Monocle runs with the user-supplied cds_dims_data file.
                  If -T not flagged, option -i is REQUIRED.
            -P 	[INT]	   Only if I is specified. 
                           Number of components to be used for dimensionality reduction by PCA to denoise.
               				Default: 50 components
   -i    [STR]    If -T is not flagged, option -i is REQUIRED.
                  IRLBA dims file or cisTopic Matrix file for weighting.
   -d    [STR]    If -T is not flagged, option -d is REQUIRED.
                  Dimensionality Reduction style of dims file (Accepts [umap|tsne])
   -D    [2|3]    Dimensions to be used for final plotting (2D or 3D plotting)
                  Default: 2 Dimensions
   -L    [STR] RGE method for branch analysis in monocle3. Must member of the following list:
   				[SimplePPT,L1graph,DDRTree] Default=SimplePPT
   -G    [FLAG]   Partition Groups. Will only build branches within the same annotation.
                  (Default=No Partition)
   -X          Retain intermediate files (Default = delete)
                  
";


if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'T'} && !defined $opt{'i'}) {die $die2};
if (defined $opt{'i'} && !defined $opt{'d'}) {die $die2};
if (defined $opt{'d'} && $opt{'d'} !~ /tsne|umap/) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = "$ARGV[0]/monocle_output"};
if (!defined $opt{'D'}) {$opt{'D'} = 2};
if (!defined $opt{'P'}) {$opt{'P'} = 50};
if (!defined $opt{'L'}) {$opt{'L'} = "SimplePPT"};

system("mkdir $opt{'O'}");

open R, ">$opt{'O'}/monocle.r";

print R "
suppressWarnings(library(monocle))
message(\"Loading Monocle3\")
# reading in matrix, annotation, and dimension data

cds_cell_data <- read.delim(\"$ARGV[0]/cds_cell_data.txt\")
cds_dims_data <- read.table(\"$ARGV[0]/cds_dims_data.txt\",header=F,row.names=1)
cds_site_data <- read.delim(\"$ARGV[0]/cds_site_data.txt\")
cds_counts_matrix <- read.table(\"$ARGV[0]/cds_counts_matrix.txt\")
message(\"Read in CDS Files.\")

rownames(cds_cell_data)==rownames(cds_dims_data)
cds_cell_data<-cbind(cds_cell_data,cds_dims_data)

# generating cell data set
feature_data <- new(\"AnnotatedDataFrame\", data = cds_site_data)
sample_data <- new(\"AnnotatedDataFrame\", data = cds_cell_data)

dimensions_data <- new(\"AnnotatedDataFrame\", data = cds_dims_data)

cds <- suppressWarnings(newCellDataSet(as.matrix(cds_counts_matrix), phenoData = sample_data, featureData = feature_data))
set.seed(2017)
pData(cds)\$cells <- NULL 
cds <- detectGenes(cds)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

";

if (defined $opt{'T'}) {
print R "
message(\"No dims file given, using Monocle3 for dimensionality reduction\")
cds<-preprocessCDS(cds,num_dim=$opt{'P'},use_tf_idf=TRUE,verbose=T)
message(\"Writing out PCA components data frame\")
write.table(as.data.frame(cds\@normalized_data_projection), file=\"$opt{'O'}/$opt{'P'}_pca.dims\",quote=F,sep=\"\\t\",row.names=T,col.names=T)

cds <- reduceDimension(cds, max_components = $opt{'D'},reduction_method = \'UMAP\',metric=\"cosine\",verbose = T)
write.table(as.data.frame(t(reducedDimS(cds))), file=\"$opt{'O'}/trajectory_$opt{'D'}D.dims\",quote=F,sep=\"\\t\",row.names=T,col.names=F)
cds <- partitionCells(cds)

";
} else {
print R "
FM <- exprs(cds)
cds\@auxOrderingData\$normalize_expr_data <- FM
irlba_pca_res<-t(read.table(\"$opt{'i'}\",header=T,row.names=1))
#irlba_pca_res = irlba_pca_res[match(rownames(irlba_pca_res),colnames(cds)),]
#added in matching, better than match 
matchnum<-rownames(irlba_pca_res) %in% colnames(cds)
irlba_pca_res = irlba_pca_res[matchnum,]
num_dim<-ncol(irlba_pca_res)
cds\@normalized_data_projection <- as.matrix(irlba_pca_res)

";

if ($opt{'d'} =~ /umap/) {
print R "
dim_no<-ncol(cds_dims_data)
S <- t(cds_dims_data)
Y <- S
W <- as.matrix(t(irlba_pca_res))
A <- S
colnames(A) <- colnames(cds)
colnames(S) <- colnames(cds)
colnames(Y) <- colnames(cds)
reducedDimA(cds) <- A
cds\@reducedDimW<- W
cds\@reducedDimS <- as.matrix(Y)
cds\@reducedDimK <- S
cds\@dim_reduce_type <- \"UMAP\"
pData(cds)\$umap_1 = reducedDimA(cds)[1, ]
pData(cds)\$umap_2 = reducedDimA(cds)[2, ]
if(dim_no==3){
pData(cds)\$umap_3 = reducedDimA(cds)[3, ]
}
";
} 

if ($opt{'d'} =~ /tsne/) {
print R "
dim_no<-ncol(cds_dims_data)
reducedDimA(cds) <- t(cds_dims_data)
cds\@auxClusteringData[[\"tSNE\"]]\$pca_components_used <- num_dim
cds\@dim_reduce_type <- \"tSNE\"
pData(cds)\$tsne_1 = reducedDimA(cds)[1, ]
pData(cds)\$tsne_2 = reducedDimA(cds)[2, ]
if(dim_no==3){
pData(cds)\$tsne_3 = reducedDimA(cds)[3, ]
}
";

}
}

if ($opt{'D'}==3) {
print R "
#3D Plotting
message(\"Generating 3D Plots\")
cds <- clusterCells(cds,verbose = T,cores=10)
cds <- suppressWarnings(partitionCells(cds)) 
";

   if (defined $opt{'G'}) {
   print R "
   cds <- learnGraph(cds, max_components = 3, RGE_method = \"$opt{'L'}\", partition_group=paste(colnames(pData(cds))[1]),do_partition=T,verbose = T)
   ";
   } else {
   print R "
   cds <- learnGraph(cds, max_components = 3, RGE_method = \"$opt{'L'}\", verbose = T)
   "; 
   }

print R "
#Writing out full CDS file
saveRDS(cds,file=\"$opt{'O'}/monocle.CDS.rds\")

#Save branch point coordinates, in two formats
dp_mst <- minSpanningTree(cds)
reduced_dim_coords <- reducedDimK(cds)
ica_space_df <- data.frame(Matrix::t(reduced_dim_coords[1:3,]))
colnames(ica_space_df) <- c(\"prin_graph_dim_1\", \"prin_graph_dim_2\",\"prin_graph_dim_3\")
ica_space_df\$sample_name <- row.names(ica_space_df)
ica_space_df\$sample_state <- row.names(ica_space_df)
edge_list <- as.data.frame(get.edgelist(dp_mst))
colnames(edge_list) <- c(\"source\", \"target\")

out_edge_list=data.frame()
i=1
source<-as.character(edge_list[i,1])
target<-as.character(edge_list[i,2])
out_edge_list<-rbind(c(ica_space_df[ica_space_df\$sample_name==source,],\"line_segment\"=i),c(ica_space_df[ica_space_df\$sample_name==target,],\"line_segment\"=i))
for (i in 2:nrow(edge_list)){
source<-as.character(edge_list[i,1])
target<-as.character(edge_list[i,2])
out_edge_list<-rbind(out_edge_list,c(ica_space_df[ica_space_df\$sample_name==source,],\"line_segment\"=i),c(ica_space_df[ica_space_df\$sample_name==target,],\"line_segment\"=i))
}
write.table(out_edge_list,file=\"$opt{'O'}/monocle3_segmentvectors.txt\",col.names=TRUE,row.names=FALSE,sep=\"\\t\",quote=FALSE)


edge_df <- merge(ica_space_df, edge_list, by.x = \"sample_name\",by.y = \"source\", all = TRUE)
edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = \"source_prin_graph_dim_1\",prin_graph_dim_2 = \"source_prin_graph_dim_2\", prin_graph_dim_3 = \"source_prin_graph_dim_3\"))
edge_df <- merge(edge_df, ica_space_df[, c(\"sample_name\",\"prin_graph_dim_1\", \"prin_graph_dim_2\", \"prin_graph_dim_3\")],by.x = \"target\", by.y = \"sample_name\", all = TRUE)
edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = \"target_prin_graph_dim_1\",prin_graph_dim_2 = \"target_prin_graph_dim_2\", prin_graph_dim_3 = \"target_prin_graph_dim_3\"))
write.table(as.matrix(edge_df),file=\"$opt{'O'}/monocle3_branchpoints.txt\",col.names=TRUE,row.names=FALSE,sep=\"\\t\",quote=FALSE)

#Plot out 3D version
plot_3d_cell_trajectory(cds,color_by=paste(colnames(pData(cds))[1]),webGL_filename=\"$opt{'O'}/trajectory_3D.html\",image_filename=\"$opt{'O'}/trajectory_3D.gif\",show_backbone=TRUE,backbone_segment_color=\"#000000\")

";

} else {
print R "
message(\"Generating Plots\")
cds <- clusterCells(cds,verbose = T,cores=10)
cds <- suppressWarnings(partitionCells(cds))
";
   if (defined $opt{'G'}) {
   print R "
   cds <- learnGraph(cds, max_components = 3, RGE_method = \"$opt{'L'}\", partition_group=paste(colnames(pData(cds))[1]),do_partition=T,verbose = T)
   ";
   } else {
   print R "
   cds <- learnGraph(cds, max_components = 3, RGE_method = \"$opt{'L'}\", verbose = T)
   "; 
   }

print R "
#Writing out full CDS file
saveRDS(cds,file=\"$opt{'O'}/monocle.CDS.rds\")

#Save branch point coordinates
dp_mst <- minSpanningTree(cds)
reduced_dim_coords <- reducedDimK(cds)
ica_space_df <- data.frame(Matrix::t(reduced_dim_coords[1:2,]))
colnames(ica_space_df) <- c(\"prin_graph_dim_1\", \"prin_graph_dim_2\")
ica_space_df\$sample_name <- row.names(ica_space_df)
ica_space_df\$sample_state <- row.names(ica_space_df)
edge_list <- as.data.frame(get.edgelist(dp_mst))
colnames(edge_list) <- c(\"source\", \"target\")
edge_df <- merge(ica_space_df, edge_list, by.x = \"sample_name\",by.y = \"source\", all = TRUE)
edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = \"source_prin_graph_dim_1\",prin_graph_dim_2 = \"source_prin_graph_dim_2\"))
edge_df <- merge(edge_df, ica_space_df[, c(\"sample_name\",\"prin_graph_dim_1\", \"prin_graph_dim_2\")],by.x = \"target\", by.y = \"sample_name\", all = TRUE)
edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = \"target_prin_graph_dim_1\",prin_graph_dim_2 = \"target_prin_graph_dim_2\"))
write.table(as.matrix(edge_df),file=\"$opt{'O'}/monocle3_branchpoints.txt\",col.names=TRUE,row.names=FALSE,sep=\"\\t\",quote=FALSE)

p<-plot_cell_trajectory(cds,cell_size=0.1,color_by = \"Cluster\",backbone_color=\"#000000\",show_state_name=T,alpha=0.3)
ggsave(plot=p,filename=\"$opt{'O'}/monocle3.cluster.plot.png\",width=5,height=4,dpi=900)
ggsave(plot=p,filename=\"$opt{'O'}/monocle3.cluster.plot.pdf\",width=5,height=4)

p<-plot_cell_trajectory(cds,cell_size=0.1,color_by = paste(colnames(pData(cds))[1]),backbone_color=\"#000000\",show_state_name=T,alpha=0.3)
ggsave(plot=p,filename=\"$opt{'O'}/monocle3.annot.plot.png\",width=5,height=4,dpi=900)
ggsave(plot=p,filename=\"$opt{'O'}/monocle3.annot.plot.pdf\",width=5,height=4)

p<-ggplot(edge_df,aes(x=source_prin_graph_dim_1,y=source_prin_graph_dim_2,xend=target_prin_graph_dim_1,yend=target_prin_graph_dim_2))+geom_segment()+geom_text(color=\"red\",label=edge_df\$sample_name,size=2,check_overlap=TRUE)
ggsave(plot=p,filename=\"$opt{'O'}/monocle3.branchname.plot.png\",width=5,height=4,dpi=900)
ggsave(plot=p,filename=\"$opt{'O'}/monocle3.branchname.plot.pdf\",width=5,height=4)

";
}

print R "
#Determine the root state.
cds<-orderCells(cds)
pr_graph_test <- principalGraphTest(cds, k=3, cores=10)
diff_access <- dplyr::add_rownames(pr_graph_test) %>% dplyr::arrange(plyr::desc(morans_test_statistic), plyr::desc(-qval))
write.table(as.matrix(diff_access),file=\"$opt{'O'}/monocle3_diffaccess.txt\",col.names=TRUE,row.names=FALSE,sep=\"\\t\",quote=FALSE)

#Overwrite monocle.CDS file with final analysis
saveRDS(cds,file=\"$opt{'O'}/monocle.CDS.rds\")

# writing out final data frames
write.table(as.matrix(pData(cds)),file=\"$opt{'O'}/monocle3_cells.txt\",col.names=TRUE,row.names=TRUE,sep=\"\\t\",quote=FALSE)
write.table(as.matrix(Biobase::exprs(cds)),file=\"$opt{'O'}/monocle3_aggragated_cells_count.txt\",col.names=TRUE,row.names=TRUE,sep=\"\\t\",quote=FALSE)
write.table(as.matrix(fData(cds)),file=\"$opt{'O'}/monocle3_features.txt\",col.names=TRUE,row.names=FALSE,sep=\"\\t\",quote=FALSE)

";
close R;
system("$Rscript $opt{'O'}/monocle.r");
if (!defined $opt{'X'}) {
    system("rm -f $opt{'O'}/monocle.r");
}
}
1;
