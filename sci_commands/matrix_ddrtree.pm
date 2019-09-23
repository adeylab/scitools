package sci_commands::matrix_ddrtree;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_ddrtree");

sub matrix_ddrtree {

@ARGV = @_;
# Defaults

getopts("O:XR:C:c:", \%opt);

$die2 = "
scitools matrix-ddrtree [options] [input matrix] [annotation file] [dims file]
   or    ddrtree-matrix
 Applies ddrt to provided matrix 

Options:
   -O   [STR]   Output prefix (default is [input].ddrt.dims)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)
   -C   [STR]   Color coding file (annot (tab) #hexColor)
   -c   [STR]   Color coding string
                  Annot=#hexColor,Annot2=#hexColor
                  
Note: Requires monocle2 and cicero R packages so dependencies need to look for those
      This works specifically with monocle2. Will be upgraded once monocle 3 is more stable. Test
";


if (!defined $ARGV[0]) {die $die2};
if (!defined $ARGV[1]) {die $die2};
if (!defined $ARGV[2]) {die $die2};

if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'C'} && defined $opt{'c'}) {die "\nSpecify either a color string (-c) or a color coding file (-C), not both!\n$die2"};
if (defined $opt{'C'}) {read_color_file($opt{'C'})};
if (defined $opt{'c'}) {read_color_string($opt{'c'})};


open SITE_DATA, ">./cds_site_data.txt";
open CELL_DATA, ">./cds_cell_data.txt";
open DIMS_DATA, ">./cds_dims_data.txt";
open COUNTS, ">./cds_counts_matrix.txt";
open LAMBDA_DIST, ">./Lambda_dist.txt";

read_annot($ARGV[1]);
read_dims($ARGV[2]);



#make sure annot matches matrix
open IN, "$ARGV[0]";
$h = <IN>; chomp $h; @H = split(/\t/, $h);
$h_out = "";
for ($i = 0; $i < @H; $i++) {
	if ((defined $CELLID_annot{$H[$i]}) && defined ($CELLID_DIMS{$H[$i]})) {
		$h_out .= "$H[$i]\t";
	}
}
#finish printing out header by removing last \t
$h_out =~ s/\t$//;
print COUNTS "$h_out\n";
#print out header

print SITE_DATA "site_name\tchr\tbp1\tbp2\tnum_cells_expressed\tsite_length\n";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$site = shift(@P);
	($chr,$bp1,$bp2) = split(/[_]/, $site);
	$siteName = "$chr\_$bp1\_$bp2";
	$siteOut = "$siteName";
	$siteOut2 = "$siteName";
	$chrID = $chr; $chrID =~ s/chr//;
	my $SITENAME_maxSignal=0;
	for ($i = 0; $i < @P; $i++) {
		if ((defined $CELLID_annot{$H[$i]}) && defined ($CELLID_DIMS{$H[$i]})) {
			$SITENAME_totalSignal{$siteName}+=$P[$i];
			
			if ($P[$i] > $SITENAME_maxSignal{$siteName}){
				$SITENAME_maxSignal{$siteName}=$P[$i];
				#print $siteName."\t$P[$i]\t".$SITENAME_maxSignal{$siteName}."\n";
			}
			
			if ($P[$i]>0) {
				$SITENAME_expressed{$siteName}++;
			}
			
			
			
			if ($P[$i]>0) {
				$CELLID_expressed{$H[$i]}++;
			}
			$siteOut .= "\t$P[$i]";
		}
	}
	print COUNTS "$siteOut\n";
	
	
	for ($i = 0; $i < @P; $i++) 
		{
			if ((defined $CELLID_annot{$H[$i]}) && defined ($CELLID_DIMS{$H[$i]})) 
			{
		
			if ($P[$i]>0) {
				$tempval = $P[$i]/$SITENAME_maxSignal{$siteName};
				#print $siteName."\t$P[$i]\t$SITENAME_maxSignal{$siteName}\t".$tempval."\n";
				$siteOut2 .= "\t$tempval";
			}
			else
			{
				$siteOut2 .= "\t0";
			}
			}
		
		}
	
	#I THINK YOU NEED THE NUMBER OF CELLS EXPRESSED PER SITE HERE NOT THE TOTAL SIGNAL
	#OLD: Have double checked! it works!
	#print SITE_DATA "$siteName\t$siteName\t$chrID\t$bp1\t$bp2\t$SITENAME_totalSignal{$siteName}\t".($bp2-$bp1)."\n";
	print SITE_DATA "$siteName\t$siteName\t$chrID\t$bp1\t$bp2\t$SITENAME_expressed{$siteName}\t".($bp2-$bp1)."\n";
	
} close IN;
close COUNTS; close SITE_DATA; 

print CELL_DATA "cells\ttimepoint\tnum_genes_expressed\n";
print DIMS_DATA "dimension_1\tdimension_2\n";

for ($i = 0; $i < @H; $i++) {
	if ((defined $CELLID_annot{$H[$i]}) && defined ($CELLID_DIMS{$H[$i]})) {
		print DIMS_DATA  join("\t",@{$CELLID_DIMS{$H[$i]}})."\n";
		print CELL_DATA "$H[$i]\t$H[$i]\t$CELLID_annot{$H[$i]}\t$CELLID_expressed{$H[$i]}\n";
	}
} close CELL_DATA; close DIMS_DATA;



open R, ">$opt{'O'}.ddrt.r";

print R "
library(monocle)
library(cicero)
library(igraph)
# reading in matrix, annotation, and dimension data
cds_cell_data <- read.delim(\"./cds_cell_data.txt\")
cds_dims_data <- read.delim(\"./cds_dims_data.txt\")
cds_site_data <- read.delim(\"./cds_site_data.txt\")
cds_counts_matrix <- read.table(\"./cds_counts_matrix.txt\")

#make it binary this is just for cicero
#threshold <- 1
#bcdata <- ifelse(cds_counts_matrix < threshold, 0, 1)
#cds_counts_matrix<-bcdata

#check if these two are in the same order and combine for mean

rownames(cds_cell_data)==rownames(cds_dims_data)

cds_cell_data<-cbind(cds_cell_data,cds_dims_data)



# generating cell data set
feature_data <- new(\"AnnotatedDataFrame\", data = cds_site_data)
sample_data <- new(\"AnnotatedDataFrame\", data = cds_cell_data)
dimensions_data <- new(\"AnnotatedDataFrame\", data = cds_dims_data)
input_cds <- newCellDataSet(as.matrix(cds_counts_matrix), phenoData = sample_data, featureData = feature_data)
input_cds\@expressionFamily <- binomialff()
input_cds\@expressionFamily\@vfamily <- \"binomialff\"



set.seed(2017)
";


print R "# add cell data
	

pData(input_cds)\$cells <- NULL	
agg_cds <- aggregate_nearby_peaks(input_cds, distance = 10000)
agg_cds <- detectGenes(agg_cds)
agg_cds <- estimateSizeFactors(agg_cds)
agg_cds <- estimateDispersions(agg_cds)

#from here it is good
dimA<-t(as.matrix(cds_dims_data,rownames=F))
rownames(dimA)<-NULL
reducedDimA(agg_cds)<-dimA


#instead of this we can potentially add in the aggregate group kmers. Ok for now
agg_cds <- clusterCells(agg_cds, verbose = F,cores=10)

clustering_DA_sites <- differentialGeneTest(agg_cds,fullModelFormulaStr = '~Cluster')

#might not need to use this
# This takes a few minutes to run
#diff_timepoint <- differentialGeneTest(agg_cds,
#                  fullModelFormulaStr=\"~timepoint + num_genes_expressed\")

ordering_sites <- row.names(clustering_DA_sites)[order(clustering_DA_sites\$qval)][1:10000]

agg_cds <- setOrderingFilter(agg_cds, ordering_sites)

agg_cds <- reduceDimension(agg_cds, max_components = 2,
          residualModelFormulaStr=\"~num_genes_expressed\",
          reduction_method = \'DDRTree\')
agg_cds <- orderCells(agg_cds)





p<-plot_cell_trajectory(agg_cds, color_by = \"timepoint\")



";

if ($color_mapping !~ /none/i) {
	print R "
	p<-p+scale_colour_manual(values = c($color_mapping))";
}



print R "

ggsave(plot=p,filename=\"$opt{'O'}.annot_plot.png\",width=5,height=4,dpi=900)
ggsave(plot=p,filename=\"$opt{'O'}.annot_plot.pdf\",width=5,height=4);

p<-plot_cell_trajectory(agg_cds, color_by = \"State\")
ggsave(plot=p,filename=\"$opt{'O'}.state_plot.png\",width=5,height=4,dpi=900)
ggsave(plot=p,filename=\"$opt{'O'}.state_plot.pdf\",width=5,height=4);


cds<-agg_cds
#Overwrite monocle.CDS file with final analysis
saveRDS(cds,file=\"monocle.CDS.rds\")


#Save branch point coordinates
dp_mst <- minSpanningTree(cds)
reduced_dim_coords <- reducedDimK(cds)
ica_space_df <- data.frame(Matrix::t(reduced_dim_coords[1:2,]))
colnames(ica_space_df) <- c(\"prin_graph_dim_1\", \"prin_graph_dim_2\")
ica_space_df\$sample_name <- row.names(ica_space_df)
ica_space_df\$sample_state <- row.names(ica_space_df)
#edge_list <- as.data.frame(get.edgelist(dp_mst))
#colnames(edge_list) <- c(\"source\", \"target\")
#edge_df <- merge(ica_space_df, edge_list, by.x = \"sample_name\",by.y = \"source\", all = TRUE)
#edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = \"source_prin_graph_dim_1\",prin_graph_dim_2 = \"source_prin_graph_dim_2\"))
#edge_df <- merge(edge_df, ica_space_df[, c(\"sample_name\",\"prin_graph_dim_1\", \"prin_graph_dim_2\")],by.x = \"target\", by.y = \"sample_name\", all = TRUE)
#edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = \"target_prin_graph_dim_1\",prin_graph_dim_2 = \"target_prin_graph_dim_2\"))
#write.table(as.matrix(edge_df),file=\"$opt{'O'}/monocle_branchpoints.txt\",col.names=TRUE,row.names=FALSE,sep=\"\\t\",quote=FALSE)

p<-plot_cell_trajectory(agg_cds, color_by = \"timepoint\")
ggsave(plot=p,filename=\"$opt{'O'}.timepoint_plot.png\",width=5,height=4,dpi=900)
ggsave(plot=p,filename=\"$opt{'O'}.timepoint_plot.pdf\",width=5,height=4);


#this is to save the dims



plot_cell_trajectory2<-function (cds, x = 1, y = 2, color_by = \"State\", show_tree = TRUE, 
    show_backbone = TRUE, backbone_color = \"black\", markers = NULL, 
    use_color_gradient = FALSE, markers_linear = FALSE, show_cell_names = FALSE, 
    show_state_number = FALSE, cell_size = 1.5, cell_link_size = 0.75, 
    cell_name_size = 2, state_number_size = 2.9, show_branch_points = TRUE, 
    theta = 0, ...) 
{
    gene_short_name <- NA
    sample_name <- NA
    sample_state <- pData(cds)\$State
    data_dim_1 <- NA
    data_dim_2 <- NA
    lib_info_with_pseudo <- pData(cds)
    if (is.null(cds\@dim_reduce_type)) {
        stop(\"Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.\")
    }
    if (cds\@dim_reduce_type == \"ICA\") {
        reduced_dim_coords <- reducedDimS(cds)
    }
    else if (cds\@dim_reduce_type \%in\% c(\"simplePPT\", \"DDRTree\")) {
        reduced_dim_coords <- reducedDimK(cds)
    }
    else {
        stop(\"Error: unrecognized dimensionality reduction method.\")
    }
    ica_space_df <- data.frame(Matrix::t(reduced_dim_coords[c(x, 
        y), ]))
    colnames(ica_space_df) <- c(\"prin_graph_dim_1\", \"prin_graph_dim_2\")
    ica_space_df\$sample_name <- row.names(ica_space_df)
    ica_space_df\$sample_state <- row.names(ica_space_df)
    dp_mst <- minSpanningTree(cds)
    if (is.null(dp_mst)) {
        stop(\"You must first call orderCells() before using this function\")
    }
    edge_list <- as.data.frame(get.edgelist(dp_mst))
    colnames(edge_list) <- c(\"source\", \"target\")
    edge_df <- merge(ica_space_df, edge_list, by.x = \"sample_name\", 
        by.y = \"source\", all = TRUE)
    edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = \"source_prin_graph_dim_1\", 
        prin_graph_dim_2 = \"source_prin_graph_dim_2\"))
    edge_df <- merge(edge_df, ica_space_df[, c(\"sample_name\", 
        \"prin_graph_dim_1\", \"prin_graph_dim_2\")], by.x = \"target\", 
        by.y = \"sample_name\", all = TRUE)
    edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = \"target_prin_graph_dim_1\", 
        prin_graph_dim_2 = \"target_prin_graph_dim_2\"))
    S_matrix <- reducedDimS(cds)
    data_df <- data.frame(t(S_matrix[c(x, y), ]))
    data_df <- cbind(data_df, sample_state)
    colnames(data_df) <- c(\"data_dim_1\", \"data_dim_2\")
    data_df\$sample_name <- row.names(data_df)
    data_df <- merge(data_df, lib_info_with_pseudo, by.x = \"sample_name\", 
        by.y = \"row.names\")
    return_rotation_mat <- function(theta) {
        theta <- theta/180 * pi
        matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 
            nrow = 2)
    }
    tmp <- return_rotation_mat(theta) \%\*\% t(as.matrix(data_df[, 
        c(2, 3)]))
    data_df\$data_dim_1 <- tmp[1, ]
    data_df\$data_dim_2 <- tmp[2, ]
    tmp <- return_rotation_mat(theta = theta) \%\*\% t(as.matrix(edge_df[, 
        c(\"source_prin_graph_dim_1\", \"source_prin_graph_dim_2\")]))
    edge_df\$source_prin_graph_dim_1 <- tmp[1, ]
    edge_df\$source_prin_graph_dim_2 <- tmp[2, ]
    tmp <- return_rotation_mat(theta) \%\*\% t(as.matrix(edge_df[, 
        c(\"target_prin_graph_dim_1\", \"target_prin_graph_dim_2\")]))
    edge_df\$target_prin_graph_dim_1 <- tmp[1, ]
    edge_df\$target_prin_graph_dim_2 <- tmp[2, ]
    markers_exprs <- NULL
    if (is.null(markers) == FALSE) {
        markers_fData <- subset(fData(cds), gene_short_name %in% 
            markers)
        if (nrow(markers_fData) >= 1) {
            markers_exprs <- reshape2::melt(as.matrix(exprs(cds[row.names(markers_fData), 
                ])))
            colnames(markers_exprs)[1:2] <- c(\"feature_id\", \"cell_id\")
            markers_exprs <- merge(markers_exprs, markers_fData, 
                by.x = \"feature_id\", by.y = \"row.names\")
            markers_exprs\$feature_label <- as.character(markers_exprs\$gene_short_name)
            markers_exprs\$feature_label[is.na(markers_exprs\$feature_label)] <- markers_exprs\$Var1
        }
    }
    if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 
        0) {
        data_df <- merge(data_df, markers_exprs, by.x = \"sample_name\", 
            by.y = \"cell_id\")
        if (use_color_gradient) {
            if (markers_linear) {
                g <- ggplot(data = data_df, aes(x = data_dim_1, 
                  y = data_dim_2)) + geom_point(aes(color = value), 
                  size = I(cell_size), na.rm = TRUE) + scale_color_viridis(name = paste0(\"value\"), 
                  ...) + facet_wrap(~feature_label)
            }
            else {
                g <- ggplot(data = data_df, aes(x = data_dim_1, 
                  y = data_dim_2)) + geom_point(aes(color = log10(value + 
                  0.1)), size = I(cell_size), na.rm = TRUE) + 
                  scale_color_viridis(name = paste0(\"log10(value + 0.1)\"), 
                    ...) + facet_wrap(~feature_label)
            }
        }
        else {
            if (markers_linear) {
                g <- ggplot(data = data_df, aes(x = data_dim_1, 
                  y = data_dim_2, size = (value * 0.1))) + facet_wrap(~feature_label)
            }
            else {
                g <- ggplot(data = data_df, aes(x = data_dim_1, 
                  y = data_dim_2, size = log10(value + 0.1))) + 
                  facet_wrap(~feature_label)
            }
        }
    }
    else {
        g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
    }
    if (show_tree) {
        g <- g + geom_segment(aes_string(x = \"source_prin_graph_dim_1\", 
            y = \"source_prin_graph_dim_2\", xend = \"target_prin_graph_dim_1\", 
            yend = \"target_prin_graph_dim_2\"), size = cell_link_size, 
            linetype = \"solid\", na.rm = TRUE, data = edge_df)
    }
    if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 
        0) {
        if (use_color_gradient) {
        }
        else {
            g <- g + geom_point(aes_string(color = color_by), 
                na.rm = TRUE)
        }
    }
    else {
        if (use_color_gradient) {
        }
        else {
            g <- g + geom_point(aes_string(color = color_by), 
                size = I(cell_size), na.rm = TRUE)
        }
    }
    if (show_branch_points && cds\@dim_reduce_type == \"DDRTree\") {
        mst_branch_nodes <- cds\@auxOrderingData[[cds\@dim_reduce_type]]\$branch_points
        branch_point_df <- subset(edge_df, sample_name \%in\% mst_branch_nodes)[,c(\"sample_name\", \"source_prin_graph_dim_1\", \"source_prin_graph_dim_2\")]
        branch_point_df\$branch_point_idx <- match(branch_point_df\$sample_name,mst_branch_nodes)
        branch_point_df <- branch_point_df[!duplicated(branch_point_df\$branch_point_idx),]
        g <- g + geom_point(aes_string(x = \"source_prin_graph_dim_1\", 
            y = \"source_prin_graph_dim_2\"), size = 5, na.rm = TRUE, 
            data = branch_point_df) + geom_text(aes_string(x = \"source_prin_graph_dim_1\", 
            y = \"source_prin_graph_dim_2\", label = \"branch_point_idx\"), 
            size = 4, color = \"white\", na.rm = TRUE, data = branch_point_df)
    }
    if (show_cell_names) {
        g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
    }
    if (show_state_number) {
        g <- g + geom_text(aes(label = sample_state), size = state_number_size)
    }
    g <- g + xlab(paste(\"Component\", x)) + 
        ylab(paste(\"Component\", y)) + theme(legend.position = \"top\", 
        legend.key.height = grid::unit(0.35, \"in\")) + theme(legend.key = element_blank()) + 
        theme(panel.background = element_rect(fill = \"white\"))
    g
    write.table(data_df,file=\"save.dims.txt\",col.names=T,row.names=F,quote=F,sep=\"\\t\")
    write.table(branch_point_df,file=\"save.bp.txt\",col.names=T,row.names=F,quote=F,sep=\"\\t\")
    write.table(edge_df,file=\"save.edgelist.txt\",col.names=T,row.names=F,quote=F,sep=\"\\t\")
    writoutdims<-data.frame(annot=data_df\$sample_name,data_dim1=data_df\$data_dim1,data_dim2=data_df\$data_dim2)
    write.table(writoutdims,file=\"ddrtee.dims\",col.names=F,row.names=F,quote=F,sep=\"\\t\")
}





#add stuff here if you want to reroot agg_cds <- orderCells(agg_cds, root_state = \"D\")

p<-plot_cell_trajectory(agg_cds, color_by = \"Pseudotime\")
ggsave(plot=p,filename=\"$opt{'O'}.lambda_plot.png\",width=5,height=4,dpi=900)
ggsave(plot=p,filename=\"$opt{'O'}.lambda_plot.pdf\",width=5,height=4);

pData(input_cds)\$Pseudotime <- pData(agg_cds)[colnames(input_cds),]\$Pseudotime
pData(input_cds)\$State <- pData(agg_cds)[colnames(input_cds),]\$State
#write out
plot_cell_trajectory2(agg_cds)

# writing out pseudotime etc so you can recreate everything
write.table(as.matrix(pData(agg_cds)),file=\"./$opt{'O'}_ddrt_aggragated_cells.txt\",col.names=TRUE,row.names=TRUE,sep=\"\\t\",quote=FALSE)
write.table(as.matrix(Biobase::exprs(agg_cds)),file=\"./$opt{'O'}_ddrt_aggragated_cells_norm_count.txt\",col.names=TRUE,row.names=TRUE,sep=\"\\t\",quote=FALSE)
write.table(as.matrix(fData(agg_cds)),file=\"./$opt{'O'}_ddrt_aggragated_features.txt\",col.names=TRUE,row.names=FALSE,sep=\"\\t\",quote=FALSE)
";






close R;

system("$Rscript $opt{'O'}.ddrt.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.ddrt.r");
}
}
1;
