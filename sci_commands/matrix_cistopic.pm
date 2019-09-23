package sci_commands::matrix_cistopic;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_cistopic");

sub matrix_cistopic {

@ARGV = @_;
getopts("O:A:c:C:S:G:r:n:T:XR:P:L:W:", \%opt);

$thrP = 0.975;

$die2 = "
scitools matrix-cistopic [options] [input matrix (or rds for topic bed output only)]
   or    cistopic

[input matrix] =	A counts matrix generated from atac_count. 
					Can be filtered prior to this command call, or 
					can supply the filtering options during cisTopic processing.

cisTopic serves as an alternative textmining algorithm to TFIDF-LSI processing.
It is to be run on a sciATAC counts matrix. For more information see:
github.com/aertslab/cisTopic/

Outputs a matrix file similar to matrix-irlba function call. To be processed through
[matrix_tsne|matrix_umap|matrix_PCA|matrix_SWNE]

cisTopic consists of 4 main steps: 
(1) generation of a binary accessibility matrix as input for LDA; 
(2) LDA and model selection; 
(3) cell state identification using the topic-cell distributions from LDA and 
(4) exploration of the region-topic distributions.

Options:
   -O   [STR]  Output prefix (default is [input].cistopic.dims)
   -c   [INT]  Number of nonZero sites per column (cell) to retain (def = 1)
   -r   [INT]  Number of nonZero sites per row (peak) to retain (def = 1)
   -n   [INT] 	Number of cores for parallel processing. (def=1)
   -A   [STR]  Annotation file is useful. cisTopic will provide a PCA with the influence of and the heatmap of Topics
   -C   [STR]  color file 
   -S   [STR]  color string
   -L   [STR]  Columns file for sparseMatrix (will try to auto-detect), currently unsupported
   -W   [STR]  Rows file for sparseMatrix (will try to auto-detect), currently unsupported
   -T   [INT]  User defined number of Topics to use. 
                  If unspecified: will generate 15,20,25,30,50,65,100 Topics,
                  and use log-liklihood estimators to select the best.
                  Specification can be a single number of a comma separated list.
                  Will use a core for each number supplied (DO NOT EXCEED A LIST LENGTH OF 10)
   -P   [FLT]  ThrP (def = $thrP)
   -G   [STR]  Genome and genes to use to annotate regions	
   -X          Retain intermediate files (def = delete)
   -R   [STR]  Rscript call (def = $Rscript)

Note: Requires cisTopic R package

";

if (!defined $ARGV[0]) {die $die2};

$prefix = $ARGV[0];
$prefix =~ s/\.gz$//;
$prefix =~ s/\.matrix$//;
$prefix =~ s/\.values$//;
$prefix =~ s/\.sparseMatrix$//;
if (!defined $opt{'O'}) {
	$opt{'O'} = $prefix;
}

if ($ARGV[0] =~ /sparseMatrix/i) {
	$sparse = 1;
	if (defined $opt{'L'}) {
		$col_file = $opt{'L'};
	} else {
		if (-e "$prefix.sparseMatrix.cols") {
			$col_file = "$prefix.sparseMatrix.cols";
		} elsif (-e "$prefix.sparseMatrix.cols.gz") {
			$col_file = "$prefix.sparseMatrix.cols.gz";
		} else {
			die "ERROR: Cannot detect cols file (e.g. $prefix.sparseMatrix.cols), please provide as -C\n";
		}
	}
	if (defined $opt{'W'}) {
		$row_file = $opt{'W'};
	} else {
		if (-e "$prefix.sparseMatrix.rows") {
			$row_file = "$prefix.sparseMatrix.rows";
		} elsif (-e "$prefix.sparseMatrix.rows.gz") {
			$row_file = "$prefix.sparseMatrix.rows.gz";
		} else {
			die "ERROR: Cannot detect rows file (e.g. $prefix.sparseMatrix.rows), please provide as -C\n";
		}
	}
} else {$sparse = 0};
if (!defined $opt{'c'}) {$opt{'c'} = 1};
if (!defined $opt{'r'}) {$opt{'r'} = 1};
if (!defined $opt{'n'}) {$opt{'n'} = 1};
if (defined $opt{'P'}) {$thrP = $opt{'P'}};
if (!defined $opt{'T'}) {$opt{'T'} = "15,20,25,30,50,65,100"};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'C'}) {read_color_file($opt{'C'})};
if (defined $opt{'S'}) {read_color_string($opt{'S'})};
if (!defined $opt{'G'}) {$opt{'G'} = "hg38"};




if ($ARGV[0] =~ /\.rds$/i) {

if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.rds$//i};

open R, ">$opt{'O'}.output_topics.r";
print R "
library(plyr)
library(cisTopic)

readRDS(\"$ARGV[0]\")
cisTopicObject <- getRegionsScores(cisTopicObject, method='Z-score', scale=TRUE)
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=$thrP, plot=FALSE)
getBedFiles(cisTopicObject, path='$opt{'O'}.topics')
";

close R;

system("$Rscript $opt{'O'}.output_topics.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.output_topics.r");
}

} else {

open R, ">$opt{'O'}.cistopic.r";
print R "
library(plyr)
library(cisTopic)
library(Matrix)
";

if ($ARGV[0] =~ /sparseMatrix/i) {
	print R "IN<-as.matrix(read.table(\"$ARGV[0]\"))
IN<-sparseMatrix(i=IN[,1],j=IN[,2],x=IN[,3])
COLS<-read.table(\"$col_file\")
colnames(IN)<-COLS\$V1
ROWS<-read.table(\"$row_file\")
row.names(IN)<-ROWS\$V1\n";
} else {
	print R "IN<-as.matrix(read.table(\"$ARGV[0]\"))\n";
}

print R "
row.names(IN)<-sub(\"_\",\"-\",sub(\"_\",\":\",row.names(IN)))

#Set up filtered binarized counts matrix
cisTopicObject <- createcisTopicObject(IN,min.cells=$opt{'r'},min.regions=$opt{'c'}, keepCountsMatrix=FALSE)
";
if (defined $opt{'A'})
{
print R "
annot<-read.table(file=\"$opt{'A'}\",header=F)
#match rows of annot
annot<-annot[match(colnames(IN), annot\$V1),]
row.names(annot)<-annot\$V1
names(annot)<-c(\"cellname\",\"LineType\")
cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = annot)
";

}



print R "
cisTopicObject <- runModels(cisTopicObject, topic=c($opt{'T'}), seed=2018, nCores=$opt{'n'}, burnin = 250, iterations = 300)

if (length(c($opt{'T'}))>1){
#Future update:Plot model log likelihood (P(D|T)) at the last iteration
pdf(file=\"$opt{'O'}.cistopic.modelselection.pdf\")
cisTopicObject <- selectModel(cisTopicObject)
logLikelihoodByIter(cisTopicObject)
dev.off()  
modelMat <- scale(cisTopicObject\@selected.model\$document_expects, center = TRUE, scale = TRUE)
} else {
modelMat<-scale(cisTopicObject\@models\$document_expects,center=TRUE,scale=TRUE)
}

#Print out cisTopic Matrix#
tModelmat<-as.data.frame(t(modelMat))
Modeldf<-as.data.frame(modelMat)
rownames(tModelmat)<-cisTopicObject\@cell.names
colnames(Modeldf)<-cisTopicObject\@cell.names
row.names(Modeldf)<-paste0(\"Topic_\",row.names(Modeldf))
write.table(Modeldf,file=\"$opt{'O'}.cistopic.matrix\",col.names=T,row.names=T,quote=F,sep=\"\\t\")
saveRDS(cisTopicObject,\"$opt{'O'}.cistopicObject.rds\")
#adding part where the contribution matrix is calculated, we use a binarization method to select for peaks that contribute to each topic
cisTopicObject <- getRegionsScores(cisTopicObject, method='Z-score', scale=TRUE)
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=$thrP, plot=FALSE)
getBedFiles(cisTopicObject, path='$opt{'O'}.topics')

saveRDS(cisTopicObject,\"$opt{'O'}.cistopicObject.rds\")
";

if (defined $opt{'A'}){
   print R "
   #temporarily disabled
   #cisTopicObject <- runPCA(cisTopicObject)
   #coordinates <- cisTopicObject\@dr[[\'PCA\']]\$ind.coord
   #write.table(coordinates,file=\"$opt{'O'}.PCA.internal.dims\",col.names=T,row.names=T,quote=F,sep=\"\\t\")
   #png(file=\"$opt{'O'}.PCA_cistopic.png\",width=12,height=12,units=\"in\",res=300)
   #plotCellStates(cisTopicObject, method=\'Biplot\', topic_contr=\'Z-score\',topics=\'all\', colorBy=c(\'LineType\'))
   #dev.off()
   #pdf(file=\"$opt{'O'}.PCA_cistopic.pdf\",width=12,height=12)
   #plotCellStates(cisTopicObject, method=\'Biplot\', topic_contr=\'Z-score\', topics=\'all\', colorBy=c(\'LineType\'))
   #dev.off()
   ";
    if (defined $opt{'S'} || defined $opt{'C'})  
   {
   print R "
   color_ch<-list(Type=c($color_mapping))
   ha_col<-HeatmapAnnotation(Type = annot\$LineType,col=color_ch)
   ";
   }
   else
   {
   print R "
   ha_col<-HeatmapAnnotation(Type = annot\$LineType)
   ";
   }
print R "
   png(\"$opt{'O'}.Heatmap_prob_cistopic.png\",width=12,height=12,units=\"in\",res=600)
   cellTopicHeatmap(cisTopicObject, method=\'Probability\',bottom_annotation=ha_col)
   dev.off()
   pdf(\"$opt{'O'}.Heatmap_prob_cistopic.pdf\",width=12,height=12)
   cellTopicHeatmap(cisTopicObject, method=\'Probability\',bottom_annotation=ha_col)
   dev.off()
   png(\"$opt{'O'}.Heatmap_zscore_cistopic.png\",width=12,height=12,units=\"in\",res=600)
   cellTopicHeatmap(cisTopicObject, method=\'Z-score\',bottom_annotation=ha_col)
   dev.off()
   pdf(\"$opt{'O'}.Heatmap_zscore_cistopic.pdf\",width=12,height=12)
   cellTopicHeatmap(cisTopicObject, method=\'Z-score\',bottom_annotation=ha_col)
   dev.off()
   ";
}
else 
{
print R "
   png(\"$opt{'O'}.Heatmap_prob_cistopic.png\",width=12,height=12,units=\"in\",res=600)
   cellTopicHeatmap(cisTopicObject, method=\'Probability\')
   dev.off()
   pdf(\"$opt{'O'}.Heatmap_prob_cistopic.pdf\",width=12,height=12)
   cellTopicHeatmap(cisTopicObject, method=\'Probability\')
   dev.off()
   png(\"$opt{'O'}.Heatmap_zscore_cistopic.png\",width=12,height=12,units=\"in\",res=600)
   cellTopicHeatmap(cisTopicObject, method=\'Z-score\')
   dev.off()
   pdf(\"$opt{'O'}.Heatmap_zscore_cistopic.pdf\",width=12,height=12)
   cellTopicHeatmap(cisTopicObject, method=\'Z-score\')
   dev.off()
";

}
if ($opt{'G'} eq "hg38") {
print R"
#adding annotation to regions
library(\"org.Hs.eg.db\")
library(\"TxDb.Hsapiens.UCSC.hg38.knownGene\")
cisTopicObject <- annotateRegions(cisTopicObject, txdb=TxDb.Hsapiens.UCSC.hg38.knownGene,annoDb='org.Hs.eg.db')
#will add more
#signaturesHeatmap(cisTopicObject, selected.signatures = 'annotation')
write.table(cisTopicObject\@region.data,file=\"$opt{'O'}.annotationregionsandscores.matrix\",col.names=T,row.names=T,quote=F,sep=\"\\t\")
saveRDS(cisTopicObject,\"$opt{'O'}.cistopicObject.rds\")
";
}
else
{
#other genomes coming later
}


close R;

system("$Rscript $opt{'O'}.cistopic.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.cistopic.r");
}

}

}
1;
