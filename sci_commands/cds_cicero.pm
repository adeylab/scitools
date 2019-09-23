package sci_commands::cds_cicero;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("cds_cicero");

sub cds_cicero {

@ARGV = @_;
# Defaults

getopts("O:R:G:k:P:p:X", \%opt);

$die2 = "
scitools cds-cicero [options] [directory containing cds files]
   or    cicero-cds

Prior to running, ensure you have ran scitools matrix-makecds. 
That function will convert matrix files into the CDS format that Cicero requires.
Runs the current version of Cicero within a directory containing files for CDS format input. 
Outputs:
Tables describing
 the coaccessibility of peaks:                          cicero_conns.txt
 the membership peaks within CCANs:                     cicero.CCANS.txt
RDS files for
  the coaccessibility of peaks:                         cicero.conns.rds
  the gene activity score matrix (in dgC matrix):       cicero.gene_activity_matrix.rds
  the working cds file (after aggregation):             cicero.CDS.rds


Options:
   -O    [STR]   Output Directory (default is [current working directory]/cicero_output)
   -R    [STR]   Rscript call (def = $Rscript)
   -G    [STR]   Genome to be used. Must be of list:
                  [human.hg38.genome,mouse.mm10.genome]
                  Default: human.hg38.genome
   -P 	 [INT] 	Distance (in bp) between Peaks for aggregation. (Default: 10000) 
   -p     [INT]   Give CCAN values to print seperatefd by commas (e.g.: 1,2,3,4), if given value \"ALL\" all will be printed out      
   -k 	 [INT] 	Number of cells to aggregate per bin. (Default: 50) 	
   -X           Retain intermediate files (Default = delete)
                  
";


if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = "$ARGV[0]/cicero_output"};
if (!defined $opt{'G'}) {$opt{'G'} = "human.hg38.genome"};
if (!defined $opt{'P'}) {$opt{'P'} = "10000"};
if (!defined $opt{'k'}) {$opt{'k'} = "50"};

system("mkdir $opt{'O'}");

open R, ">$opt{'O'}/cicero.r";

print R "
suppressWarnings(library(cicero))
message(\"Loading Cicero\")
# reading in matrix, annotation, and dimension data

cds_cell_data <- read.delim(\"$ARGV[0]/cds_cell_data.txt\")
cds_dims_data <- read.delim(\"$ARGV[0]/cds_dims_data.txt\",header=F)
cds_site_data <- read.delim(\"$ARGV[0]/cds_site_data.txt\")
cds_counts_matrix <- read.table(\"$ARGV[0]/cds_counts_matrix.txt\")
message(\"Read in CDS Files.\")

# generating cell data set
feature_data <- new(\"AnnotatedDataFrame\", data = cds_site_data)
sample_data <- new(\"AnnotatedDataFrame\", data = cds_cell_data)

if(ncol(cds_dims_data)==4){
colnames(cds_dims_data)<-c(\"cellID\",\"Dimension_1\",\"Dimension_2\",\"Dimension_3\")
row.names(cds_dims_data)<-cds_dims_data\$cellID
dimension_reduction<-cds_dims_data[2:ncol(cds_dims_data)]
}else{
colnames(cds_dims_data)<-c(\"cellID\",\"Dimension_1\",\"Dimension_2\")
row.names(cds_dims_data)<-cds_dims_data\$cellID
dimension_reduction<-cds_dims_data[2:ncol(cds_dims_data)]
}

message(\"Setting up CDS matrix, binarized for Cicero.\")


cds <- suppressWarnings(newCellDataSet(as.matrix(cds_counts_matrix), phenoData = sample_data, featureData = feature_data,expressionFamily = VGAM::negbinomial.size(),lowerDetectionLimit = 0))

cds\@expressionFamily <- VGAM::binomialff()
cds\@expressionFamily\@vfamily <- \"binomialff\"

pData(cds)\$temp <- NULL
fData(cds)\$chr <- as.numeric(as.character(fData(cds)\$chr))
fData(cds)\$bp1 <- as.numeric(as.character(fData(cds)\$bp1))
fData(cds)\$bp2 <- as.numeric(as.character(fData(cds)\$bp2))
cds <- cds[order(fData(cds)\$chr, fData(cds)\$bp1),]

#Running necessary preprocessing steps
set.seed(2017)
pData(cds)\$cells <- NULL
cds <- aggregate_nearby_peaks(cds, distance = $opt{'P'})
cds <- detectGenes(cds)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

#Spliting cds to a pre-aggregation data set
input_cds<-cds

message(\"Using Supplied Dims File.\")
cds <- make_cicero_cds(cds, reduced_coordinates = dimension_reduction, k=$opt{'k'})

message(\"Running Cicero\")
";

if ($opt{'G'} eq "human.hg38.genome") {
print R "

human.hg38.genome<-read.table(\"/home/groups/oroaklab/refs/hg38/hg38.bedtools.genome\",header=F)
conns <- run_cicero(cds, human.hg38.genome)
gene_annotation<-read.table(\"/home/groups/oroaklab/refs/hg38/gencode.v29.Gene.annotation.bed\",header=F)
colnames(gene_annotation)<-c(\"chr\",\"start\",\"end\",\"gene\")
";
} else {
print R "

data(\"$opt{'G'}\")
conns <- run_cicero(cds, $opt{'G'})
gene_annotation<-read.table(\"/home/groups/oroaklab/refs/mm10/mm10.RefGene.names.bed\",header=F)
colnames(gene_annotation)<-c(\"chr\",\"start\",\"end\",\"gene\")
";
}

print R "
message(\"Sample Cicero Output:\")
head(conns)

message(\"Generating CCANS.\")
CCAN_assigns <- generate_ccans(conns)

message(\"Sample CCANs Output:\")
head(CCAN_assigns)

saveRDS(cds,\"$opt{'O'}/cicero.CDS.rds\")
write.table(as.data.frame(conns), file=\"$opt{'O'}/cicero.output.txt\",quote=F,sep=\"\\t\",row.names=F,col.names=T)
write.table(as.data.frame(CCAN_assigns), file=\"$opt{'O'}/cicero.CCANS.txt\",quote=F,sep=\"\\t\",row.names=F,col.names=F)

message(\"Sample CCANs Output:\")
input_cds <- annotate_cds_by_site(input_cds, gene_annotation)
# generate unnormalized gene activity matrix
unnorm_ga <- build_gene_activity_matrix(input_cds, conns)

#set up peak count per cell
num_genes<-pData(input_cds)$num_genes_expressed
names(num_genes) <- row.names(pData(input_cds))

cicero_gene_activities<-normalize_gene_activities(unnorm_ga,num_genes)
write.table(as.data.frame(as.matrix(cicero_gene_activities)),file=\"$opt{'O'}/cicero.gene_activity_matrix.txt\",row.names=T,col.names=T,sep=\"\\t\",quote=F)

";

close R;
system("$Rscript $opt{'O'}/cicero.r");

if (defined $opt{'p'}) {

   #this is the part that that does the plotting when defined
  if($opt{'p'} eq "ALL")
  {
   open IN, "$opt{'O'}/cicero.CCANS.txt";
   while ($l = <IN>) {
      chomp $l;
      @P = split(/\t/, $l);
      $CCAN_include{$P[1]}++;
  } close IN; } else {
   chomp $opt{'p'};
   @P = split(/,/, $opt{'p'});
   foreach $CCAN (@P){
      $CCAN_include{$CCAN}++;
   }
  }
  if ($opt{'G'} eq "human.hg38.genome"){$plot_genome="hg38"}elsif($opt{'G'} eq "human.hg19.genome"){$plot_genome="hg19"}elsif($opt{'G'} eq "mouse.mm10.genome"){$plot_genome="mm10"};
   foreach $CCAN_incl (sort keys %CCAN_include)
   {
      system("awk \'{if(\$2==$CCAN_incl){print}}\' $opt{'O'}/cicero.CCANS.txt > $opt{'O'}/cicero.CCANS_$CCAN_incl.tmp");
      system("scitools cicero_plot -G $plot_genome $opt{'O'} $opt{'O'}/cicero.CCANS_$CCAN_incl.tmp");
   }

}

if (!defined $opt{'X'}) {
    system("rm -f $opt{'O'}/cicero.r");
    system("rm -f $opt{'O'}/*.tmp");
}
}
1;
