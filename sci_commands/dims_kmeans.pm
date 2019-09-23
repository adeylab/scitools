package sci_commands::dims_kmeans;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("dims_kmeans");

sub dims_kmeans {

@ARGV = @_;
# Defaults
$range_default = "1-15";
$xdim = 1;
$ydim = 2;

getopts("O:SXR:K:D:s:P:x:y:pt", \%opt);

$die2 = "
scitools dims-kmeans [options] [input dims]
   or    kmeans

Performs kmeans clustering on a dims file or matrix file.
If using a matrix file, use the -t flag.
User may want to first run kmeans clustering with the -S flag.
  This will generate a silhouette plot for a more empirical expected cluster count.
They may then rerun this script without the -S flag, and the -K flag specified.

Options:
   -K   [INT]   K-value (number of clusters, Required) if \"S\" defined script does silhouette analysis with the K value as max k value 
   -S           If Flagged, script does silhouette analysis
   -O   [STR]   Output prefix (default is [input].Kmeans_K[K].annot)
   -D   [RNG]   Range of dimensions to include (range format)
                (e.g. 1-5,6,8; def = $range_default)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)
   -p           Plot the output using specified dims file
   -P   [DIMS]  Plot the resulting Kmeans clustering using this
                specified dims file rather than argument dims file.
   -x   [INT]   X-dimension to plot (def = $xdim)
   -y   [INT]   Y-dimension to plot (def = $ydim)
   -s   [STR]   scitools call (def = $scitools)
   -t           Transpose matrix before KMEANS
                (flag if [input dims] is a cistopic or irlba matrix)
";

if (!defined $ARGV[0] | !defined $opt{'K'}) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'x'}) {$xdim = $opt{'x'}};
if (defined $opt{'y'}) {$ydim = $opt{'y'}};
if (!defined $opt{'D'}) {$opt{'D'} = $range_default};
if (defined $opt{'p'} && !defined $opt{'P'}) {$opt{'P'} = $ARGV[0]};
read_ranges($opt{'D'});

if (defined !$opt{'t'}) {
open IN, "$ARGV[0]";
$dim_line = <IN>; close IN;
chomp $dim_line; @DL = split(/\t/, $dim_line);
if (@DL < $RANGE_VALUES[@RANGE_VALUES-1]) {
        die "ERROR: The dimension ranges (max = $RANGE_VALUES[@RANGE_VALUES-1]) specified are beyond the number of dimensions in $ARGV[0] (".@DL.")!\n";
}; 
} else {
$dim_line=`wc -l < $ARGV[0]`;
chomp $dim_line;
if ($dim_line < $RANGE_VALUES[@RANGE_VALUES-1]) {
  die "ERROR: The dimension ranges (max = $RANGE_VALUES[@RANGE_VALUES-1]) specified are beyond the number of dimensions in $ARGV[0] ($dim_line)!\n";
};
};

if (defined $opt{'S'})
{
  print "Silhouette analysis\n";
#if silhouette analysis
open R, ">$opt{'O'}.silhouette_maxK$opt{'K'}.r";
print R "
library(\"cluster\")";

if (defined $opt{'t'}) {
print R "
df<-as.matrix(read.table(\"$ARGV[0]\",row.names=1)[$range_R_set,])
df<-as.data.frame(t(df))
";
} else {
print R "
df<-as.matrix(read.table(\"$ARGV[0]\",row.names=1)[,$range_R_set])
";
};

print R "
silhouette_score <- function(k){
  km <- kmeans(df, centers = k, nstart=50)
  ss <- silhouette(km\$cluster, dist(df))
  mean(ss[, 3])
}
k <- 2:$opt{'K'}
avg_sil <- sapply(k, silhouette_score)
png(
  \"$opt{'O'}.silhouette_maxK$opt{'K'}.png\",
  width     = 3.25,
  height    = 3.25,
  units     = \"in\",
  res       = 1200,
  pointsize = 4
)
par(
  mar      = c(5, 5, 2, 2),
  xaxs     = \"i\",
  yaxs     = \"i\",
  cex.axis = 2,
  cex.lab  = 2
)
plot(k, type=\'b\', avg_sil, xlab=\'Number of clusters\', ylab=\'Average Silhouette Scores\', frame=FALSE)
dev.off()

"; close R;

system("$Rscript $opt{'O'}.silhouette_maxK$opt{'K'}.r");

if (!defined $opt{'X'}) {
  system("rm -f $opt{'O'}.silhouette_maxK$opt{'K'}.r");}
} else {
print "Kmeans analysis\n";
open R, ">$opt{'O'}.Kmeans_K$opt{'K'}.r";

if (defined $opt{'t'}) {
print R "
D<-as.matrix(read.table(\"$ARGV[0]\",row.names=1)[$range_R_set,])
D<-as.data.frame(t(D))
";
} else {
print R "
D<-as.matrix(read.table(\"$ARGV[0]\",row.names=1)[,$range_R_set])
";
};

print R "
K<-kmeans(D,$opt{'K'})
ANN<-as.matrix(K\$cluster)
ANN[,1]<-sub(\"\^\", \"K\", ANN[,1])
write.table(ANN,file=\"$opt{'O'}.Kmeans_K$opt{'K'}.annot\",col.names=FALSE,row.names=TRUE,quote=FALSE,sep=\"\\t\")
"; close R;

system("$Rscript $opt{'O'}.Kmeans_K$opt{'K'}.r");

if (defined $opt{'P'}) {
	print STDERR "PLOTTING: $scitools plot-dims -O $opt{'O'}.Kmeans_K$opt{'K'} -A $opt{'O'}.Kmeans_K$opt{'K'}.annot -x $xdim -y $ydim $opt{'P'}\n";
	system("$scitools plot-dims -O $opt{'O'}.Kmeans_K$opt{'K'} -A $opt{'O'}.Kmeans_K$opt{'K'}.annot -x $xdim -y $ydim $opt{'P'}");
}

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.Kmeans_K$opt{'K'}.r");
}
}
}
1;
