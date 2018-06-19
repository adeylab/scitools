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

getopts("O:XR:K:D:s:P:x:y:p", \%opt);

$die2 = "
scitools dims-kmeans [options] [input dims]
   or    kmeans

Performs kmeans clustering on a dims file

Options:
   -K   [INT]   K-value (number of clusters, Required)
   -O   [STR]   Output prefix (default is [input].Kmeans_K[K].annot)
   -D   [RNG]   Range of dimensions to include (range format)
                (e.g. 1-5,6,8; def = $range_default)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)
   -p           Plot the output using specified dims file
   -P   [DIMS]  Plot the resulting Kmeans clustering using the
                specified dims file (for more plot options use
                scitools plot-dims)
   -x   [INT]   X-dimension to plot (def = $xdim)
   -y   [INT]   Y-dimension to plot (def = $ydim)
   -s   [STR]   scitools call (def = $scitools)

";

if (!defined $ARGV[0] | !defined $opt{'K'}) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'x'}) {$xdim = $opt{'x'}};
if (defined $opt{'y'}) {$ydim = $opt{'y'}};
if (!defined $opt{'D'}) {$opt{'D'} = $range_default};
if (defined $opt{'p'} && !defined $opt{'P'}) {$opt{'P'} = $ARGV[0]};
read_ranges($opt{'D'});

open IN, "$ARGV[0]";
$dim_line = <IN>; close IN;
chomp $dim_line; @DL = split(/\t/, $dim_line);
if (@DL < $RANGE_VALUES[@RANGE_VALUES-1]) {
	die "ERROR: The dimension ranges (max = $RANGE_VALUES[@RANGE_VALUES-1]) specified are beyond the number of dimensions in $ARGV[0] (".@DL.")!\n";
}

open R, ">$opt{'O'}.Kmeans_K$opt{'K'}.r";
print R "
D<-read.table(\"$ARGV[0]\",row.names=1)[,$range_R_set]
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
1;
