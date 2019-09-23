package sci_commands::dims_dbscan;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("dims_dbscan");

sub dims_dbscan {

@ARGV = @_;
# Defaults
$range_default = "1-2";
$xdim = 1;
$ydim = 2;
$epsilon = 2;
$minPts = 10;

getopts("O:XR:D:s:P:x:y:e:m:pt", \%opt);

$die2 = "
scitools dims-dbscan [options] [input dims]
   or    dbscan

Performs dbscan clustering on a dims file

Options:
   -O   [STR]   Output prefix (default is [input].dbscan.annot)
   -e   [FLT]   Epsilon value, akin to distance (def = $epsilon)
   -m   [INT]   Minimum points to seed a cluster (def = $minPts)
   -D   [RNG]   Range of dimensions to include (range format)
                (e.g. 1-5,6,8; def = $range_default)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)
   -p           Plot on the dims file (-P will toggle)
   -P   [DIMS]  Plot the resulting dbscan clustering using the
                specified dims file (for more plot options use
                scitools plot-dims)
   -x   [INT]   X-dimension to plot (def = $xdim)
   -y   [INT]   Y-dimension to plot (def = $ydim)
   -s   [STR]   scitools call (def = $scitools)
   -t           Transpose matrix before DBSCAN
                (flag if [input dims] is a cistopic or irlba matrix)

Requires the dbscan R package.

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'x'}) {$xdim = $opt{'x'}};
if (defined $opt{'y'}) {$ydim = $opt{'y'}};
if (defined $opt{'e'}) {$epsilon = $opt{'e'}};
if (defined $opt{'m'}) {$minPts = $opt{'m'}};
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

open R, ">$opt{'O'}.dbscan.r";

if (defined $opt{'t'}) {
print R "
library(dbscan)
DIMS<-as.matrix(read.table(\"$ARGV[0]\",row.names=1)[$range_R_set,])
DIMS<-as.data.frame(t(DIMS))
";
} else {
print R "
library(dbscan)
DIMS<-as.matrix(read.table(\"$ARGV[0]\",row.names=1)[,$range_R_set])
";
};

print R "
DIST<-dist(DIMS)
DBS<-dbscan(DIST,$epsilon,minPts = $minPts)
ANN<-as.matrix(DBS\$cluster)
rownames(ANN)<-rownames(DIMS)
ANN[,1]<-sub(\"\^\", \"D\", ANN[,1])
write.table(ANN,file=\"$opt{'O'}.dbscan.annot\",col.names=FALSE,row.names=TRUE,sep=\"\\t\",quote=FALSE)
"; close R;

system("$Rscript $opt{'O'}.dbscan.r");

if (defined $opt{'P'}) {
        print STDERR "PLOTTING: $scitools plot-dims -O $opt{'O'}.dbscan -A $opt{'O'}.dbscan.annot -x $xdim -y $ydim $opt{'P'}\n";
        system("$scitools plot-dims -O $opt{'O'}.dbscan -A $opt{'O'}.dbscan.annot -x $xdim -y $ydim $opt{'P'}");
}

if (!defined $opt{'X'}) {
        system("rm -f $opt{'O'}.dbscan.r");
}

}
1;
