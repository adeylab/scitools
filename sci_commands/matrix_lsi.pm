package sci_commands::matrix_lsi;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_lsi");

sub matrix_lsi {

@ARGV = @_;
# Defaults
$range_default = "1-15";
$block1_size = 100000;

getopts("O:D:d:XR:", \%opt);

$die2 = "
scitools matrix-lsi [options] [tfidf matrix]
   or    lsi

Options:
   -O   [STR]   Output prefix (def = [input].tfidf.LSI.matrix)
   -D   [RNG]   Range of dimensions to include (range format)
                 (e.g. 1-5,6,8; def = $range_default)
   -X           Retain intermediate files (def = delete)
   -B           Big matrix (toggle if matrix exceeds R limit of 2^31 elements)
   -N   [INT]   Number of rows to output a sfirst block (~half the matrix
                 should not exceed ~120k for a ~15k cell matrix, second block
                 is the rest of the row; def = $block1_size)
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires 'svd' R package

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
if (!defined $opt{'D'}) {$opt{'D'} = $range_default};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'N'}) {$block1_size = $opt{'N'}};
read_ranges($opt{'D'});

open PARAMS, ">$opt{'O'}.LSI.params";
print PARAMS "scitools matrix-lsi ON $ARGV[0]
Dimensions: $range_R_set\n";
close PARAMS;

open R, ">$opt{'O'}.LSI.r";
print R "
library(svd)

IN<-read.table(\"$ARGV[0]\")

SVD<-svd(as.matrix(IN))

D<-diag(SVD\$d[$range_R_set])
U<-as.matrix(SVD\$u[,$range_R_set])
V<-as.matrix(SVD\$v[,$range_R_set])

options(digits=6)
LSI<-t(scale(V %*% D %*% t(U)))
colnames(LSI)<-colnames(IN)
rownames(LSI)<-rownames(IN)
";

if (!defined $opt{'B'}) {
print R "
write.table(LSI,file=\"$opt{'O'}.LSI.matrix\",sep=\"\\t\",row.names=TRUE,col.names=TRUE,quote=FALSE)
";
close R;
system("$Rscript $opt{'O'}.LSI.r");
} else {
print R "
write.table(head(LSI,$block1_size),file=\"$opt{'O'}.LSI.matrix\",sep=\"\\t\",row.names=TRUE,col.names=TRUE,quote=FALSE)
write.table(tail(LSI,(nrow(LSI)-$block1_size)),file=\"$opt{'O'}.LSI.2.matrix\",sep=\"\\t\",row.names=TRUE,col.names=FALSE,quote=FALSE)
";
close R;
system("$Rscript $opt{'O'}.LSI.r");
system("cat $opt{'O'}.LSI.2.matrix >> $opt{'O'}.LSI.matrix && rm -f $opt{'O'}.LSI.2.matrix");
}

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.LSI.r");
}

}
1;
