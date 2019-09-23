package sci_commands::counts2tsne;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("counts2tsne");

sub counts2tsne {

@ARGV = @_;
# Deafults
$irlba_dims = 50;
$tsne_dims = 2;
$perp = 30;

getopts("O:XD:R:P:d:", \%opt);

$die2 = "
scitools counts2tsne [options] [input counts matrix]

Takes in a counts matrix (filtered) and will run tfidf,
irlba, and then rtsne.

Options:
   -O   [STR]   Output prefix (default is [input].tSNE)
   -D   [INT]   Dimensions to use from irlba (def = $irlba_dims)
   -d   [INT]   Dimensions to embed tSNE (def = $tsne_dims)
   -P   [INT]   Perplexity (def = $perp)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires irlba and Rtsne R packages

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'D'}) {$irlba_dims = $opt{'D'}};
if (defined $opt{'d'}) {$tsne_dims = $opt{'d'}};
if (defined $opt{'P'}) {$perp = $opt{'P'}};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

# IFIDF code from matrix-tfidf
read_matrix_stats($ARGV[0]);
open IN, "$ARGV[0]";
open OUT, ">$opt{'O'}.tfidf";
$h = <IN>; print OUT "$h";
while ($l = <IN>) {
	chomp $l;
	@P = split (/\t/, $l);
	$rowID = shift(@P);
	print OUT "$rowID";
	for ($i = 0; $i < @P; $i++) {
		$tf = ($P[$i]/$MATRIX_CellID_signal{$MATRIX_COLNAMES[$i]});
		$idf = (log(1+($matrix_colNum/($MATRIX_feature_signal{$rowID}+1))));
		$score = sprintf("%.6f", $tf*$idf);
		print OUT "\t$score";
	}
	print OUT "\n";
} close IN; close OUT;

open R, ">$opt{'O'}.counts2tsne.r";

print R "
library(irlba)
library(Rtsne)

# load in tfidf matrix
TFIDF<-as.matrix(read.table(\"$opt{'O'}.tfidf\"))

# irlba
IRLBA<-irlba(TFIDF,$irlba_dims)

# tSNE
TSNE<-Rtsne(IRLBA\$v,dims=$tsne_dims,perplexity=$perp,check_duplicates=FALSE,pca=FALSE)

# output
rownames(TSNE\$Y)<-colnames(TFIDF)
write.table(TSNE\$Y,file=\"$opt{'O'}.counts2tsne.dims\",sep=\"\\t\",col.names=FALSE,row.names=TRUE,quote=FALSE)

";

close R;

system("$Rscript $opt{'O'}.counts2tsne.r 2>/dev/null");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.counts2tsne.r $opt{'O'}.tfidf");
}

}
1;