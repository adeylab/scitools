package sci_commands::matrix_tsne;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_tsne");

sub matrix_tsne {

@ARGV = @_;
# Deafults
$dims = 2;
$perp = 30;

getopts("O:XD:R:P:d", \%opt);

$die2 = "
scitools matrix-tsne [options] [input matrix/dims]
   or    dims-tsne
         tsne

Produces a dims file with tSNE coordinates, will auto
detect if it is a matrix or dims input file.
Default is matrix.

Options:
   -O   [STR]   Output prefix (default is [input].tSNE)
   -D   [INT]   Dimensions to embed tSNE in (def = $dims)
   -d           Input is a dims file (if file suffix is not .dims)
   -P   [INT]   Perplexity (def = $perp)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires Rtsne R package

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//; $opt{'O'} =~ s/\.dims$//};
if (defined $opt{'D'}) {$dims = $opt{'D'}};
if (defined $opt{'P'}) {$perp = $opt{'P'}};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

open R, ">$opt{'O'}.p$perp.tSNE.r";

print R "
library(Rtsne)";
if (defined $opt{'d'} || $ARGV[0] =~ /dims$/) {
	print R "
tIN<-read.table(\"$ARGV[0]\",row.names=1)
TSNE<-Rtsne(tIN,dims=$dims,perplexity=$perp,check_duplicates=FALSE,pca=FALSE)";
} else {
	print R "
tIN<-t(read.table(\"$ARGV[0]\"))
TSNE<-Rtsne(tIN,dims=$dims,perplexity=$perp,check_duplicates=FALSE)";
}
print R "
DIMS<-TSNE\$Y
rownames(DIMS)<-rownames(tIN)
write.table(DIMS,file=\"$opt{'O'}.p$perp.tSNE.dims\",sep=\"\\t\",col.names=FALSE,row.names=TRUE,quote=FALSE)
"; close R;

system("$Rscript $opt{'O'}.p$perp.tSNE.r 2>/dev/null");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.p$perp.tSNE.r");
}

}
1;
