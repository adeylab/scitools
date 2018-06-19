package sci_commands::matrix_pca;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_pca");

sub matrix_pca {

@ARGV = @_;
getopts("O:XR:", \%opt);

$die2 = "
scitools matrix-pca [options] [input matrix]
   or    pca

Produces a dims file with PCA coordinates

Options:
   -O   [STR]   Output prefix (default is [input].PCA.dims)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires methods R package

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

open R, ">$opt{'O'}.PCA.r";
print R "
library(methods)
tM<-t(read.table(\"$ARGV[0]\"))
D<-dist(tM,method=\"euclidean\",diag=FALSE,upper=FALSE)
PCA<-prcomp(D,center=TRUE)
write.table(PCA\$x,file=\"$opt{'O'}.PCA.dims\",sep=\"\\t\",row.names=TRUE,col.names=FALSE,quote=FALSE)
";
close R;

system("$Rscript $opt{'O'}.PCA.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.PCA.r");
}

}
1;
