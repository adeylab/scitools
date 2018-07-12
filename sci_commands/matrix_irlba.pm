package sci_commands::matrix_irlba;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_irlba");

sub matrix_irlba {

$dims = 50;

@ARGV = @_;
getopts("O:XR:D:", \%opt);

$die2 = "
scitools matrix-irlba [options] [input matrix]
   or    irlba

Produces a dims file with irlba coordinates

Options:
   -O   [STR]   Output prefix (default is [input].irlba.dims)
   -D   [INT]   Max dimensions to compute (def = $dims)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires irlba R package

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

open R, ">$opt{'O'}.irlba.r";
print R "
library(irlba)
IN<-as.matrix(read.table(\"$ARGV[0]\"))
I<-irlba(IN, $dims)
rownames(I\$v)<-colnames(IN)
write.table(I\$v, file=\"$opt{'O'}.irlba_$dims.dims\",col.names=FALSE,row.names=TRUE,quote=FALSE,sep=\"\\t\")
";
close R;

system("$Rscript $opt{'O'}.irlba.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.irlba.r");
}

}
1;
