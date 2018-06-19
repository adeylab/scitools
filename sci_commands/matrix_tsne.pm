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

getopts("O:XD:R:P:", \%opt);

$die2 = "
scitools matrix-tsne [options] [input matrix]
   or    tsne

Produces a dims file with tSNE coordinates

Options:
   -O   [STR]   Output prefix (default is [input].tSNE)
   -D   [INT]   Dimensions to embed tSNE in (def = $dims)
   -P   [INT]   Perplexity (def = $perp)
   -X           Retain intermediate files (def = delete)
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires Rtsne R package

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix$//};
if (defined $opt{'D'}) {$dims = $opt{'D'}};
if (defined $opt{'P'}) {$perp = $opt{'P'}};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};

open R, ">$opt{'O'}.p$perp.tSNE.r";
print R "
library(Rtsne)
tIN<-t(read.table(\"$ARGV[0]\"))
TSNE<-Rtsne(tIN,dims=$dims,perplexity=$perp,check_duplicates=FALSE)
write.table(TSNE\$Y,file=\"$opt{'O'}.p$perp.tSNE.raw\",col.names=FALSE,row.names=FALSE)
"; close R;

system("$Rscript $opt{'O'}.p$perp.tSNE.r 2>/dev/null");

open IN, "$ARGV[0]";
$h = <IN>; close IN;
chomp $h; $h =~ s/\r//;
@CELLID_list = split(/\t/, $h);

open IN, "$opt{'O'}.p$perp.tSNE.raw";
open OUT, ">$opt{'O'}.p$perp.tSNE.dims";
while ($l = <IN>) {
	chomp $l;
	$cellID = shift(@CELLID_list);
	print OUT "$cellID\t$l\n";
} close IN; close OUT;

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.p$perp.tSNE.r $opt{'O'}.p$perp.tSNE.raw");
}

}
1;
