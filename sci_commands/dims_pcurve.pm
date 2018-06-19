package sci_commands::dims_pcurve;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("dims_pcurve");

sub dims_pcurve {

@ARGV = @_;
# Defaults
$range_default = "1-2";

getopts("O:D:XR:", \%opt);

$die2 = "
scitools dims-pcurve [options] [input dims]
   or    pcurve

Projects a principle curve through the dimensions specified

Options:
   -O   [STR]   Output prefix (default is [input].pcurve)
   -D   [RNG]   Range of dimensions to include (range format)
                (e.g. 1-5,6,8; def = $range_default)
   -R   [STR]   Rscript call (def = $Rscript)
   -X           Retain intermediate files (def = delete)

Requires the princurve R package.

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.dims$//; $opt{'O'} =~ s/\.txt$//};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (!defined $opt{'D'}) {$opt{'D'} = $range_default};
read_ranges($opt{'D'});

open R, ">$opt{'O'}.pcurve.r";
print R "
library(princurve)
IN <- as.matrix(read.table(\"$ARGV[0]\",row.names=1))
PCURVE <- principal.curve(IN[,$range_R_set])
write.table(PCURVE\$s,file=\"$opt{'O'}.pcurve.proj.dims\",sep=\"\\t\",col.names=FALSE,row.names=TRUE,quote=FALSE)
write.table(IN[,$range_R_set],file=\"$opt{'O'}.pcurve.orig.dims\",sep=\"\\t\",col.names=FALSE,row.names=TRUE,quote=FALSE)
LAMBDA <- cbind(PCURVE\$lambda)
rownames(LAMBDA) <- rownames(PCURVE\$s)
write.table(LAMBDA,file=\"$opt{'O'}.pcurve.lambda\",sep=\"\\t\",col.names=FALSE,row.names=TRUE,quote=FALSE)
"; close R;

system("$Rscript $opt{'O'}.pcurve.r");

if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.pcurve.r");
}

}
1;
