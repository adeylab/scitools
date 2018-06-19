package sci_commands::atac_mergepeak;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("atac_mergepeak");

sub atac_mergepeak {

@ARGV = @_;
getopts("b:O:", \%opt);

$die2 = "
scitools atac-mergepeak [options] [peak bed file] [peak bed file 2] ...
   or    mergepeak(s)

Options:
   -O   [STR]   Output prefix (default is stdout, adds .bed)
   -b   [STR]   Bedtools call (def = $bedtools)

";

if (!defined $ARGV[1]) {die $die2};
if (defined $opt{'O'}) {$opt{'O'} =~ s/\.bed//; open OUT, ">$opt{'O'}.bed"};

$bed_list = ""; for ($i = 0; $i < @ARGV; $i++) {$bed_list .= "$ARGV[$i] "};
open IN, "cat $bed_list | $bedtools sort -i - 2>/dev/null | $bedtools merge -i - 2>/dev/null |";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	if (defined $opt{'O'}) {
		print OUT "$P[0]\t$P[1]\t$P[2]\t$P[0]_$P[1]_$P[2]\n";
	} else {
		print "$P[0]\t$P[1]\t$P[2]\t$P[0]_$P[1]_$P[2]\n";
	}
} close IN;

if (defined $opt{'O'}) {close OUT};

}
1;
