package sci_commands::group_specific;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("group_specific");

sub group_specific {

@ARGV = @_;

getopts("b:", \%opt);

$die2 = "

scitools group-specific [options] [output prefix] [bed1] [bed2] ... (bedN)

Will take speaks called for each group and intersect to pull those specific to each group.
Will also produce a merged file with all peaks and peaks in every group specific set will
have the smae peaks coordinates as the peaks in the merged.

Options:
   -b   [STR]   bedtools call (def = $bedtools)
   -X           retain intermediate files (def = no)

";

if (!defined $ARGV[2]) {die $die2};
if (defined $opt{'b'}) {$bedtools = $opt{'b'}};

$files = "";
for ($i = 1; $i < @ARGV; $i++) {$files .= "$ARGV[$i] "};
open OUT, ">$ARGV[0].group_specific.merged.bed";
open IN, "cat $files | $bedtools sort -i - | $bedtools merge -i - |";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	print OUT "$P[0]\t$P[1]\t$P[2]\t$P[0]_$P[1]_$P[2]\n";
} close IN; close OUT;

open OUT, "| $bedtools sort - > $ARGV[0].group_specific.bed";
for ($i = 1; $i < @ARGV; $i++) {
	$files = "";
	for ($j = 1; $j < @ARGV; $j++) {
		if ($i != $j) {
			$files .= "$ARGV[$j] ";
		}
	}
	$name = $ARGV[$i]; $name =~ s/\.bed$//;
	system("cat $files | $bedtools sort -i - | $bedtools merge -i - > $ARGV[0].$name.inv_tmp.bed");
	system("$bedtools intersect -v -b $ARGV[0].$name.inv_tmp.bed -a $ARGV[$i] > $ARGV[0].$name.spc_tmp.bed");
	open IN, "$bedtools intersect -b $ARGV[0].merged.bed -a $ARGV[0].$name.spc_tmp.bed |";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		print OUT "$P[0]\t$P[1]\t$P[2]\t$name\_specific\n";
	} close IN;
	
	if (!defined $opt{'X'}) {
		system("rm -f $ARGV[0].*_tmp.bed");
	}
} close OUT;


}
1;
