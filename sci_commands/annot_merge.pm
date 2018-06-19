package sci_commands::annot_merge;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("annot_merge");

sub annot_merge {

$joiner = "_";

@ARGV = @_;
getopts("O:xj:", \%opt);

$die2 = "
scitools annot-merge [options] [annot1] [annot2] (optional annot N etc...)
   or    merge-annot

Will merge annotaitons and include only the intersect

Options:
   -O   [STR]   Output file name / prefix (def = annot1 prefix w/ merge)
   -x           Include cellIDs not present in all files
   -j   [STR]   Character to join annotaitons (def = $joiner)

";

if (!defined $ARGV[1]) {die $die2};
if (defined $opt{'j'}) {$joiner = $opt{'j'}};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.annot$//; $opt{'O'} .= ".merged"};
$opt{'O'} =~ s/\.annot$//;

for ($annotID = 0; $annotID < @ARGV; $annotID++) {
	read_annot($ARGV[$annotID]);
	%{$ANNOTID_CELLID_annot{$annotID}} = %CELLID_annot;
	%CELLID_annot = (); %ANNOT_count = ();
	open IN, "$ARGV[$annotID]";
	while ($l = <IN>) {
		chomp $l; ($cellID,$annot) = split(/\t/, $l);
		$CELLID_annotCT{$cellID}++;
	}
}

open OUT, ">$opt{'O'}.annot";
foreach $cellID (keys %CELLID_annotCT) {
	if ($CELLID_annotCT{$cellID} == @ARGV || defined $opt{'x'}) {
		$new_annot = "";
		for ($annotID = 0; $annotID < @ARGV; $annotID++) {
			if (defined $ANNOTID_CELLID_annot{$annotID}{$cellID}) {
				$new_annot .= $ANNOTID_CELLID_annot{$annotID}{$cellID}.$joiner;
			}
		}
		$new_annot =~ s/$joiner$//;
		print OUT "$cellID\t$new_annot\n";
	}
}
close OUT;

}
1;
