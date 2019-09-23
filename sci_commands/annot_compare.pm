package sci_commands::annot_compare;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("annot_compare");

sub annot_compare {

@ARGV = @_;
getopts("O:", \%opt);

$die2 = "
scitools annot-compare [options] [annot1] [annot2]
   or    compare-annot

Will compare the proportions in each of the annotaitons across files.
CellIDs must be the same (will take the intersect only)

Options:
   -O   [STR]   Output prefix (def = annot1.annot2)
                Adds .compare.txt suffix

";

if (!defined $ARGV[1]) {die $die2};

if (!defined $opt{'O'}) {
	$pfx1 = $ARGV[0]; $pfx1 =~ s/\.annot$//;
	$pfx2 = $ARGV[1]; $pfx2 =~ s/\.annot$//;
	$opt{'O'} = "$pfx1.$pfx2";
}
$opt{'O'} =~ s/\.txt$//; $opt{'O'} =~ s/\.compare$//;

open IN, "$ARGV[0]";
while ($l = <IN>) {
	chomp $l;
	($cellID,$annot) = split(/\t/, $l);
	$CELLID_annot1{$cellID} = $annot;
} close IN;

$cell_match_ct = 0;
open IN, "$ARGV[1]";
while ($l = <IN>) {
	chomp $l;
	($cellID,$annot2) = split(/\t/, $l);
	if (defined $CELLID_annot1{$cellID}) {
		$annot1 = $CELLID_annot1{$cellID};
		$cell_match_ct++;
		$ANNOT1_ANNOT2_ct{$annot1}{$annot2}++;
		$ANNOT1_ANNOT2_ct{$annot2}{$annot1}++;
	}
} close IN;

open OUT, ">$opt{'O'}.compare.txt";
print OUT "Annot1\tAnnot2\tN_A2_in_A1\tPct_A2_in_A1\tPct_total\n";
foreach $annot1 (keys %ANNOT1_ANNOT2_ct) {
	$annot_ct = 0;
	foreach $annot2 (keys %{$ANNOT1_ANNOT2_ct{$annot1}}) {
		$annot_ct += $ANNOT1_ANNOT2_ct{$annot1}{$annot2};
	}
	foreach $annot2 (keys %{$ANNOT1_ANNOT2_ct{$annot1}}) {
		$A1_pct = sprintf("%.3f", ($ANNOT1_ANNOT2_ct{$annot1}{$annot2}/$annot_ct)*100);
		$G_pct = sprintf("%.3f", ($ANNOT1_ANNOT2_ct{$annot1}{$annot2}/$cell_match_ct)*100);
		print OUT "$annot1\t$annot2\t$ANNOT1_ANNOT2_ct{$annot1}{$annot2}\t$A1_pct\t$G_pct\n";
	}
} close OUT;

}
1;