package sci_commands::annot_proportions;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("annot_proportions");

sub annot_proportions {

@ARGV = @_;
getopts("O:A:", \%opt);

$die2 = "
scitools annot-proportions [options] -A [annot file] [check file]

Produces proportions of annotations in respective file.
Supports matrix, values, annot, list files (incl sparseMatrix cols)

Options:
   -O   [STR]   Output prefix (default is check file prefix)
   -A   [STR]   Annotation file (required)
   
";

if (!defined $ARGV[0]) {die "Provide a check file!\n$die2"};
if (!defined $opt{'A'}) {die "Provide an annotation file as -A option!\n$die2"};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.\D+$//};
read_annot($opt{'A'});

if ($ARGV[0] =~ /\.gz$/) {
	open IN, "$zcat $ARGV[0] |";
} else {
	open IN, "$ARGV[0]";
}

if ($ARGV[0] =~ /\.matrix/) {
	$l = <IN>; close IN; chomp $l;
	@P = split(/\t/, $l);
	foreach $cellID (@P) {
		$CELLIDS{$cellID}++;
	}
} else {
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		$CELLIDS{$P[0]}++;
	}
	close IN;
}

foreach $cellID (keys %CELLIDS) {
	$cellCT++;
	if (!defined $CELLID_annot{$cellID}) {
		$CELLID_annot{$cellID} = "Undefined";
	}
	$annot = $CELLID_annot{$cellID};
	$ANNOT_cellCT{$annot}++;
}

open OUT, ">$opt{'O'}.proportions.txt";
foreach $annot (sort {$b<=>$a} keys %ANNOT_cellCT) {
	$frac = sprintf("%.3f", $ANNOT_cellCT{$annot}/$cellCT);
	print OUT "$annot\t$ANNOT_cellCT{$annot}\t$frac\n";
} close OUT;

}
1;
