package sci_commands::pcurve_center;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("pcurve_center");

sub pcurve_center {

@ARGV = @_;
# Defaults
$max = 1;
$type = "centered";

getopts("O:A:a:M:P:r", \%opt);

$die2 = "
scitools pcurve-center [options] [input pcurve.lambda]
   or    center-pcurve
         lambda-center
         center-lambda

Centers a pcurve lambda and adjusts the scale.
Default is to center on all cells and extend lambda out.

Options:
   -O   [STR]   Output prefix (default is [prefix].centered.lambda)
   -A   [STR]   Annotation file (optional)
   -a   [STR]   Annotation that should be used to center the pcurve
                  (requires -A to be specified)
   -M   [FLT]   Maximum absolute value to adjust lambda to (either side of 0)
                  (def = $max)
   -P   [STR]   Progression type (centered/cen, progressive/pro)
                  (def = $type, if -a it will only center on that annotation)
   -r           Reverse order (defaults to current order)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.lambda$//;
if (defined $opt{'M'}) {$max = $opt{'M'}};
if (defined $opt{'P'}) {$type = $opt{'P'}};
if ($type !~ /^(cen|pro)/i) {
	die "\nERROR: Cannot determine preogresstion type ($type) - specify as centered/cen or progressive/pro\n";
}
if (defined $opt{'A'}) {read_annot($opt{'A'})};
if (defined $opt{'r'}) {$r = -1} else {$r = 1};

read_values($ARGV[0]);

if (defined $opt{'a'}) {
	$type = "centered";
	if (defined $ANNOT_count{$opt{'a'}}) {
		$annot_value_sum = 0; $annot_cell_count = 0; $center_value = 0;
		foreach $cellID (keys %CELLID_value) {
			if ($CELLID_annot{$cellID} eq $opt{'a'}) {
				$annot_cell_count++;
				$annot_value_sum += $CELLID_value{$cellID};
			}
		}
		if ($annot_cell_count>0) {
			$center_value = $annot_value_sum/$annot_cell_count;
			if (abs($value_max - $center_value)>abs($value_min - $center_value)) {
				$progression_span = abs($value_max-$center_value);
			} else {
				$progression_span = abs($value_min - $center_value);
			}
		} else {
			die "\nERROR: No cells in labda file $ARGV[0] have the annotation $opt{'a'} that was provided.\n";
		}
	} else {
		die "\nERROR: Provided annotation $opt{'a'} cannot be found in the annotaiton file provided: $opt{'A'}, or the annotation file was not specified (-A)\n";
	}
} elsif ($type =~ /^cen/i) {
	$center_value = $value_mean;
	if (abs($value_max - $center_value)>abs($value_min - $center_value)) {
		$progression_span = abs($value_max-$center_value);
	} else {
		$progression_span = abs($value_min - $center_value);
	}
} else {
	$center_value = $value_min;
	$progression_span = $value_range;
}

open OUT, ">$opt{'O'}.centered.lambda";
foreach $cellID (sort {$CELLID_value{$a}<=>$CELLID_value{$b}} keys %CELLID_value) {
	$lambda = ((($CELLID_value{$cellID}-$center_value)/$progression_span)*$max*$r);
	print OUT "$cellID\t$lambda\n";
} close OUT;

}
1;
