package sci_commands::pcurve_prune;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("pcurve_prune");

sub pcurve_prune {

@ARGV = @_;
# Defaults
$prune = 0.1;

getopts("O:o:p:l:e:", \%opt);

$die2 = "
scitools pcurve-prune [options] [pcurve prefix]
   or    prune-pcurve

Prunes cells distant from the pcurve. It will print a new file for
the orig and proj dims, and lambda.

Will search for [prefix].orig.dims
                [prefix].proj.dims
                [prefix].lambda (or [prefix].centered.lambda)

Options:
   -O   [STR]   Output prefix (default is [prefix].pruned.lambda)
   -e   [FLT]   Fraction of cells to exclude by distance (def = $prune)
   -o   [STR]   Original dims file (def = auto find [prefix].orig.dims)
   -p   [STR]   Pcurve dims file (def = auto find [prefix].proj.dims)
   -l   [STR]   Lambda file (def = auto find [prefix].lambda)
                (if [prefix].centered.lambda is found, it will use it)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.lambda$//;

if (!defined $opt{'o'}) {
	if (-e "$ARGV[0].orig.dims") {
		$opt{'o'} = "$ARGV[0].orig.dims";
	} else {
		die "\n\nERROR: Cannot find $ARGV[0].orig.dims, check your prefix ($ARGV[0])\n";
	}
}
if (!defined $opt{'p'}) {
	if (-e "$ARGV[0].proj.dims") {
		$opt{'p'} = "$ARGV[0].proj.dims";
	} else {
		die "\n\nERROR: Cannot find $ARGV[0].proj.dims, check your prefix ($ARGV[0])\n";
	}
}
if (!defined $opt{'l'}) {
	if (-e "$ARGV[0].centered.lambda") {
		$opt{'l'} = "$ARGV[0].centered.lambda";
	} elsif (-e "$ARGV[0].lambda") {
		$opt{'l'} = "$ARGV[0].lambda";
	} else {
		die "\n\nERROR: Cannot find $ARGV[0].centered.lambda OR $ARGV[0].lambda, check your prefix ($ARGV[0])\n";
	}
}

read_dims($opt{'o'});
read_pcurve_dims($opt{'p'});
read_values($opt{'l'});

# calculate the distances
%CELLID_pcurveDist = ();
$cells_in_set = 0;
foreach $cellID (keys %CELLID_DIMS) {
	if (defined $CELLID_PCURVE_DIMS{$cellID}[0]) {
		$dist_sum = 0;
		for ($dim = 1; $dim < @{$CELLID_DIMS{$cellID}}; $dim++) {
			if (defined $CELLID_PCURVE_DIMS{$cellID}[$dim]) {
				$dist_sum += ($CELLID_PCURVE_DIMS{$cellID}[$dim] - $CELLID_DIMS{$cellID}[$dim])**2;
			}
		}
		$CELLID_pcurveDist{$cellID} = sqrt($dist_sum);
		$cells_in_set++;
	}
}

# figure out which ones to keep & print to new lambda file
system("cp $opt{'o'} $opt{'O'}.pruned.orig.dims"); # this stays the same
open LAMBDA, ">$opt{'O'}.pruned.lambda";
open PROJ, ">$opt{'O'}.pruned.proj.dims";
$excluded_ct = 0; $checked_ct = 0;
foreach $cellID (sort {$CELLID_pcurveDist{$b}<=>$CELLID_pcurveDist{$a}} keys %CELLID_pcurveDist) {
	$checked_ct++;
	if (($checked_ct/$cells_in_set)<=$prune) {
		$excluded_ct++;
		$threshold = $CELLID_pcurveDist{$cellID};
	} else {
		print LAMBDA "$cellID\t$CELLID_value{$cellID}\n";
		print PROJ "$cellID";
		for ($dim = 1; $dim < @{$CELLID_PCURVE_DIMS{$cellID}}; $dim++) {
			print PROJ "\t$CELLID_PCURVE_DIMS{$cellID}[$dim]";
		} print PROJ "\n";
	}
} close LAMBDA; close PROJ;

print STDERR "INFO: $excluded_ct cells were excluded out of $cells_in_set, with a distance of $threshold or greater from the pcurve.\n";

}
1;
