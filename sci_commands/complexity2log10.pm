package sci_commands::complexity2log10;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("complexity2log10");

sub complexity2log10 {

@ARGV = @_;
getopts("O:F:N:c:L:", \%opt);

# defaults
$complexity_range = "0,100";
$read_thresh = 1000;

$die2 = "
scitools complexity2log10 [options] [complexity file]

Will make a log10 unique reads values file

Options:
   -O   [STR]   Output prefix (default is complexity file prefix)
   -L   [STR]   List of CellIDs to retain
   -F   [STR]   Filter-bam log file (overrides -c and -N)
   -N   [STR]   Min unique reads (def = $read_thresh)
   -c   [MIN,MAX] Mina nd max complexity percent to retain (def = $complexity_range)

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {
	$opt{'O'} = $ARGV[0];
	$opt{'O'} =~ s/\.complexity\.txt$//;
}
if (defined $opt{'L'}) {
	open IN, "$opt{'L'}";
	while ($l = <IN>) {
		chomp $l; $l =~ s/\s.+$//;
		$CELLID_inList{$l} = 1;
	} close IN;
}
# opt F to override
if (defined $opt{'F'}) {
	open FILT, "$opt{'F'}";
	$null = <FILT>;
	while ($l = <FILT>) {
		chomp $l; $l =~ s/^\s+//;
		@P = split(/\s+/, $l);
		if ($P[0] eq "c") {
			$opt{'c'} = $P[1];
			print STDERR "INFO: -c detected in filt file as $opt{'c'}\n";
		} elsif ($P[0] eq "N") {
			$opt{'N'} = $P[1];
			print STDERR "INFO: -N detected in filt file as $opt{'N'}\n";
		}
	} close FILT;
}
if (defined $opt{'c'}) {
	$complexity_range = $opt{'c'};
}
($min_comp,$max_comp) = split(/,/, $complexity_range);
if (defined $opt{'N'}) {
	$read_thresh = $opt{'N'};
}

open IN, "$ARGV[0]";
open OUT, ">$opt{'O'}.log10_unique_reads.values";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	if (!defined $opt{'L'} || defined $CELLID_inList{$P[0]}) {
		if ($P[3] >= $read_thresh && $P[4] <= $max_comp && $P[4] >= $min_comp) {
			$log10ur = sprintf("%.4f", log($P[4])/log(10));
			print OUT "$P[1]\t$log10ur\n";
		}
	}
} close IN; close OUT;

}
1;