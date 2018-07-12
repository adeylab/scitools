package sci_commands::atac_enrich;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("atac_enrich");

sub atac_enrich {

@ARGV = @_;
getopts("O:b:", \%opt);

$die2 = "
scitools atac-enrich [options] [peaks.bed] [peak_sets.bed] [feature_sets.bed]
   or    atac-enrichment

This tools will perform a test for the enrichment of features (eg. motifs) within
a set, or sets of peaks. For example, peak_sets.bed may be cicero-linked CCANs.

Each bed must have 4 columns: chr, start, end, annotation; and the peak_sets
must be peaks that exactly match up with peaks in the peaks.bed file.

Output is columns = features and rows = peak sets, p-val for each.

Options:
   -O   [STR]   Output prefix (default is [peak_sets].[feature_sets])
   -b   [STR]   Bedtools call (def = $bedtools)
   
";

if (!defined $ARGV[2]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[1]; $opt{'O'} =~ s/\.bed$//; $opt{'O'} .= ".".$ARGV[2]; $opt{'O'} =~ s/\.bed$//};

open LOG, ">$opt{'O'}.atac_enrich.log";
$ts = localtime(time);
print LOG "$ts\tProgram called: $ARGV[0] $ARGV[1] $ARGV[2]\n";

# count peaks
open IN, "$ARGV[0]";
while ($l = <IN>) {$peakCT++};
close IN;

$ts = localtime(time);
print LOG "$ts\tTotal peaks in $ARGV[0] = $peakCT\n";

# annotate all peaks for features
open IN, "$bedtools intersect -a $ARGV[0] -b $ARGV[2] -wa -wb |";
while ($l = <IN>) {
	chomp $l;
	($peak_chr,$peak_start,$peak_end,$peakID,$feature_chr,$feature_start,$feature_end,$featureID) = split(/\t/, $l);
	$peakID = $peak_chr."_".$peak_start."_".$peak_end; # ensure that peakID format is this
	if (!defined $FEATURE_PEAK_hit{$featureID}{$peakID}) {
		push @{$PEAKID_features{$peakID}}, $featureID;
		if (!defined $FEATUREID_peak_count{$featureID}) {
			$FEATUREID_peak_count{$featureID}=1;
			$featureCT++;
		} else {
			$FEATUREID_peak_count{$featureID}++;
		}
		$FEATURE_PEAK_hit{$featureID}{$peakID} = 1;
	}
	$peak_feature_intersect_count++;
} close IN;

$ts = localtime(time);
print LOG "$ts\tFeature count = $featureCT\n";
print LOG "$ts\tPeak - Feature intersect count = $peak_feature_intersect_count\n";

# load in peak sets
open IN, "$ARGV[1]";
while ($l = <IN>) {
	chomp $l;
	($peak_chr,$peak_start,$peak_end,$setID) = split(/\t/, $l);
	$peakID = $peak_chr."_".$peak_start."_".$peak_end;
	$SETID_peakCT{$setID}++;
	push @{$SETID_peakList{$setID}}, $peakID;
	$peakSetCT++;
} close IN;

$ts = localtime(time);
print LOG "$ts\tPeak set count = $peakSetCT\n";

# precompute logfact
$ts = localtime(time);
print LOG "$ts\tBuilding logfact hash (max=$peakCT) ...";
build_logfact($peakCT);
$ts = localtime(time);
print LOG "done. ($ts)\n";

# perform all by all
open OUT, ">$opt{'O'}.pvals.matrix";
open COUNTS, ">$opt{'O'}.set_counts.matrix";
$header = "";
foreach $featureID (sort keys %FEATUREID_peak_count) {
	$header .= "$featureID\t";
} $header =~ s/\t$//;
print OUT "$header\n";
print COUNTS "$header\n";

$ts = localtime(time);
print LOG "$ts\tCalculating p-values for each factor for each set ...";

foreach $setID (sort keys %SETID_peakCT) {
	print OUT "$setID"; print COUNTS "$setID";
	foreach $featureID (sort keys %FEATUREID_peak_count) {
		$peaks_with_feature = $FEATUREID_peak_count{$featureID};
		$peaks_without_feature = $peakCT - $FEATUREID_peak_count{$featureID};
		$set_peaks_with_feature = 0;
		foreach $peakID (@{$SETID_peakList{$setID}}) {
			$added_peak = 0;
			foreach $check_feature (@{$PEAKID_features{$peakID}}) {
				if ($check_feature eq $featureID && $added_peak < 1) {$set_peaks_with_feature++; $added_peak=1};
			}
		}
		$set_peaks_total = $SETID_peakCT{$setID};
		if ($set_peaks_with_feature > 0) {
			$hypergeom_cumulative_pval = hypergeom($peaks_with_feature,$peaks_without_feature,$set_peaks_total,$set_peaks_with_feature);
		} else {
			$hypergeom_cumulative_pval = 1;
		}
		print COUNTS "\t$peaks_with_feature,$peaks_without_feature,$set_peaks_total,$set_peaks_with_feature";
		print OUT "\t$hypergeom_cumulative_pval";
		push @ALL_PVALS, $hypergeom_cumulative_pval;
	}
	print OUT "\n"; print COUNTS "\n";
}

print LOG " done.\n";
close OUT; close COUNTS; close LOG;

}

sub hypergeom {
    ($n, $m, $N, $i) = @_;
	
	$max_possible = ($n + $m + $N + $i);
	if (!defined $logfact[ $max_possible ]) {add_logfact($max_possible)};
	
    $loghyp1 = $logfact[ $m ] + $logfact[ $n ]
             + $logfact[ $N ] + $logfact[ $m + $n - $N ];
    $loghyp2 = $logfact[ $i ] + $logfact[ $n - $i ] 
             + $logfact[ $m + $i - $N ] + $logfact[ $N - $i ] 
             + $logfact[ $m + $n ];

    return exp($loghyp1 - $loghyp2);
}

sub build_logfact {
	$log_max = $_[0];
	$value = 0;
	for ($current = 1; $current <= $log_max; $current++) {
		$value += log($current);
		$logfact[$current] = $value;
	}
}

sub add_logfact {
	$new_max = $_[0];
	$value = $logfact[$log_max];
	for ($current = $log_max+1; $current <= $new_max; $current++) {
		$value += log($current);
		$logfact[$current] = $value;
	}
	$log_max = $new_max;
}

1;