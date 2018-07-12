package sci_commands::signal_make;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("signal_make");

sub signal_make {

@ARGV = @_;

# Defaults
$span_size = 5000;
$sub_size = 500;
$slide_size = 250;
$matching_annot = "TRUE";

getopts("O:A:a:W:w:l:b:s:Xn:rv", \%opt);

$die2 = "
scitools signal-make [options] [bam file] [bed file of features]
   or    make-signal

This tool will generate a specialized matrix file with multiple entries
over each window provided. Columns are the annotaiton and the sub-window
and rows are the bed file features. It also produces a z-scored singal
matrix where the z-score is produced for all sub windows within the
annotation.

Options:
   -O   [STR]   Output prefix (default is [bam prefix].[bed prefix].signal)
   -A   [STR]   Annotation file (will aggregate signal over the annotations)
                  Note: if not provided, will assume signal for every cellID
                  wich is not recommended - for plotting individual cells, use
                  'scitools plot-reads'
   -a   [STR]   Comma separated list of annotations to include (requires -A)
                  Will position annotations in the specified order
   -r           Order rows by the order of the input bed file (overrides -n)
   -n   [STR]   Sort by the signal for this annotation.
                  Default is to sort by all (if no -A), by subsets of matching
                  bed and cell annotations (ie bed names are the same as annot
                  names), or by the first annotation listed.
   -W   [INT]   Window size - the span of the signal view (def = $span_size)
   -w   [INT]   Sub-window size for signal (def = $sub_size)
   -l   [INT]   Sub-window slide size (def = $slide_size)
   -b   [STR]   Bedtools call (def = $bedtools)
   -s   [STR]   Samtools call (def = $samtools)
   -v           Verbose
   -X           Do not delete intermediate files (def = delete)

";

if (!defined $ARGV[1]) {die $die2};
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotaiton file (-A) if specifying annotations to filter (-a)!\n$die2"};
if (!defined $opt{'O'}) {
	$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/bam$//;
	$opt{'O'} .= $ARGV[1]; $opt{'O'} =~ s/\.bed$//;
}
if (defined $opt{'W'}) {$span_size = $opt{'W'}};
if (defined $opt{'w'}) {$sub_size = $opt{'w'}};
if (defined $opt{'l'}) {$slide_size = $opt{'l'}};

if (!defined $opt{'A'}) {
	print STDERR "
WARNING: No annotation specified - will treat as a single annotation.\n";
	open IN, "$samtools view $ARGV[0] |";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		($cellID,$other) = split(/:/, $P[0]);
		$CELLID_annot{$cellID} = $cellID;
		$ANNOT_read_count{'all'}++;
	} close IN;
} else {
	read_annot($opt{'A'});
	open IN, "$samtools view $ARGV[0] |";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		$cellID = $P[0]; $cellID =~ s/:.+$//;
		if (defined $CELLID_annot{$cellID}) {
			$ANNOT_read_count{$CELLID_annot{$cellID}}++;
		}
	} close IN;
}

if (defined $opt{'v'}) {
	foreach $annot (keys %ANNOT_read_count) {
		print STDERR "INFO: Annot = $annot, total reads = $ANNOT_read_count{$annot}\n";
	}
}

if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		if (!defined $ANNOT_read_count{$annot}) {
			die "ERROR: Annotation specified in -a ($annot) was not found in the annotaiton file $opt{'A'}\n";
		}
		$ANNOT_include{$annot} = 1;
	}
} else {
	foreach $annot (sort keys %ANNOT_read_count) {
		$ANNOT_include{$annot} = 1;
		push @ANNOT_LIST, $annot;
	}
}

if (defined $opt{'v'}) {
	for ($annot_pos = 0; $annot_pos < @ANNOT_LIST; $annot_pos++) {
		print STDERR "INFO: Annot pos = $annot_pos, Annot = $ANNOT_LIST[$annot_pos], Annot ct = $ANNOT_read_count{$ANNOT_LIST[$annot_pos]}\n";
	}
}

$flank_size = int($span_size/2);

# make temporary bed file with subwindows
open IN, "$ARGV[1]";
open OUT, "| $bedtools sort -i - > $opt{'O'}.signal.subwin.bed";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$chr = $P[0];
	$start = $P[1];
	$end = $P[2];
	$winID = "$chr\_$start\_$end";
	if (defined $P[3]) {
		$fgroup = $P[3];
		if (!defined $ANNOT_include{$P[3]}) {
			$matching_annot = "FALSE";
		}
	} else {
		$fgroup = "Feature";
	}
	$FGROUP_count{$fgroup}++;
	$WINID_fgroup{$winID} = $fgroup;
	$FGROUP_WINID_list{$fgroup}{$winID} = 1;
	$WINID_list{$winID} = 1;
	$midPt = int(($start+$end)/2);
	if (($midPt-$flank_size)>0) {
		$subNum = 0;
		for ($i = ($midPt-$flank_size); $i <= ($midPt+$flank_size)-$sub_size; $i += $slide_size) {
			print OUT "$chr\t$i\t".($i+$sub_size)."\t$winID\t$subNum\n";
			$subNum++;
		}
		$totwin = $subNum-1;
	}
	$featureCT++;
} close IN; close OUT;

if (defined $opt{'v'}) {
	print STDERR "INFO: Total windows: $totwin, Features = $featureCT\n";
}

# intitalize counts
foreach $winID (keys %WINID_list) {
	for ($i = 0; $i <= $totwin; $i++) {
		foreach $annot (keys %ANNOT_include) {
			$WINID_SUB_ANNOT_count{$winID}{"$annot\_window_$i"} = 0;
		}
	}
}

# intersect to get counts
open IN, "$bedtools intersect -abam $ARGV[0] -b $opt{'O'}.signal.subwin.bed -bed -wa -wb 2>/dev/null |";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$cellID = $P[3]; $cellID =~ s/:.+$//;
	if (defined $CELLID_annot{$cellID}) {
		$annot = $CELLID_annot{$cellID};
		$winID = $P[15];
		$subWin = $P[16];
		$subID = "$annot\_window_$subWin";
		$SUBID_annot{$subID} = $annot;
		$WINID_SUB_ANNOT_count{$winID}{$subID}++;
		$passing_intersections++;
	}
} close IN;

if (defined $opt{'v'}) {
	print STDERR "INFO: Passing Intersections of windows and bam: $passing_intersections\n";
}

# print the signal matrix
open OUT, ">$opt{'O'}.raw.signal";
open NRM, ">$opt{'O'}.norm.signal";

$header = ""; @COLNAMES = ();
for ($annot_pos = 0; $annot_pos < @ANNOT_LIST; $annot_pos++) {
	$annot = $ANNOT_LIST[$annot_pos];
	for ($subWin = 0; $subWin <= $totwin; $subWin++) {
		$header .= "$annot\_window_$subWin\t";
		push @COLNAMES, "$annot\_window_$subWin";
	}
} $header =~ s/\t$//;
print OUT "$header\n";
print NRM "$header\n";

# order rows
$midWin = int($totwin/2); @ROW_ORDER = ();
if (defined $opt{'r'}) {
	open IN, "$ARGV[1]";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		$winID = "$P[0]\_$P[1]\_$P[2]";
		if (defined $WINID_list{$winID}) {
			push @ROW_ORDER, $winID;
		}
	} close IN;
} elsif (!defined $opt{'A'}) {
	$midWinID = "all_window_$midWin";
	foreach $fgroup (sort keys %FGROUP_WINID_list) {
		foreach $winID (sort {$WINID_SUB_ANNOT_count{$b}{$midWinID}<=>$WINID_SUB_ANNOT_count{$a}{$midWinID}} keys %{$FGROUP_WINID_list{$fgroup}}) {
			push @ROW_ORDER, $winID;
		}
	}
} else {
	if (defined $opt{'n'}) {
		$midWinID = "$opt{'n'}\_window_$midWin";
		foreach $winID (sort {$WINID_SUB_ANNOT_count{$b}{$midWinID}<=>$WINID_SUB_ANNOT_count{$a}{$midWinID}} keys %WINID_SUB_ANNOT_count) {
			push @ROW_ORDER, $winID;
		}
	} elsif ($matching_annot = "TRUE") {
		for ($annot_pos = 0; $annot_pos < @ANNOT_LIST; $annot_pos++) {
			$annot = $ANNOT_LIST[$annot_pos];
			$midWinID = "$annot\_window_$midWin";
			foreach $winID (sort {$WINID_SUB_ANNOT_count{$b}{$midWinID}<=>$WINID_SUB_ANNOT_count{$a}{$midWinID}} keys %{$FGROUP_WINID_list{$annot}}) {
				push @ROW_ORDER, $winID;
			}
		}
	} else {
		($annot,$null) = split(/_window_/, $COLNAMES[0]);
		print STDERR "INFO: Annotaiton for cells provided, but annotaitons of peaks do not all match cell annotations.
      If this is incorrect - double check that annotaitons and bed feature names match.
      Proceeding with sorting on the first listed annotation: $annot\n";
		$midWinID = "$annot\_window_$midWin";
		foreach $fgroup (sort keys %FGROUP_WINID_list) {
			foreach $winID (sort {$WINID_SUB_ANNOT_count{$b}{$midWinID}<=>$WINID_SUB_ANNOT_count{$a}{$midWinID}} keys %{$FGROUP_WINID_list{$fgroup}}) {
				push @ROW_ORDER, $winID;
			}
		}
	}
}

if (defined $opt{'v'}) {
	print STDERR "INFO: Total rows that have been ordered is: ".@ROW_ORDER."\n";
}

# print raw and norm
for ($row_pos = 0; $row_pos < @ROW_ORDER; $row_pos++) {
	$winID = $ROW_ORDER[$row_pos];
	print OUT "$winID"; print NRM "$winID";
	for ($colNum = 0; $colNum < @COLNAMES; $colNum++) {
		$subID = $COLNAMES[$colNum];
		$annot = $SUBID_annot{$subID};
		if (!defined $WINID_SUB_ANNOT_count{$winID}{$subID}) {
			$WINID_SUB_ANNOT_count{$winID}{$subID} = 0;
		}
		print OUT "\t$WINID_SUB_ANNOT_count{$winID}{$subID}";
		$nrm = sprintf("%.4f", ($WINID_SUB_ANNOT_count{$winID}{$subID}/$ANNOT_read_count{$annot})*($featureCT*$totwin));
		print NRM "\t$nrm";
		$ANNOT_posCT{$annot}++;
		$ANNOT_sum{$annot}+=$WINID_SUB_ANNOT_count{$winID}{$subID};
		push @{$ANNOT_values{$annot}}, $WINID_SUB_ANNOT_count{$winID}{$subID};
	}
	print OUT "\n"; print NRM "\n";
}
close OUT; close NRM;

# calculate each annot mean
foreach $annot (keys %ANNOT_posCT) {
	$ANNOT_mean{$annot} = $ANNOT_sum{$annot}/$ANNOT_posCT{$annot};
	$stdev_sum = 0;
	foreach $value (@{$ANNOT_values{$annot}}) {
		$stdev_sum += ($ANNOT_mean{$annot} - $value)**2;
	}
	$ANNOT_stdev{$annot} = sqrt($stdev_sum/$ANNOT_posCT{$annot});
}

open OUT, ">$opt{'O'}.zscore.signal";
print OUT "$header\n";

for ($row_pos = 0; $row_pos < @ROW_ORDER; $row_pos++) {
	$winID = $ROW_ORDER[$row_pos];
	print OUT "$winID";
	for ($colNum = 0; $colNum < @COLNAMES; $colNum++) {
		$subID = $COLNAMES[$colNum];
		$annot = $SUBID_annot{$subID};
		$zscore = sprintf("%.4f", ($WINID_SUB_ANNOT_count{$winID}{$subID}-$ANNOT_mean{$annot})/$ANNOT_stdev{$annot});
		print OUT "\t$zscore";
	}
	print OUT "\n";
}
close OUT;


if (!defined $opt{'X'}) {
	system("rm -f $opt{'O'}.signal.subwin.bed");
}

}
1;
