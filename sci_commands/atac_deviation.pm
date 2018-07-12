package sci_commands::atac_deviation;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("atac_deviation");

sub atac_deviation {

@ARGV = @_;
$args = join("\t", @ARGV);

# Defaults:
$permCT = 100;
$binCT = 100;
$minTFCT = 10;
$TSS_flanking = 20000;

getopts("O:X:b:P:B:F:G:", \%opt);

$die2 = "
scitools atac-deviation [options] [counts matrix, may be unfiltered] [feature bed / gene list]
   or    deviation

Options:
   -O   [STR]   Output prefix directory (default is matrix, then bed prefix; adds .dev)
   -I   [STR]   File of feature-peak intersections (col1 = peak ID, col2 = comma sep
                 list of features assigned to the ID) - does not require feature bed
                 to be specified when used.
   -P   [INT]   Permutation count (def = $permCT)
   -B   [INT]   Bin count (def = $binCT)
   -F   [INT]   Min peaks for a feature to include it (def = $minTFCT)
   -G   [STR]   Gene info (refGene.txt formats) - required if using a gene list
                 Shortcut eg: hg38, hg19, mm10
                 Gene list must be annotated, ie: geneName (tab) GeneSetName
                 OR fasta-like: >annotation, then subsequent lines as the
                    genes associated with it (can be comma-sep)
   -S   [INT]   Flanking size (out from TSS in bp, def = $TSS_flanking)
   -b   [STR]   Bedtools call (def = $bedtools)
   -X           Retain intermediate files (def = remove)

";

if (!defined $ARGV[0] || (!defined $opt{'I'} && !defined $ARGV[1])) {die $die2};
if (!defined $opt{'O'}) {
	$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.matrix//;
	$opt{'O'} .= ".".$ARGV[1]; $opt{'O'} =~ s/\.bed//;
} else {$opt{'O'} =~ s/\.dev//};
if (defined $opt{'F'}) {$minTFCT = $opt{'F'}};
if (defined $opt{'P'}) {$permCT = $opt{'P'}};
if (defined $opt{'B'}) {$binCT = $opt{'B'}};

system("mkdir $opt{'O'}.dev");

open LOG, ">$opt{'O'}.dev/atac_deviation.log";
$ts = localtime(time);
print LOG "$ts\tatac-deviation called\n\t\t\t\t$args\n";

open IN, "$ARGV[0]";
$h = <IN>; chomp $h; @H = split(/\t/, $h);
for ($i = 0; $i < @H; $i++) {
	$CELLID_totalCT{$H[$i]} = 0;
}
open OUT, ">$opt{'O'}.dev/peaks.bed";
$siteCT = 0; $sum_all_sites_all_cells = 0;
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$siteID = shift(@P);
	($chr,$start,$end) = split(/[:-_]/, $siteID);
	print OUT "$chr\t$start\t$end\t$siteID\n";
	$siteCT++;
	$SITEID_totalCT{$siteID} = 0;
	for ($i = 0; $i < @H; $i++) {
		$SITEID_totalCT{$siteID}+=$P[$i];
		$CELLID_totalCT{$H[$i]}+=$P[$i];
		$CELLID_SITEID_ct{$H[$i]}{$siteID} = $P[$i];
		$sum_all_sites_all_cells+=$P[$i];
	}
} close IN; close OUT;

$ts = localtime(time);
print LOG "$ts\tMatrix read:
\t\t\t\t$siteCT sites
\t\t\t\t$sum_all_sites_all_cells total signal\n";

if (defined $opt{'G'}) {
	$genes_found = 0; $geneCT = 0; $genes_missing = 0;
	$ts = localtime(time);
	print LOG "$ts\tGene file specified - reading in refgene file.\n";
	if (defined $REF{$opt{'G'}}) {
		$ref_file = $REF{$opt{'G'}};
		$opt{'G'} = $ref_file;
		$opt{'G'} =~ s/\.fa$/\.refGene.txt/;
	}
	read_refgene($opt{'G'});
	
	print LOG "\t\t\t\tMatching genes to coordinates and builing annotated bed.\n";
	open IN, "$ARGV[1]";
	open OUT, ">$opt{'O'}.dev/genes_TSS_$TSS_flanking.bed";
	$gene_spec_mode = 0;
	while ($l = <IN>) {
		chomp $l;
		if ($l =~ /^>/) {$gene_spec_mode=1};
		if ($gene_spec_mode<1) {
			($listed_gene,$listed_annot) = split(/\t/, $l);
			$geneCT++;
			process_gene($listed_gene,$listed_annot);
		} else {
			if ($l =~ /^>/) {
				$listed_annot = $l;
				$listed_annot =~ s/^>//;
			} else {
				$l =~ s/\s//g;
				@GENES = split(/,/, $l);
				foreach $listed_gene (@GENES) {
					$geneCT++;
					process_gene($listed_gene,$listed_annot);
				}
			}
		}
	} close IN; close OUT;
	if ($genes_found<1) {die "ERROR: Gene mode was specified but no genes provided could be found in the refgene file: $opt{'G'}\n"};
}

sub process_gene {
	$gene = $_[0];
	$annot = $_[1];
	if (defined $GENENAME_geneID{$gene}) {
		$geneID = $GENENAME_geneID{$gene};
	} else {$geneID = $gene};
	if (defined $GENEID_coords{$geneID}) {
		$genes_found++;
		($chr,$start,$end) = split(/[:-_]/, $GENEID_coords{$geneID});
		if ($GENEID_strand{$geneID} =~ /\+/) {
			$TSS_pos = $start;
		} else {
			$TSS_pos = $end;
		}
		$flank_start = ($TSS_pos-$TSS_flanking); if ($flank_start<1) {$flank_start=1};
		$flank_end = ($TSS_pos+$TSS_flanking);
		print OUT "$chr\t$flank_start\t$flank_end\t$annot\n";
	} else {
		$genes_missing++;
	}
}

# if this file is not defined, make it
if (!defined $opt{'I'}) {
	
	if (defined $opt{'G'}) {
		open IN, "$opt{'O'}.dev/genes_TSS_$TSS_flanking.bed";
	} else {
		open IN, "$ARGV[1]";
	}
	
	while ($l = <IN>) {
		chomp $l;
		($chr,$start,$end,$TF) = split(/\t/, $l);
		if (!defined $TF || $TF eq "") {
			die "ERROR: Line = $l, cannot identify the feature set annotation. If it is a gene file, specify -G.\n";
		}
		$TF_count{$TF}++;
	} close IN;
	
	open IN, "$bedtools intersect -a $opt{'O'}.dev/peaks.bed -b $ARGV[1] -wa -wb |";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		if (!defined $SITEID_TF_ct{$P[3]}{$P[7]}) {
			$SITEID_TF_ct{$P[3]}{$P[7]} = 1;
			$TF_siteCT{$P[7]}++;
		} else {
			$SITEID_TF_ct{$P[3]}{$P[7]}++;
		}
		if (!defined $SITEID_TFct{$P[3]}) {
			$SITEID_TFct{$P[3]} = 1;
		} else {
			$SITEID_TFct{$P[3]}++;
		}
	} close IN;

	open OUT, ">$opt{'O'}.dev/site_intersects.txt";
	foreach $siteID (keys %SITEID_TF_ct) {
		$TF_list = "";
		foreach $TF (sort keys %{$SITEID_TF_ct{$siteID}}) {
			$TF_list .= "$TF,";
		} $TF_list =~ s/,$//;
		print OUT "$siteID\t$TF_list\n";
	} close OUT;

	open OUT, ">$opt{'O'}.dev/feature_peak_counts.txt";
	print OUT "#TF\tIntersectSites\tFracTFinPeaks\tFracPeaksInTF\n";
	foreach $TF (sort {$TF_siteCT{$b}<=>$TF_siteCT{$a}} keys %TF_siteCT) {
		$frac1 = sprintf("%.2f", ($TF_siteCT{$TF}/$siteCT)*100);
		$frac2 = sprintf("%.2f", ($TF_siteCT{$TF}/$TF_count{$TF})*100);
		print OUT "$TF\t$TF_siteCT{$TF}\t$frac1\t$frac2\n";
	} close OUT;
	%TF_siteCT = (); %SITEID_TF_ct = ();
}

# stratify the peaks into bins
$fraction_increment = 1/$binCT;
$peaks_per_bin = $siteCT/$binCT;
$ts = localtime(time);
print LOG "$ts\tStratifying peaks into $binCT like-signal-bins, with $peaks_per_bin peaks per bin.\n";

# OPTION 2: stratify by read DENSITY of peak, to account for the size of the peak and not purely counts
foreach $siteID (keys %SITEID_totalCT) {
	($chr,$start,$end) = split(/[:-_]/, $siteID);
	$SITEID_density{$siteID} = $SITEID_totalCT{$siteID}/($end-$start);
}

$site_signal_rank = 0;
$binID = 1;
@{$BINID_siteArray{$binID}} = ();
$current_cut = $fraction_increment;
open OUT, ">$opt{'O'}.dev/site_binIDs.txt";
foreach $siteID (sort {$SITEID_density{$a}<=>$SITEID_density{$b}} keys %SITEID_density) { # peak read density stratification
#foreach $siteID (sort {$SITEID_totalCT{$a}<=>$SITEID_totalCT{$b}} keys %SITEID_totalCT) { # read count stratification
	$site_signal_rank++;
	if (($site_signal_rank/$siteCT)>$current_cut) {
		$current_cut+=$fraction_increment;
		$binID++;
		@{$BINID_siteArray{$binID}} = ();
	}
	$SITEID_binID{$siteID} = $binID;
	push @{$BINID_siteArray{$binID}}, $siteID;
	print OUT "$siteID\t$binID\n";
} close OUT;

# shuffle peaks in the bins
foreach $binID (keys %BINID_siteArray) {
	@shuffled = shuffle(@{$BINID_siteArray{$binID}});
	@{$BINID_siteArray{$binID}} = @shuffled;
}

# load in the peak IDs with TF mappings
if (defined $opt{'I'}) {open IN, "$opt{'I'}"} else {open IN, "$opt{'O'}.dev/site_intersects.txt"};
while ($l = <IN>) {
	chomp $l;
	($siteID, $TF_list) = split(/\t/, $l);
	@TF_SET = split(/,/, $TF_list);
	foreach $TF (@TF_SET) {
		if (!defined $TF_siteCT{$TF}) {$TFCT++};
		$TF_siteCT{$TF}++;
		$TF_BINID_count{$TF}{$SITEID_binID{$siteID}}++;
		$TF_totalCT{$TF}+=$SITEID_totalCT{$siteID};
		foreach $cellID (keys %CELLID_totalCT) {
			$TF_CELLID_observed{$TF}{$cellID}+=$CELLID_SITEID_ct{$cellID}{$siteID};
		}
	}
}

$ts = localtime(time);
print LOG "$ts\tFiltering for passing TFs and reporting bins:\n";

$passTFCT = 0;
foreach $TF (keys %TF_siteCT) {
	if ($TF_siteCT{$TF}>=$minTFCT) {
		$TF_passing{$TF} = 1;
		$passTFCT++;
		$TF_binList = "";
		for ($binID = 1; $binID <= $binCT; $binID++) {
			if (defined $TF_BINID_count{$TF}{$binID}) {
				$TF_binList .= "$TF_BINID_count{$TF}{$binID},";
			} else {
				$TF_binList .= "0,";
			}
		} $TF_binList =~ s/,$//;
		print LOG "\t\t\t\t$TF\t$TF_binList\n";
	}
}
$ts = localtime(time);
print LOG "$ts\tLoaded in TF mappings - $TFCT total TFs, $passTFCT passing by having $minTFCT peaks.\n";
print LOG "$ts\tStarting calculations...\n";

# calculate expected counts for each cell for each TF
# also generate background permutations and calculate the deviation

open DEV, ">$opt{'O'}.dev/deviations.txt";
$header = "";
foreach $cellID (sort keys %CELLID_totalCT) {
	$header .= "$cellID\t";
} $header =~ s/\t$//;
print DEV "$header\n";

open OUT, ">$opt{'O'}.dev/feature_stats.txt";
print OUT "#TF\tCellID\tObserved\tExpected\tDeviation\n";

open VAR, ">$opt{'O'}.dev/feature_variability.txt";

foreach $TF (keys %TF_passing) {
	
	$ts = localtime(time);
	print LOG "\t\t\t\t$ts\tPerforming calculations for TF: $TF\n";
	
	# calculate expected for each cell
	foreach $cellID (keys %CELLID_totalCT) {
		$TF_CELLID_expected{$TF}{$cellID} = $TF_totalCT{$TF} * ( $CELLID_totalCT{$cellID} / $sum_all_sites_all_cells );
	}
	
	# generate background prime expected "randExpected" for each TF prime is 0, all others are permutations
	for ($permID = 0; $permID <= $permCT; $permID++) {
		print LOG "\t\t\t\t\t$ts\tPermutation $permID ...\n";
		foreach $binID (keys %{$TF_BINID_count{$TF}}) {
			for ($add_peak = 0; $add_peak < $TF_BINID_count{$TF}{$binID}; $add_peak++) {
				$siteID = $BINID_siteArray{$binID}[$add_peak];
				$pTF_PERMID_totalCT{$TF}{$permID}+=$SITEID_totalCT{$siteID};
				foreach $cellID (keys %CELLID_totalCT) {
					$pTF_PERMID_CELLID_observed{$TF}{$permID}{$cellID}+=$CELLID_SITEID_ct{$cellID}{$siteID};
				}

			}
		}
		foreach $cellID (keys %CELLID_totalCT) {
			$pTF_PERMID_CELLID_expected{$TF}{$permID}{$cellID} = $pTF_PERMID_totalCT{$TF}{$permID} * ( $CELLID_totalCT{$cellID} / $sum_all_sites_all_cells );
			if ($permID>0) {
				$TF_CELLID_deviation_sum_of_squares{$TF}{$cellID} += ( $pTF_PERMID_CELLID_observed{$TF}{$permID}{$cellID} - $pTF_PERMID_CELLID_expected{$TF}{$permID}{$cellID} )**2; # using background TF expected
				$TF_CELLID_variability_sum_of_squares{$TF}{$cellID} += ( $pTF_PERMID_CELLID_observed{$TF}{$permID}{$cellID} - $pTF_PERMID_CELLID_expected{$TF}{$permID}{$cellID} )**2; # using background TF expected
				#$TF_CELLID_deviation_sum_of_squares{$TF}{$cellID} += ( $pTF_PERMID_CELLID_observed{$TF}{$permID}{$cellID} - $TF_CELLID_expected{$TF}{$cellID} )**2; # using TF expected versus background TF expected
				#$TF_CELLID_variability_sum_of_squares{$TF}{$cellID} += ( $pTF_PERMID_CELLID_observed{$TF}{$permID}{$cellID} - $TF_CELLID_expected{$TF}{$cellID} )**2; # using TF expected versus background TF expected
			}
		}
		
		# reshuffle the site bins
		foreach $binID (keys %BINID_siteArray) {
			@shuffled = shuffle(@{$BINID_siteArray{$binID}});
			@{$BINID_siteArray{$binID}} = @shuffled;
		}
		
	}
	
	# calculate deviation & print out info / calculations for the TF
	print DEV "$TF";
	$TF_variability_numerator_sum_of_squares{$TF} = 0;
	$TF_variability_denominator_sum_of_squares{$TF} = 0;
	$pTF_variability_numerator_sum_of_squares{$TF} = 0;
	foreach $cellID (sort keys %CELLID_totalCT) {
		
		if ($TF_CELLID_deviation_sum_of_squares{$TF}{$cellID} != 0) {
			
			$TF_CELLID_deviation{$TF}{$cellID} = ( $TF_CELLID_observed{$TF}{$cellID} - $TF_CELLID_expected{$TF}{$cellID} ) / sqrt( $TF_CELLID_deviation_sum_of_squares{$TF}{$cellID} / $permCT );
			
			$pTF_CELLID_deviation{$TF}{$cellID} = ( $pTF_CELLID_observed{$TF}{'0'}{$cellID} - $pTF_PERMID_CELLID_expected{$TF}{'0'}{$cellID} ) / sqrt( $TF_CELLID_deviation_sum_of_squares{$TF}{$cellID} / $permCT );
			#$pTF_CELLID_deviation{$TF}{$cellID} = ( $pTF_CELLID_observed{$TF}{'0'}{$cellID} - $TF_CELLID_expected{$TF}{$cellID} ) / sqrt( $TF_CELLID_deviation_sum_of_squares{$TF}{$cellID} / $permCT );
			
			$TF_variability_numerator_sum_of_squares{$TF} += ( $TF_CELLID_observed{$TF}{$cellID} - $TF_CELLID_expected{$TF}{$cellID} )**2;
			
			$pTF_variability_numerator_sum_of_squares{$TF} += ( $pTF_CELLID_observed{$TF}{'0'}{$cellID} - $pTF_PERMID_CELLID_expected{$TF}{'0'}{$cellID} )**2;
			#$pTF_variability_numerator_sum_of_squares{$TF} += ( $pTF_CELLID_observed{$TF}{'0'}{$cellID} - $TF_CELLID_expected{$TF}{$cellID} )**2;
			
			$TF_variability_denominator_sum_of_squares{$TF} += ( $TF_CELLID_variability_sum_of_squares{$TF}{$cellID} / $permCT );
			
			print OUT "$TF\t$cellID\t$TF_CELLID_observed{$TF}{$cellID}\t$TF_CELLID_expected{$TF}{$cellID}\t$TF_CELLID_deviation{$TF}{$cellID}\n";
			print DEV "\t$TF_CELLID_deviation{$TF}{$cellID}";
		} else {
			print DEV "\t0";
		}
		
	} print DEV "\n";

	$TF_variability{$TF} = sqrt( $TF_variability_numerator_sum_of_squares{$TF} / $TF_variability_denominator_sum_of_squares{$TF} );
	
	print VAR "$TF\t$TF_variability{$TF}\n"
	
}
close OUT;
close DEV;
close VAR;

open SRT, ">$opt{'O'}.dev/feature_variability.ordered.txt";
$rank = 0;
foreach $TF (sort {$TF_variability{$b}<=>$TF_variability{$a}} keys %TF_variability) {
	print SRT "$rank\t$TF\t$TF_variability{$TF}\n";
	$rank++;
}
close SRT;

}

sub shuffle {
	
	@initial = @_;
	$i = $#initial+1;
	while ( --$i ) {
		my $j = int rand( $i+1 );
		@initial[$i,$j] = @initial[$j,$i];
	}
	return @initial;
	
}

1;