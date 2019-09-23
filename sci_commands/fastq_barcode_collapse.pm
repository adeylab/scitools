package sci_commands::fastq_barcode_collapse;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("fastq_barcode_collapse");

sub fastq_barcode_collapse {

# defaults
$hdist = 1;
$jdist = 2;
$excl_list = "AAAAA,CCCCC,GGGGG,TTTTTT,NN";

@ARGV = @_;
getopts("O:H:x:v:N:oJ:", \%opt);

$die2 = "
scitools fastq-barcode-collapse [options] [index1_read_fastq],(index2_read_fastq) [read1_fastq] (read2_fastq)
   or:   fastq-collapse
         barcode-collapse
		 
Designed for microwell-ATAC-seq barcodes.

Will take in fastq files and collapse by barcodes into a sci-format fastq.
If two index files are needed for the cell barcode, provide the second in
the same argument separated by a comma. Output will be gzipped.

Note: Will greedily assign to barcodes matching by hamming distance.
(eg: barc1->barc2 & barc2->barc3; so barc1->barc3)
Each hamming distance match is a 'jump', and only a certain number of jumps are allowed
before reverting to the original barcode sequence.

Failed reads (ie excluded by tag patterns) will not be reported at all.

Read names are:
@[barcode with -N prefix]:[read number].[jumps](.ORIG=original barcode)#0/1(|2)

Options:
   -O   [STR]   Output prefix (def = index1 prefix . collapsed)
   -H   [STR]   Hamming distance, max=2, 1 is strongly recommended (def = $hdist)
   -J   [STR]   Maximum jump distance to final barcode (def = $jdist)
   -x   [STR]   Comma separated exclusion list of sequences
                 def = $excl_list
   -N   [STR]   Add in this pattern to the cell barcode names
                 (eg. RUN1; ideally no special characters other than _)
   -o           Include original barcode in the read name (def = no)

";

if (!defined $ARGV[1]) {die $die2};
if (!defined $opt{'O'}) {
	$opt{'O'} = $ARGV[0];
	$opt{'O'} =~ s/,.+$//;
	$opt{'O'} =~ s/\.gz$//;
	$opt{'O'} =~ s/\.fq$//;
	$opt{'O'} =~ s/\.fastq//;
}
if (defined $opt{'H'}) {$hdist = $opt{'H'}};
if (defined $opt{'J'}) {$jdist = $opt{'J'}};
if (defined $opt{'x'}) {$excl_list = $opt{'x'}};
if (defined $opt{'N'}) {
	if ($opt{'N'} !~ /_$/) {$opt{'N'} .= "_"};
	$opt{'N'} =~ s/(:|\"|\'|,|\.|\@|\>|\<|\@|\#)/_/;
}

@EXCLUSION_PATTERNS = split(/,/, $excl_list);

open LOG, ">$opt{'O'}.collapse.log";
$ts = localtime(time);
print LOG "INFO: $ts, program called.\n\tOptions:\n";
foreach $option (keys %opt) {
	print LOG "\t\t$option\t$opt{$option}\n";
}
print LOG "\tIndex Files: $ARGV[0]\n\tRead1 File: $ARGV[1]\n";
if (defined $ARGV[2]) {print LOG "\tRead2 File: $ARGV[2]\nINFO: Processing index fastq files...\n"};

# read in all barcodes

if ($ARGV[0] =~ /,/) {
	($index1_file,$index2_file) = split(/,/, $ARGV[0]);
	if ($index1_file =~ /\.gz$/) {
		open IX1, "$zcat $index1_file |";
	} else {
		open IX1, "$index1_file";
	}
	if ($index2_file =~ /\.gz$/) {
		open IX2, "$zcat $index2_file |";
	} else {
		open IX2, "$index2_file";
	}
} else {
	if ($ARGV[0] =~ /\.gz$/) {
		open IX1, "$zcat $ARGV[0] |";
	} else {
		open IX1, "$ARGV[0]";
	}
}

$tag_exclusion_ct = 0; $raw_barc_total = 0; $processed_reads = 0; $increment = 1000000; $report = $increment;
while ($tag = <IX1>) {
	chomp $tag; $seq = <IX1>; chomp $seq; $null = <IX1>; $qual = <IX1>; chomp $qual;
	if (defined $index2_file) {
		$tag2 = <IX2>; chomp $tag2; $seq2 = <IX2>; chomp $seq2; $null = <IX2>; $qual2 = <IX2>; chomp $qual2;
		$seq .= $seq2;
		$qual .= $qual2;
	}
	$tag =~ s/^\@//; $tag =~ s/\s.+$//;
	$exclude = 0;
	foreach $pattern (@EXCLUSION_PATTERNS) {
		if ($seq =~ /$pattern/i) {
			$exclude++;
		}
	}
	if ($exclude == 0) {
		if (!defined $RAW_BARC_ct{$seq}) {$raw_barc_total++};
		$RAW_BARC_ct{$seq}++;
		$TAG_barc{$tag} = $seq;
	} else {
		$TAG_exclude{$tag}++;
		$tag_exclusion_ct++;
	}
	$processed_reads++;
	if ($processed_reads>=$report) {
		$ts = localtime(time);
		print LOG "\t$ts, $report reads processed. $tag_exclusion_ct excluded due to tag sequence.\n";
		$report += $increment;
	}
}
close IX1; if (defined $index2_file) {close IX2};

$ts = localtime(time);
print LOG "INFO: $ts, n=$raw_barc_total total raw barcodes found, n=$tag_exclusion_ct reads were excluded based on exclusion sequences.\nINFO: Processing from least to most abundant ...\n";

# go from lowest count to highest
$increment = 0.05; $report = $increment;
@BASES = ("A","C","G","T","N");
foreach $barc (sort {$RAW_BARC_ct{$a}<=>$RAW_BARC_ct{$b}} keys %RAW_BARC_ct) {
	$match = "null"; $match_count = $RAW_BARC_ct{$barc};
	@SEQ = split(//, $barc);
	for ($pos = 0; $pos < @SEQ; $pos++) {
		foreach $subBase (@BASES) {
			if ($SEQ[$pos] ne $subBase) {
				@SUB_SEQ = @SEQ;
				$SUB_SEQ[$pos] = $subBase;
				$test_barc = join("", @SUB_SEQ);
				if (defined $RAW_BARC_ct{$test_barc}) {
					if ($RAW_BARC_ct{$test_barc} > $match_count) {
						$match = $test_barc;
						$match_count = $RAW_BARC_ct{$test_barc};
#						print STDERR "DEBUG: $barc matches $test_barc w/ hd1; match count is $match_count\n";
					}
				}
				if ($hdist>1) {
					for ($pos2 = 0; $pos2 < @SEQ; $pos2++) {
						foreach $subBase2 (@BASES) {
							@SUB_SEQ2 = @SUB_SEQ;
							$SUB_SEQ2[$pos2] = $subBase2;
							$test_barc2 = join("", @SUB_SEQ2);
							if (defined $RAW_BARC_ct{$test_barc}) {
								if ($RAW_BARC_ct{$test_barc2} > $match_count) {
									$match = $test_barc2;
									$match_count = $RAW_BARC_ct{$test_barc2};
								}
							}
						}
					}
				}
			}
		}
		if ($match ne "null") {
			$RAW_BARC_match{$barc} = $match;
		}
	}
	$barcs_processed++;
	if (($barcs_processed/$raw_barc_total) >= $report) {
		$ts = localtime(time);
		print LOG "\t$ts, $report fraction processed (n=$barcs_processed)\n";
		$report += $increment;
	}
}
$ts = localtime(time);
print LOG "INFO: $ts, Matching complete - generating output.\n";

# now go through reads and rename them in sci format
if ($ARGV[1] =~ /\.gz$/) {
	open R1, "$zcat $ARGV[1] |";
} else {
	open R1, "$ARGV[1]";
}
open O1, "| gzip > $opt{'O'}.collapsed.1.fq.gz";

if (defined $ARGV[2]) {
	if ($ARGV[2] =~ /\.gz$/) {
		open R2, "$zcat $ARGV[2] |";
	} else {
		open R2, "$ARGV[2]";
	}
	open O2, "| gzip > $opt{'O'}.collapsed.2.fq.gz";
}
$jump_excluded = 0;
while ($tag = <R1>) {
	chomp $tag; $seq = <R1>; chomp $seq; $null = <R1>; $qual = <R1>; chomp $qual;
	if (defined $ARGV[2]) {
		$tag2 = <R>; chomp $tag2; $seq2 = <R2>; chomp $seq2; $null = <R2>; $qual2 = <R2>; chomp $qual2;
	}
	$tag =~ s/^\@//; $tag =~ s/\s.+$//;
	if (!defined $TAG_barc{$tag}) {
		if (!defined $TAG_exclude{$tag}) {
			print LOG "WARNING: $tag in read1 was not found in the index fastq file!\n";
		}
	} else {
		$barc = $TAG_barc{$tag};
		$barc_jumps = 0;
		while (defined $RAW_BARC_match{$barc}) {
			$barc_jumps++;
			$new = $RAW_BARC_match{$barc};
			$barc = $new;
		}
		$JUMP_HIST{$barc_jumps}++;
		
#		print STDERR "DEBUG: $TAG_barc{$tag} jumped $barc_jumps to match $barc\n";
		
		if ($barc_jumps > $jdist) {
			$jump_excluded++;
			$barc = $TAG_barc{$tag};
			$barc_jumps = 0;
		}
		
		if (defined $opt{'N'}) {$barc = $opt{'N'}.$barc};
		$BARC_collapsed_ct{$barc}++;
		
		if (defined $opt{'o'}) {
			$new_tag = "\@$barc:$BARC_collapsed_ct{$barc}.$barc_jumps.ORIG=$TAG_barc{$tag}";
		} else {
			$new_tag = "\@$barc:$BARC_collapsed_ct{$barc}.$barc_jumps";
		}
		
		print O1 "$new_tag#0/1\n$seq\n\+\n$qual\n";
		if (defined $ARGV[2]) {
			print O2 "$new_tag#0/2\n$seq2\n\+\n$qual2\n";
		}
	}
}

close O1; if (defined $ARGV[2]) {close O2};

print LOG "INFO: Complete, n=$jump_excluded reads were uncollapsed due to > $jdist jumps.\nINFO: Barcode jump histogram:\n\tJUMPS\tCOUNT\n";
foreach $jump_count (sort {$a<=>$b} keys %JUMP_HIST) {
	print LOG "\t$jump_count\t$JUMP_HIST{$jump_count}\n";
}

open REPORT, ">$opt{'O'}.collapse_counts.values";
foreach $barc (sort {$BARC_collapsed_ct{$b}<=>$BARC_collapsed_ct{$a}} keys %BARC_collapsed_ct) {
	print REPORT "$barc\t$BARC_collapsed_ct{$barc}\n";
} close REPORT;
	
}
1;
