package sci_commands::index_performance;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("index_performance");

sub index_performance {

@ARGV = @_;

getopts("O:I:R:s:t:A:G:b:NT:", \%opt);

# DEFAULTS
@LETTERS = ("0", "A", "B", "C", "D", "E", "F", "G", "H");
%LETTER_NUM = ("A"=>"1", "B"=>"2", "C"=>"3", "D"=>"4", "E"=>"5", "F"=>"6", "G"=>"7", "H"=>"8");
$gradient_def = "BuYlRd";
$bias = 0.65;
$threshold = 1000;
$minUniq = 1000;

$die2 = "
scitools index-perform [options] [fastq, bam, annot, complexity, or list]

Will generate the read counts from each well of each tier of indexing.
Should be on a pre-filetered since it will plot all barcode combos present.

Requires barcode cell IDs.

Options:
   -O   [STR]   Output prefix, will create a folder (def = input prefix)
   -I   [STR]   Index file
         (default = $VAR{'SCI_index_file'})
         (Index names must be in form of: [Tier]_[set]_[i5/i7]_[A-H/1-12])
   -A   [STR]   Annotation file (only include cell IDs in the annot file)
   -t   [INT]   Threshold of reads for a plate to include (def = $threshold)
   -T   [INT]   Min unique reads (if complexity file, def = $minUniq)
   -N           Do not log scale (eg for cell counts)
   -G   [GRD]   Color gradient for plots (def = $gradient_def)
                For all available gradients, run 'scitools gradient'
   -b   [FLT]   Gradient bias (def = $bias)
   -R   [STR]   Rscript call (def = $Rscript)
   -s   [STR]   samtools call (if bam; def = $samtools)

Note: Requires ggplot2 R package for plotting

";

if (!defined $ARGV[0]) {die $die2};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.gz$//; $opt{'O'} =~ s/\.fq$//; $opt{'O'} =~ s/\.bam$//;
if (-e "$opt{'O'}.index_performance") {die "$opt{'O'}.index_performance Directory already exists! Exiting!\n\n"};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (defined $opt{'I'}) {$VAR{'SCI_index_file'} = $opt{'I'}};
if (defined $opt{'t'}) {$threshold = $opt{'t'}};
if (defined $opt{'T'}) {$minUniq = $opt{'T'}};
if (defined $opt{'A'}) {read_annot($opt{'A'})};
if (!defined $opt{'G'}) {$opt{'G'} = $gradient_def};
if (defined $opt{'b'}) {$bias = $opt{'b'}};
$gradient_function = get_gradient($opt{'G'});
$gradient_function =~ s/\)$//; $gradient_function .= ",bias=$bias)";

read_indexes($VAR{'SCI_index_file'});

%CELLID_count = ();
if ($ARGV[0] =~ /\.bam$/) {
	open IN, "$samtools view $ARGV[0] 2>/dev/null |";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		$cellID = $P[0]; $cellID =~ s/:.+$//;
		$CELLID_count{$cellID}++;
	} close IN;
} else {
	if ($ARGV[0] =~ /fastq|fq/) {
		if ($ARGV[0] =~ /\.gz$/) {
			open IN, "$zcat $ARGV[0] |";
		} elsif ($ARGV[0] =~ /\.fq$/) {
			open IN, "$ARGV[0]";
		}
		while ($cellID = <IN>) {
			chomp $cellID; $null = <IN>; $null = <IN>; $null = <IN>;
			$cellID =~ s/:.+$//; $cellID =~ s/^\@//;
			$CELLID_count{$cellID}++;
		} close IN;
	} elsif ($ARGV[0] =~ /complexity/) {
		open IN, "$ARGV[0]";
		while ($l = <IN>) {
			chomp $l;
			@P = split(/\t/, $l);
			$CELLID_count{$P[1]}++;
		} close IN;
	} elsif ($ARGV[0] =~ /annot|annotation|list|txt/) {
		open IN, "$ARGV[0]";
		while ($l = <IN>) {
			chomp $l;
			@P = split(/\s+/, $l);
			$CELLID_count{$P[0]}++;
		} close IN;
	} else {
		die "ERROR: Cannot determine file type!\n";
	}
}

$NEX_set_count = 0; $PCR_set_count = 0;
foreach $cellID (keys %CELLID_count) {
	if (!defined $opt{'A'} || defined $CELLID_annot{$cellID}) {
		
		# get indexes
		$nex_i7 = substr($cellID,0,8);
		$pcr_i7 = substr($cellID,8,10);
		$nex_i5 = substr($cellID,18,8);
		$pcr_i5 = substr($cellID,26,10);
		
		if (!defined $INDEX_POS_SEQ_id{'1'}{$nex_i7} ||
			!defined $INDEX_POS_SEQ_id{'2'}{$pcr_i7} ||
			!defined $INDEX_POS_SEQ_id{'3'}{$nex_i5} ||
			!defined $INDEX_POS_SEQ_id{'4'}{$pcr_i5}) {
			print STDERR "WARNING: Barcode combo: $cellID ($nex_i7)($pcr_i7)($nex_i5)($pcr_i5) has a barcode not present in the index file!\n";
		}
		
		#split index name to map to coordinates
		$nex_i7_set = $INDEX_POS_SEQ_id{'1'}{$nex_i7};
		$pcr_i7_set = $INDEX_POS_SEQ_id{'2'}{$pcr_i7};
		$nex_i5_set = $INDEX_POS_SEQ_id{'3'}{$nex_i5};
		$pcr_i5_set = $INDEX_POS_SEQ_id{'4'}{$pcr_i5};
				
		$nex_i7_col = $INDEX_POS_SEQ_well{'1'}{$nex_i7};
		$pcr_i7_col = $INDEX_POS_SEQ_well{'2'}{$pcr_i7};
		$nex_i5_row = $INDEX_POS_SEQ_well{'3'}{$nex_i5};
		$pcr_i5_row = $INDEX_POS_SEQ_well{'4'}{$pcr_i5};
		
		#make and store the wellID and counts
		$nex_id = $LETTER_NUM{$nex_i5_row}.",".$nex_i7_col;
		$pcr_id = $LETTER_NUM{$pcr_i5_row}.",".$pcr_i7_col;
		$NEX_SET_count{$nex_i5_set.$nex_i7_set}+=$CELLID_count{$cellID};
		$PCR_SET_count{$pcr_i5_set.$pcr_i7_set}+=$CELLID_count{$cellID};
		$NEX_SET_WELLID_count{$nex_i5_set.$nex_i7_set}{$nex_id}+=$CELLID_count{$cellID};
		$PCR_SET_WELLID_count{$pcr_i5_set.$pcr_i7_set}{$pcr_id}+=$CELLID_count{$cellID};
		
		#get stats on broader sets
		#add in later ##########################
		
	}
}

system("mkdir $opt{'O'}.index_performance");

@PLOT_FILES = ();

foreach $nexSet (keys %NEX_SET_WELLID_count) {
	if ($NEX_SET_count{$nexSet} >= $threshold) {
		open OUT, ">$opt{'O'}.index_performance/nex_$nexSet.wellID.counts";
		print OUT "#row\tcol\tcount\n";
		for ($row = 1; $row <= 8; $row++) {
			for ($col = 1; $col <= 12; $col++) {
				$wellID = $row.",".$col;
				if (defined $NEX_SET_WELLID_count{$nexSet}{$wellID}) {
					print OUT "$row\t$col\t$NEX_SET_WELLID_count{$nexSet}{$wellID}\n";
				} else {
					print OUT "$row\t$col\t0\n";
				}
			}
		}
		close OUT;
		push @PLOT_FILES, "$opt{'O'}.index_performance/nex_$nexSet.wellID.counts";
		push @PLOT_TITLE, "Transposase Index Combination: $nexSet";
	} else {
		print STDERR "NEX Set: $nexSet has $NEX_SET_count{$nexSet} reads which is < threshold of $threshold.\n";
	}
}

foreach $pcrSet (keys %PCR_SET_WELLID_count) {
	if ($PCR_SET_count{$pcrSet} >= $threshold) {
		open OUT, ">$opt{'O'}.index_performance/pcr_$pcrSet.wellID.counts";
		for ($row = 1; $row <= 8; $row++) {
			for ($col = 1; $col <= 12; $col++) {
				$wellID = $row.",".$col;
				if (defined $PCR_SET_WELLID_count{$pcrSet}{$wellID}) {
					print OUT "$row\t$col\t$PCR_SET_WELLID_count{$pcrSet}{$wellID}\n";
				} else {
					print OUT "$row\t$col\t0\n";
				}
			}
		}
		close OUT;
		push @PLOT_FILES, "$opt{'O'}.index_performance/pcr_$pcrSet.wellID.counts";
		push @PLOT_TITLE, "PCR Index Combination: $nexSet";
	} else {
		print STDERR "PCR Set: $pcrSet has $PCR_SET_count{$pcrSet} reads which is < threshold of $threshold.\n";
	}
}

open R, ">$opt{'O'}.index_performance/plot_plates.r";
print R "library(ggplot2)
$gradient_function\n";

for ($i = 0; $i<@PLOT_FILES; $i++) {
	$plotFile = $PLOT_FILES[$i];
	$plotTitle = $PLOT_TITLE[$i];
	print R "
# Plot File = $plotFile
IN<-read.table(\"$plotFile\")
colnames(IN)<-c(\"row\",\"col\",\"value\")
PLT<-ggplot(data=IN) + theme_bw() + ggtitle(\"$plotTitle\") +
	xlab(\"Column\") + ylab(\"Row\") +";
if (!defined $opt{'N'}) {
	print R "	geom_tile(aes(col,row,fill=log10(value))) +
	labs(fill=\"Log10\nReads\") +";
} else {
print R "	geom_tile(aes(col,row,value)) +
	labs(fill=\"Count\") +";
} print R "
	scale_y_reverse(breaks=c(1,2,3,4,5,6,7,8)) +
	scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12)) +
	scale_fill_gradientn(colours=gradient_funct(99))
ggsave(plot=PLT,filename=\"$plotFile.png\",height=4,width=6)
ggsave(plot=PLT,filename=\"$plotFile.pdf\",height=4,width=6)
";
}

close R;

system("$Rscript $opt{'O'}.index_performance/plot_plates.r");

}
1;
