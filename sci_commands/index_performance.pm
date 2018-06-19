package sci_commands::index_performance;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("index_performance");

sub index_performance {

@ARGV = @_;

# Defaults
$gradient_def = "BuYlRd";

getopts("O:I:R:s:t:A:xG:", \%opt);

# DEFAULTS
@LETTERS = ("0", "A", "B", "C", "D", "E", "F", "G", "H");
%LETTER_NUM = ("A"=>"1", "B"=>"2", "C"=>"3", "D"=>"4", "E"=>"5", "F"=>"6", "G"=>"7", "H"=>"8");
%WELL_xy = ();
foreach $letter (keys %LETTER_NUM) {
	for ($number = 1; $number <= 12; $number++) {
		$well = $letter.$number;
		$WELL_xy{$well} = "$LETTER_NUM{$letter}\t$number";
	}
}
$threshold = 1000000;

$die2 = "
scitools index-perform [options] [fastq or bam]

Will generate the read counts from each well of each tier of indexing.
Should be on a pre-filetered since it will plot all barcode combos present.

Options:
   -O   [STR]   Output prefix, will create a folder (def = input prefix)
   -I   [STR]   Index file
         (default = $VAR{'SCI_index_file'})
         (Index names must be in form of: [Tier]_[set]_[i5/i7]_[A-H/1-12])
   -A   [STR]   Annotation file (only include cell IDs in the annot file)
   -t   [INT]   Threshold of reads for a plate to include (def = $threshold)
   -x           Do not report plate-plate-well stats (just plate-well)
   -G   [GRD]   Color gradient for plots (def = $gradient_def, biased 0.65)
                For all available gradients, run 'scitools gradient'
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
if (defined $opt{'A'}) {read_annot($opt{'A'})};
if (!defined $opt{'G'}) {$opt{'G'} = $gradient_def};
$gradient_function = get_gradient($opt{'G'});
$gradient_function =~ s/\)$//; $gradient_function .= ",bias=0.65)";

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
	if ($ARGV[0] =~ /\.gz$/) {
		open IN, "$zcat $ARGV[0] |";
	} elsif ($ARGV[0] =~ /\.fq$/) {
		open IN, "$ARGV[0]";
	} else {die "\n\nCannot determine file input type! Provide either a fastq (can be gzipped) or a bam file!\n\n$die2"};
	while ($cellID = <IN>) {
		chomp $cellID; $null = <IN>; $null = <IN>; $null = <IN>;
		$cellID =~ s/:.+$//; $cellID =~ s/^\@//;
		$CELLID_count{$cellID}++;
	} close IN;
}

$NEX_set_count = 0; $PCR_set_count = 0;
foreach $cellID (keys %CELLID_count) {
	if (!defined $opt{'A'} || defined $CELLID_annot{$cellID}) {

		$ix1 = substr($cellID,0,8);
		$ix2 = substr($cellID,8,10);
		$ix3 = substr($cellID,18,8);
		$ix4 = substr($cellID,26,10);
		
		if (!defined $INDEX_POS_SEQ_id{'1'}{$ix1} ||
			!defined $INDEX_POS_SEQ_id{'2'}{$ix2} ||
			!defined $INDEX_POS_SEQ_id{'3'}{$ix3} ||
			!defined $INDEX_POS_SEQ_id{'4'}{$ix4}) {
			print STDERR "WARNING: Barcode combo: $cellID ($ix1)($ix2)($ix3)($ix4) has a barcode not present in the index file!\n";
		}
		
		$NEX_set = $INDEX_POS_SEQ_id{'3'}{$ix3}.$INDEX_POS_SEQ_id{'1'}{$ix1};
		$PCR_set = $INDEX_POS_SEQ_id{'4'}{$ix4}.$INDEX_POS_SEQ_id{'2'}{$ix2};
		
		$NEX_well = $INDEX_POS_SEQ_well{'3'}{$ix3}.$INDEX_POS_SEQ_well{'1'}{$ix1};
		$PCR_well = $INDEX_POS_SEQ_well{'4'}{$ix4}.$INDEX_POS_SEQ_well{'2'}{$ix2};
		
		if (!defined $NEX_SET_total{$NEX_set}) {$NEX_set_count++};
		if (!defined $PCR_SET_total{$PCR_set}) {$PCR_set_count++};
		
		$NEX_SET_total{$NEX_set}+=$CELLID_count{$cellID};
		$PCR_SET_total{$PCR_set}+=$CELLID_count{$cellID};
		
		$NEX_SET_WELL_total{$NEX_set}{$NEX_well}+=$CELLID_count{$cellID};
		$PCR_SET_WELL_total{$PCR_set}{$PCR_well}+=$CELLID_count{$cellID};
		
		$NEX_SET_PCR_SET_total{$NEX_set}{$PCR_set}+=$CELLID_count{$cellID};
		
		$NEX_SET_PCR_SET_NEX_WELL_total{$NEX_set}{$PCR_set}{$NEX_well}+=$CELLID_count{$cellID};
		$NEX_SET_PCR_SET_PCR_WELL_total{$NEX_set}{$PCR_set}{$PCR_well}+=$CELLID_count{$cellID};
		
	}
}

system("mkdir $opt{'O'}.index_performance");

open SUMMARY, ">$opt{'O'}.index_performance/summary.txt";
$ts = localtime(time);
print SUMMARY "$ts scitools index-performance on $ARGV[0]
Transposase-based index totals:\n";
foreach $NEX_set (keys %NEX_SET_total) {
	print SUMMARY "  $NEX_set\t$NEX_SET_total{$NEX_set}\n";
}
print SUMMARY "PCR-based index totals:\n";
foreach $PCR_set (keys %PCR_SET_total) {
	print SUMMARY "  $PCR_set\t$PCR_SET_total{$PCR_set}\n";
}

$plates_to_plot = 0;
open OUT, ">$opt{'O'}.index_performance/plate_performance.txt";
foreach $NEX_set (keys %NEX_SET_WELL_total) {
	if ($NEX_SET_total{$NEX_set}>=$threshold) {
		foreach $NEX_well (sort keys %WELL_xy) {
			$coord = $WELL_xy{$NEX_well};
			if (defined $NEX_SET_WELL_total{$NEX_set}{$NEX_well}) {
				print OUT "NEX\t$NEX_set\tALL\t$NEX_well\t$NEX_SET_WELL_total{$NEX_set}{$NEX_well}\t$coord\n";
			} else {
				print OUT "NEX\t$NEX_set\tALL\t$NEX_well\t0\t$coord\n";
			}
		}
		$plates_to_plot++;
	}
}
foreach $PCR_set (keys %PCR_SET_WELL_total) {
	if ($PCR_SET_total{$PCR_set}>=$threshold) {
		foreach $PCR_well (sort keys %WELL_xy) {
			$coord = $WELL_xy{$PCR_well};
			if (defined $PCR_SET_WELL_total{$PCR_set}{$PCR_well}) {
				print OUT "PCR\tALL\t$PCR_set\t$PCR_well\t$PCR_SET_WELL_total{$PCR_set}{$PCR_well}\t$coord\n";
			} else {
				print OUT "PCR\tALL\t$PCR_set\t$PCR_well\t0\t$coord\n";
			}
		}
		$plates_to_plot++;
	}
}
if ($PCR_set_count > 1 && $NEX_set_count > 1 && !defined $opt{'x'}) {
	foreach $NEX_set (keys %NEX_SET_WELL_total) {
		foreach $PCR_set (keys %PCR_SET_WELL_total) {
			if ($NEX_SET_PCR_SET_total{$NEX_set}{$PCR_set}>=$threshold) {
				foreach $NEX_well (sort keys %WELL_xy) {
					$coord = $WELL_xy{$NEX_well};
					if (defined $NEX_SET_PCR_SET_NEX_WELL_total{$NEX_set}{$PCR_set}{$NEX_well}) {
						print OUT "NEX\t$NEX_set\t$PCR_set\t$NEX_well\t$NEX_SET_PCR_SET_NEX_WELL_total{$NEX_set}{$PCR_set}{$NEX_well}\t$coord\n";
					} else {
						print OUT "NEX\t$NEX_set\t$PCR_set\t$NEX_well\t0\t$coord\n";
					}
				}
				$plates_to_plot++;
			}
		}
	}
	foreach $NEX_set (keys %NEX_SET_WELL_total) {
		foreach $PCR_set (keys %PCR_SET_WELL_total) {
			if ($NEX_SET_PCR_SET_total{$NEX_set}{$PCR_set}>=$threshold) {
				foreach $PCR_well (sort keys %WELL_xy) {
					$coord = $WELL_xy{$PCR_well};
					if (defined $NEX_SET_PCR_SET_PCR_WELL_total{$NEX_set}{$PCR_set}{$PCR_well}) {
						print OUT "PCR\t$NEX_set\t$PCR_set\t$PCR_well\t$NEX_SET_PCR_SET_PCR_WELL_total{$NEX_set}{$PCR_set}{$PCR_well}\t$coord\n";
					} else {
						print OUT "PCR\t$NEX_set\t$PCR_set\t$PCR_well\t0\t$coord\n";
					}
				}
				$plates_to_plot++;
			}
		}
	}
}
close OUT;

if ($plates_to_plot<2) {
	die "\nThe number of plates passing filters to plot is zero. Try reducing the threshold to plot. (currently $threshold)\n";
}

$plot_height = 0.5+(int(($plates_to_plot/3)+1)*2.5);

open R, ">$opt{'O'}.index_performance/plot_plates.r";
print R "
library(ggplot2)
$gradient_function
IN<-read.table(\"$opt{'O'}.index_performance/plate_performance.txt\")
colnames(IN)<-c(\"Plate_type\",\"NEX_set\",\"PCR_set\",\"Plate_well\",\"Well_count\",\"Row\",\"Column\")
PLT<-ggplot(data=IN) + theme_bw() +
	geom_tile(aes(Column,Row,fill=log10(Well_count+1))) +
	scale_y_reverse(breaks=c(1,2,3,4,5,6,7,8)) +
	xlab(\"Column\") + ylab(\"Row\") +
	facet_wrap(c(\"Plate_type\",\"NEX_set\",\"PCR_set\"),ncol=3,labeller=label_wrap_gen(multi_line=FALSE)) +
	theme(strip.background=element_rect(fill=\"transparent\")) +
	scale_fill_gradientn(colours=gradient_funct(99)) +
	scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11,12)) +
	labs(fill=\"Log10\nReads\") +
	theme(strip.background=element_blank(),
		panel.grid=element_blank(),
		axis.line=element_blank(),
		axis.ticks=element_blank(),
		panel.background=element_blank(),
		plot.background=element_blank())
ggsave(plot=PLT,filename=\"$opt{'O'}.index_performance/plate_performance.png\",height=$plot_height,width=12)
ggsave(plot=PLT,filename=\"$opt{'O'}.index_performance/plate_performance.pdf\",height=$plot_height,width=12)
"; close R;

system("$Rscript $opt{'O'}.index_performance/plot_plates.r");

}
1;
