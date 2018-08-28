package sci_commands::fastq_from_geo;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("fastq_from_geo");

sub fastq_from_geo {

@ARGV = @_;

$die2 = "
scitools fastq-from-geo [output prefix] [fastq1.gz] [fastq2.gz] [barcode_fastq(3).gz]
   or    geo-fastq
         fastq-geo
         geo2fastq
   sra can be substituted for geo as well for calling the command.

This tool takes in three fastq files where one fastq is the complete
demultiplexed cell barcode and outputs paired fastq files with the
read name set tot he barcode (standard scitools format)

";

if (!defined $ARGV[3]) {die $die2};

if ($ARGV[1] =~ /.gz$/) {open R1, "$zcat $ARGV[1] |"} else {open R1, "$ARGV[1]"};
if ($ARGV[2] =~ /.gz$/) {open R2, "$zcat $ARGV[2] |"} else {open R2, "$ARGV[2]"};
if ($ARGV[3] =~ /.gz$/) {open IX, "$zcat $ARGV[3] |"} else {open IX, "$ARGV[3]"};

open O1, "| $gzip > $ARGV[0].1.fq.gz";
open O2, "| $gzip > $ARGV[0].2.fq.gz";

$read_number = 0;

while ($null = <R1>) {
	$seq1 = <R1>; chomp $seq1; $null = <R1>; $qual1 = <R1>; chomp $qual1;
	$null = <R2>; $seq2 = <R2>; chomp $seq2; $null = <R2>; $qual2 = <R2>; chomp $qual2;
	$null = <IX>; $barc = <IX>; chomp $barc; $null = <IX>; $null = <IX>;
	print O1 "\@$barc:$read_number#0/1\n$seq1\n\+\n$qual1\n";
	print O2 "\@$barc:$read_number#0/2\n$seq2\n\+\n$qual2\n";
	$read_number++;
}

close R1; close R2; close IX; close O1; close O2;

}

1;