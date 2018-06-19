package sci_commands::fastq_merge;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("fastq_merge");

sub fastq_merge {

@ARGV = @_;

getopts("O:r:", \%opt);

$die2 = "
scitools fastq-merge [options] [output prefix] [read1 fq comma sep] [read2 fq comma sep]
   or    merge-fastq

Will merge fastq files that have the same set of cellIDs present
Adds in a run identifier (after read number) based on the order of fastqs.

Options:
   -r   [STR]   Comma separated run IDs (def = 1,2,3... etc)

";

if (defined $opt{'O'}) {unshift @ARGV, $opt{'O'}};
if (!defined $ARGV[2]) {die $die2};
$ARGV[0] =~ s/\.gz$//; $ARGV[0] =~ s/\.fq$//; $ARGV[0] =~ s/\.$//; 
@R1_FQ = split(/,/, $ARGV[1]);
@R2_FQ = split(/,/, $ARGV[2]);
if (length(@R1_FQ) != length(@R2_FQ)) {die "\n\nERROR: must specify the same number of read 1 and read 2 fastq files!\n"};
if (defined $opt{'r'}) {
	@RUNIDS = split(/,/, $opt{'r'});
	if (length(@RUNIDS) != length(@R1_FQ)) {
		die "\n\nERROR: When specifying run IDs you must specify the same number as input fastq files\n";
	}
}
open O1, "| $gzip > $ARGV[0].1.fq.gz";
open O2, "| $gzip > $ARGV[0].2.fq.gz";
for ($id = 0; $id <= @R1_FQ; $id++) {
	if ($R1_FQ[$id] =~ /\.gz$/) {open R1, "$zcat $R1_FQ[$id] |"} else {open R1, "$R1_FQ[$id]"};
	if ($R2_FQ[$id] =~ /\.gz$/) {open R2, "$zcat $R2_FQ[$id] |"} else {open R2, "$R2_FQ[$id]"};
	if (defined $opt{'r'}) {$runID = $RUNIDS[$id]} else {$runID = $id};
	while ($r1tag = <R1>) {
		chomp $r1tag; $r1seq = <R1>; $null = <R1>; $r1qual = <R1>;
		$r2tag = <R2>; chomp $r2tag; $r2seq = <R2>; $null = <R2>; $r2qual = <R2>;
		$r1tag =~ s/\#.+$//; $r2tag =~ s/\#.+$//;
		print O1 "$r1tag.$runID#0/1\n$r1seq\+\n$r1qual";
		print O2 "$r2tag.$runID#0/2\n$r2seq\+\n$r2qual";
	} close R1; close R2;
} close O1; close O2;

}
1;
