package sci_commands::fastq_align;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("fastq_align");

sub fastq_align {

@ARGV = @_;

# Defaults
$threads = 1;
$memory = "2G";

getopts("t:b:s:A:L:rO:m:", \%opt);

$die2 = "
scitools fastq-align [options] [bwa reference] [output_prefix] [read1.fq] (read2.fq)
   or    align-fastq
   or    align

Produces a sorted bam file. Read 2 is optional.

Options:
   -t   [INT]   Threads for alignment (def = $threads)
   -b   [STR]   Bwa call (def = $bwa)
   -s   [STR]   Samtools call (def = $samtools)
   -m   [MEM]   Samtools sort mex memory per thread, K/M/G (def = $memory)

Bwa reference shortcuts:
$ref_shortcuts

";

if (!defined $ARGV[2]) {die $die2};
if (defined $REF{$ARGV[0]}) {
	$ref_file = $REF{$ARGV[0]};
} else {
	$ref_file = $ARGV[0];
}

if (defined $opt{'O'}) {
	$out_prefix = $opt{'O'};
} else {
	$out_prefix = $ARGV[1];	
}
$out_prefix =~ s/\.bam$//;
$out_prefix =~ s/\.$//;

if (defined $opt{'b'}) {$bwa = $opt{'b'}};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (!defined $opt{'t'}) {$opt{'t'} = $threads};
if (!defined $opt{'m'}) {$opt{'m'} = $memory};

if (defined $ARGV[3]) {
	$align_command = "$bwa mem -t $opt{'t'} $ref_file $ARGV[2] $ARGV[3] 2>> $out_prefix.align.log | $samtools view -bSu - 2>> $out_prefix.align.log | $samtools sort -m $opt{'m'} -T $out_prefix.TMP - > $out_prefix.bam 2>> $out_prefix.align.log";
} else { # single ended
	$align_command = "$bwa mem -t $opt{'t'} $ref_file $ARGV[2] 2>> $out_prefix.align.log | $samtools view -bSu - 2>> $out_prefix.align.log | $samtools sort -m $opt{'m'} -T $out_prefix.TMP - > $out_prefix.bam 2>> $out_prefix.align.log";
}

#print "Running: $align_command\n";
system($align_command);

}
1;
