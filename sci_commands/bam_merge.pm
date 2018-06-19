package sci_commands::bam_merge;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_merge");

sub bam_merge {

@ARGV = @_;
getopts("s:", \%opt);

$die2 = "
scitools bam-merge [options] [output bam] [input bam1] [input bam2] ...
   or    merge-bam

Simple wrapper for samtools merge
(use samtools merge for more complicated operations)

Options:
   -s   [STR]   Samtools call (def = $samtools)

";

if (defined $opt{'O'}) {unshift @ARGV, $opt{'O'}};
if (!defined $ARGV[1]) {die $die2};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};

$ARGV[0] =~ s/\.bam$//;
$args = "";
for ($i = 1; $i < @ARGV; $i++) {$args .= "$ARGV[$i] "};

system("$samtools merge $ARGV[0].bam $args");

}
1;
