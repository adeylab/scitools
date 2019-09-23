package sci_commands::met_fastqtrim;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("met_fastqtrim");

sub met_fastqtrim {

@ARGV = @_;

getopts("O:n:p:", \%opt);
$die2 = "
scitools met_fastqtrim [options] [FastQ File]

sciMET/sciNOMe Trimming of Reads to remove adaptor content.
[FastQ File] 			= 		Generated through fastq-dump script.



Options:

   -o 	[STR]	Output Prefix. 
   				(Default: all text within last \"/\" and first \"\.\"\ of [Input File])
   -O 	[STR] 	Output Directory
   				(Default: Current Directory)	
";

if (!defined $ARGV[0]) {die $die2};

if (!defined $opt{'O'}) {$opt{'O'} = getcwd()};
if (!defined $opt{'o'}) {$opt{'o'}=$ARGV[0]; my @o = split(/\//,$opt{'o'}); $opt{'o'}=$o[-1]; @o = split(/\./,$opt{'o'}); $opt{'o'}=$o[0]};

system("/home/users/mulqueen/tools/trim_galore_zip/trim_galore --gzip -a AGATCGGAAGAGC -O . $ARGV[0]"); 
}
1;
}