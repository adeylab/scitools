package sci_commands::dependencies;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("dependencies");

sub dependencies {

@ARGV = @_;
getopts("R:F:O:o:I:1:2:A:B:P:W:K:", \%opt);

$die2 = "
scitools dependencies [options] [report output file]
  or     depend

This will run a check to make sure the required dependencies
are present and command-line callable. If defaults are not
used (shown below), then specify the executables using
options. It will then verify R packages are installed and
output all software & R package status in a report file.

Options:
   -s   [STR]   Samtools call (def = $samtools)
   -b   [STR]   Bedtools call (def = $bedtools)
   -B   [STR]   Bwa call (def = $bwa)
   -m   [STR]   Macs2 call (def = $macs2)
   -R   [STR]   Rscript call (def = $Rscript)
   -S   [STR]   Scitools call (def = $scitools)
                (for self-calling)
   -P   [STR]   Python call (def = $Pscript)
   -W   [STR]   Bowtie2 call (def = $bowtie2)
   -K   [STR]   Bismark call (def = $bismark)

Hardcoded:      Gzip (def = $gzip)
                Zcat (def = $zcat)

";

if (!defined $ARGV[0]) {die $die2};
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (defined $opt{'b'}) {$bedtools = $opt{'b'}};
if (defined $opt{'B'}) {$bwa = $opt{'B'}};
if (defined $opt{'m'}) {$macs2 = $opt{'m'}};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'S'}) {$scitools = $opt{'S'}};
if (defined $opt{'P'}) {$Pscript = $opt{'P'}};
if (defined $opt{'W'}) {$bowtie2 = $opt{'W'}};
if (defined $opt{'K'}) {$bismark = $opt{'K'}};

$executable_warnings = 0;
$R_package_warnings = 0;

open OUT, ">$ARGV[0]";
$ts = localtime(time);
print OUT "$ts scitools dependencies check.

Checking command-line executables:\n";

@DEPENDENCIES = ($samtools,$bedtools,$bwa,$macs2,$Rscript,$scitools,$gzip,$zcat,$Pscript,$bowtie2,$bismark);
foreach $dependency (@DEPENDENCIES) {
	$path = "";
	open WHICH, "which $dependency 2>/dev/null |";
	$path = <WHICH>; close WHICH;
	if ($path =~ /\//) {
		print OUT "  \"$dependency\" found. Path = $path";
	} else {
		print OUT "  WARNING: \"$dependency\" is not command-line callable!\n    Ensure that the executable is correct. If the tool is not installed, scitools functions that use this tool will not be useable.\n";
		$executable_warnings++;
	}
}

print OUT "\nChecking R packages:\n"; close OUT;

@R_PACKAGES = ("ggplot2","svd","Rtsne","methods","dbscan","chromVAR","chromVARmotifs","irlba","princurve");
open R, ">$ARGV[0].r";
foreach $package (@R_PACKAGES) {
	print R "if (!require($package)) {print(\"  R: WARNING: $package cannot be loaded! Some scitools functions will not be useable unless it is installed.\",quote=FALSE)} else {print(\"  $package found and loadable!\",quote=FALSE)}\n";
} close R;
system ("$Rscript $ARGV[0].r >> $ARGV[0] 2>/dev/null && rm -f $ARGV[0].r");

open IN, "$ARGV[0]";
while ($l = <IN>) {if ($l =~ /R: WARNING:/) {$R_package_warnings++}};
close IN;

open OUT, ">>$ARGV[0]";
print OUT "\nDependency check completed:
  $executable_warnings failed command-line executables
  $R_package_warnings failed R packages
"; close OUT;

print "\nDependency check completed with: $executable_warnings failed command-line executables, and $R_package_warnings failed R packages.\n\n";

}
1;
