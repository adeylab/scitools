package sci_commands::bam_filter;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_filter");

sub bam_filter {

@ARGV = @_;
getopts("s:O:A:a:N:C:c:L:v", \%opt);

$die2 = "
scitools bam-filter [options] [bam file]
   or    filter-bam

Filters a SCI bam file by a variety of parameters

Options:
   -O   [STR]   Output prefix (default is bam file prefix)
   -A   [STR]   Annotation file
   -a   [STR]   Comma separated list of annotations to include (requires -A)
   -N   [INT]   Minimum read count for each barcode to keep
                (performs much faster if -C is provided)
   -C   [STR]   Complexity file (can be comma separated)
   -c   [STR]   min,max percent complexity to include (requires -C)
                if just one number, will assume max
   -L   [STR]   File listing indexes to include
   -v           Keep the reciprocal of specified filters
   -s   [STR]   Samtools call (def = $samtools)
   
";

if (!defined $ARGV[0]) {die $die2};
$filter_params = "";
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.bam$//; $filter_params .= "O=$opt{'O'}_"};
if (defined $opt{'c'} && !defined $opt{'C'}) {die "\nMust provide a complexity file (-C) if specifying min and max complexity to filter (-c)!\n$die2"};
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotaiton file (-A) if specifying annotations to filter (-a)!\n$die2"};
if (defined $opt{'c'}) {
	if ($opt{'c'} =~ /,/) {
		($min_compl,$max_compl) = split(/,/, $opt{'c'});
	} else {
		$min_compl = 0;
		$max_compl = $opt{'c'};
	}
	$filter_params .= "c=$opt{'c'}_";
}
if (defined $opt{'s'}) {$samtools = $opt{'s'}};
if (defined $opt{'A'}) {read_annot($opt{'A'}); $filter_params .= "A=$opt{'A'}_"};
if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
	$filter_params .= "a=$opt{'a'}_";
}

if (defined $opt{'C'}) {read_complexity($opt{'C'}); $filter_params .= "C=$opt{'C'}_"};
if (defined $opt{'L'}) {
	open LIST, "$opt{'L'}";
	while ($list_line = <LIST>) {chomp $list_line; $CELLID_inList{$list_line} = 1};
	close LIST;
	$filter_params .= "L=$opt{'L'}_";
}

$filter_params =~ s/_$//;

open LOG, ">$opt{'O'}.filt.log";
$ts = localtime(time);
print LOG "$ts\tscitools bam-filter called:\n";
foreach $option (keys %opt) {print LOG "   $option   $opt{$option}\n"};

$included_reads = 0; $total_reads = 0;
if (defined $opt{'N'} && !defined $opt{'C'}) {
	%CELLID_uniq_reads = ();
	open IN, "$samtools view $ARGV[0] 2>/dev/null |";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		if ($P[2] !~ /(M|Y|L|K|G|Un|Random|Alt|_)/i) {
			$cellID = $P[0];
			$cellID =~ s/:.*$//;
			if (!defined $CELLID_uniq_reads{$cellID}) {
				$CELLID_uniq_reads{$cellID} = 1;
			} else {
				$CELLID_uniq_reads{$cellID}++;
			}
		}
	} close IN;
}

$out_header = "";
open IN, "$samtools view -h $ARGV[0] 2>/dev/null |";
open OUT, "| $samtools view -bS - > $opt{'O'}.filt.bam 2>/dev/null";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$lineType = "READ";
	if ($P[0] =~ /^\@/) {
		if ($P[0] =~ /\@RG/) {
			$lineType = "RG";
		} else {
			$out_header .= "$l\n";
			$lineType = "HEADER";
		}
	}
	
	if ($lineType eq "READ" || $lineType eq "RG") {
	
		if ($lineType eq "READ" && $out_header ne "done") {
			print OUT "$out_header\@PG\tID:scitools_bam-filter_$filter_params\tVN:$version\n";
			$out_header = "done";
		}
	
		if ($lineType eq "READ") {
			$cellID = $P[0];
			$cellID =~ s/:.*$//;
		} else {
			($null,$cellID) = split(/:/, $P[1]);
		}
		
		@INCLUDE_FLAGS = ();
		
		if (defined $opt{'a'}) {
			if (defined $CELLID_annot{$cellID}) {
				$annot = $CELLID_annot{$cellID};
				if (defined $ANNOT_include{$annot}) {
					push @INCLUDE_FLAGS, 1;
				} else {push @INCLUDE_FLAGS, 0};
			} else {push @INCLUDE_FLAGS, 0};
		}
		
		if (defined $opt{'N'}) {
			if ($CELLID_uniq_reads{$cellID}>=$opt{'N'}) {
				push @INCLUDE_FLAGS, 1;
			} else {push @INCLUDE_FLAGS, 0};
		}
		
		if (defined $opt{'L'}) {
			if (defined $CELLID_inList{$cellID}) {
				push @INCLUDE_FLAGS, 1;
			} else {push @INCLUDE_FLAGS, 0};
		}
		
		if (defined $opt{'c'}) {
			if ($CELLID_complexity{$cellID}<=$max_compl&&$CELLID_complexity{$cellID}>=$min_compl) {
				push @INCLUDE_FLAGS, 1;
			} else {push @INCLUDE_FLAGS, 0};
		}
		
		if (defined $opt{'v'}) {
			$print = 1;
			foreach $flag (@INCLUDE_FLAGS) {if ($flag>0) {$print = 0}};
			if ($print>0) {
				if ($lineType eq "READ") {
					print OUT "$l\n";
					$included_reads++;
				} else {
					$out_header .= "$l\n";
				}
			}
		} else {
			$print = 1;
			foreach $flag (@INCLUDE_FLAGS) {if ($flag<1) {$print = 0}};
			if ($print>0) {
				if ($lineType eq "READ") {
					print OUT "$l\n";
					$included_reads++;
				} else {
					$out_header .= "$l\n";
				}
			}
		}
		if ($lineType eq "READ") {$total_reads++};
	}
} close IN; close OUT;

$ts = localtime(time);
print LOG "$ts\tTotal reads: $total_reads, Included reads: $included_reads\n";
close LOG;

}
1;
