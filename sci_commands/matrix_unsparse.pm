package sci_commands::matrix_unsparse;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_unsparse");

sub matrix_unsparse {

@ARGV = @_;
getopts("O:z", \%opt);

$die2 = "

scitools matrix-unsparse (options) [sparse_matrix] [cellIDs.txt] [rowIDs.txt or bed]
   or    unsparse
         make-unsparse

Options:
   -O   [STR]   Output prefix (def = sparse_matrix prefix)
   -z           Gzip output
   
";

if (!defined $ARGV[2]) {die $die2};

if (!defined $opt{'O'}) {
	$opt{'O'} = $ARGV[0];
	$opt{'O'} =~ s/\.gz//;
	$opt{'O'} =~ s/\.sparse_matrix//;
}

if ($ARGV[1] =~ /\.gz$/) {
	open CID, "$zcat $ARGV[1] |";
} else {
	open CID, "$ARGV[1]";
}

@CELLIDs = (); push @CELLIDs, "null";
$header = "";
while ($l = <CID>) {
	chomp $l;
	push @CELLIDs, $l;
	$header .= "$l\t";
} close CID;
$header =~ s/\t$//;

if ($ARGV[2] =~ /\.gz$/) {
	open RID, "$zcat $ARGV[2] |";
} else {
	open RID, "$ARGV[2]";
}

@ROWIDs = (); push @ROWIDs, "null";
while ($l = <RID>) {
	chomp $l;
	if ($ARGV[2] =~ /bed/) {
		@P = split(/\t/, $l);
		$rowName = $P[0]."_".$P[1]."_".$P[2];
		push @ROWIDs, $rowName;
	} else {
		push @ROWIDs, $l;
	}
} close RID;

if ($ARGV[0] =~ /\.gz$/) {
	open MTX, "$zcat $ARGV[0] |";
} else {
	open MTX, "$ARGV[0]";
}
$null = <MTX>;
while ($l = <MTX>) {
	chomp $l;
	($rowID,$cellID,$value) = split(/\s/, $l);
	$COUNTS{$rowID}{$cellID} = $value;
} close MTX;


if (defined $opt{'z'}) {
	open OUT, "| $gzip > $opt{'O'}.matrix";
} else {
	open OUT, ">$opt{'O'}.matrix";
}

print OUT "$header\n";
for ($rowID = 1; $rowID < @ROWIDs; $rowID++) {
	print OUT "$ROWIDs[$rowID]";
	for ($cellID = 1; $cellID < @CELLIDs; $cellID++) {
		if (defined $COUNTS{$rowID}{$cellID}) {
			print OUT "\t$COUNTS{$rowID}{$cellID}";
		} else {
			print OUT "\t0";
		}
	}
	print OUT "\n";
} close OUT;

}
1;
