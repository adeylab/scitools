package sci_commands::matrix_tfidf;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("matrix_tfidf");

sub matrix_tfidf {

@ARGV = @_;

getopts("O:L:Tz", \%opt);

$die2 = "
scitools matrix-tfidf [options] [counts matrix or sparseMatrix values file]
   or    tfidf

If a regular or a sparseMatrix is provided, the optput will be in the same
format as the input.

Options:
   -O   [STR]   Output prefix (default is [input].tfidf)
   -T           Use the natural log transformed TF value.
   -L   [BASE]  Log norm with base specified.
                N for natural, def = no additional log norm
   -z           Gzip output

";

if (!defined $ARGV[0]) {die $die2};

if ($ARGV[0] =~ /\.matrix/ || $ARGV[0] =~ /\.matrix\.gz/) {

	if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.gz$//; $opt{'O'} =~ s/\.matrix$//};

	read_matrix_stats($ARGV[0]);

	if ($ARGV[0] =~ /\.gz$/) {
		open IN, "$zcat $ARGV[0] |";
	} else {
		open IN, "$ARGV[0]";
	}
	
	if (defined $opt{'z'}) {
		open OUT, "| $gzip > $opt{'O'}.tfidf.gz";
	} else {
		open OUT, ">$opt{'O'}.tfidf";
	}
	
	$h = <IN>; print OUT "$h";
	while ($l = <IN>) {
		chomp $l;
		@P = split (/\t/, $l);
		$rowID = shift(@P);
		print OUT "$rowID";
		for ($i = 0; $i < @P; $i++) {
			if ($P[$i]>0) {
				$tf = ($P[$i]/$MATRIX_CellID_signal{$MATRIX_COLNAMES[$i]});
				$idf = (log(1+($matrix_colNum/($MATRIX_feature_signal{$rowID}+1))));
				if (defined $opt{'T'}) {
					$raw_score = log($tf * 100000)*$idf;
                } else {
                    $raw_score = $tf*$idf;
                }
				if (defined $opt{'L'}) {
					if ($opt{'L'} =~ /N/i) {
						$raw_score = log($raw_score);
					} else {
						$raw_score = log($raw_score)/log($opt{'L'})
					}
				}
				$score = sprintf("%.6f", $raw_score);
			} else {
				$score = "0";
			}
			print OUT "\t$score";
		}
		print OUT "\n";
	} close IN; close OUT;

} elsif ($ARGV[0] =~ /(sparseMatrix|values)/) {
	
	if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]; $opt{'O'} =~ s/\.gz$//; $opt{'O'} =~ s/\.values$//; $opt{'O'} =~ s/\.sparseMatrix$//};
	
	if ($ARGV[0] =~ /\.gz$/) {
		open IN, "$zcat $ARGV[0] |";
	} else {
		open IN, "$ARGV[0]";
	}
	while ($l = <IN>) {
		chomp $l;
		($row,$col,$val) = split(/\s/, $l);
		if (!defined $COLID_sum{$col}) {$colNum++};
		$ROWID_sum{$row}+=$val;
		$COLID_sum{$col}+=$val;
	} close IN;
	
	if (defined $opt{'z'}) {
		open OUT, "| $gzip > $opt{'O'}.sparseMatrix.tfidf.gz";
	} else {
		open OUT, ">$opt{'O'}.sparseMatrix.tfidf";
	}
	
	if ($ARGV[0] =~ /\.gz$/) {
		open IN, "$zcat $ARGV[0] |";
	} else {
		open IN, "$ARGV[0]";
	}
	while ($l = <IN>) {
		chomp $l;
		($row,$col,$val) = split(/\s/, $l);
		$tf = ($val/$COLID_sum{$col});
		$idf = (log(1+($colNum/($ROWID_sum{$row}+1))));
			if (defined $opt{'T'}) {
                $raw_score = log($tf * 100000)*$idf;
            } else {
        		$raw_score = $tf*$idf;
            }
		if (defined $opt{'L'}) {
			if ($opt{'L'} =~ /N/i) {
				$raw_score = log($raw_score);
			} else {
				$raw_score = log($raw_score)/log($opt{'L'})
			}
		}
		$score = sprintf("%.6f", $raw_score);
		print OUT "$row\t$col\t$score\n";
	} close IN; close OUT;
	
}

}
1;
