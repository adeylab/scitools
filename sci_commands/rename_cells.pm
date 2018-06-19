package sci_commands::rename_cells;

use sci_utils::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("rename_cells");

sub rename_cells {

@ARGV = @_;

# Defaults
$naming_scheme = "Cell_[number]";

getopts("O:A:N:R:Dx", \%opt);

$die2 = "
scitools rename-cells [options] [input file] (additional input file) etc...
   
Will rename cells according to an annotation, or create an annotation.

Options:
   -O   [STR]   Output prefix (default is input file without file
                extension and 'renamed')
   -R   [STR]   Annotation file for renaming (uses the provided names)
                  Note: if a cell is not in this file, it will be excluded!

   -N   [STR]   New naming schema, includes specific variables:
                  text       = includes in name
                  [number]   = number of the barcode (required)
                  [annot]    = annotation from option -A
                  [orig]     = the original cell ID
                  Example:     MySample_[annot]_[number]
                  Default:     $naming_scheme
   -A   [STR]   Annotation file to include as [annot] in new names
   -x           Exclude cells not in the annot file (-A)
                  (def = annot: Cell)
   -D           Do not create new files, just an annotaiton file to later
                  be used as the -R option

";

if (!defined $ARGV[0]) {die $die2};

if (defined $opt{'O'}) {
	$out = $opt{'O'};
} else {
	$out = $ARGV[0];
	$out =~ s/(\.fq\.gz|\.fq|\.fastq|\.fastq\.gz|\.bam|\.sam|\.matrix|\.tfidf|\.tf|\.values|\.annot|\.annotation|\.dims|\.LSI)$//i;
}

if (defined $opt{'R'} && defined $opt{'A'}) {die "ERROR: Both -R and -A cannot be defined. -R should be used as the only option (other than -O if used).\n"};

%ORIGINAL_newID = ();
if (defined $opt{'R'}) {
	print STDERR "\nUsing existing renaming annotaiton file: $opt{'R'}\n";
	read_annot($opt{'R'});
	foreach $cellID (keys %CELLID_annot) {
		$ORIGINAL_newID{$cellID} = $CELLID_annot{$cellID};
	}
} else {
	# parse scheme
	if (defined $opt{'N'}) {$naming_scheme = $opt{'N'}};
	if ($naming_scheme !~ /\[number\]/) {die "ERROR: Naming scheme must contain the [number] field.\n"};
	if ($naming_scheme =~ /\[annot\]/ && !defined $opt{'A'}) {die "\nERROR: When specifying [annot] in naming scheme, must provide an annot file as -A.\n"};
	if (defined $opt{'A'}) {read_annot($opt{'A'})};
}

$newID_number = 0;
for ($in_file = 0; $in_file < @ARGV; $in_file++) {

	# figure out file input
	if ($ARGV[$in_file] =~ /(\.fq\.gz|\.fq|\.fastq|\.fastq\.gz)$/i) { # fastq
		if ($ARGV[$in_file] =~ /\.gz$/) {
			open IN, "$zcat $ARGV[$in_file] |";
		} else {
			open IN, "$ARGV[$in_file]";
		}
		if (!defined $opt{'D'}) {
			open OUT, "| $gzip > $out.renamed.fq.gz";
		}
		while ($tag = <IN>) {
			chomp $tag; $seq = <IN>; chomp $seq; $null = <IN>; $qual = <IN>; chomp $qual;
			$tag =~ s/^\@//; ($origID,$tail_info) = split(/:/, $tag);
			if (!defined $ORIGINAL_newID{$origID}) {
				$ORIGINAL_newID{$origID} = rename_cell($origID);
			}
			if (!defined $opt{'D'} && $ORIGINAL_newID{$origID} !~ /00EXCL00/) {
				$newID = $ORIGINAL_newID{$origID};
				print OUT "\@$newID:$tail_info\n$seq\n\+\n$qual\n";
			}
		} close IN;
		if (!defined $opt{'D'}) {close OUT};
	} elsif ($ARGV[$in_file] =~ /(\.bam|\.sam)$/i) { # bam/sam
		if ($ARGV[$in_file] =~ /\.bam$/) {
			open IN, "$samtools view -h $ARGV[$in_file] |";
		} else {
			open IN, "$ARGV[$in_file]";
		}
		if (!defined $opt{'D'}) {
			open OUT, "| $samtools view -bS - > $out.renamed.bam 2>/dev/null";
		}
		while ($l = <IN>) {
			chomp $l;
			@P = split(/\t/, $l);
			if ($P[0] =~ /^\@/) {
				if ($P[0] =~ /\@RG/) {
					$RG_lines = "TRUE";
					$origID = $P[1]; $origID =~ s/^ID://;
					$ORIGINAL_newID{$origID} = rename_cell($origID);
					$newID = $ORIGINAL_newID{$origID};
					$out_line .= "\@RG\tID:$newID\tSM:$newID\tLB:$newID\tPL:SCI\n";
				} else {
					$out_line = $l;
				}
			} else {
				($origID,$tail_info) = split(/:/, $P[0]);
				if ($RG_lines eq "TRUE" && !defined $ORIGINAL_newID{$origID}) {
					$out_line = "00EXCL00";
				} else {
					if (!defined $ORIGINAL_newID{$origID}) {
						$ORIGINAL_newID{$origID} = rename_cell($origID);
					}
					$newID = $ORIGINAL_newID{$origID};
					$P[0] = "$newID:$tail_info";
					for ($i = 10; $i < @P; $i++) {
						if ($P[$i] =~ /^RG:Z:/) {
							$P[$i] = "RG:Z:$newID";
						}
					}
					$out_line = join("\t", @P);
				}
			}
			if (!defined $opt{'D'} && $out_line !~ /00EXCL00/) {
				print OUT "$out_line\n";
			}
		} close IN;
		if (!defined $opt{'D'}) {close OUT};
	} elsif ($ARGV[$in_file] =~ /(\.matrix|\.tfidf|\.tf|\.LSI)$/i) { # matrix
		@FP = split(/\./, $ARGV[$in_file]); $extension = pop(@FP);
		open IN, "$ARGV[$in_file]";
		if (!defined $opt{'D'}) {
			open OUT, ">$out.renamed.$extension";
		}
		$header = <IN>; chomp $header;
		@OH = split(/\t/, $header);
		$new_header = ""; @NH = ();
		for ($i = 0; $i < @OH; $i++) {
			$origID = $OH[$i];
			if (!defined $ORIGINAL_newID{$origID}) {
				$ORIGINAL_newID{$origID} = rename_cell($origID);
			}
			$newID = $ORIGINAL_newID{$origID};
			if ($newID !~ /00EXCL00/) {
				$new_header .= "$newID\t";
				push @NH, $newID;
			} else {
				$OH[$i] = "00EXCL00";
			}
		} $new_header =~ s/\t$//;
		if (!defined $opt{'D'}) {
			print OUT "$new_header\n";
			while ($l = <IN>) {
				chomp $l;
				@P = split(/\t/, $l);
				$siteID = shift(@P);
				$out_line = "$siteID";
				for ($i = 0; $i < @OH; $i++) {
					if ($OH[$i] !~ /00EXCL00/) {
						$out_line .= "\t$P[$i]";
					}
				}
				print OUT "$new_header\n";
			}
			close OUT;
		} close IN;
	} elsif ($ARGV[$in_file] =~ /(\.dims|\.values|\.annot|\.annotation)$/i) { # dims/values/annot
		@FP = split(/\./, $ARGV[$in_file]); $extension = pop(@FP);
		open IN, "$ARGV[$in_file]";
		if (!defined $opt{'D'}) {
			open OUT, ">$out.renamed.$extension";
		}
		while ($l = <IN>) {
			chomp $l;
			@P = split(/\t/, $l);
			$origID = $P[0];
			if (!defined $ORIGINAL_newID{$origID}) {
				$ORIGINAL_newID{$origID} = rename_cell($origID);
			}
			$newID = $ORIGINAL_newID{$origID};
			if (!defined $opt{'D'} && $newID !~ /00EXCL00/) {
				$P[0] = $newID;
				$out_line = join("\t", @P);
				print OUT "$out_line\n";
			}
		} close IN;
		if (!defined $opt{'D'}) {close OUT};
	} else {
		print STDERR "ERROR: Cannot determine input file type of $ARGV[$in_file] based ont he file extension.\n";
	}
}

if (!defined $opt{'R'}) { # print out cross-mapping annot files
	open DIR, ">$out.rename.annot";
	open REV, ">$out.rename_reverse.annot";
	foreach $origID (keys %ORIGINAL_newID) {
		if ($ORIGINAL_newID{$origID} !~ /00EXCL00/) {
			print DIR "$origID\t$ORIGINAL_newID{$origID}\n";
			print REV "$ORIGINAL_newID{$origID}\t$origID\n";
		}
	}
	close DIR; close REV;
}

}

# subroutine specific for this command:
sub rename_cell { # note - include an error message if -R is specified and a cellID is found that is not in the -R annot
	if (!defined $opt{'R'}) {
		$newID_number++;
		$rename_origID = $_[0];
		if (defined $opt{'x'}) {
			$rename_annot = "00EXCL00";
		} else {
			$rename_annot = "Cell";
		}
		if (defined $CELLID_annot{$rename_origID}) {
			$rename_annot = $CELLID_annot{$rename_origID};
		}
		$rename_newID = $naming_scheme;
		$rename_newID =~ s/\[annot\]/$rename_annot/;
		$rename_newID =~ s/\[orig\]/$rename_origID/;
		$rename_newID =~ s/\[number\]/$newID_number/;
		return $rename_newID;
	} else {
		return "00EXCL00";
	}
}

1;
