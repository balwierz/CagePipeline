package promoteromeLib;
use strict;
use warnings;
#use Math::GSL::SF qw(:all);
our $VERSION = '1.00';
use base 'Exporter';
our @EXPORT = qw(readCumulative regressPowerLaw h2counts gnuplotRevcumLine readRawSummary
	normalize addSample writeSummary processSAM processBed corr txt libSamBarStructureFromFile
	writeSummaryNormAndRaw cmpID cumulative_rev_unnorm readFit writeNormRaw cmpKawajiKey );

sub readCumulative 
{
	my($file) = @_;
	my @data;
	open(F, $file) or die "Cannot open $file $!\n";
	while(<F>)
	{
		chomp;
		my($x, $y) = split(/\t/);
		push @data, [($x, $y)];
	}
	close F;
	#print STDERR "Read " . scalar @data . " lines\n";
	return \@data;
}

sub cmpID
{
	my($hRef, $a, $b) = @_;
	if($hRef->{$a}->{chr} ne $hRef->{$b}->{chr})
	{
		return $hRef->{$a}->{chr} cmp $hRef->{$b}->{chr};
	}
	elsif($hRef->{$a}->{str} ne $hRef->{$b}->{str})
	{
		return $hRef->{$a}->{str} cmp $hRef->{$b}->{str};
	}
	return $hRef->{$a}->{pos} <=> $hRef->{$b}->{pos};
}

sub regressPowerLaw
{
	# works per one sample
	# takes a reference to an array of two element (x,y) arrays
	# slope needs to be a string because of rounding
	my ($dataRef, $region_startX, $region_startY, $slope) = @_;
	
	my %zetaTable = 
	(
		"1.25", 4.59511182584295,
		"1.15", 7.25469458506813,
		"1.16", 6.83874082424536,
		"1.2",  5.59158,
		"2.15", 1.524,
		"2.25", 1.46,
		"2.16", 1.517,
		"2.2",  1.491,
		"2",    1.645
	);
	my $zetaSlope = $zetaTable{$slope};
	die "Do not have zeta in the table $.\n" if not defined $zetaSlope;
		#Math::GSL::SF::gsl_sf_zeta($slope);
	my $mux   = 0;
	my $muy   = 0;
	my $varx  = 0;
	my $vary  = 0;
	my $covar = 0;
	my $n     = 0;
	for(my $i=0; $i<@{$$dataRef}; $i++)
	{
		#print $$dataRef->[$i]->[0], "\t";
		my $x = $$dataRef->[$i]->[0];
		my $y = $$dataRef->[$i]->[1];
		if ($x > $region_startX && $y > $region_startY)
		{
			#print $x, "\t";
			$x = log($x);
			$y = log($y);
			$mux   += $x;
			$muy   += $y;
			$varx  += $x * $x;
			$vary  += $y * $y;
			$covar += $x * $y;
			++$n;
		}
	}
	# we demand at least 20 points to fit the line robustly
	#print $n, "\n";
	if($n < 14)
	{
		print "Not enough points\n";
		return(0, 0, 0, 0);
	}
	$mux   /= $n;
	$muy   /= $n;
	$varx  /= $n;
	$vary  /= $n;
	$covar /= $n;
	$varx  -= $mux * $mux;
	$vary  -= $muy * $muy;
	$covar -= $mux * $muy;
	# we fit y = b + a*x  // a is negative in our case.
	# the code below for a is actually calculating principal
	# direction that is going down.
	
	my $eigVal1 = 0.5 * ($varx + $vary + sqrt(($varx-$vary)**2 + 4*$covar**2));
	my $eigVal2 = 0.5 * ($varx + $vary - sqrt(($varx-$vary)**2 + 4*$covar**2));
		
	my $a =
		($vary - $varx) / (2.0 * $covar) -
		sqrt(1 +
				 (($vary - $varx) / (2.0 * $covar)) *
				 (($vary - $varx) / (2.0 * $covar)));
	my $b    = $muy - $a * $mux;
	my $beta = $slope / (-$a);	# slope is the desired revcum slope and positive (*=-1 applied); a is revcum slope and negative.
	my $lambda =
		exp($b / (-$a) + log($slope * $zetaSlope) / (-$a) - log(1000000.0) / (-$a));
		#This  constant is zeta(1.25) = Sum[n^(-1.25), {n, 1, Infinity}]
	my $newten     = exp(log(10. / $lambda) / $beta);
	my $newone     = exp(log(1.0 / $lambda) / $beta);
	my $newhundred = exp(log(100. / $lambda) / $beta);
	#print join(" ", txt($mux), txt($muy), txt($a), txt($b), txt($beta), txt($lambda), 
	#	txt($newone), txt($newten), txt($newhundred)), "\n";
	my $txt = "Eig1: $eigVal1 Eig2: $eigVal2 Ratio: " . $eigVal1/$eigVal2 . " a: $a b: $b beta: $beta lambda: $lambda" ;
	return($a, $b, $beta, $lambda, $txt);
}

sub txt
{
	return sprintf("%8g", $_[0]);
}

sub gnuplotRevcumLine
{
	my($file, $outF, $a, $b) = @_;
	open(G, "|GDFONTPATH=/usr/share/fonts/liberation GNUPLOT_DEFAULT_GDFONT=LiberationSans-Regular gnuplot") or die;
	print G "set terminal png\n";
	print G "set output '$outF'\n";
	print G "set log\n";
	print G "plot [0.9:*][1:*] '$file' w l notitle, exp($b)*x**$a\n";
	close G;
}

sub readRawSummary
{
	my($file) = @_;
	my @result;
	open(F, $file) or die "$! $file $.";
	while(<F>)
	{
		next if /^#/;
		chomp;
		my($id, $chr, $str, $pos, $raw) = split(/\t/);
		push @result, [($id, $raw)];
	}
	close F;
	return \@result;
}


sub normalize
{
	my ($exprRef, $beta, $lambda) = @_;
	print "$beta $lambda\n";
	$beta   = 1 / $beta;
	$lambda = 1 / $lambda;
	foreach my $id (keys %$exprRef)
	{
		if ($exprRef->{$id}->{count} == 0 )
		{
			$exprRef->{$id}->{norm} = 0;
		}
		else
		{
			$exprRef->{$id}->{norm} = exp(log($exprRef->{$id}->{count} * $lambda) * $beta);
		}
	}
}

#sub addSample
#{
	#my ($allExprRef, $exprRef) = @_;
#}

sub writeNormRaw
{
	my($h, $file, $dir, $sample) = @_;
	open(my $fh, ">", $file) or die;
	# find out all the chromosomes
	my %allChrStr = ();
	foreach my $key (keys %$h)
	{
		$allChrStr{$h->{$key}->{chr}.$h->{$key}->{str}} = 1;
	}
	foreach my $chrStr (keys %allChrStr)
	{
		mkdir $dir.$chrStr;
		undef $allChrStr{$chrStr};
		open($allChrStr{$chrStr}, ">", $dir.$chrStr.'/'.$sample) or die "cannot open " .$dir.$chrStr.'/'.$sample, "\n";
	}
	print STDERR "Sorting keys...";
	
	my @keys = sort {cmpKawajiKey($a, $b)} keys %$h;
	print STDERR " [DONE]\n";
	foreach my $id (@keys)
	{
		my $fid = $allChrStr{$h->{$id}->{chr}.$h->{$id}->{str}};
		my $line = join ("\t", $h->{$id}->{chr}, $h->{$id}->{str}, $h->{$id}->{pos},
			$h->{$id}->{count}, $h->{$id}->{norm}) . "\n";
		print $fid $line;
		print $fh $line;
	}
	close $fh;

	foreach my $chrStr (keys %allChrStr)
	{
		close($allChrStr{$chrStr});
	}
}

# this is old
sub writeSummaryNormAndRaw
{
	my ($file, $exprRef, $allBarcodes, $totals) = @_;
	# exprRef has structure ->{id}->{barcode}->[0 1] raw and normalized
	
	open(T, ">$file.totals") or die $!;
	foreach my $t (@$totals)
	{
		print T $t,"\n";
	}
	close T;
	
	open(F, ">", $file) or die $!;
	print F "#id\tchr\tstr\tpos";
	foreach my $b (@$allBarcodes)
	{
		print F "\t$b.norm";
	}
	foreach my $b (@$allBarcodes)
	{
		print F "\t$b.raw";
	}
	print F "\n";
	foreach my $id (sort keys %$exprRef )
	{
		my($chr, $str, $pos) = split("_", $id);
		my @line = ();
		foreach my $barcode (@$allBarcodes)
		{
			push @line, exists $exprRef->{$id}->{$barcode} ? $exprRef->{$id}->{$barcode}->[1] : "0";  #norm
		}
		foreach my $barcode (@$allBarcodes)
		{
			push @line, exists $exprRef->{$id}->{$barcode} ? $exprRef->{$id}->{$barcode}->[0] : "0";  #raw
		}
		print F join("\t", $id, $chr, $str, $pos, @line), "\n";
	}
	close F;	
}

sub writeSummary
{
	my ($file, $exprRef) = @_;
	open(F, ">", $file) or die $!;
	for(my $i=0; $i<@$exprRef; $i++)
	{
		my($chr, $str, $pos) = split("_", $exprRef->[$i]->[0]);
		print F join("\t", $exprRef->[$i]->[0], $chr, $str, $pos, $exprRef->[$i]->[1], $exprRef->[$i]->[2]), "\n";
	}
	close F;
}

sub libSamBarStructureFromFile
{
	my ($dir, $rawData) = @_;
	my @fileList;
	opendir(D, $dir) or die "$! $dir";
	while (my $f = readdir D)
	{
	  if ($f =~ /\.sam$/ || $f =~ /\.sam.gz$/)
	  {
	    my ($library, $sample, $barcode) = split(/\./, $f);
	    push @{$rawData->{$library}->{$sample}}, $barcode;
	    push @fileList, $f;
	  }
	}
	closedir D;
	@fileList; #yes, it is by value
}

sub processBed
{
	# bed uses 0-based coordinates. and indexes places in between nucleotides
	# (-p- bonds!)
	# so the bed coordinates a,b go to a+1,b in the standard (nucleotide) notation
	# It should be always that a+1==b, but we make an IF just in case it changes
	# in the future.
	# from Kawaji-san
	#It is formatted as BED ( http://genome.ucsc.edu/FAQ/FAQformat.html#format1 ), as indicated in suffix, which adopt zero-start. Accordingly, it should be interpreted as followings:
	#
	#    1    2    3     1-start coordinate (conventional)
	# 0    1    2    3   0-start cooridnate (BED or indication on UCSC browser)
	#----A----T----G----
	#         T----G----  (reads on forward strand)
	#----T----A           (reads on revers strand)
	#
	#When a read start 5'-TG... on forward strand, I put chr:1..2,+.
	#When a read start 5'-AT... on reverse strand, I put chr:1..2,-. 
	my($hRef, $f) = @_;
	open(F, "zcat '$f' |") or die;
	my %ret;
	#$ret{tagCount} = 0;
	$ret{mappedCount} = 0;
	#$ret{genome}   = undef;
	while(<F>)
	{
		chomp;
		my($chr, $beg, $end, $id, $count, $str) = split(/\s+/);
		my $tssPosition = $str eq '+' ? $beg+1 : $end;
		$hRef->{$id}->{count} = $count;
		$hRef->{$id}->{chr} = $chr;
		$hRef->{$id}->{str} = $str;
		$hRef->{$id}->{pos} = $tssPosition;
		$ret{mappedCount} += $count;
	}
	close F;
	print STDERR "Read ". $ret{mappedCount} . " tags\n";
	return \%ret;
}

sub h2counts
{
	my($h) = @_;
	my @ret;
	foreach my $v (values %$h)
	{
		#print $v; die;
		push @ret, $v->{count};
	}
	print STDERR scalar @ret . " values exported\n";
	return \@ret;
}

sub processSAM
{
	# Reads SAM file, returns number of tags from the header
	# stores data in $hRef in a format
	# $hRef->{$id}->{count}
	#	$hRef->{$id}->{chr} = $chr;
	#	$hRef->{$id}->{str} = $strand;
	#	$hRef->{$id}->{pos} = $tssPosition;
	#######################################################
	my($hRef, $f) = @_;
	open(F, "zcat $f |") or die;
	my %ret;
	$ret{tagCount} = 0;
	$ret{mappedCount} = 0;
	$ret{genome}   = undef;
	
	while (<F>)
	{
		chomp;
		if ($_ =~ /^\@/)    #comment
		{
			if($_ =~ /^\@CO\s+TagCount:(\d+)$/)
			{
				$ret{tagCount} = $1;
			}
			if($_ =~ /^\@SQ\s+(.*)/)
			{
				my(@fields) = split(/\s+/);
				foreach my $pair (@fields)
				{
					my($key, $value) = split(/:/,$pair);
					if($key eq "AS")  #genome
					{
						if(defined $ret{genome} && $ret{genome} ne $value)
						{
							die "File $f contains multiple genomes\n";
						}
						$ret{genome} = $value;
					}
				}
			}
			next;
		}
		#last; #uncomment if you want just a structure
		my ($seq, $flag, $chr, $leftpos, $quality, $cigar, $mrm, $mpos, $isize, $genomicSeq, $queryQuality, @rest) =
			split(/\t/);

		### get expected number of tags mapping to this place
		my ($xc, $xw);
		foreach my $entry (@rest)
		{
			my ($tag, $dataType, $value) = split(/\:/, $entry);

			#print $tag, " ", $dataType, " ", $value;		die;
			if ($tag eq "XC")    #global tag count
			{
				$xc = $value;      # we assume dataType == f float
			}
			elsif ($tag eq "XW")    # local mapping weight
			{
				$xw = $value;         # we assume dataType == f float
			}
		}
		my $expectedNum = $xc * $xw;
		$ret{mappedCount} += $expectedNum;
		next
			if not $expectedNum;    # sometimes we have a weight of zero...

		### get strand from the flag
		my $strand = $flag & 0x10;
		if ($flag - $strand)
		{
			die "Some unknown flag in SAM file exists, exitting\n";
		}
		$strand = $strand ? "-" : "+";

		### get tss position
		my $tssPosition;
		if ($cigar =~ /^(\d+)M$/)
		{
			$tssPosition = $strand eq "+" ? $leftpos : $leftpos - $1 + 1;
		}
		else
		{
			die "Sequence $seq in SAM file is not a perfect match. Unsupported\n";
		}

		### save it in the hash
		my $id = join("_", $chr, $strand, $tssPosition);
		if (not exists $hRef->{$id})
		{
			$hRef->{$id}->{count} = 0;
			$hRef->{$id}->{chr} = $chr;
			$hRef->{$id}->{str} = $strand;
			$hRef->{$id}->{pos} = $tssPosition;
		}
		$hRef->{$id}->{count} += $expectedNum;
	}
	close F;	
	return \%ret;
}


sub corr
{
  my ($x, $y) = @_;
  my $N             = scalar @{$x};
  my $sum_sq_x      = 0;
  my $sum_sq_y      = 0;
  my $sum_coproduct = 0;
  my $mean_x        = $$x[0];
  my $mean_y        = $$y[0];
  for (my $i = 1 ; $i < $N ; $i++)
  {
    my $sweep   = ($i) / ($i + 1);
    my $delta_x = $$x[$i] - $mean_x;
    my $delta_y = $$y[$i] - $mean_y;
    $sum_sq_x      += $delta_x * $delta_x * $sweep;
    $sum_sq_y      += $delta_y * $delta_y * $sweep;
    $sum_coproduct += $delta_x * $delta_y * $sweep;
    $mean_x        += $delta_x / ($i + 1);
    $mean_y        += $delta_y / ($i + 1);
  }
  my $pop_sd_x    = sqrt($sum_sq_x / $N);
  my $pop_sd_y    = sqrt($sum_sq_y / $N);
  my $cov_x_y     = $sum_coproduct / $N;
  my $correlation = $cov_x_y / ($pop_sd_x * $pop_sd_y);
  return $correlation;
}

sub cumulative_rev_unnorm
{
	# unnormalized reverse cumulative of the counts in the array
	# counts should be positive numbers.
	my ($tab) = @_;
	my $ret = ""; # returned text
	my $n = scalar @$tab;
	my @tab = sort {$a <=> $b} @$tab;
	my $lastval = -3.52923756295872957883;#nothing magic. just a negative value.
	my $lastplateau = $n;
	for(my $i=0; $i<$n; $i++)
	{
		if($tab[$i] != $lastval)
		{
			if($i)
			{
				$ret .= $tab[$i] . "\t" . $lastplateau . "\n";
			}
		  $ret .= $tab[$i] . "\t" . ($n-$i) . "\n";
		  $lastplateau = $n-$i;
		  $lastval = $tab[$i];
		}
	}
	return $ret;
}

sub readFit
{
	my($f) = @_;
	open(my $fh, $f) or die $!;
	$_ = <$fh>;
	close $fh;
	return split(/\t/);
}

sub cmpKawajiKey
{
	my $tmp = shift @_;
	my $tmp2 = shift @_;
	($tmp) =~ /^chr(\S+)\:\d+\.\.(\d+)\,([\+\-])$/ or die "$. unparsable $tmp2 $tmp\n";
	my $c1 = $1;
	my $s1 = $3;
	my $p1 = $2;
	($tmp2) =~ /^chr(\S+)\:\d+\.\.(\d+)\,([\+\-])$/ or die "$. unparsable $tmp2 $tmp\n";
	my $c2 = $1;
	my $s2 = $3;
	my $p2 = $2;
	return (($c1 cmp $c2) or ($s1 cmp $s2) or ($p1 <=> $p2));
}

1;
