#!/import/bc2/soft/app/perl/5.10.1/Linux/bin/perl -w
use strict;
use promoteromeLib;
use List::Util qw(first max min sum);
use YAML::Any;
$| = 1;

######################### CONFIG ######################

my @dendriticReplicates = qw
(
	CNhs10855
	CNhs11062
);

# each of these has 13M tags, 8.6M TSS together, quite a lot of outliers in one direction
my @mouseGoodReplicates = qw
(
	CNhs10623
	CNhs11072
);
# really bad sigma 0.20
my @mouse2 = qw
(
	CNhs11208
	CNhs11048
);
#3139-67G5	3139-67G5	Mouse E17.5 Internal control	Mus musculus	standard whole cell	standard whole cell	3139-67G5.rna	3139-67G5.rna	total RNA	comment:Internal control RNA for FANTOM5; rna_sample_type:total RNA; extract_method:standard whole cell; RIN:8.5; A260A280:2; A260A230:2.1	OP-HELISCOPE-CAGE-v3.12	747	OP-HELISCOPE-CAGE-v3.12	CNhs11208	CNhs11208	OP-HELISCOPE-sequencing-v1.0	OP-HELISCOPE-sequencing-v1.0	NA	NA	fc2.ch13	/sequencedata/heliscope/srf/2010_10_21_R2_30Q_1100FOV	/sequencedata/heliscope/srf/2010_10_21_R2_30Q_1100FOV
#3139-67G5	3139-67G5	Mouse E17.5 Internal control	Mus musculus	standard whole cell	standard whole cell	3139-67G5.rna	3139-67G5.rna	total RNA	comment:Internal control RNA for FANTOM5; rna_sample_type:total RNA; extract_method:standard whole cell; RIN:8.5; A260A280:2; A260A230:2.1	OP-HELISCOPE-CAGE-v3.12	722	OP-HELISCOPE-CAGE-v3.12	CNhs11048	CNhs11048	OP-HELISCOPE-sequencing-v1.0	OP-HELISCOPE-sequencing-v1.0	NA	NA	fc1.ch13	/sequencedata/heliscope/srf/2010_10_08_R4_30Q_1100FOV	/sequencedata/heliscope/srf/2010_10_08_R4_30Q_1100FOV

my @mouse3 = qw
(
	CNhs12230
	CNhs11262
);

my @kawaji = qw
(
	CNhs10479
	CNhs10749
	CNhs11778
);

my $conf_file = shift @ARGV or die "specify a config file!\n";
my $config = YAML::Any::LoadFile($conf_file);
my $dir = $config->{dir};
my $tmpDir = $config->{tmpDir};
my $cumDir = $config->{cumDir};
my $mappedDir = $config->{mappedDir};
my $normDir = $config->{normDir};

my $promselector = \@mouse3; #\@dendriticReplicates;
#range of 20000 tags
my $W = 1.0 / log(20000);
my $rawCutoff = 5; # tags
my $rawCutoffUpper = 99999999999999;
my $doScatter = shift @ARGV;
my $differExprFract = 0.000; # put it to something like 0.005 if you believe these are NOT full replicates

####################### CODE ##########################
print STDERR "make scatter file = " . ((defined $doScatter && $doScatter) ? "T" : "F") . "\n";
my $twopi = 3.14159265 * 2;
my $expr;

# read totals:
my @tot = ();
foreach my $i (0..1)
{
	my $s = $promselector->[$i];
	# read totals:
	open(F, $tmpDir.$mappedDir.$s.".mappedCount") or die $!;
	$_ = <F>; chomp;
	push @tot, $_;
	close(F);

	#we read both the files with absolute numbers and with tpm
	print STDERR "Reading expression $s\n";
	open(my $fh, $tmpDir.$normDir.$s.".nr") or die $!;
	while(<$fh>)
	{
		chomp;
		my($chr, $str, $pos, $raw, $norm) = split(/\t/);
		my $key = join(" ", $chr, $str, $pos);
		$expr->{raw}->{$key}->[$i] = $raw;
		$expr->{norm}->{$key}->[$i] = $norm;
	}
	close $fh;
}

my $totratio = $tot[0] / $tot[1];

my @sd_norm = (); # filled by the loop below
my @sd_tpm  = ();
my @sam = ();
if( defined $doScatter && $doScatter)
{
	open(SCATTER, ">", "scatter") or die;
	open(VAR, ">", "variances") or die;
}

print STDERR "Selecting pairs from ". (scalar keys %{$expr->{raw}}) . " TSSs\n";

# this is for speedup
my @IDs = keys %{$expr->{raw}};
my @bins;
my @bins2;
if( defined $doScatter && $doScatter)
{
	print STDERR "Sorting keys...";
	#@IDs = sort{cmpKey($a, $b)} keys %{$expr->{raw}};
	print STDERR "  [done]\n";
	#for(my $i=0; $i<2000; $i++)
	#{
#		$bins[$i] = 0;
#	}
}
else
{
	#@IDs = keys %{$expr->{raw}};
}

foreach my $key (@IDs)
{
	if( defined $doScatter && $doScatter)
	{
		my $r0 = defined $expr->{raw}->{$key}->[0] ? $expr->{raw}->{$key}->[0] : 0 ;
		my $r1 = defined $expr->{raw}->{$key}->[1] ? $expr->{raw}->{$key}->[1] : 0 ;
		
		my $n0 = defined $expr->{norm}->{$key}->[0] ? $expr->{norm}->{$key}->[0] : 0 ;
		my $n1 = defined $expr->{norm}->{$key}->[1] ? $expr->{norm}->{$key}->[1] : 0 ;
		print SCATTER $r0, "\t", $r1 , "\t", $key, "\n";
		# variance binned
		#if($r0 > 0  && $r1 > 0)
		{
			my $mean = log( ($n1 + $n0) / 2 );
			#(log($n0+0.01) + log($n1+0.01)) / 2 ;
			my $bin = int($mean * 10 + 0.5) + 300;
			#print VAR $bin, "\n";
			push @{$bins[$bin]}, ((log($n0+0.1) - log($n1+0.1))/2)**2;
			push @{$bins2[$bin]},((log($n0+0.1) - log($n1+0.1))/2);
		}
	}
	if ( defined $expr->{raw}->{$key}->[0] && $expr->{raw}->{$key}->[0] >= $rawCutoff &&
		 defined $expr->{raw}->{$key}->[1] && $expr->{raw}->{$key}->[1] >= $rawCutoff &&
		 $expr->{raw}->{$key}->[0] <= $rawCutoffUpper &&
		 $expr->{raw}->{$key}->[1] <= $rawCutoffUpper)
	{
		#squared log deviation
		my $d = log($expr->{norm}->{$key}->[1] / $expr->{norm}->{$key}->[0]);
		$d *= $d;
		push @sd_norm,  $d;
		
		$d = log(($expr->{raw}->{$key}->[1] / $expr->{raw}->{$key}->[0]) * $totratio);
		$d *= $d;
		push @sd_tpm,	$d;
		
		# 1/n + 1/m
		my $samp = (1.0 / $expr->{raw}->{$key}->[0]) + (1.0 / $expr->{raw}->{$key}->[1]);
		push @sam, $samp;
	}
}



if( defined $doScatter && $doScatter)
{
	for(my $i=0; $i<@bins; $i++)
	{
		if( defined $bins[$i])
		{
			print VAR exp(($i-300)/10), "\t",
			sum(@{$bins[$i]})/(scalar @{$bins[$i]}) - (sum(@{$bins2[$i]})/(scalar @{$bins2[$i]}))**2,
			"\t", (scalar @{$bins[$i]}), "\n";
		}
		else
		{
			print VAR exp(($i-300)/10), "\t", 0, "\t", 0, "\n";
		}
	}
	close SCATTER;
	close VAR;
}
#die if defined $doScatter && $doScatter;
my $num = @sam;
print STDERR "Selected $num points\n";
foreach my $sdref ((\@sd_norm, \@sd_tpm))
{
	#print STDERR "Fitting\n";
	my $sqmax = 10.0;
	my $sqmin = 0.0000;
	my $sq = undef; # two sigma squared
	my $pi = undef;
	while (abs(2.0 * ($sqmax - $sqmin) / ($sqmax + $sqmin)) > 0.005)
	{
		$sq = 0.5 * ($sqmax + $sqmin);
		#print "Current sq $sq\n";

		#calculate
		my $dev = 1.0;
		$pi  = $differExprFract;
		my $dsig = undef;
		while ($dev > 0.005)
		{

			#calculate new pi
			my $tot  = 0;
			$dsig = 0;
			my $totLogLikelihood = 0;
			for (my $i = 0 ; $i < $num ; ++$i)
			{
				my $Ldif = $pi * $W;
				my $v = 1.0 / ($sq + $sam[$i]);
				my $Lsame =
					(1.0 - $pi) * exp(-0.5 * $$sdref[$i] * $v) * sqrt($v * 1/$twopi)
					;    #real number is 1/sqrt(2pi) 0.39894228
					# shouldn't it be 1/(2pi) ??
				$tot += $Ldif / ($Lsame + $Ldif);
				$dsig += 0.5 * $Lsame * $v * ($$sdref[$i] * $v - 1) / ($Lsame + $Ldif);
				# the numerator looks good. This is d(lsame)/d(sq)
				# the denominator comes from d(log(L))/d(sq) = d(lsame)/d(sq) / (L)
				$totLogLikelihood += log($Lsame + $Ldif);
			}
			my $oldpi = $pi;
			$pi    = $tot / $num;
			if ($pi < 0.0000001)
			{
				$pi = 0;
			}
			if ($oldpi > 0 || $pi > 0)
			{
				$dev = 0.5 * abs($pi - $oldpi) / ($pi + $oldpi);
			}
			else
			{
				$dev = 0;
			}
			print STDERR "2sigma^2 $sq totLogLikelihood $totLogLikelihood dsig $dsig pi $pi dev $dev\n";
		}
		if ($dsig < 0)
		{
			$sqmax = $sq;
		}
		else
		{
			$sqmin = $sq;
		}
	}
	print "sigma=",  txt(sqrt(0.5*$sq)), " sigma^2=".txt($sq/2) ." pi=" . txt($pi). "\n";
}

sub cmpKey
{
	my($a, $b) = @_;
	my($c1, $s1, $p1) = split(/ /, $a);
	my($c2, $s2, $p2) = split(/ /, $b);
	return (($c1 cmp $c2) || ($s1 cmp $s2) || ($p1 <=> $p2));
}

#now calculate for large and small sig
#$sq   = 0.0001;
#$tot  = 0;
#$dsig = 0;
#for ($i = 0 ; $i < $num ; ++$i)
#{
    #$Ldif  = $pi * $W;
    #$v     = 1.0 / ($sq + $sam[$i]);
    #$Lsame = (1.0 - $pi) * exp(-0.5 * $sd[$i] * $v) * sqrt($v * 0.39894228);
    #$tot += $Ldif / ($Lsame + $Ldif);
    #$dsig += 0.5 * $Lsame * $v * ($sd[$i] * $v - 1) / ($Lsame + $Ldif);
#}
#$oldpi = $pi;
#$pi    = $tot / $num;
#print STDERR "sq 0.0001 new pi is $pi dsig is $dsig\n";

#$sq   = 4.0;
#$tot  = 0;
#$dsig = 0;
#for ($i = 0 ; $i < $num ; ++$i)
#{
    #$Ldif  = $pi * $W;
    #$v     = 1.0 / ($sq + $sam[$i]);
    #$Lsame = (1.0 - $pi) * exp(-0.5 * $sd[$i] * $v) * sqrt($v * 0.39894228);
    #$tot += $Ldif / ($Lsame + $Ldif);
    #$dsig += 0.5 * $Lsame * $v * ($sd[$i] * $v - 1) / ($Lsame + $Ldif);
#}
#$oldpi = $pi;
#$pi    = $tot / $num;
#print STDERR "sq 4.0 new pi is $pi dsig is $dsig\n";

__END__
for allcage: first normalized, than tpm
2008-02-17
2*sigma^2 = 0.0969328588867188, sigma = 0.220150924239167
2*sigma^2 = 0.0974211352539062, sigma = 0.220704706852738
2*sigma^2 = 0.233406103515625, sigma = 0.34161828369953
2*sigma^2 = 0.225593681640625, sigma = 0.335852409281685
2*sigma^2 = 0.116708051757813, sigma = 0.241565779610661
2*sigma^2 = 0.102059760742187, sigma = 0.225897942379061

Redone 16.02.2009
pi=0.00258953294277914, 2*sigma^2 = 0.09674072265625, sigma = 0.219932629066551
pi=0.002543236144098, 2*sigma^2 = 0.09735107421875, sigma = 0.220625331975671
pi=0.00139435589446938, 2*sigma^2 = 0.2325439453125, sigma = 0.340986763168675
pi=0.00129141045096554, 2*sigma^2 = 0.2252197265625, sigma = 0.33557393117054
pi=0.00221328502217111, 2*sigma^2 = 0.11688232421875, sigma = 0.241746069480716
pi=0.00211657555154217, 2*sigma^2 = 0.10223388671875, sigma = 0.226090564507622

allcage, 0h comparison between libraries (pi is not apropriate though)
pi=0.00492157101442631, 2*sigma^2 = 0.65185546875, sigma = 0.57090080957641
pi=0.00480391996551706, 2*sigma^2 = 0.67138671875, sigma = 0.579390506804349
pi=0.00671753280807771, 2*sigma^2 = 0.70556640625, sigma = 0.593955556523382
pi=0.00556543501941274, 2*sigma^2 = 0.71533203125, sigma = 0.598051850281395
pi=0.00392978940408358, 2*sigma^2 = 0.66650390625, sigma = 0.577279787559724
pi=0.00273307786609931, 2*sigma^2 = 0.66650390625, sigma = 0.577279787559724



for pma
2008-04-22 
2*sigma^2 = 0.274421318359375, sigma = 0.370419571809708
2*sigma^2 = 0.0981535498046875, sigma = 0.221532785163605
2*sigma^2 = 0.334967587890625, sigma = 0.409247839267738
2*sigma^2 = 0.225593681640625, sigma = 0.335852409281685
2*sigma^2 = 0.120614262695313, sigma = 0.245575103273227
2*sigma^2 = 0.102059760742187, sigma = 0.225897942379061
