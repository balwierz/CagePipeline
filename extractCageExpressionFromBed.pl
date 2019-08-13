#!/usr/bin/perl -w
use strict;
use PerlIO::gzip;
use YAML::Any;
$|++;

######################## config
my $conf_file = shift @ARGV or die "specify a config file!\n";
my $config = YAML::Any::LoadFile($conf_file);

# my $dir = $config->{dir}; # unused!
my $tmpDir = $config->{tmpDir};
my $cumDir = $config->{cumDir};
my $mappedDir = $config->{mappedDir};
my $normDir = $config->{normDir};
my $normDirSplit = $config->{normDirSplit};

##################### CODE

# a hash which keps the positions from Kawaji's file
# {position} --> clusterID
# the chromosome and strand are fixed since we run this script
# once for each combination
my %pos2clus;


my $f = shift @ARGV or die;
my $outF = shift @ARGV or die;
my $sampleF = shift @ARGV or die;
my $taskID = shift @ARGV or die;
my $ignoreSingletons = shift @ARGV or 0;
$taskID --; #because they start at 1 instead of 0

# we need to translate taskID to a chrstr combination
open(K, "ls $tmpDir"."$normDirSplit |") or die;
my @allChrStr = <K>;
close K;
@allChrStr = sort @allChrStr; #just to make sure
my $thisChrStr = $allChrStr[$taskID];
chomp $thisChrStr;



# a global variable keeping all expression levels
#my $expr; #key {str}{pos}{library} value is normalized expression
#my $exprRaw;

my @list; #of samples
open(my $fH, $sampleF) or die;
while(<$fH>)
{
	chomp;
	next if /^#/;
	push @list, $_;
}
close $fH;

my $nSamples = scalar @list;
my %sample2index;
my @index2sample;
for(my $i=0; $i<@list; $i++)
{
	my @v = split(/\t/, $list[$i]);
	my $sName = "";
	if(@v == 1)
		{$sName = $v[0];}
	elsif(@v == 2)
		{$sName = $v[1];}
	else
	{
		print "Input sample file in a wrong format ", $list[$i]; die;
	}
	$list[$i] = $sName;
	$sample2index{$sName} = $i;
	$index2sample[$i] = $sName;
}

#print join "\t", @list, "\n";

# read the BED file
if($f =~ /\.gz$/)
{
	open($fH, "<:gzip", $f) or die;
}
else
{
	open($fH, $f) or die;
}
my @allClus;
while(<$fH>)
{
	chomp;
	next if $_ =~ /^\s*$/;
	my($chr, $beg, $end, $clus, $foo, $str) = split(/\t/);
	
	# ignore the clusters from other chr.str:
	#print $chr, $str, " ", $thisChrStr, "\n";
	next if $chr.$str ne $thisChrStr;
	
	
	push @allClus, $clus;
	# update the hash:
	for(my $i=$beg+1; $i<=$end; $i++)
	{
		$pos2clus{$i} = $clus;
	}
}
close $fH;

# global varialbles holding expression
my $exprNorm;
my $exprRaw;
my $sumNorm;
my $maxNorm;
my $sumRaw;
my $maxRaw;
# now we go to our normalized data and update the expression:
my $d = "$tmpDir/$normDirSplit/$thisChrStr/";
opendir(my $eD, $d) or die;
my $i = 0;
#while(my $f = readdir $eD)
#{
foreach my $f (@list)   #if($f !~ /^\./ && exists $sample2index{$f}) # don't read other samples
{
	open(my $fH, $d.$f) or die "Missing data: $d $f\n";
	print STDERR "Reading $thisChrStr $f\n";
	while(<$fH>)
	{
		chomp;
		my($chr, $str, $pos, $raw, $norm) = split(/\t/);
		if(exists $pos2clus{$pos})
		{
			my $clus = $pos2clus{$pos};
			if(! $ignoreSingletons || $raw > 1)
			{
				$exprNorm->{$clus}{$f} += $norm;
				$exprRaw ->{$clus}{$f} += $raw;
				$sumNorm->{$clus} += $norm;
				$sumRaw ->{$clus} += $raw;
				if((! defined $maxNorm->{$clus}) || $norm > $maxNorm->{$clus})
				{
					$maxNorm->{$clus} = $norm;
				}
				if((! defined $maxRaw->{$clus}) || $raw > $maxRaw->{$clus})
				{
					$maxRaw->{$clus} = $raw;
				}
			}
		}
	}
	close $fH;
	#$allNormalized{$f} = 1;
	$i++;
}
	# print the size 
	#if(!($i % 10))
	#{
		#print STDERR "Size: ".`ps -p $$ h -o rss`;
	#}
#}
#closedir $eD;

# optionally: sort the clusters
#@allClus = sort {cmpID($a, $b)} @allClus;


# and finally print the result:

open(O, ">$outF.normalized") or die;
open(O2, ">$outF.raw") or die;
foreach my $clus (@allClus)
{
        print O $clus;
        #print O "\t", exists $maxNorm->{$clus} ? $maxNorm->{$clus} : 0, "\t", exists $sumNorm->{$clus} ? $sumNorm->{$clus} : 0;
        print O2 $clus;
        #print O2 "\t", exists $maxRaw->{$clus} ? $maxRaw->{$clus} : 0, "\t", exists $sumRaw->{$clus} ? $sumRaw->{$clus} : 0;
        foreach my $sample (@list)
        {
                if(exists $exprNorm->{$clus}->{$sample})
                {
                        print O "\t", $exprNorm->{$clus}->{$sample};
                        print O2 "\t", $exprRaw->{$clus}->{$sample};
                }
                else
                {
                        print O "\t", "0";
                        print O2 "\t", "0";
                }
        }
        print O "\n";
        print O2 "\n";
}
close O;
close O2;


sub cmpID
{
	my($id1, $id2) = @_;
	$id1 =~ /(chr\S+):(\d+)\.\.(\d+),([\-\+])/ or die;
	my($beg1, $end1) = ($2, $3);
	$id2 =~ /(chr\S+):(\d+)\.\.(\d+),([\-\+])/ or die;
	my($beg2, $end2) = ($2, $3);	
	return $beg1 <=> $beg2; # they cannot overlap (I hope!)
}
	

