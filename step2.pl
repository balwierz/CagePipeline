#!/usr/bin/perl -w
use strict;
use lib '.';
use promoteromeLib;
use YAML::Any;
$| = 1;

# This file was written for clustering of data from Fantom 5.
# This script fits the power-laws to the cumulatives produced in step1.pl



######################## config
my $conf_file = shift @ARGV or die "specify a config file!\n";
my $theOnlyFile = shift @ARGV;  # second parameter if not empty is just one file to run on: for parallelizing
my $config = YAML::Any::LoadFile($conf_file);

my $dir = $config->{dir};
my $tmpDir = $config->{tmpDir};
my $cumDir = $config->{cumDir};
my $mappedDir = $config->{mappedDir};
my $normDir = $config->{normDir};
my $fitDir = $config->{fitDir};
my $fitRegionF	= $config->{fitRegionF};

##################### CODE
mkdir $tmpDir.$fitDir;

my $fitRegion = readRegions($fitRegionF);

if(defined $theOnlyFile)
{
	my($desc, $sample) = split(/\./, $theOnlyFile);
	process($sample);
}
else
{
	opendir(my $D, $tmpDir.$cumDir) or die;
	while(my $f = readdir $D)
	{
		if(-f $tmpDir.$cumDir.$f) #&& $f =~ /\.revcum$/
		{
			if(! -e $tmpDir.$fitDir.$f )
			{
				process($f);
			}
		}
	}
	close $D;
}

sub process
{
	my ($f) = @_;
	my ($sampleID) = split(/\./, $f);
	print STDERR $sampleID, "\n";
	my($startX, $startY);
	if(exists $fitRegion->{$sampleID})
	{
	 	($startX, $startY) = @{$fitRegion->{$sampleID}};
	}
	else
	{
		$startX = 5;
		$startY = 1000;
	}
	my $d = readCumulative($tmpDir.$cumDir.$f);
	my ($a, $b, $beta, $lambda, $txt) = regressPowerLaw(\$d, $startX, $startY, "1.15");
	print $txt, "\n";
	gnuplotRevcumLine($tmpDir.$cumDir.$f, $tmpDir.$cumDir.$f.'.revcum.png', $a, $b);
	
	open(my $outH, ">", $tmpDir.$fitDir.$f) or die;
	print $outH join("\t", $f, $a, $b, $beta, $lambda), "\n";
	close $outH;
}

sub readRegions
{
	my($f) = @_;
	open (my $fH, $f) or die;
	my %ret;
	while(<$fH>)
	{
		chomp;
		next if /^\s*#/;
		my($lib, $startX, $endX, $startY) = split(/\t/);
		my @tmp = ($startX, $startY);
		$ret{$lib} = \@tmp;
	}
	close $fH;
	return \%ret;
}
