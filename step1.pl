#!/usr/bin/perl -w
#$ -S /import/bc2/soft/app/perl/5.10.1/Linux/bin/perl
#$ -P project_zavolan
#$ -q fs_very_long_hm
#$ -l mem_total=10000M
#$ -o step1.out
#$ -e step1.err
#$ -j n
#$ -N step1
#$ -cwd

# This file was written for clustering of data from Fantom 5.
# This script produces cumulatives in the tmpDir directory.

# TODO: it should also make a file with sample id => description, depth, and so on.


# necessery to exec the RIGHT perl in the system call
# $ENV{PATH} = '/import/bc2/soft/app/perl/5.10.1/Linux/bin:' . $ENV{PATH};

use strict;
use lib '.';
use promoteromeLib;
use YAML::Any;
$| = 1;


######################## config
my $conf_file = shift @ARGV or die "specify a config file!\n";
my $theOnlyFile = shift @ARGV;  # second parameter if not empty is just one file to run on: for parallelizing
my $config = YAML::Any::LoadFile($conf_file);

my $dir = $config->{dir};
my $tmpDir = $config->{tmpDir};
my $cumDir = $config->{cumDir};
my $mappedDir = $config->{mappedDir};
my $normDir = $config->{normDir};

##################### CODE

# make sure the dires are created
mkdir $tmpDir.$cumDir;
mkdir $tmpDir.$mappedDir;
mkdir $tmpDir.$normDir;

if(defined  $theOnlyFile)
{
	processFile($dir, $theOnlyFile);
}
else
{
	opendir(my $D, $dir) or die;
	while(my $f = readdir $D)
	{
		if(-f $dir.$f || -l $dir.$f)
		{
			processFile($dir, $f);
		}
	}
	close $D;
}

sub processFile
{
	my ($dir, $f) = @_;
	my ($sampleDesc, $sampleID) = split(/\./, $f);
	my $outF = $tmpDir.$cumDir.$sampleID;
	if(! -e $outF)
	{
		print STDERR $f, "\n";
		my %h;
		my $sampleInfo = processBed(\%h, $dir.$f);
		my $cumTxt = cumulative_rev_unnorm(h2counts(\%h));
		open(A, '>', $outF ) or die "cannot write to $outF\n";
		print A $cumTxt;
		close A;

		open(B, '>', $tmpDir.$mappedDir.$sampleID) or die;  #.".mappedCount"
		print B $sampleInfo->{mappedCount}, "\n";
		close B;
	}
	else
	{
		print STDERR "Sample $sampleID is already processed. Quitting.\n";
	}
}

