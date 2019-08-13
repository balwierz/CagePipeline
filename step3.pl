#!/import/bc2/soft/app/perl/5.10.1/Linux/bin/perl -w
use strict;
use lib '.';
use promoteromeLib;
use YAML::Any;
$| = 1;

# This file was written for clustering of data from Fantom 5.
# This script normalizes the tag counts using cumulative fits from step2.pl


######################## config
my $conf_file = shift @ARGV or die "specify a config file!\n";
my $theOnlyFile = shift @ARGV;  # second parameter if not empty is just one file to run on: for parallelizing
my $config = YAML::Any::LoadFile($conf_file);

my $dir = $config->{dir};
my $tmpDir = $config->{tmpDir};
my $cumDir = $config->{cumDir};
my $mappedDir = $config->{mappedDir};
my $normDir = $config->{normDir};
my $normDirSplit = $config->{normDirSplit};
my $fitDir = $config->{fitDir};

##################### CODE

mkdir $tmpDir.$normDir;
mkdir $tmpDir.$normDirSplit;

if(defined $theOnlyFile)
{
	process($theOnlyFile)
}
else
{
	opendir(my $D, $dir) or die;
	while(my $f = readdir $D)
	{
		if(-f $dir.$f || -l $dir.$f)
		{
			process($f);
		}
	}
	close $D;
}


sub process
{
	my($f) = @_;
	my ($sampleDesc, $sampleID) = split(/\./, $f);
	#my $outF = $tmpDir.$normDir.$sampleID;
	
	#if(-f $dir.$f && ! -e $outF)
	{
		print STDERR $sampleID, "\n";
		my %h;
		my($f2, $a, $b, $beta, $lambda) = readFit($tmpDir.$fitDir.$sampleID);
		print STDERR "Reading RAW expression\n";
		processBed(\%h, $dir.$f);
		print STDERR "Normalizing\n";
		normalize(\%h, $beta, $lambda); #adds {norm} keys to the hash
		print STDERR "Writing norm and raw\n";
		writeNormRaw(\%h, $tmpDir.$normDir.$sampleID, $tmpDir.$normDirSplit, $sampleID);
		system('gzip', $tmpDir.$normDir.$sampleID);
	}
}
