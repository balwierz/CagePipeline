#!/mnt/biggles/opt/perl/bin/perl -w
use strict;
use YAML::Any;

my $conf_file = shift @ARGV or die "Privide a conf file!\n";
my $doOnlyBadSamples = shift @ARGV || 0;
my $config = YAML::Any::LoadFile($conf_file);
my $dir = $config->{dir};
my $tmpDir = $config->{tmpDir};
my $cumDir = $config->{cumDir};

my %badSamples;
open(A, $config->{fitRegionF}) or die;
while(<A>)
{
	chomp;
	s/#.*$//;
	next if not $_;
	my ($sample) = split(/\s+/);
	$badSamples{$sample} = 1;
}
close(A);

opendir(D, $dir) or die;
while(my $file = readdir D)
{
	if( -f $dir.$file && $file =~ /\S+\.(CNhs\d+)\.\S+\.ctss\./ && (!$doOnlyBadSamples || defined $badSamples{$1}))
	{
		my ($sampleDesc, $sampleID) = split(/\./, $file);
		if(! -e $tmpDir.$cumDir.$sampleID)
		{
			print STDERR "Processing: ", $file,"\n";
			open(F, ">submit.pl") or die;
			print F
			'#!/mnt/biggles/opt/perl/bin/perl'."\n".
			"#SBATCH -e slurm-%J.err\n".
			"#SBATCH -o slurm-%J.out\n".
			"#SBATCH --mem 5G\n".
			"system('/usr/bin/perl', 'step1.pl', '$conf_file', '$file');"."\n".
			"system('/usr/bin/perl', 'step2.pl', '$conf_file', '$file');"."\n".
			"system('/usr/bin/perl', 'step3.pl', '$conf_file', '$file');"."\n";
			close F;
			system('sbatch', 'submit.pl');
		}
		else
		{
			print STDERR "Already processed: $sampleID\n";
		}
	}
	else
	{
		print "Ignoring $file\n"
	}
}	
closedir D;



