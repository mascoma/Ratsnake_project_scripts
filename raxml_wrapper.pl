#!/usr/bin/perl -w

=head1 NAME

raxml_wrapper.pl

=head1 SYNOPSIS
This Perl script batches jobs to RAxML.

  raxml_wrapper.pl -- this program will batch jobs to raxml

Options:

 --help        Show brief help and exit
 --indir       Input directory with Phylip alignments to be analyzed
 --config      Contains the RAxML binary name and the common run parameters
 --outdir      Output directory

Your flavor of RAxML must be in your path

The config must specify parameters for the raxml batch run, such as

RAXML=raxml
PARAMETERS=-T 8 -m PROTGAMMAJTTF -f d

Usage example:

perl raxml_wrapper.pl -c config.raxml -i input_dir/ -o output_dir/

=head1 DESCRIPTION

Run RAxML on a set of alignments

Usage example:

perl raxml_wrapper.pl -c config.raxml -i input_dir/ -o output_dir/

=head1 AUTHOR

Apurva Narechania
anarechania *a|t* amnh.org

=head1 COPYRIGHT

Copyright (c) 2008-2009 American Museum of Natural History

This library is free software;  
you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

# ----------------------------------------------------

#####
#####     

#####SETUP#####

use strict;
use Getopt::Long;
use Pod::Usage;

my ($help, $indir, $config, $outdir);
GetOptions(
    'h|help'          => \$help,
    'i|indir=s'       => \$indir,
    'o|outdir=s'      => \$outdir,
    'c|config=s'      => \$config,
    ) or pod2usage;

pod2usage if $help;

for my $option ($indir, $outdir, $config){
    (warn ("Missing a required option\n") and pod2usage)
	unless ($option);
}

`mkdir -p $outdir/results`;
`mkdir -p $outdir/stderr`;

#####MAIN#####

# get the config
my $conf = parse_config ($config);

# do the raxml run for every aln
opendir (I, "$indir");
my @alns = sort (readdir (I));
shift @alns;
shift @alns;
closedir (I);

foreach my $aln (@alns){
    warn "RAXML $aln\n";
    raxml_run ($aln, $indir, $conf, $outdir);
}


#####SUBS#####

sub raxml_run{
    my $aln  = shift;
    my $in   = shift;
    my $conf = shift;
    my $out  = shift;

    # raxml
    my $raxml = "$conf->{RAXML}";
    ($raxml .= " $conf->{'PARAMETERS'}") if ($conf->{'PARAMETERS'});
    $raxml .= " -s $in/$aln";
    $raxml .= " -n $aln";
    $raxml .= " -w $out/results";
    `$raxml &> $out/stderr/$aln.stderr`;
}

sub parse_config{
    my $file = shift;

    my %config;
    open (F, "$file");
    while (my $line = <F>){
        chomp $line;
        my ($key, $value) = split (/\=/, $line, 2);
        $config{$key} = $value;

    }
    close (F);

    return (\%config);
}
