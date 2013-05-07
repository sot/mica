#!/usr/bin/env perl

use strict;
use warnings;
use JSON;
use Data::Dumper;

use lib '/proj/sot/ska/jeanproj/git/mica/mica/starcheck/lib';
use StarcheckParser;


my $starcheck_file = $ARGV[0];
my $starcheck = StarcheckParser->new_starcheck($starcheck_file);
my %header = $starcheck->get_header_lines();
my $obs_data= $starcheck->get_all();
my %out = ('header' => \%header,
           'obsdata' => $obs_data);

my $json = to_json(\%out, {pretty => 1});
print $json;


