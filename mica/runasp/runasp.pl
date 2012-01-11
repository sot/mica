#!/usr/bin/env perl
# a quick script to go through job JSON files and attempt to run
# pipes.  Intended to become a lightweight "server" if needed.
# No error checking at present...

use strict;
use warnings;
use JSON;
use IO::File;
use Data::Dumper;
use App::Env;
use Shell::GetEnv;
use Carp;

my $job_dir = $ARGV[0];

my $ascds = App::Env->new('ASCDS');
my $ocat_auth_file = "/proj/sot/ska/data/aspect_authorization/set_ascds_ocat_vars.csh";
my $ocat = Shell::GetEnv->new( 'tcsh', "source $ocat_auth_file");
for my $my_ocat ('ASCDS_OCAT_UNAME', 'ASCDS_OCAT_SERVER', 'ASCDS_OCAT_PWORD'){
    $ascds->setenv( $my_ocat, $ocat->envs()->{$my_ocat});
}



sub get_job{
    my $jobfile = shift;
    open(my $fh,'<',"$jobfile") or carp $!;
    my $jobtext = do { local $/;<$fh>; };
    close($fh);
    my $job = decode_json($jobtext);
    return $job;
}

use File::chdir;

sub run_pipe{
    my ($job, $cmd_str) = @_;
    local $CWD = $job->{'dir'};
    print `pwd`;
    my $jstr = "flt_run_pipe -i $job->{indir} -o $job->{outdir} "
             . "-r $job->{root} -t $job->{pipe_ped} "
             . "-a \"INTERVAL_START\"=$job->{istart} "
             . "-a \"INTERVAL_STOP\"=$job->{istop} "
             . "-a obiroot=$job->{obiroot} "
             . "-a revision=1 " ;
    if (defined $cmd_str){
        $jstr .= $cmd_str;
    }
    print $jstr, "\n";
    $ascds->qexec("$jstr");
}

#use File::Temp qw(tempfile);
use File::Copy;
use File::Spec::Functions qw(catdir);
sub cut_stars{
    my $job = shift;
    use Data::Dumper;
    print Dumper $job;
    my $outdir = catdir($job->{dir}, $job->{outdir});
    print "$outdir \n";
    my @starfiles = glob("${outdir}/*stars.txt");
    if (not scalar(@starfiles)){
        croak("No star.txt files found matching ${outdir}/*stars.txt");
    }
    copy("$starfiles[0]", "$starfiles[0].orig") 
        or die "Backup Copy failed: $!";    
    open(my $star_fh, $starfiles[0]);
    my @starlines = <$star_fh>;
    my @newlist = @starlines;
    for my $skip_slot (@{$job->{skip_slot}}){
        @newlist = grep(!/\s*${skip_slot}\s+1/, @newlist);
    }
    open(my $star_wfh, '>', $starfiles[0]);
    print $star_wfh @newlist;
    close($star_wfh);
}

sub run_job{
    my $job = shift;
#    print Dumper $job;
    if (scalar($job->{'skip_slot'})){
        run_pipe($job, '-S check_star_data');
        cut_stars($job);
        run_pipe($job, '-s check_star_data');
    }
    else{
        run_pipe($job);
    }
}

my @jobs = glob("${job_dir}/*job");

#print Dumper \@jobs;

JOB:
for my $jobfile (@jobs){
    if ((-e "${jobfile}.done") or (-e "${jobfile}.lock")){
        next JOB;
    }
    else{
        system("touch ${jobfile}.lock");
        my $job = get_job($jobfile);
        run_job($job);
        system("touch ${jobfile}.done");
    }
    
}
