#!/usr/bin/env perl

# DO NOT EDIT the /cluster/bin/scripts copy of this file --
# edit ~/kent/src/hg/utils/automation/doBlastzChainNet.pl instead.

# $Id: doBlastzChainNet.pl,v 1.33 2010/04/12 16:33:12 hiram Exp $

# to-do items:
# - lots of testing
# - better logging: right now it just passes stdout and stderr,
#   leaving redirection to a logfile up to the user
# - -swapBlastz, -loadBlastz
# - -tDb, -qDb
# - -tUnmasked, -qUnmasked
# - -axtBlastz
# - another Gill wish list item: save a lav header (involves run-blastz-ucsc)
# - 2bit / multi-sequence support when abridging?
# - reciprocal best?
# - hgLoadSeq of query instead of assuming there's a $qDb database?

use Getopt::Long;
use warnings;
use strict;
use FindBin qw($Bin);
use lib "$Bin";
use HgAutomate;
use HgRemoteScript;
use HgStepManager;
use File::Basename;

# Hardcoded paths/command sequences:
my $getFileServer = '/cluster/bin/scripts/fileServer';
my $blastzRunUcsc = "$Bin/blastz-run-ucsc";
my $partition = "$Bin/partitionSequence.pl";
my $clusterLocal = '/scratch/hg';
my $clusterSortaLocal = '/iscratch/i';
my @clusterNAS = ('/cluster/bluearc', '/san/sanvol1');
my $clusterNAS = join('/... or ', @clusterNAS) . '/...';
my @clusterNoNo = ('/cluster/home', '/projects');
my @fileServerNoNo = ('kkhome', 'kks00');
my @fileServerNoLogin = ('kkusr01', '10.1.1.3', '10.1.10.11',
			 'sanhead1', 'sanhead2', 'sanhead3', 'sanhead4',
			 'sanhead5', 'sanhead6', 'sanhead7', 'sanhead8');

# Option variable names, both common and peculiar to doBlastz:
use vars @HgAutomate::commonOptionVars;
use vars @HgStepManager::optionVars;
use vars qw/
    $opt_blastzOutRoot
    $opt_swap
    $opt_chainMinScore
    $opt_chainLinearGap
    $opt_tRepeats
    $opt_qRepeats
    $opt_readmeOnly
    $opt_ignoreSelf
    $opt_syntenicNet
    $opt_noDbNameCheck
    $opt_inclHap
    $opt_noLoadChainSplit
    $opt_loadChainSplit
    $opt_swapDir
    $opt_asmId
    $opt_skipDownload
    $opt_trackHub
    /;
###DIANA COMBINED ABOVE $OPT_TASMID AND $OPT_QASMID INTO $OPT_ASMID

# Specify the steps supported with -continue / -stop:
my $stepper = new HgStepManager(
    [ { name => 'partition',  func => \&doPartition },
      { name => 'blastz',     func => \&doBlastzClusterRun },
      { name => 'cat',        func => \&doCatRun },
      { name => 'chainRun',   func => \&doChainRun },
      { name => 'chainMerge', func => \&doChainMerge },
      { name => 'net',        func => \&netChains },
      { name => 'load',       func => \&loadUp },
      { name => 'download',   func => \&doDownloads },
      { name => 'cleanup',    func => \&cleanup },
      { name => 'syntenicNet',func => \&doSyntenicNet }
    ]
			       );

# Option defaults:
### DINA CHANGED 'KU' TO LOCALHOST BELOW AND HGWDEV TO LOCAL HOST AS WELL
my $bigClusterHub = 'localhost';
# my $smallClusterHub = 'encodek';
my $smallClusterHub = 'localhost';
my $dbHost = 'localhost';
my $workhorse = 'localhost';
my $defaultChainLinearGap = "loose";
my $defaultChainMinScore = "1000";	# from axtChain itself
my $defaultTRepeats = "";		# for netClass option tRepeats
my $defaultQRepeats = "";		# for netClass option qRepeats
my $defaultSeq1Limit = 30;
my $defaultSeq2Limit = 100;

sub usage {
  # Usage / help / self-documentation:
  my ($status, $detailed) = @_;
  my $base = $0;
  $base =~ s/^(.*\/)?//;
  # Basic help (for incorrect usage):
  print STDERR "
usage: $base DEF
options:
";
  print STDERR $stepper->getOptionHelp();
print STDERR <<_EOF_
    -blastzOutRoot dir    Directory path where outputs of the blastz cluster
                          run will be stored.  By default, they will be
                          stored in the $HgAutomate::clusterData build directory , but
                          this option can specify something more cluster-
                          friendly: $clusterNAS .
                          If dir does not already exist it will be created.
                          Blastz outputs are removed in the cleanup step.
    -swap                 DEF has already been used to create chains; swap
                          those chains (target for query), then net etc. in
                          a new directory:
                          $HgAutomate::clusterData/\$qDb/$HgAutomate::trackBuild/blastz.\$tDb.swap/
    -chainMinScore n      Add -minScore=n (default: $defaultChainMinScore) to the
                                  axtChain command.
    -chainLinearGap type  Add -linearGap=<loose|medium|filename> to the
                                  axtChain command.  (default: loose)
    -tRepeats table       Add -tRepeats=table to netClass (default: rmsk)
    -qRepeats table       Add -qRepeats=table to netClass (default: rmsk)
    -ignoreSelf           Do not assume self alignments even if tDb == qDb
    -syntenicNet          Perform optional syntenicNet step
    -noDbNameCheck        ignore Db name format
    -inclHap              include haplotypes *_hap* in chain/net, default not
    -loadChainSplit       load split chain tables, default is not split tables
    -swapDir path         directory to work in for swap, default:
                          /hive/data/genomes/qDb/bed/blastz.tDb.swap/
    -asmId assemblyHubId  full name for assembly hub,
                          e.g.: GCF_007474595.1_mLynCan4_v1.p
    -skipDownload         do not construct the downloads directory
    -trackHub             construct big* files for track hub
_EOF_
  ;
###DIANA AGAIN REMOVED TASMID AND QASMID AND
###COMBINED INTO ASMID
print STDERR &HgAutomate::getCommonOptionHelp('dbHost' => $dbHost,
				      'workhorse' => $workhorse,
				      'fileServer' => '',
				      'bigClusterHub' => $bigClusterHub,
				      'smallClusterHub' => $smallClusterHub);
print STDERR "
Automates UCSC's blastz/chain/net pipeline:
    1. Big cluster run of blastz.
    2. Small cluster consolidation of blastz result files.
    3. Small cluster chaining run.
    4. Sorting and netting of chains on the fileserver
       (no nets for self-alignments).
    5. Generation of liftOver-suitable chains from nets+chains on fileserver
       (not done for self-alignments).
    6. Generation of axtNet and mafNet files on the fileserver (not for self).
    7. Addition of gap/repeat info to nets on hgwdev (not for self).
    8. Loading of chain and net tables on hgwdev (no nets for self).
    9. Setup of download directory on hgwdev.
    10.Optional (-syntenicNet flag): Generation of syntenic mafNet files.
DEF is a Scott Schwartz-style bash script containing blastz parameters.
This script makes a lot of assumptions about conventional placements of
certain files, and what will be in the DEF vars.  Stick to the conventions
described in the -help output, pray to the cluster gods, and all will go
well.  :)

To use this script outside the UCSC infrastructure, use options:
    -dbHost=localhost (when there is no local genome database to load results)

    -smallClusterHub=localhost -bigClusterHub=localhost -fileServer=localhost
    This assumes the process is performed on your parasol hub machine, and
    thus all the references to 'localhost' are this parasol hub machine.
    Verify your .ssh keys are correct: 'ssh localhost' should function OK.

    -swapDir=/some/path/blastz.targetDb.swap/ work directory for -swap
    -skipDownload - leaves all constructed files in the working directory
    -trackHub - constructs bigChain and bigMaf files to use in a track hub
";
  # Detailed help (-help):
  print STDERR "
Assumptions:
1. $HgAutomate::clusterData/\$db/ is the main directory for database/assembly \$db.
   $HgAutomate::clusterData/\$tDb/$HgAutomate::trackBuild/blastz.\$qDb.\$date/ will be the directory
   created for this run, where \$tDb is the target/reference db and
   \$qDb is the query.  (Can be overridden, see #10 below.)
   $dbHost:$HgAutomate::goldenPath/\$tDb/vs\$QDb/ (or vsSelf)
   is the directory where downloadable files need to go.
   LiftOver chains (not applicable for self-alignments) go in this file:
   $HgAutomate::clusterData/\$tDb/$HgAutomate::trackBuild/liftOver/\$tDbTo\$QDb.over.chain.gz
   a copy is kept here (in case the liftOver/ copy is overwritten):
   $HgAutomate::clusterData/\$tDb/$HgAutomate::trackBuild/blastz.\$qDb.\$date/\$tDb.\$qDb.over.chain.gz
   and symbolic links to the liftOver/ file are put here:
   $dbHost:$HgAutomate::goldenPath/\$tDb/liftOver/\$tDbTo\$QDb.over.chain.gz
   $dbHost:$HgAutomate::gbdb/\$tDb/liftOver/\$tDbTo\$QDb.over.chain.gz
2. DEF's SEQ1* variables describe the target/reference assembly.
   DEF's SEQ2* variables describe the query assembly.
   If those are the same assembly, then we're doing self-alignments and
   will drop aligned blocks that cross the diagonal.
3. DEF's SEQ1_DIR is either a directory containing one nib file per
   target sequence (usually chromosome), OR a complete path to a
   single .2bit file containing all target sequences.  This directory
   should be in $clusterLocal or $clusterSortaLocal .
   SEQ2_DIR: ditto for query.
4. DEF's SEQ1_LEN is a tab-separated dump of the target database table
   chromInfo -- or at least a file that contains all sequence names
   in the first column, and corresponding sizes in the second column.
   Normally this will be $HgAutomate::clusterData/\$tDb/chrom.sizes, but for a
   scaffold-based assembly, it is a good idea to put it in $clusterSortaLocal
   or $clusterNAS
   because it will be a large file and it is read by blastz-run-ucsc
   (big cluster script).
   SEQ2_LEN: ditto for query.
5. DEF's SEQ1_CHUNK and SEQ1_LAP determine the step size and overlap size
   of chunks into which large target sequences are to be split before
   alignment.  SEQ2_CHUNK and SEQ2_LAP: ditto for query.
6. DEF's SEQ1_LIMIT and SEQ2_LIMIT decide what the maximum number of
   sequences should be for any partitioned file (the files created in the
   tParts and qParts directories).  This limit only effects SEQ1 or SEQ2
   when they are 2bit files.  Some 2bit files have too many contigs.  This
   reduces the number of blastz hippos (jobs taking forever compared to
   the other jobs).  SEQ1_LIMIT defaults to $defaultSeq1Limit and SEQ2_LIMIT defaults to $defaultSeq2Limit.
7. DEF's BLASTZ_ABRIDGE_REPEATS should be set to something nonzero if
   abridging of lineage-specific repeats is to be performed.  If so, the
   following additional constraints apply:
   a. Both target and query assemblies must be structured as one nib file
      per sequence in SEQ*_DIR (sorry, this rules out scaffold-based
      assemblies).
   b. SEQ1_SMSK must be set to a directory containing one file per target
      sequence, with the name pattern \$seq.out.spec.  This file must be
      a RepeatMasker .out file (usually filtered by DateRepeats).  The
      directory should be under $clusterLocal or $clusterSortaLocal .
      SEQ2_SMSK: ditto for query.
8. DEF's BLASTZ_[A-Z] variables will be translated into blastz command line
   options (e.g. BLASTZ_H=foo --> H=foo, BLASTZ_Q=foo --> Q=foo).
   For human-mouse evolutionary distance/sensitivity, none of these are
   necessary (blastz-run-ucsc defaults will be used).  Here's what we have
   used for human-fugu and other very-distant pairs:
BLASTZ_H=2000
BLASTZ_Y=3400
BLASTZ_L=6000
BLASTZ_K=2200
BLASTZ_Q=$HgAutomate::clusterData/blastz/HoxD55.q
   Blastz parameter tuning is somewhat of an art and is beyond the scope
   here.  Webb Miller and Jim can provide guidance on how to set these for
   a new pair of organisms.
9. DEF's PATH variable, if set, must specify a path that contains programs
   necessary for blastz to run: blastz, and if BLASTZ_ABRIDGE_REPEATS is set,
   then also fasta-subseq, strip_rpts, restore_rpts, and revcomp.
   If DEF does not contain a PATH, blastz-run-ucsc will use its own default.
10. DEF's BLASTZ variable can specify an alternate path for blastz.
11. DEF's BASE variable can specify the blastz/chain/net build directory
    (defaults to $HgAutomate::clusterData/\$tDb/$HgAutomate::trackBuild/blastz.\$qDb.\$date/).
12. SEQ?_CTGDIR specifies sequence source with the contents of full chrom
    sequences and the contig randoms and chrUn.  This keeps the contigs
    separate during the blastz and chaining so that chains won't go through
    across multiple contigs on the randoms.
13. SEQ?_CTGLEN specifies a length file to be used in conjunction with the
    special SEQ?_CTGDIR file specified above which contains the random contigs.
14. SEQ?_LIFT specifies a lift file to lift sequences in the SEQ?_CTGDIR
    to their random and chrUn positions.  This is useful for a 2bit file that
    has both full chrom sequences and the contigs for the randoms.
15. SEQ2_SELF=1 specifies the SEQ2 is already specially split for self
    alignments and to use SEQ2 sequence for self alignment, not just a
    copy of SEQ1
16. TMPDIR - specifies directory on cluster node to keep temporary files
    Typically TMPDIR=/scratch/tmp
17. All other variables in DEF will be ignored!

" if ($detailed);
  exit $status;
}


# Globals:
my %defVars = ();
my ($DEF, $tDb, $qDb, $QDb, $isSelf, $selfSplit, $buildDir, $fileServer);
my ($swapDir, $asmId, $splitRef, $inclHap, $secondsStart, $secondsEnd, $dbExists, $qDbExists);
###AGAIN,DIANA SWITCHED QUERY AND TARGET TO JUST ID, MIGHT BE TO SPLIT ALL GENOMES THAN QUERY OUT OF THIS SCRIPT
sub isInDirList {
  # Return TRUE if $dir is under (begins with) something in dirList.
  my ($dir, @dirList) = @_;
  my $pat = '^(' . join('|', @dirList) . ')(/.*)?$';
  return ($dir =~ m@$pat@);
}

sub enforceClusterNoNo {
  # Die right away if user is trying to put cluster output somewhere
  # off-limits.
  my ($dir, $desc) = @_;
  if (&isInDirList($dir, @clusterNoNo)) {
    die "\ncluster outputs are forbidden to go to " .
      join (' or ', @clusterNoNo) . " so please choose a different " .
      "$desc instead of $dir .\n\n";
  }
  # use this only if it exists, this is UCSC infrastructure:
  if ( -e $getFileServer ) {
    my $testFileServer = `$getFileServer $dir/`;
    if (scalar(grep /^$testFileServer$/, @fileServerNoNo)) {
      die "\ncluster outputs are forbidden to go to fileservers " .
        join (' or ', @fileServerNoNo) . " so please choose a different " .
        "$desc instead of $dir (which is hosted on $testFileServer).\n\n";
    }
  }
}

#DIANA COMBINING T AND Q INTO A AGAIN IN SUB BELOW
sub checkOptions {
  # Make sure command line options are valid/supported.
  my $ok = GetOptions(@HgStepManager::optionSpec,
		      @HgAutomate::commonOptionSpec,
		      "blastzOutRoot=s",
		      "swap",
		      "chainMinScore=i",
		      "chainLinearGap=s",
		      "tRepeats=s",
		      "qRepeats=s",
		      "readmeOnly",
		      "ignoreSelf",
                      "syntenicNet",
                      "noDbNameCheck",
                      "inclHap",
                      "noLoadChainSplit",
                      "loadChainSplit",
                      "swapDir=s",
                      "asmId=s",
                      "skipDownload",
                      "trackHub"
		     );
  &usage(1) if (!$ok);
  &usage(0, 1) if ($opt_help);
  &HgAutomate::processCommonOptions();
  my $err = $stepper->processOptions();
  usage(1) if ($err);
  $dbHost = $opt_dbHost if ($opt_dbHost);
  if ($opt_swap) {
    if ($opt_continue) {
      if ($stepper->stepPrecedes($opt_continue, 'net')) {
	warn "\nIf -swap is specified, then -continue must specify a step ".
	  "of \"net\" or later.\n";
	&usage(1);
      }
    } else {
      # If -swap is given but -continue is not, force -continue and tell
      # $stepper to reevaluate options:
      $opt_continue = 'chainMerge';
      $err = $stepper->processOptions();
      usage(1) if ($err);
    }
    if ($opt_stop) {
      if ($stepper->stepPrecedes($opt_stop, 'chainMerge')) {
	warn "\nIf -swap is specified, then -stop must specify a step ".
	"of \"chainMerge\" or later.\n";
	&usage(1);
      }
    }
  }
  if ($opt_blastzOutRoot) {
    if ($opt_blastzOutRoot !~ m@^/\S+/\S+@) {
      warn "\n-blastzOutRoot must specify a full path.\n";
      &usage(1);
    }
    &enforceClusterNoNo($opt_blastzOutRoot, '-blastzOutRoot');
    if (! &isInDirList($opt_blastzOutRoot, @clusterNAS)) {
      warn "\n-blastzOutRoot is intended to specify something on " .
	"$clusterNAS, but I'll trust your judgment " .
	"and use $opt_blastzOutRoot\n\n";
    }
  }
  $workhorse = $opt_workhorse if ($opt_workhorse);
  $bigClusterHub = $opt_bigClusterHub if ($opt_bigClusterHub);
  $smallClusterHub = $opt_smallClusterHub if ($opt_smallClusterHub);
}

#########################################################################
# The following routines were taken almost verbatim from blastz-run-ucsc,
# so may be good candidates for libification!  unless that would slow down
# blastz-run-ucsc...
# nfsNoodge() was removed from loadDef() and loadSeqSizes() -- since this
# script will not be run on the cluster, we should fully expect files to
# be immediately visible.

sub loadDef {
  # Read parameters from a bash script with Scott's param variable names:
  my ($def) = @_;
  my $fh = &HgAutomate::mustOpen("$def");
  while (<$fh>) {
    s/^\s*export\s+//;
    next if (/^\s*#/ || /^\s*$/);
    if (/(\w+)\s*=\s*(.*)/) {
      my ($var, $val) = ($1, $2);
      while ($val =~ /\$(\w+)/) {
	my $subst = $defVars{$1};
	if (defined $subst) {
	  $val =~ s/\$$1/$subst/;
	} else {
	  die "Can't find value to substitute for \$$1 in $DEF var $var.\n";
	}
      }
      $defVars{$var} = $val;
    }
  }
  close($fh);
}

sub loadSeqSizes {
  # Load up sequence -> size mapping from $sizeFile into $hashRef.
  my ($sizeFile, $hashRef) = @_;
  my $fh = &HgAutomate::mustOpen("$sizeFile");
  while (<$fh>) {
    chomp;
    my ($seq, $size) = split;
    $hashRef->{$seq} = $size;
  }
  close($fh);
}

# end shared stuff from blastz-run-ucsc
#########################################################################

sub requireVar {
  my ($var) = @_;
  die "Error: $DEF is missing variable $var\n" if (! defined $defVars{$var});
}

sub requirePath {
  my ($var) = @_;
  my $val = $defVars{$var};
  die "Error: $DEF $var=$val must specify a complete path\n"
    if ($val !~ m@^/\S+/\S+@);
  if ( -d $val ) {
    my $fileCount = `find $val -maxdepth 1 -type f | wc -l`;
    chomp $fileCount;
    if ($fileCount < 1) {
	die "Error: $DEF variable: $var=$val specifies an empty directory.\n";
    }
  } elsif ( ! -s $val ) {
    die "Error: $DEF variable: $var=$val is not a file or directory.\n";
  }
}

sub requireNum {
  my ($var) = @_;
  my $val = $defVars{$var};
  die "Error: $DEF variable $var=$val must specify a number.\n"
    if ($val !~ /^\d+$/);
}

my $oldDbFormat = '[a-z][a-z](\d+)?';
my $newDbFormat = '[a-z][a-z][a-z][A-Z][a-z][a-z0-9](\d+)?';
my $patchDbFormat = 'grc[A-Z][0-9]+P[0-9]+';
sub getDbFromPath {
  # Require that $val is a full path that contains a recognizable db as
  # one of its elements (possibly the last one).
  my ($var) = @_;
  my $val = $defVars{$var};
  my $db;
  my $dbFromName = basename($val);
  $dbFromName =~ s/.2bit//;
  if (! $opt_noDbNameCheck) {
    if ( $val =~ m@^/\S+/($oldDbFormat|$newDbFormat|$patchDbFormat)((\.2bit)|(/(\S+)?))?$@) {
      $db = $1;
    } else {
      die "Error: $DEF variable $var=$val must be a full path with " .
        "a recognizable database as one of its elements.\n"
    }
  }
  if ($opt_noDbNameCheck) {
    $db = $dbFromName;
  } else {
    if (! defined($db)) {
      if ($val =~ m#^/hive/data/genomes/#) {
	$val =~ s#^/hive/data/genomes/##;
	$val =~ s#/.*##;
	$db = $val;
	warn "Warning: assuming database $db from /hive/data/genomes/<db>/ path\n";
      } elsif ($val =~ m#^/scratch/data/#) {
	$val =~ s#^/scratch/data/##;
	$val =~ s#/.*##;
	$db = $val;
	warn "Warning: assuming database $db from /scratch/data/<db>/ path\n";
      }
    }
  }
return $db;
}

sub checkDef {
  # Make sure %defVars contains what we need and looks consistent with
  # our assumptions.
  foreach my $s ('SEQ1_', 'SEQ2_') {
    foreach my $req ('DIR', 'LEN', 'CHUNK', 'LAP') {
      &requireVar("$s$req");
    }
    &requirePath($s . 'DIR');
    &requirePath($s . 'LEN');
    &requireNum($s . 'CHUNK');
    &requireNum($s . 'LAP');
  }
  $tDb = &getDbFromPath('SEQ1_DIR');
  $qDb = &getDbFromPath('SEQ2_DIR');
  $isSelf = $opt_ignoreSelf ? 0 : ($tDb eq $qDb);
  # special split on SEQ2 for Self alignments
  $selfSplit = $defVars{'SEQ2_SELF'} || 0;
  $QDb = $isSelf ? 'Self' : ucfirst($qDb);
  if ($isSelf && $opt_swap) {
    die "-swap is not supported for self-alignments\n" .
        "($DEF has $tDb as both target and query).\n";
  }
  HgAutomate::verbose(1, "$DEF looks OK!\n" .
	  "\ttDb=$tDb\n\tqDb=$qDb\n\ts1d=$defVars{SEQ1_DIR}\n" .
	  "\tisSelf=$isSelf\n");
  if ($defVars{'SEQ1_SMSK'} || $defVars{'SEQ2_SMSK'} ||
      $defVars{'BLASTZ_ABRIDGE_REPEATS'}) {
    &requireVar('BLASTZ_ABRIDGE_REPEATS');
    foreach my $s ('SEQ1_', 'SEQ2_') {
      my $var = $s. 'SMSK';
      &requireVar($var);
      &requirePath($var);
    }
    HgAutomate::verbose(1, "Abridging repeats!\n");
  }
}

sub doPartition {
  # Partition the sequence up before blastz.
  my $paraHub = $opt_blastzOutRoot ? $bigClusterHub : $workhorse;
  my $runDir = "$buildDir/run.blastz";
  my $targetList = "$tDb.lst";
  my $queryList = $isSelf ? $targetList :
	($opt_ignoreSelf ? "$qDb.ignoreSelf.lst" : "$qDb.lst");
  if ($selfSplit) {
    $queryList = "$qDb.selfSplit.lst"
  }
  my $tPartDir = '-lstDir tParts';
  my $qPartDir = '-lstDir qParts';
  my $outRoot = $opt_blastzOutRoot ? "$opt_blastzOutRoot/psl" : '../psl';

  my $seq1Dir = $defVars{'SEQ1_CTGDIR'} || $defVars{'SEQ1_DIR'};
  my $seq2Dir = $defVars{'SEQ2_CTGDIR'} || $defVars{'SEQ2_DIR'};
  my $seq1Len = $defVars{'SEQ1_CTGLEN'} || $defVars{'SEQ1_LEN'};
  my $seq2Len = $defVars{'SEQ2_CTGLEN'} || $defVars{'SEQ2_LEN'};
  my $seq1Limit = (defined $defVars{'SEQ1_LIMIT'}) ? $defVars{'SEQ1_LIMIT'} :
    $defaultSeq1Limit;
  my $seq2Limit = (defined $defVars{'SEQ2_LIMIT'}) ? $defVars{'SEQ2_LIMIT'} :
    $defaultSeq2Limit;
  my $seq2MaxLength = `awk '{print \$2}' $seq2Len | sort -rn | head -1`;
  chomp $seq2MaxLength;
  my $bundleParts = 0;
  # OK to bundle parts list bits into 2bit files when not abridging
  $bundleParts = 1 if ( ! $defVars{'BLASTZ_ABRIDGE_REPEATS'} );

  my $partitionTargetCmd =
    ("$partition $defVars{SEQ1_CHUNK} $defVars{SEQ1_LAP} " .
     "$seq1Dir $seq1Len -xdir xdir.sh -rawDir $outRoot $seq1Limit " .
     "$tPartDir > $targetList");
  my $partitionQueryCmd =
    (($isSelf && (! $selfSplit)) ?
     '# Self-alignment ==> use target partition for both.' :
     "$partition $defVars{SEQ2_CHUNK} $defVars{SEQ2_LAP} " .
     "$seq2Dir $seq2Len $seq2Limit " .
     "$qPartDir > $queryList");
  &HgAutomate::mustMkdir($runDir);
  my $whatItDoes =
"It computes partitions of target and query sequences into chunks of the
specified size for the blastz cluster run.  The actual splitting of
sequence is not performed here, but later on by blastz cluster jobs.";
  my $bossScript = newBash HgRemoteScript("$runDir/doPartition.bash", $paraHub,
				      $runDir, $whatItDoes, $DEF);
  $bossScript->add(<<_EOF_
$partitionTargetCmd
export L1=`wc -l < $targetList`
$partitionQueryCmd
export L2=`wc -l < $queryList`
export L=`echo \$L1 \$L2 | awk '{print \$1*\$2}'`
echo "cluster batch jobList size: \$L = \$L1 * \$L2"
_EOF_
    );
  if ($bundleParts) {
  $bossScript->add(<<_EOF_
if [ -d tParts ]; then
  echo 'constructing tParts/*.2bit files'
  ls tParts/*.lst | sed -e 's#tParts/##; s#.lst##;' | while read tPart
  do
    sed -e 's#.*.2bit:##;' tParts/\$tPart.lst \\
      | twoBitToFa -seqList=stdin $seq1Dir stdout \\
        | faToTwoBit stdin tParts/\$tPart.2bit
  done
fi
if [ -d qParts ]; then
  echo 'constructing qParts/*.2bit files'
  ls qParts/*.lst | sed -e 's#qParts/##; s#.lst##;' | while read qPart
  do
    sed -e 's#.*.2bit:##;' qParts/\$qPart.lst \\
      | twoBitToFa -seqList=stdin $seq2Dir stdout \\
        | faToTwoBit stdin qParts/\$qPart.2bit
  done
fi
_EOF_
    );
  }
  $bossScript->execute();
  my $mkOutRootHost = $opt_blastzOutRoot ? $paraHub : $fileServer;
  my $mkOutRoot =     $opt_blastzOutRoot ? "mkdir -p $opt_blastzOutRoot;" : "";
  &HgAutomate::run("$HgAutomate::runSSH $mkOutRootHost " .
		   "'(cd $runDir; $mkOutRoot sh -ef xdir.sh)'");
}

sub doBlastzClusterRun {
  # Set up and perform the big-cluster blastz run.
  my $paraHub = $bigClusterHub;
  my $runDir = "$buildDir/run.blastz";
  my $targetList = "$tDb.lst";
  my $outRoot = $opt_blastzOutRoot ? "$opt_blastzOutRoot/psl" : '../psl';
  my $queryList = $isSelf ? $targetList :
	($opt_ignoreSelf ? "$qDb.ignoreSelf.lst" : "$qDb.lst");
  if ($selfSplit) {
    $queryList = "$qDb.selfSplit.lst"
  }
###DIANA REMOVED THE FIRST PART OF THIS THAT KILLS THE SCRIPT IF WE ALREADY HAVE FILES
###PARTITIONED
###ASSUMING AS WE ARE AVOIDING USING PARASOL AND  THIS WILL BE
###BROKEN UP INTO CHUNKS THAT WE DO NOT WANT IT TO END BECAUSE WE HAVE FILES
  # First, make sure we got through the partitioning already
  if (! -e "$runDir/$targetList" && ! $opt_debug) {
    die "doBlastzClusterRun: there's no target list file " .
        "so start over without the -continue align.\n";
  }
  if (! -e "$runDir/$queryList" && ! $opt_debug) {
    die "doBlastzClusterRun: there's no query list file" .
        "so start over without the -continue align.\n";
  }
  my $templateCmd = ("$blastzRunUcsc -outFormat psl " .
		     ($isSelf ? '-dropSelf ' : '') .
		     '$(path1) $(path2) ../DEF ' .
		     #'{check out exists ' .
		     $outRoot . '/$(file1)/$(file1)_$(file2).psl');
  &HgAutomate::makeGsub($runDir, $templateCmd);
  `touch "$runDir/para_hub_$paraHub"`;
  my $whatItDoes = "It sets up and performs the big cluster blastz run.";
  #HAS A NEW COUNT VARIABLE THAT MIGHT BE IN THE SCRIPT COUNTING FILES/PARTS
  my $countvariable = "\$COUNT";
  my $bossScript = new HgRemoteScript("$runDir/doClusterRun.csh", $paraHub,
				      $runDir, $whatItDoes, $DEF);
  ###DIANA GETS RID OF THIS LINE BECAUSE PARA WONT WORK# my $paraRun = &HgAutomate::paraRun();
  my $gensub2 = &HgAutomate::gensub2();
  $bossScript->add(<<_EOF_
$gensub2 $targetList $queryList gsub jobList
###AS WELL AS THIS ONE SO IT DOESN'T TRY  AND SUBMIT#pararun#
cp /gpfs/scratch/withomas/project_GenomeAnnot/scripts/step1_alignment/Diana_scripts/scriptSetup.sh . #copy the slurm script for 1-1990 jobs
### BUT COPIES A STARTER SLURM SCRIPT TO THE CURRENT DIRECTORY, AND THEN RUNS IT BELOW
sh scriptSetup.sh 
_EOF_
    );
  $bossScript->execute();
}	#	sub doBlastzClusterRun {}

sub doCatRun {
  # Do a small cluster run to concatenate the lowest level of chunk result
  # files from the big cluster blastz run.  This brings results up to the
  # next level: per-target-chunk results, which may still need to be
  # concatenated into per-target-sequence in the next step after this one --
  # chaining.
  my $paraHub = $smallClusterHub;
  my $runDir = "$buildDir/run.cat";
###DIANA REMOVING THESE 9 LINES THAT MAKES SURE NOTHING IS IN THE FOLDERS BECAUSE BREAKING
###UP PARASOL PIPELINE
###
###
###  
  # First, make sure previous stage was successful
  my $successFile = "$buildDir/run.blastz/failed_jobs.txt";
  if (-e $successFile) {
    die "doCatRun: looks like doBlastzChainNet was not completely successful" .
      " Check failed_jobs.txt to see which job failed.  Run again the failed jobs" .
	  "and Either run with -continue chainRun or some later stage,  and run again.\n";
    }
### DIANA CHANGED FILE TO SUCCESS FILE TO FAILED JOB TO MAKE SURE THERE WERE NO ISSUE AND ALL FILES THERE INSTEAD OF JUST RUN TIME BECAUSE PARASOL WONT BE USED  
  &HgAutomate::mustMkdir($runDir);
  &HgAutomate::makeGsub($runDir,
      "./cat.csh \$(path1) ../pslParts/\$(file1).psl.gz");
  `touch "$runDir/para_hub_$paraHub"`;

  my $outRoot = $opt_blastzOutRoot ? "$opt_blastzOutRoot/psl" : '../psl';

  my $fh = &HgAutomate::mustOpen(">$runDir/cat.csh");
  print $fh <<_EOF_
#!/bin/csh -ef
find $outRoot/\$1/ -name "*.psl" | xargs cat | gzip -c > \$2
_EOF_
  ;
  close($fh);

  my $whatItDoes =
"It sets up and performs a small cluster run to concatenate all files in
each subdirectory of $outRoot into a per-target-chunk file.";
  my $countvariable = "\$COUNT";
  my $bossScript = new HgRemoteScript("$runDir/doCatRun.csh", $paraHub,
				      $runDir, $whatItDoes, $DEF);
  ###DIANA MAKING SURE IT DOES NOT RUN ON PARASOL AND ADDING IN COUNT VARIABLE#my $paraRun = &HgAutomate::paraRun();
  my $gensub2 = &HgAutomate::gensub2();
  $bossScript->add(<<_EOF_
(cd $outRoot; find . -maxdepth 1 -type d | grep '^./') \\
        | sed -e 's#/\$##; s#^./##' > tParts.lst
chmod a+x cat.csh
$gensub2 tParts.lst single gsub jobList
mkdir -p ../pslParts
# DIANA GETS RID OF PARARUN AND COPYS AND RUNS SLURM SCRIPT BELOW
cp /gpfs/scratch/withomas/project_GenomeAnnot/scripts/step1_alignment/Diana_scripts/short_scriptSetup.sh . #copy the slurm script for 1-1990 jobs
sh short_scriptSetup.sh 
_EOF_
    );
  $bossScript->execute();
}	#	sub doCatRun {}

sub makePslPartsLst {
  # Create a pslParts.lst file the subdirectories of pslParts; if some
  # are for subsequences of the same sequence, make a single .lst line
  # for the sequence (single chaining job with subseqs' alignments
  # catted together).  Otherwise (i.e. subdirs that contain small
  # target seqs glommed together by partitionSequences) make one .lst
  # line per partition.
  return if ($opt_debug);
  opendir(P, "$buildDir/pslParts")
    || die "Couldn't open directory $buildDir/pslParts for reading: $!\n";
  my @parts = readdir(P);
  closedir(P);
  my $partsLst = "$buildDir/axtChain/run/pslParts.lst";
  my $fh = &HgAutomate::mustOpen(">$partsLst");
  my %seqs = ();
  my $count = 0;
  foreach my $p (@parts) {
    $p =~ s@^/.*/@@;  $p =~ s@/$@@;
    $p =~ s/\.psl\.gz//;
    next if ($p eq '.' || $p eq '..');
    if ($p =~ m@^(\S+:\S+):\d+-\d+$@) {
      # Collapse subsequences (subranges of a sequence) down to one entry
      # per sequence:
      $seqs{$1} = 1;
    } else {
      print $fh "$p\n";
      $count++;
    }
  }
  foreach my $p (keys %seqs) {
    print $fh "$p:\n";
    $count++;
  }
  close($fh);
  if ($count < 1) {
    die "makePslPartsLst: didn't find any pslParts/ items.";
  }
}


sub doChainRun {
  # Do a small cluster run to chain alignments to each target sequence.
  my $paraHub = $smallClusterHub;
  my $runDir = "$buildDir/axtChain/run";
###DIANA REMOVING CLEAN STEP AGAIN FOR REASONS ABOVE
###
###
###
###
###
###
### AGAIN BELOW IS  MAKING SURE  THAT A FAILED JOBS TEXT DOES NOT EXIST BEFORE MOVING ON, NOT RUN  TIME NOT EXISTING
  # First, make sure previous stage was successful
  my $successFile = "$buildDir/run.cat/failed_jobs.txt";
  if (-e $successFile) {
    die "doChainRun: looks like doCatRun was not completely successful" .
      " Check failed_jobs.txt to see which job failed.  Run again the failed jobs" .
	  "and Either run with -continue chainRun or some later stage,  and run again.\n";
    }
  &HgAutomate::mustMkdir($runDir);
  &HgAutomate::makeGsub($runDir,
	       "sh chain.csh \$(file1) chain/\$(file1).chain");
  `touch "$runDir/para_hub_$paraHub"`;

  my $seq1Dir = $defVars{'SEQ1_CTGDIR'} || $defVars{'SEQ1_DIR'};
  my $seq2Dir = $defVars{'SEQ2_CTGDIR'} || $defVars{'SEQ2_DIR'};
  my $matrix = $defVars{'BLASTZ_Q'} ? "-scoreScheme=$defVars{BLASTZ_Q} " : "";
  my $minScore = $opt_chainMinScore ? "-minScore=$opt_chainMinScore" : "";
  my $linearGap = $opt_chainLinearGap ? "-linearGap=$opt_chainLinearGap" :
	"-linearGap=$defaultChainLinearGap";
  my $fh = &HgAutomate::mustOpen(">$runDir/chain.csh");
  print $fh  <<_EOF_
#!/bin/csh -ef
zcat ../../pslParts/\$1*.psl.gz \\
| axtChain -psl -verbose=0 $matrix $minScore $linearGap stdin \\
    $seq1Dir \\
    $seq2Dir \\
    stdout \\
| chainAntiRepeat $seq1Dir \\
    $seq2Dir \\
    stdin \$2
_EOF_
    ;
  if (exists($defVars{'SEQ1_LIFT'})) {
  print $fh <<_EOF_
set c=\$2:t:r
echo "lifting \$2 to \${c}.lifted.chain"
liftUp liftedChain/\${c}.lifted.chain \\
    $defVars{'SEQ1_LIFT'} carry \$2
rm \$2
mv liftedChain/\${c}.lifted.chain \$2
_EOF_
    ;
  }
  if (exists($defVars{'SEQ2_LIFT'})) {
  print $fh <<_EOF_
set c=\$2:t:r
echo "lifting \$2 to \${c}.lifted.chain"
liftUp -chainQ liftedChain/\${c}.lifted.chain \\
    $defVars{'SEQ2_LIFT'} carry \$2
rm \$2
mv liftedChain/\${c}.lifted.chain \$2
_EOF_
    ;
  }
  close($fh);

  &makePslPartsLst();

  my $whatItDoes =
"It sets up and performs a small cluster run to chain all alignments
to each target sequence.";
  my $bossScript = new HgRemoteScript("$runDir/doChainRun.csh", $paraHub,
				      $runDir, $whatItDoes, $DEF);
  my $countvariable = "\$COUNT";
  my $gensub2 = &HgAutomate::gensub2();
  $bossScript->add(<<_EOF_
chmod a+x chain.csh
$gensub2 pslParts.lst single gsub jobList
mkdir -p chain liftedChain
#HERE DIANA IS COPYING OVER HER  SLURM SETUP FILE, AND REMOVING PARARUN LINE
cp /gpfs/scratch/withomas/project_GenomeAnnot/scripts/step1_alignment/Diana_scripts/short_scriptSetup.sh . #copy the slurm script for 1-1990 jobs
#WHILE ALSO NOT PURGING THE DIRECTORY AFTER BY REMOVING RMDIR LIFTEDCHAIN
sh short_scriptSetup.sh 
_EOF_
  );
  $bossScript->execute();
}	#	sub doChainRun {}

sub postProcessChains {
  # chainMergeSort etc.
  my $runDir = "$buildDir/axtChain";
  my $chain = "$tDb.$qDb.all.chain.gz";
### AGAIN WE ARE N OT STARTING CLEAN HERE REMOVING LINES
###
###
###
###
###
###
###
###
### AGAIN CHANGING SUCCESS FILE
  # First, make sure previous stage was successful
  my $successFile = "$buildDir/axtChain/run/failed_jobs.txt";
  if (-e $successFile) {
    die "postProcessChains: looks like doChainRun was not completely successful" .
      " Check failed_jobs.txt to see which job failed.  Run again the failed jobs" .
	  "and Either run with -continue chainRun or some later stage,  and run again.\n";
    }
  my $cmd="$HgAutomate::runSSH $workhorse nice ";
  $cmd .= "'find $runDir/run/chain -name \"*.chain\" ";
  $cmd .= "| chainMergeSort -inputList=stdin ";
  $cmd .= "| nice gzip -c > $runDir/$chain'";
  &HgAutomate::run($cmd);
  if ($splitRef) {
    &HgAutomate::run("$HgAutomate::runSSH $fileServer nice " .
	 "chainSplit $runDir/chain $runDir/$chain");
  }
  &HgAutomate::nfsNoodge("$runDir/$chain");
}	#	sub postProcessChains {}


sub getAllChain {
  # Find the most likely candidate for all.chain from a previous run/step.
  my ($runDir) = @_;
  my $chain;
  if (-e "$runDir/$tDb.$qDb.all.chain.gz") {
    $chain = "$tDb.$qDb.all.chain.gz";
  } elsif (-e "$runDir/$tDb.$qDb.all.chain") {
    $chain = "$tDb.$qDb.all.chain";
  } elsif (-e "$runDir/all.chain.gz") {
    $chain = "all.chain.gz";
  } elsif (-e "$runDir/all.chain") {
    $chain = "all.chain";
  } elsif ($opt_debug) {
    $chain = "$tDb.$qDb.all.chain.gz";
  }
  return $chain;
}

###  DIANA REMOVED 2 DIFFERENT FUNCTIONS
### SWAP CHAINS
### AND
### SWAPGLOBALS
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###
###


sub doChainMerge {
  # If -swap, swap chains from other org;  otherwise, merge the results
  # from the chainRun step.
  if ($opt_swap) {
    &swapChains();
    &swapGlobals();
  } else {
    &postProcessChains();
  }
}

###DIANA REMOVED ABOUT 200 LINES HERE
###ALMOSTS ALL HAVING TO DEAL WITH NETTING
###PERHAPS SHE SPPLIT THE CODE TO GET ALL CHAINS FIRST BEFORE NETTING OR
### NETS SEPARATLY

sub getBlastzParams {
  my %vars;
  # Return parameters in BLASTZ_Q file, or defaults, for README.txt.
  my $matrix =
"           A    C    G    T
      A   91 -114  -31 -123
      C -114  100 -125  -31
      G  -31 -125  100 -114
      T -123  -31 -114   91";
  if ($defVars{'BLASTZ_Q'}) {
    my $readLineLimit = 100;  # safety valve to get out if reading nonsense
    my $linesRead = 0;
    my $fh = &HgAutomate::mustOpen($defVars{'BLASTZ_Q'});
    my $line;
    my $matrixFound = 0;
    while (!$matrixFound && ($linesRead < $readLineLimit) && ($line = <$fh>)) {
      ++$linesRead;
      next if (($line =~ m/^#/) || ($line =~ m/^$/));
      if ($line =~ m/^\s*A\s+C\s+G\s+T\s*$/) {
        $matrixFound = 1;
      } else {
         chomp $line;
         $line =~ s/\s+//g;
         $line =~ s/#.*//;
         die "can not find tag=value in $defVars{BLASTZ_Q}" if ($line !~ /=/);
         my ($tag, $value) = split('=',$line);
         # ignore O E gap_open_penalty gap_extend_penalty
         next if ($tag eq "O" || $tag eq "E"
               || $tag eq "gap_open_penalty" || $tag eq "gap_extend_penalty");
         $vars{$tag} = $value;
      }
    }
    die "can not find score matrix in $defVars{BLASTZ_Q}" if (!$matrixFound);
    $line =~ s/^   // if (length($line) > 22);
    $matrix = '        ' . $line;
    foreach my $base ('A', 'C', 'G', 'T') {
      $line = <$fh>;
      die "Too few lines of $defVars{BLASTZ_Q}" if (! $line);
      if ($line !~ /^[ACGT]?\s*-?\d+\s+-?\d+\s+-?\d+\s+-?\d+\s*$/) {
	die "Can't parse this line of $defVars{BLASTZ_Q}:\n$line";
      }
      $line =~ s/^[ACGT] //;
      $matrix .= "      $base " . $line;
    }
    chomp $matrix;
    $line = <$fh>;
    if ($line && $line =~ /\S/) {
      warn "\nWarning: BLASTZ_Q matrix file $defVars{BLASTZ_Q} has " .
           "additional contents after the matrix -- those are ignored " .
	   "by blastz.\n\n";
    }
    close($fh);
  }
  my $o = $defVars{'BLASTZ_O'} || 400;
  my $e = $defVars{'BLASTZ_E'} || 30;
  my $k = $defVars{'BLASTZ_K'} || 3000;
  my $l = $defVars{'BLASTZ_L'} || 3000;
  my $h = $defVars{'BLASTZ_H'} || 2000;
  my $blastzOther = '';
  foreach my $var (sort keys %defVars) {
    if ($var =~ /^BLASTZ_(\w)$/) {
      my $p = $1;
      if ($p ne 'K' && $p ne 'L' && $p ne 'H' && $p ne 'Q') {
	if ($blastzOther eq '') {
	  $blastzOther = 'Other lastz
parameters specifically set for this species pair:';
	}
	$blastzOther .= "\n    $p=$defVars{$var}";
      }
    }
  }
  return ($matrix, $o, $e, $k, $l, $h, $blastzOther);
}

sub commafy {
  # Assuming $num is a number, add commas where appropriate.
  my ($num) = @_;
  $num =~ s/(\d)(\d\d\d)$/$1,$2/;
  $num =~ s/(\d)(\d\d\d),/$1,$2,/g;
  return($num);
}

sub describeOverlapping {
  # Return some text describing how large sequences were split.
  my $lap;
  my $chunkPlusLap1 = $defVars{'SEQ1_CHUNK'} + $defVars{'SEQ1_LAP'};
  my $chunkPlusLap2 = $defVars{'SEQ2_CHUNK'} + $defVars{'SEQ2_LAP'};
  if ($chunkPlusLap1 == $chunkPlusLap2) {
    $lap .= "Any sequences larger\n" .
"than " . &commafy($chunkPlusLap1) . " bases were split into chunks of " .
&commafy($chunkPlusLap1) . " bases
overlapping by " . &commafy($defVars{SEQ1_LAP}) . " bases for alignment.";
  } else {
    $lap .= "Any $tDb sequences larger\n" .
"than " . &commafy($chunkPlusLap1) . " bases were split into chunks of " .
&commafy($chunkPlusLap1) . " bases overlapping
by " . &commafy($defVars{SEQ1_LAP}) . " bases for alignment.  " .
"A similar process was followed for $qDb,
with chunks of " . &commafy($chunkPlusLap2) . " overlapping by " .
&commafy($defVars{SEQ2_LAP}) . ".";
  }
  $lap .= "  Following alignment, the
coordinates of the chunk alignments were corrected by the
blastz-normalizeLav script written by Scott Schwartz of Penn State.";
  return $lap;
}

###DIANA REMOVING READ ME BEING MADE, PROBABLY UNNECESARY TO TRY AND MODIFY A README SCRIPT
### WHEN THINGS ABOVE HAVE CHANGED SO MUCH, ALSO GOT  RID  OF REFENCES AS IS HER SCRIPT
### YUP LOOKED BELOW AND DIANA CREATED READ ME ON HER OWN
###SO FROM HERE GETS RID OF 
###sub dumpDownloadReadme
###sub installDownloads
### NOT DOWNLOADING ANY OF THE GENOMES IN THE SCRIPT
###sub doDownloads
###sub cleanup
###NO NEED TO CLEAN UP DOWNLOADED FILES IF NOT DOWNLOADING
###sub doSyntenicNet
###NOT NETTING IN THIS SCRIPT
###
###
###
###
###


#########################################################################
#
# -- main --

# Prevent "Suspended (tty input)" hanging:
&HgAutomate::closeStdin();

#$opt_debug = 1;

&checkOptions();

&usage(1) if (scalar(@ARGV) != 1);
$secondsStart = `date "+%s"`;
chomp $secondsStart;
($DEF) = @ARGV;

$inclHap = "";
$inclHap = "-inclHap" if ($opt_inclHap);
&loadDef($DEF);
&checkDef();

my $seq1IsSplit = (`wc -l < $defVars{SEQ1_LEN}` <=
		   $HgAutomate::splitThreshold);
my $seq2IsSplit = (`wc -l < $defVars{SEQ2_LEN}` <=
		   $HgAutomate::splitThreshold);

# might be an assembly hub build
$asmId = $opt_asmId ? $opt_asmId : "";

# Undocumented option for quickly generating a README from DEF:
if ($opt_readmeOnly) {
  $splitRef = $opt_swap ? $seq2IsSplit : $seq1IsSplit;
  &swapGlobals() if $opt_swap;
  &dumpDownloadReadme("/tmp/README.txt");
  exit 0;
}

my $date = `date +%Y-%m-%d`;
chomp $date;
$buildDir = $defVars{'BASE'} ||
  "$HgAutomate::clusterData/$tDb/$HgAutomate::trackBuild/blastz.$qDb.$date";

if ($opt_swap) {
  my $inChain = &getAllChain("$buildDir/axtChain");
  if (! defined $inChain) {
    die "-swap: Can't find $buildDir/axtChain/[$tDb.$qDb.]all.chain[.gz]\n" .
        "which is required for -swap.\n";
  }
  if ($opt_swapDir) {
    $swapDir = $opt_swapDir;
  } else {
    $swapDir = "$HgAutomate::clusterData/$qDb/$HgAutomate::trackBuild/blastz.$tDb.swap";
  }
  &HgAutomate::mustMkdir("$swapDir/axtChain");
  $splitRef = $seq2IsSplit;
  &HgAutomate::verbose(1, "Swapping from $buildDir/axtChain/$inChain\n" .
	      "to $swapDir/axtChain/$qDb.$tDb.all.chain.gz .\n");
} else {
  if (! -d $buildDir) {
    &HgAutomate::mustMkdir($buildDir);
  }
if (! $opt_blastzOutRoot &&
  $stepper->stepPrecedes($stepper->getStartStep(), 'chainRun')) {
    &enforceClusterNoNo($buildDir,
	    'blastz/chain/net build directory (or use -blastzOutRoot)');
  }
  $splitRef = $seq1IsSplit;
  &HgAutomate::verbose(1, "Building in $buildDir\n");
}

if (! -e "$buildDir/DEF") {
  &HgAutomate::run("cp $DEF $buildDir/DEF");
}

$fileServer = &HgAutomate::chooseFileServer($opt_swap ? $swapDir : $buildDir);

# may be working on a 2bit file that does not have a database browser
$dbExists = 0;
$dbExists = 1 if (&HgAutomate::databaseExists($dbHost, $tDb));
# may be working with a query that does not have a database
$qDbExists = 0;
$qDbExists = 1 if (&HgAutomate::databaseExists($dbHost, $qDb));

# When running -swap, swapGlobals() happens at the end of the chainMerge step.
# However, if we also use -continue with some step later than chainMerge, we
# need to call swapGlobals before executing the remaining steps.
if ($opt_swap &&
    $stepper->stepPrecedes('chainMerge', $stepper->getStartStep())) {
  &swapGlobals();
}

$stepper->execute();

$secondsEnd = `date "+%s"`;
chomp $secondsEnd;
my $elapsedSeconds = $secondsEnd - $secondsStart;
my $elapsedMinutes = int($elapsedSeconds/60);
$elapsedSeconds -= $elapsedMinutes * 60;

HgAutomate::verbose(1,
	"\n *** All done !  Elapsed time: ${elapsedMinutes}m${elapsedSeconds}s\n");
HgAutomate::verbose(1,
	" *** Make sure that goldenPath/$tDb/vs$QDb/README.txt is accurate.\n")
  if ($stepper->stepPrecedes('load', $stepper->getStopStep()));
HgAutomate::verbose(1,
	" *** Add {chain,net}$QDb tracks to trackDb.ra if necessary.\n")
  if ($stepper->stepPrecedes('net', $stepper->getStopStep()));
HgAutomate::verbose(1,
	"\n\n");

