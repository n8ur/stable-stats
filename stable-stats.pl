#!/usr/bin/perl -w

# stable-stat.pl
# version 0.99a -- generate statistics from file containing MJD and phase data
# Dec. 15, 2015
#
# Copyright 2005-2015 by John R. Ackermann  N8UR (jra@febo.com)
# Licensed under the GPL version 2 or later; see the file COPYING
# included with this distribution.  I request, but do not require, that
# any modifications that correct bugs or errors, or increase the program's
# functionality, be sent via email to the author at the address above.

use strict;
#use bignum;
use POSIX qw(tmpnam);
use IO::File;
use Getopt::Std;
use Statistics::OLS;
use n8ur qw(lower_case trim squash collapse round);
use n8ur_gpib qw(checkSRQ serviceSRQ logline);

###############################################################################
# local paths
###############################################################################

my $grace = "/usr/bin/grace";
my $work_dir = "/tmp/stable-stats/";
my $data_dir = "/data/";
my $phase_template = "/var/local/grace-templates/phase.template";
my $sigma_template = "/var/local/grace-templates/sigma.template";
my $freq_template = "/var/local/grace-templates/freq.template";

# this is relative to web directory
my $image_dir = "images/";

###############################################################################

my $opt;
my $data_file;
my $df;
my $html_file;
my $phase_file;
my $flat_phase_file;
my $freq_file;
my $sigma_file;
my $comment_file;
my $adev_table_file;
my $i;
my $counter;
my $sum;
my $reading;
my @tmpxarray;
my @tmpyarray;
my @xarray;
my @yarray;
my @flat_array;
my @freq_xarray;
my @freq_yarray;
my $y_axis_label;
my $days;
my $seconds;
my $total_offset;
my $linear_offset;
my $offset_per_sample;
my $freq_offset;
my $aging_day;
my $aging_week;
my $aging_month;
my $aging_year;
my $sample_size;
my $new_index;
my $ls;
my $ls_freq;
my $intercept;
my $slope;
my $R_squared;
my $freq_slope;
my $freq_intercept;
my $freq_R_squared;
my $xmin;
my $xmax;
my $ymin;
my $ymax;
my $range;
my $flat_range;
my $avY;
my $avX;
my $scale;
my @xsigma;
my @ysigma;
my @xsigma_all;
my @ysigma_all;
my $adev;
my $time_diff;
my $tau;

#----------
# usage and option initialization
my $opt_string = 'hn:a:t:b:x:y:g:s:f:';
sub usage() {
print STDERR << "EOF";

usage: $0 [-h] [-n top|center|bottom] [-a samples ] [-t <seconds>]
	[-b <basename>] [-x pixels -y pixels] [-g title] [-s subtitle]
	-f filename

This program reads a data file (or STDIN) with lines containing MJD/phase
data pairs.  If no basename is specified, it will print useful statistics
to STDOUT and terminate.

If a basename is specified, it will silently create
basename.html, containing tables with those statistics; basename-phase.png,
a PNG graph of the phase data; and basename-sigma.png, a PNG graph of the
Allan Deviation data.

If a file basename-comments.txt exists in the working
directory, its contents will be included at the top of the HTML page.  The
comments file can contain most valid HTML tags.

-h	: this (help) message

-a	: average datafile by n samples

-n	: normalize phase data, setting top|center|bottom as zero

-t 	: tau0 (sample interval); if not specified, calculated from MJD
	  If averaging, specify tau0 after average taken

-b	: basename for output files; they'll be put in $work_dir

-g	: optional graph title (1 line; applied to both -phase and -sigma)

-s	: optional graph subtitle (1 line; applied to both -phase and -sigma)
		-- will have \"Phase via\" or \"Overlapping Allan Deviation via\" prepended

-x,-y	: output dimension in pixels; default 550x450 

-f	: data file name; if "-" read from STDIN

EOF
}

#----------------------
getopts( "$opt_string", \my %opt ) or usage() and exit;
usage() and exit if $opt{h};
usage() and exit if !$opt{f};

my $normalize = 0;
if ($opt{n}) {
	$normalize = lower_case($opt{n});
	}

my $average = 0;
if ($opt{a}) {
	$average = $opt{a};
	}

my $tau0 = 0;
if ($opt{t}) {
	$tau0 = $opt{t};
	}

my $basename = "";
if ($opt{b}) {
	$basename = $opt{b};
	$html_file = $basename . ".html";
	$phase_file = $basename . "-phase.png";
	$flat_phase_file = $basename . "-flat_phase.png";
	$sigma_file = $basename . "-sigma.png";
	$freq_file = $basename . "-freq.png" ;
	$adev_table_file = $basename . "-adev.dat";
	$comment_file = $basename . "-comment.txt";
	}

my $title = "";
if ($opt{g}) {
	$title = $opt{g};
	}

my $subtitle = "";
if ($opt{s}) {
	$subtitle = $opt{s};
	}

my $x = 550;
if ($opt{x}) {
	$x = $opt{x};
	}

my $y = 450;
if ($opt{y}) {
	$y = $opt{y};
	}

if ($opt{f}) {
	$data_file = $opt{f};
	}
else {
	die "No data file specified!\n";
	}

#--------

# pretty printer for time values
sub timeprint() {
	my ($time) = @_;
	my $scale = 1;
	my $label = "s";
	if (abs($time) < 1e0) {$scale = 1e3;$label = "ms";}
	if (abs($time) < 1e-3) {$scale = 1e6;$label = "us";}
	if (abs($time) < 1e-6) {$scale = 1e9;$label = "ns";}
	if (abs($time) < 1e-9) {$scale = 1e12;$label = "ps";}

	$time = $time * $scale;
	my $result = sprintf("%.3f %s",$time,$label);
	return $result;
}

# ADEV calculation
sub adev_calc() {
	my ($tau0,$sample_size,$avg,@yarray) = @_;
	my $sum = 0.0;
	my $i;
	my $y;
	my $last_y;
	my $dy;
	my $samples = 1;
	my $adev;

	# requiring more than bare minimum of points to avoid whipsaw
	if ($sample_size > 3) {
		for ($i=$avg; $i<$sample_size; $i = $i + $avg) {
			$y = $yarray[$i] - $yarray[$i-$avg];
			if ($i > $avg) {
				$dy = $y - $last_y;
				$sum += $dy * $dy;
				}
			$last_y = $y;
			$samples++;
			}
	}
	return $adev = (sqrt($sum / (2.0 * ($samples - 1))))/($tau0*$avg);
}

sub oadev_calc() {

# Overlapping Allan Deviation
# Sigma^2(Tau) = 
# 1 / (2*(N-2*m)*Tau^2) * Sum(X[i+2*m]-2*X[i+m]+X[i], i=1, i=N-2*m)
# Tau is the averaging time, N is the total number of points,
# and Tau = m*Tau0, where Tau0 is the basic measurement interval       

my ($tau0,$N,$m,@phase) = @_;

my $sigma = 0;
my $numGaps = 0;
my $sum;
my $i;

my $tau = $tau0 * $m;

for($i = 0; $i < ($N-2*$m); $i++) {
	$sum = 0;
	if(($phase[$i+2*$m]==0 ||  $phase[$i+$m]==0 || $phase[$i]==0) 
		&& $i!=0 && $i!=($N-2*$m-1)) {
			$numGaps++;
		} else {
			$sum = $phase[$i+2*$m] - 2*$phase[$i+$m] + $phase[$i];
			$sigma += $sum * $sum;
		}
}
              
$sigma = $sigma / (2*($N-$numGaps-2*$m)*$tau*$tau);
$sigma = sqrt($sigma);
return $sigma;
}


sub last_mod() {
	my ($filename) = @_;
	my @Stats = stat($filename);

	my $ModTime = $Stats[9];
	my @Weekdays = ('Sunday', 'Monday', 'Tuesday', 'Wednesday',
			 'Thursday', 'Friday', 'Saturday');
	my @Months = ('Jan', 'Feb', 'Mar', 'Apr', 'May',
		 'Jun', 'Jul', 'Aug', 'Sep', 'Oct','Nov', 'Dec');
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)
		= gmtime($ModTime);
	my $Weekday = $Weekdays[$wday];
	my $Month = $Months[$mon];
	$year += 1900;
	my $last_mod = 
		sprintf("%02u:%02u:%02u %02u-%s-%4u (UTC)",
		$hour,$min,$sec,$mday,$Month,$year);
	return $last_mod;
}

# screen output routine
sub results() {
	my $missing = ($seconds/$tau0)-$sample_size;
	if ($missing >= 1) {
		$missing = sprintf("\(%u missing readings\)",$missing);
		}
	else {
		$missing = "";
		}
		
	print "\n";
	print "Data file: ",$data_file,"\n";
	if ($data_file ne "-") {
		print "Last modified: ",&last_mod($data_file),"\n";
		}
	print "\n";
	printf "%u readings over %.2f days, tau0 = %u seconds %s\n",
		$sample_size,$days,$tau0,$missing;
	print "MJD date range from ",$xarray[0]," to ",$xarray[$#xarray],"\n";
	print "\n";
	print "End-to-end phase change:\t",&timeprint($total_offset),"\n";
	print "Phase change per day:\t\t",&timeprint($total_offset/$days),"\n";
	print "Peak-to-peak phase range:\t",&timeprint(abs($ymin-$ymax)),"\n";
	print "\n";
	printf "Frequency offset (linear, R2=%3.3F):\t%4.3E\n",
		$R_squared,$linear_offset;
	printf "Frequency offset (endpoints):\t%4.3E\n",$total_offset/$seconds;
	print "\n";
	print "Time difference (REF-DUT, last sample): ",&timeprint($time_diff),"\n";
	print "\n";
	print "Tau\t\tADEV\n";
	for ($i=0;$i<=$#xsigma;$i++) {
		printf "%3.2E\t%4.3E\n",$xsigma[$i],$ysigma[$i];
		}	
}

# html output routine
sub html_results {
	my $outfile;
    	$outfile = IO::File->new($html_file, "w");

	my $missing = ($seconds/$tau0)-$sample_size;
	if ($missing >= 1) {
		$missing = sprintf("\(%u missing readings\)",$missing);
		}
	else {
		$missing = "";
		}

	# NOTE -- HEADER for www.febo.com; the next three lines won't be
	# useful for anyone else!!!
	print $outfile "<HTML>\n<HEAD>\n<TITLE>$title</TITLE>\n</HEAD>\n";
	print $outfile "<!--\#INCLUDE VIRTUAL=\"/system/header.html\" -->\n";
	print $outfile "<H1>$title</H1>\n";

	# read comment file if it exists
	my $comment;
	$comment = new IO::File $comment_file, "r";
	if (defined $comment) {
		print $outfile "</CENTER>\n<BLOCKQUOTE>\n";
		while ($reading = <$comment>) {
			print $outfile $reading;
			}
	print $outfile "</BLOCKQUOTE>\n<P>\n<HR WIDTH=75%\n<CENTER>\n<P>\n";
	}

	print $outfile "<CENTER>\n";
	print $outfile "Data file: <B>",$data_file,"</B>\n";
	if ($data_file ne "-") {
		print $outfile "<BR>\n";
		print $outfile "Last modified: <B>",
			&last_mod($data_file),"</B>\n";
		}
	print $outfile "<P>\n";
	print $outfile "MJD date range from ",$xarray[0],
		" to ",$xarray[$#xarray],"\n";

	print $outfile "<BR>\n";
	printf $outfile
		"%u readings over %.2f days, tau0 = %u seconds %s\n",
		$sample_size,$days,$tau0,$missing;

	print $outfile "<P>\n";
	print $outfile "<TABLE WIDTH = 60%>\n";
	print $outfile "<TR><TD>End-to-end phase change:</TD><TD ALIGN=RIGHT>",
		&timeprint($total_offset),"</TD></TR>\n";
	print $outfile "<TR><TD>Phase change per day:</TD><TD ALIGN=RIGHT>",
		&timeprint($total_offset/$days),"</TD></TR>\n";
	print $outfile "<TR><TD>Peak-to-peak phase range:</TD><TD ALIGN=RIGHT>",
		&timeprint(abs($range)),"</TD></TR>\n";
	print $outfile "</TABLE>\n";
	print $outfile "<P>\n";
	print $outfile "<TABLE WIDTH = 60%>\n";
	print $outfile "<TR><TD>Time difference (REF-DUT, last sample):</TD>",
		"<TD ALIGN=RIGHT>",&timeprint($time_diff),"</TD></TR>\n";
	print $outfile "</TABLE>\n";
	print $outfile "<P>\n";

	print $outfile "<TABLE WIDTH = 60%>\n";
	printf $outfile "<TR><TD>Frequency offset (linear, R<sup>2</sup>=%3.3F):</TD><TD ALIGN=RIGHT>%4.2E</TD></TR>\n", $R_squared,$linear_offset;
	printf $outfile "<TR><TD>Frequency offset (endpoints):</TD><TD ALIGN=RIGHT>%4.2E</TD></TR>\n",$total_offset/$seconds;
	print $outfile "</TABLE>\n";
	print $outfile "<P>\n";

	print $outfile "<TABLE WIDTH = 60%>\n";
	printf $outfile "<TR><TD>Drift (linear, R<sup>2</sup>=%3.3F):</TD><TD ALIGN=RIGHT>%4.2E/day</TD></TR>\n", $freq_R_squared, $aging_day;
	print $outfile "</TABLE>\n";
	print $outfile "<TABLE WIDTH = 60%>\n";
	printf $outfile "<TR><TD ALIGN=LEFT>(%4.2E/week, %4.2E/month, %4.2E/year)</TD></TR>\n", $aging_week, $aging_month, $aging_year;
	print $outfile "</TABLE>\n";
	print $outfile "<P>\n";

	print $outfile "<P>\n<HR WIDTH=75%\n<P>\n";
	print $outfile "<P>\n<H4>Raw Phase</H4>\n<P>\n";
	print $outfile "<IMG SRC=$image_dir$phase_file>\n";
	print $outfile "<P><HR><P>\n";
	print $outfile "<P>\n<H4>Offset Removed</H4>\n<P>\n";
	print $outfile "<TABLE WIDTH = 60%>\n";
	print $outfile "<TR><TD>Peak-to-peak noise:</TD><TD ALIGN=RIGHT>",
		&timeprint(abs($flat_range)),"</TD></TR>\n";
	print $outfile "</TABLE>\n";
	print $outfile "<P>\n";
	print $outfile "<IMG SRC=$image_dir$flat_phase_file>\n";
	print $outfile "<P><HR><P>\n";
	print $outfile "<P>\n<H4>Frequency</H4>\n<P>\n";
	print $outfile "<IMG SRC=$image_dir$freq_file>\n";
	print $outfile "<P><HR><P>\n";
 	print $outfile "<H4>Overlapping Allan Deviation</H4><P>\n";
	print $outfile "<TABLE WIDTH=60% BORDER=1>\n";
	print $outfile "<TR><TD ALIGN=CENTER>Tau</TD><TD ALIGN=CENTER>ADEV</TD>
		<TD>&nbsp</TD><TD ALIGN=CENTER>Tau</TD><TD ALIGN=CENTER>ADEV</TD></TR>\n";

	for ($i=0;$i<=($#xsigma/2);$i++) {
		if ($xsigma[$i]) {
			printf $outfile
		 	"<TR><TD ALIGN=CENTER>%3.2E</TD><TD ALIGN=CENTER>%4.3E</TD>",
			$xsigma[$i],
			$ysigma[$i];
			print $outfile "<TD>&nbsp&nbsp</TD>";
			}
		if ($xsigma[$i+round(0,($#xsigma+1)/2)]) {
			printf $outfile
			"<TD ALIGN=CENTER>%3.2E</TD><TD ALIGN=CENTER>%4.3E</TD></TR>\n",
			$xsigma[$i+round(0,($#xsigma+1)/2)],
			$ysigma[$i+round(0,($#ysigma+1)/2)];
			}
		}	
	print $outfile "</TABLE>\n";
	print $outfile "<P>\n";
	print $outfile "<IMG SRC=$image_dir$sigma_file>\n";
	print $outfile "<P>\n";
	print $outfile "This page was generated by <A HREF='http://www.febo.com/pages/tools/stable-stats/stable-stats.html'><TT>stable-stats</TT></A>\n";
	print $outfile "<P>\n";
	print $outfile "</CENTER>\n";
	# NOTE:  FOOTER for www.febo.com -- it won't be useful
	# for anyone else!!!
	print $outfile "<!--\#INCLUDE VIRTUAL=\"/system/footer.html\" -->\n";
}

sub phase_graph {
	my ($pngfile,$xarray,$yarray) = @_;
	my $reading;
	my $tmpfile;
	my $infile;
	my $outfile;
	my $workstring;
	my $i;

	if ($subtitle) {
		$workstring = "Phase via " . $subtitle . " (tau0 = $tau0 seconds)";
		}
	else {
		$workstring = " ";
		}
	 
	do { $tmpfile = tmpnam() }
    		until $outfile = IO::File->new($tmpfile, O_RDWR|O_CREAT|O_EXCL);
	
	$infile = new IO::File $phase_template, "r";
	if (defined $infile) {
		while ($reading = <$infile>) {
			$reading =~ s/##TITLE##/$title/;
			$reading =~ s/##SUBTITLE##/$workstring/;
			$reading =~ s/##XAXIS LABEL##/MJD/;
			$reading =~ s/##YAXIS LABEL##/$y_axis_label/;
			print $outfile $reading;
			}
		for ($i=0;$i<@$xarray;$i++) {
			print $outfile $xarray->[$i]," ",$yarray->[$i],"\n";
			}
		print $outfile "&";
	 }
	
	my @args = ($grace,"-hardcopy","-hdevice","PNG","-printfile",
		$pngfile,$tmpfile,"-pexec","autoscale");
	system(@args) == 0
             or die "system @args failed: $?";

	unlink($tmpfile) or die "Couldn't unlink $tmpfile!";
}

sub freq_graph {
	my ($pngfile,$xarray,$yarray) = @_;
	my $reading;
	my $tmpfile;
	my $infile;
	my $outfile;
	my $workstring;
	my $freq_y_axis_label = "Frequency";
	my $i;

	if ($subtitle) {
		$workstring = "Frequency via " . $subtitle . " (tau0 = $tau0 seconds)";
		}
	else {
		$workstring = " ";
		}
	 
	do { $tmpfile = tmpnam() }
    		until $outfile = IO::File->new($tmpfile, O_RDWR|O_CREAT|O_EXCL);
	
	$infile = new IO::File $freq_template, "r";
	if (defined $infile) {
		while ($reading = <$infile>) {
			$reading =~ s/##TITLE##/$title/;
			$reading =~ s/##SUBTITLE##/$workstring/;
			$reading =~ s/##XAXIS LABEL##/MJD/;
			$reading =~ s/##YAXIS LABEL##/$freq_y_axis_label/;
			print $outfile $reading;
			}

		for ($i=1;$i<@$xarray;$i++) {
			print $outfile $xarray->[$i]," ",$yarray->[$i],"\n";
			}
		print $outfile "&";
	 }
	
	my @args = ($grace,"-hardcopy","-hdevice","PNG","-printfile",
		$pngfile,$tmpfile,"-pexec","autoscale");
	system(@args) == 0
             or die "system @args failed: $?";

	unlink($tmpfile) or die "Couldn't unlink $tmpfile!";
}

sub sigma_graph {
	my $pngfile = shift;
	my $reading;
	my $tmpfile;
	my $infile;
	my $outfile;
	my $workstring;
	my $i;

	if ($subtitle) {
		$workstring = "Overlapping ADEV (many tau) via ".$subtitle;
		}
	else {
		$workstring = " ";
		}
	 
	do { $tmpfile = tmpnam() }
    		until $outfile = IO::File->new($tmpfile, O_RDWR|O_CREAT|O_EXCL);

	$infile = new IO::File $sigma_template, "r";
	if (defined $infile) {
		while ($reading = <$infile>) {
			$reading =~ s/##TITLE##/$title/;
			$reading =~ s/##SUBTITLE##/$workstring/;
			$reading =~ s/##XAXIS LABEL##/Tau/;
			$reading =~ s/##YAXIS LABEL##/Allan Deviation/;
			print $outfile $reading;
			}
		for ($i=0;$i<=$#xsigma_all;$i++) {
			print $outfile $xsigma_all[$i]," ",$ysigma_all[$i],"\n";
			}
		print $outfile "&";
	 }
	
	my @args = ($grace,"-hardcopy","-hdevice","PNG","-printfile",
		$pngfile,$tmpfile,"-pexec","autoscale");
	system(@args) == 0
             or die "system @args failed: $?";

	unlink($tmpfile) or die "Couldn't unlink $tmpfile!";
}

###############################################################################
# main program
# #############################################################################

# Add data points
my $foo;
my $bar;
my $baz;
$df = new IO::File $data_file, "r";
if (defined $df) {
	while ($reading = <$df>) {
		if (substr($reading,0,1) =~ "[0-9]") {
			($foo,$bar,$baz) = split(/\s/,$reading);
			push(@tmpxarray,$foo);
			push(@tmpyarray,$bar);
		}
	}
} else {
	exit;
	}

# make sure there's enough data to make it all worthwhile
if ($#tmpxarray <10) {
	print "$basename: Not enough data yet!\n";
	exit 0;
	}

# if we are averaging data
if ($average > 1) {
	$counter = 0;
	$sum = 0;
	$new_index = 0;
	for ($i=0; $i<=$#tmpxarray; $i++) {
                $counter++;
                $sum+=$tmpyarray[$i];
                if ($counter == $average) {
			$xarray[$new_index] = $tmpxarray[$i];
			$yarray[$new_index] = $sum/$average;
			$new_index++;
                        $counter = 0;
                        $sum = 0;
                }
        }
}
else {
	@xarray = @tmpxarray;
	@yarray = @tmpyarray;
	}

# calculate tau0 if necessary -- average from first ten MJD values
# trying to average over full file instead
if ($tau0 == 0) {
	$tau0 = (($xarray[$#xarray]-$xarray[0])/($#xarray-1))* 86400;
	if ($tau0 >= 500) {
		# long taus tend to round up, so truncate instead
		$tau0 = int($tau0);
		}
	if (($tau0 < 500) && ($tau0 >= 10)) {
		$tau0 = round(0,$tau0);
		}
	if ($tau0 < 10) {
		$tau0 = round(1,$tau0);
		}
}

# create frequency from phase
$new_index = 0;
for ($i=1; $i<=$#xarray; $i++) {
	$freq_xarray[$new_index] = $xarray[$i];
	$freq_yarray[$new_index] = ($yarray[$i] - $yarray[$i-1])/$tau0;
	$new_index++;
}

# set up linear regression for phase	
$ls = Statistics::OLS->new(); 
$ls->setData (\@xarray, \@yarray);
$ls->regress();

# set up linear regression for freq
$ls_freq = Statistics::OLS->new(); 
$ls_freq->setData (\@freq_xarray, \@freq_yarray);
$ls_freq->regress();

# get useful information
($intercept, $slope) = $ls->coefficients();
$R_squared = $ls->rsq();
($freq_intercept, $freq_slope) = $ls_freq->coefficients();
$freq_R_squared = $ls_freq->rsq();

# For some reason, the returned slope seems to be per day, not per sample.
# I *think this is because of the MJD numbering of the X axis.  So, this 
# could potentially break with differently formatted data.
$aging_day = $freq_slope;
$aging_week = $freq_slope*7;
$aging_month = $freq_slope*30;
$aging_year = $freq_slope*365 ;

$sample_size = $ls->size();    
($xmin, $xmax, $ymin, $ymax) = $ls->minMax();
($avX, $avY) = $ls->av();
$range = $ymax-$ymin;

# do some calculations
$total_offset = $yarray[$sample_size-1]-$yarray[0];
$time_diff = $yarray[$sample_size-1];
$days = $xmax-$xmin;
$seconds = $days*86400;

# For some reason, the returned slope seems to be per day, not per sample.
# I *think this is because of the MJD numbering of the X axis.  So, this 
# could potentially break with differently formatted data.
$linear_offset = $slope/86400;
$offset_per_sample = $linear_offset*$tau0;

# create array with offset removed
for ($i=0; $i<$sample_size; $i++) {
	$flat_array[$i] = $yarray[$i]-($offset_per_sample*$i);
}

# set up linear regression on flat array	
$ls = Statistics::OLS->new(); 
$ls->setData (\@xarray, \@flat_array);
$ls->regress();
# get useful information
my ($flat_xmin, $flat_xmax, $flat_ymin, $flat_ymax) = $ls->minMax();
my ($flat_avX, $flat_avY) = $ls->av();
$flat_range = $flat_ymax - $flat_ymin;

# Do ADEV calculation

for ($i=1; $i<($sample_size/2); $i = $i*2) {
	$adev = &oadev_calc($tau0,$sample_size,$i,@flat_array);
	$tau = $tau0*$i;
	push(@xsigma,$tau);
	push(@ysigma,$adev);
	}	

# do adev for all tau
open ADEVTABLEFILE, "> $adev_table_file" or die "Can't open $adev_table_file$!";
for ($i=1; $i<($sample_size/2); $i++) {
	$adev = &oadev_calc($tau0,$sample_size,$i,@flat_array);
	$tau = $tau0*$i;
	push(@xsigma_all,$tau);
	push(@ysigma_all,$adev);
	printf ADEVTABLEFILE "%u %1.3e\n",$tau,$adev;
	# make it "many tau" by trimming large values
#	if ($i > 512) { $i++ };
#	if ($i > 256) { $i++ };
#	if ($i > 128) { $i++ };
	if ($i > 48) { $i += 1 };
# for >32, previously was $i+=2; don't remember why
	if ($i > 32) { $i += 1 };
	if ($i > 16) { $i += 1 };
	}	
close ADEVTABLEFILE;

# get magnitude of y axis
$scale = 1;
$y_axis_label = "Seconds";
if ($range < 1e0) {$scale = 1e3;$y_axis_label = "Milliseconds"}
if ($range < 1e-3) {$scale = 1e6;$y_axis_label = "Microseconds"}
if ($range < 1e-6) {$scale = 1e9;$y_axis_label = "Nanoseconds"}
if ($range < 1e-9) {$scale = 1e12;$y_axis_label = "Picoseconds"}

# get magnitude of flat y axis
my $flat_scale = 1;
my $flat_y_axis_label = "Seconds";
if ($flat_range < 1e0) {$flat_scale = 1e3;$flat_y_axis_label = "Milliseconds"}
if ($flat_range < 1e-3) {$flat_scale = 1e6;$flat_y_axis_label = "Microseconds"}
if ($flat_range < 1e-6) {$flat_scale = 1e9;$flat_y_axis_label = "Nanoseconds"}
if ($flat_range < 1e-9) {$flat_scale = 1e12;$flat_y_axis_label = "Picoseconds"}

# do some calculations

# scale y values
for ($i=0; $i<$sample_size; $i++) {
	$yarray[$i] = ($yarray[$i] * $scale);
	$flat_array[$i] = ($flat_array[$i] * $flat_scale);
	}

# normalize y values
if ($normalize eq "top") {
	for ($i=0; $i<$sample_size; $i++) {
		$yarray[$i] -= $ymax * $scale;
		$flat_array[$i] -= $flat_ymax * $flat_scale;
	}
}	
if ($normalize eq "center") {
	for ($i=0; $i<$sample_size; $i++) {
		$yarray[$i] -= ($ymin + ($range/2))*$scale;
		$flat_array[$i] -= ($flat_ymin + ($flat_range/2))*$flat_scale;
		
	}
}
if ($normalize eq "bottom") {
	for ($i=0; $i<$sample_size; $i++) {
		$yarray[$i] -= $ymin * $scale;
		$flat_array[$i] -= $flat_ymin * $flat_scale;
	}
}


if ($basename) {
	&html_results();
	&phase_graph($phase_file,\@xarray,\@yarray);
	&freq_graph($freq_file,\@freq_xarray,\@freq_yarray);
	&sigma_graph($sigma_file);
	$title .= " (Offset Removed)";
	$y_axis_label = $flat_y_axis_label;
	&phase_graph($flat_phase_file,\@xarray,\@flat_array);
	}
else {
	&results();
	}


exit 0;
