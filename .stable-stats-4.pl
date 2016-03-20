#!/usr/bin/perl -w

# stable-stat.pl
# version 0.23 -- generate statistics from file containing MJD and phase data
# 12 January 2005
#
# Copyright 2005 by John R. Ackermann  N8UR (jra@febo.com)
# Licensed under the GPL version 2 or later; see the file COPYING
# included with this distribution.  I request, but do not require, that
# any modifications that correct bugs or errors, or increase the program's
# functionality, be sent via email to the author at the address above.

use strict;
use POSIX qw(tmpnam);
use IO::File;
use Getopt::Std;
use Statistics::OLS;
use n8ur qw(round lower_case);

###############################################################################
# local paths
###############################################################################

my $grace = "/usr/bin/grace";
my $work_dir = "./";
my $phase_template = "./phase.template";
my $sigma_template = "./sigma.template";

###############################################################################

my $opt;
my $data_file;
my $df;
my $html_file;
my $phase_file;
my $sigma_file;
my $i;
my $reading;
my @xarray;
my @yarray;
my $y_axis_label;
my $days;
my $seconds;
my $total_offset;
my $linear_offset;
my $sample_size;
my $ls;
my $intercept;
my $slope;
my $R_squared;
my $xmin;
my $xmax;
my $ymin;
my $ymax;
my $range;
my $avY;
my $avX;
my $scale;
my @xsigma;
my @ysigma;
my $avar;

#----------
# usage and option initialization
my $opt_string = 'hn:t:b:x:y:g:s:f:';
sub usage() {
print STDERR << "EOF";

usage: $0 [-h] [-n top|center|bottom] [-t <seconds>] [-b <basename>]
	  [-x pixels -y pixels] [-g title] [-s subtitle]  -f filename

This program reads a data file (or STDIN) with lines containing MJD/phase
data pairs.  If no basename is specified, it will print useful statistics
to STDOUT and terminate.  If a basename is specified, it will silently create
basename.html, containing tables with those statistics; basename-phase.png,
a PNG graph of the phase data; and basename-sigma.png, a PNG graph of the
Allan Deviation data.

-h	: this (help) message

-n	: normalize phase data, setting top|center|bottom as zero

-t 	: tau (sample interval); if not specified, calculated from MJD

-b	: basename for output files; they'll be put in $work_dir

-g	: optional graph title (1 line; applied to both -phase and -sigma)

-s	: optional graph subtitle (1 line; applied to both -phase and -sigma)
		-- will have \"Phase via\" or \"Allan Deviation via\" prepended

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

my $tau = 0;
if ($opt{t}) {
	$tau = $opt{t};
	}

my $basename = "";
if ($opt{b}) {
	$basename = $opt{b};
	$html_file = $basename . ".html";
	$phase_file = $basename . "-phase.png";
	$sigma_file = $basename . "-sigma.png";
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
sub adev() {
	my ($tau,$sample_size,$avg,@yarray) = @_;
	my $sum = 0.0;
	my $i;
	my $y;
	my $last_y;
	my $dy;
	my $samples = 1;
	my $adev;

	if ($sample_size > 2) {
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
	return $adev = (sqrt($sum / (2.0 * ($samples - 1))))/($tau*$avg);
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
	my $missing = ($seconds/$tau)-$sample_size;
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
	printf "%u readings over %.2f days, tau = %u seconds %s\n",
		$sample_size,$days,$tau,$missing;
	print "MJD date range from ",$xarray[0]," to ",$xarray[$#xarray],"\n";
	print "\n";
	print "End-to-End drift:\t",&timeprint($total_offset),"\n";
	print "Drift per day:\t\t",&timeprint($total_offset/$days),"\n";
	print "Peak-to-Peak drift:\t",&timeprint(abs($ymin-$ymax)),"\n";
	print "\n";
	printf "Frequency Offset (linear, R2=%3.3F):\t%4.3E\n",
		$R_squared,$linear_offset;
	printf "Frequency Offset (endpoints):\t%4.3E\n",$total_offset/$seconds;
	print "\n";

	print "Tau\t\tADEV\n";
	for ($i=0;$i<=$#xsigma;$i++) {
		printf "%3.2E\t%4.3E\n",$xsigma[$i],$ysigma[$i];
		}	
}

# html output routine
sub html_results() {
	my $outfile;
    	$outfile = IO::File->new($html_file, "w");

	my $missing = ($seconds/$tau)-$sample_size;
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
		"%u readings over %.2f days, tau = %u seconds %s\n",
		$sample_size,$days,$tau,$missing;

	print $outfile "<P>\n";
	print $outfile "<TABLE WIDTH = 60%>\n";
	print $outfile "<TR><TD>End-to-End drift:</TD><TD ALIGN=RIGHT>",
		&timeprint($total_offset),"</TD></TR>\n";
	print $outfile "<TR><TD>Drift per day:</TD><TD ALIGN=RIGHT>",
		&timeprint($total_offset/$days),"</TD></TR>\n";
	print $outfile "<TR><TD>Peak-to-Peak drift:</TD><TD ALIGN=RIGHT>",
		&timeprint(abs($ymin-$ymax)),"</TD></TR>\n";
	print $outfile "</TABLE>\n";
	print $outfile "<P>\n";
	print$outfile "<TABLE WIDTH = 60%>\n";
	printf $outfile "<TR><TD>Frequency Offset (linear, R<sup>2</sup>=%3.3F):</TD><TD ALIGN=RIGHT>%4.3E</TD></TR>\n", $R_squared,$linear_offset;
	printf $outfile "<TR><TD>Frequency Offset (endpoints):</TD><TD ALIGN=RIGHT>%4.3E</TD></TR>\n",$total_offset/$seconds;
	print $outfile "</TABLE>\n";
	print $outfile "<P>\n";
	print $outfile "<IMG SRC=$phase_file>\n";
	print $outfile "<P>\n";

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
	print $outfile "<IMG SRC=$sigma_file>\n";
	print $outfile "<P>\n";
	print $outfile "This page was generated by <A HREF='http://www.febo.com/time-freq/tools/stable-stats.html'><TT>stable-stats</TT></A>\n";
	print $outfile "<P>\n";
	print $outfile "</CENTER>\n";
	# NOTE:  FOOTER for www.febo.com -- it won't be useful
	# for anyone else!!!
	print $outfile "<!--\#INCLUDE VIRTUAL=\"/system/footer.html\" -->\n";
}

sub phase_graph() {
	my $reading;
	my $tmpfile;
	my $infile;
	my $outfile;
	my $workstring;
	my $i;

	if ($subtitle) {
		$workstring = "Phase via " . $subtitle;
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
		for ($i=0;$i<=$#xarray;$i++) {
			print $outfile $xarray[$i]," ",$yarray[$i],"\n";
			}
		print $outfile "&";
	 }
	
	my @args = ($grace,"-hardcopy","-hdevice","PNG","-printfile",
		$phase_file,$tmpfile,"-pexec","autoscale");
	system(@args) == 0
             or die "system @args failed: $?";

	unlink($tmpfile) or die "Couldn't unlink $tmpfile!";
}

sub sigma_graph() {
	my $reading;
	my $tmpfile;
	my $infile;
	my $outfile;
	my $workstring;
	my $i;

	if ($subtitle) {
		$workstring = "Allan Deviation via " . $subtitle;
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
		for ($i=0;$i<=$#xsigma;$i++) {
			print $outfile $xsigma[$i]," ",$ysigma[$i],"\n";
			}
		print $outfile "&";
	 }
	
	my @args = ($grace,"-hardcopy","-hdevice","PNG","-printfile",
		$sigma_file,$tmpfile,"-pexec","autoscale");
	system(@args) == 0
             or die "system @args failed: $?";

	unlink($tmpfile) or die "Couldn't unlink $tmpfile!";
}

###############################################################################
# main program
# #############################################################################

# Add data points
$df = new IO::File $data_file, "r";
if (defined $df) {
	while ($reading = <$df>) {
		if (substr($reading,0,1) =~ "[0-9]") {
			($x,$y) = split(/\s/,$reading);
			push(@xarray,$x);
			push(@yarray,$y);
		}
	}
} else {
	exit;
	}

# make sure there's enough data to make it all worthwhile
if ($#xarray <10) {die "Not enough data yet!\n"};

# calculate tau if necessary -- average from first ten MJD values
if ($tau == 0) {
	$tau = ($xarray[10]-$xarray[0])*8640;
	if ($tau < 10) {
		$tau = round(1,$tau);
	} else {
		$tau = round(0,$tau);
		}
}

# set up linear regression	
$ls = Statistics::OLS->new(); 
$ls->setData (\@xarray, \@yarray);
$ls->regress();

# get useful information
($intercept, $slope) = $ls->coefficients();
$R_squared = $ls->rsq();
$sample_size = $ls->size();    
($xmin, $xmax, $ymin, $ymax) = $ls->minMax();
($avX, $avY) = $ls->av();
$range = $ymax-$ymin;

# do some calculations
$total_offset = $yarray[$sample_size-1]-$yarray[0];
$days = $xmax-$xmin;
$seconds = $days*86400;

# For some reason, the returned slope seems to be per day, not per sample.
# This could potentially break with differently formatted data.
$linear_offset = $slope/86400;

# Do AVAR calculation
for ($i=1; $i<($sample_size/3); $i = $i*2) {
	$avar = &adev($tau,$sample_size,$i,@yarray);
	push(@xsigma,$tau*$i);
	push(@ysigma,$avar);
	}	

# normalize y values
if ($normalize eq "top") {
	for ($i=0; $i<$sample_size; $i++) {
		$yarray[$i] -= $ymax;
	}
}	
if ($normalize eq "center") {
	for ($i=0; $i<$sample_size; $i++) {
		$yarray[$i] -= $ymin + ($range/2);
	}
}
if ($normalize eq "bottom") {
	for ($i=0; $i<$sample_size; $i++) {
		$yarray[$i] -= $ymin;
	}
}

# get magnitude of y axis
$scale = 1;
$y_axis_label = "Seconds";
if ($range < 1e0) {$scale = 1e3;$y_axis_label = "Milliseconds"}
if ($range < 1e-3) {$scale = 1e6;$y_axis_label = "Microseconds"}
if ($range < 1e-6) {$scale = 1e9;$y_axis_label = "Nanoseconds"}
if ($range < 1e-9) {$scale = 1e12;$y_axis_label = "Picoseconds"}

# scale y values
for ($i=0; $i<$sample_size; $i++) {
	$yarray[$i] = ($yarray[$i] * $scale);
	}

if ($basename) {
	&html_results();
	&phase_graph();
	&sigma_graph();
	}
else {
	&results();
	}


exit 0;
