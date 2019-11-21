#!/usr/bin/perl

use warnings;
use strict;
use DateTime;

open my $stream, '<', "1912794.csv" or die;

my @new;

my $i = 0;
my $j = 1;
push @new, "HourlyWindDirection,HourlyWindSpeed,PrevHourlyWindDirection,PrevHourlyWindSpeed,TimeOfDay,TimeOfYear\n";
my $prev = "*,*";
my $start = DateTime->new(
	year	=> 2010,
	month	=> 1,
	day 	=> 1,
	);
foreach my $line (<$stream>) {
	if ($j == 1) {
		$j++;
		next;
	}
	chomp $line;
	$line =~ s/"//g;
	my @row = split(',',$line);
	my $date = $row[1];
	# print $date,"\n";
	$date =~ /(....)-(..)-(..)T(..):(..)/;
	my $year = $1;
	my $month = $2;
	my $day = $3;
	my $hour = $4;
	my $minute = $5;
	
	my $dt = DateTime->new(
		year		=> $year,
		month		=> $month,
		day			=> $day,
		hour 		=> $hour,
		minute		=> $minute,
		);
	my $total_seconds = ($dt->epoch() - $start->epoch())/(365*24*60*60);
	# print $total_seconds, "\n";
	if (!length $hour) {
		print "bing\n";
		print $line, "\n";
		next;
	}
	my $time_of_day = ($hour * 60 + $minute)/(60*24);
	# print $time,"\n";
	push @new, "$row[54],$row[56],$prev,$time_of_day,$total_seconds\n";
	$prev = "$row[54],$row[56]";
	# $i++;
	# if ($i == 10) {
		# last;
	# }
}
close $stream;

open my $write, '>', "new_data_1.csv" or die;
print $write @new;
close $write;