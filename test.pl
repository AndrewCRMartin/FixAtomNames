#!/usr/bin/perl
use strict;

my @angles2 = ( 0, 30,  60,  90, 120, 150, 180, 360, -30, -60, -60, -90, -120, -150, -180, 59.248);
my @angles1 = (60, 90, 120, 150, 180, 210, 240,  60,  30,   0, 360, -30,  -60,  -90, -120, -63.011);

my $nAngles = scalar(@angles1);

for(my $i=0; $i<$nAngles; $i++)
{
    printf ("%4d - %4d = %4d\n", $angles2[$i], $angles1[$i], diff($angles2[$i], $angles1[$i]));
}

sub diff
{
    my($angle2, $angle1) = @_;
    my $diff = $angle2 - $angle1;
    $diff += 360 while($diff < 0);
    $diff -= 360 while($diff > 360);
    return($diff);
}
