#!/usr/local/bin/perl

$nlines = 0;
$current_type = 800;

MAIN: {
        while (<>) {
                $nlines++;
                $type = substr($_, 9, 4);
                if ($nlines == 1 || $type == -1) {print; next;}
                if ($type =~ /^[ ]*$/ || !($type =~ /^[ ]*[0-9]+[ ]*$/)) {print; last;}
                printf("%s%4d%s", substr($_,0,9), $current_type, substr($_,13,9999));
                $current_type++;
        }
        while (<>) {print;}
}

