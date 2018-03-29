<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */

/*
 * returns the weighted mean
 */
function array_weightmean($array, $s, $w)
{
    $array = array_rmna($array);
    $average = 0;
    $weight = 0;
    $stdev = 0;
    $count = 0;
    $min = 0;
    
    $i = count($array[$s]);
    for ($x = 0; $x < $i; $x ++) {
        if (! is_numeric($array[$s][$x] * 1))
            continue;
        if (! is_numeric($array[$w][$x] * 1))
            continue;
        $average += $array[$s][$x] * $array[$w][$x];
        $weight += $array[$w][$x];
        $count ++;
    }
    if ($count != 0) {
        $average = "$average" / "$weight";
        return ("$average");
    } else {
        return (0);
    }
}

?>
