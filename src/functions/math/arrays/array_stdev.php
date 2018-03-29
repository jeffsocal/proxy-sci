<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */

/*
 * returns the standard deviation
 */
function array_stdev($array)
{
    $array = array_rmna($array);
    $average = 0;
    $stdev = 0;
    $count = 0;
    $i = count($array);
    $average = array_mean($array);
    for ($x = 0; $x < $i; $x ++) {
        // normal standard distribution
        $stdev += ($average - $array[$x]) * ($average - $array[$x]);
        $count ++;
    }
    if ($count > 1)
        $stdev = sqrt($stdev / ($count - 1));
    return ($stdev);
}

?>
