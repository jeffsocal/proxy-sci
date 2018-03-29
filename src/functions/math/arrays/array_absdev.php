<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */

/*
 * returns the absolute deviation
 */
function array_absdev($array)
{
    $array = array_rmna($array);
    $average = 0;
    $stdev = 0;
    $count = 0;
    $i = count($array);
    $average = array_mean($array);
    for ($x = 0; $x < $i; $x ++) {
        // normal standard distribution
        $stdev += abs($average - $array[$x]);
        $count ++;
    }
    if ($count > 1)
        $stdev = ($stdev / ($count - 1));
    return ($stdev);
}

?>
