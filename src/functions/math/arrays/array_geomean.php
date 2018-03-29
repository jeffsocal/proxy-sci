<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */

/*
 * returns the geometric mean
 */
function array_geomean($array)
{
    $array = array_rmna($array);
    $average = 0;
    $stdev = 0;
    $count = 0;
    $min = 0;
    $i = sizeof($array);
    foreach ($array as $x) {
        if ($x > 0) {
            $average += log($x);
            $count ++;
        }
    }
    if ($count != 0) {
        // geometric mean
        $average = exp($average / $count);
        return ($average);
    } else {
        return (0);
    }
}

?>
