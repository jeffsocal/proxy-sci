<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */

/*
 * returns the numeric frequency
 */
function array_frequency($array, $binwidth = null)
{
    $array = array_rmna($array);
    
    $min = array_min($array);
    $max = array_max($array);
    
    if (is_null($binwidth))
        $binwidth = ($max - $min) / 10;
    
    $bins = round(($max - $min) / $binwidth);
    
    $average = 0;
    $count = 0;
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
