<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
function histogram($array, $binwidth = null)
{
    $num = sizeof($array);
    $array = array_rmna($array);
    $min = array_min($array) * .9;
    $max = array_max($array) * 1.1;
    $dif = $max - $min;
    
    $seg = 32;
    if (! is_null($binwidth)) {
        $seg = ceil($dif / $binwidth);
    }
    $binwidth = ceil($dif / $seg);
    
    $table['y'] = $table['x'] = array_fill(0, $seg, 0);
    
    foreach ($table['y'] as $n => $v) {
        
        $x_val = $min + ($n * $binwidth);
        $bound_lower = $x_val - $binwidth * .5;
        $bound_upper = $x_val + $binwidth * .5;
        
        $table['x'][$n] = $x_val;
        $table['y'][$n] = sizeof(array_range($array, ">= $bound_lower && < $bound_upper"));
    }
    
    return $table;
}

function density($array, $binwidth = NULL)
{
    $min = array_min($array);
    $max = array_max($array);
    $dif = $max - $min;
    
    $seg = 48;
    
    if (! is_null($binwidth)) {
        $seg = ceil($dif / $binwidth);
    }
    $binwidth = ceil($dif / $seg);
    
    $table = histogram($array, $binwidth);
    
    // $max = array_max($table['y']);
    // foreach ($table['y'] as $n => $v) {
    // $table['y'][$n] = number_format($table['y'][$n] / $max, 6);
    // }
    
    // $table = cartesian_smoothed($table['x'], $table['y'], 6);
    
//     $sum = array_sum($table['y']);
//     foreach ($table['y'] as $n => $v) {
//         $table['y'][$n] = number_format($table['y'][$n] / $sum, 6);
//     }
    
    return $table;
}
?>
