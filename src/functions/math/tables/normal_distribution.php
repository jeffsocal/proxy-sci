<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */

/*
 * returns the normal distribution
 */
function normal_distribution($x, $stdev, $h = 1)
{
    $x_v = $x;
    $x_i = $x_v - $stdev * 5;
    $x_j = $x_v + $stdev * 5;
    
    $x_d = ($x_j - $x_i) / 1000;
    
    $e = 2.718281828;           // (Euler's number)
    $c = $stdev / 2.35482;
    $c = 2 * $c ^ 2;
    
    $i = 0;
    $table = [
        'x',
        'y'
    ];
    
    for ($x = $x_i; $x <= $x_j; $x += $x_d) {
        
        $b = pow(($x - $x_v), 2);
        $y = $h * pow($e, - ($b / $c));
        
        $table['x'][] = $x;
        $table['y'][] = $y;
    }
    return $gaussianArray;
}

?>
