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
function linear_regression($x, $y)
{
    
    /*
     * calculate number points
     */
    $n = count($x);
    
    /*
     * ensure both arrays of points are the same size
     */
    if ($n != count($y)) {
        systemError("linear_regression(): Number of elements in coordinate arrays do not match.");
    }
    /*
     * calculate sums
     */
    $x_sum = array_sum($x);
    $y_sum = array_sum($y);
    
    $xx_sum = 0;
    $xy_sum = 0;
    
    foreach ($x as $i => $v) {
        $xy_sum += ($x[$i] * $y[$i]);
        $xx_sum += ($x[$i] * $x[$i]);
    }
    
    /*
     * calculate slope
     */
    $m = (($n * $xy_sum) - ($x_sum * $y_sum)) / (($n * $xx_sum) - ($x_sum * $x_sum));
    
    /*
     * calculate intercept
     */
    $b = ($y_sum - ($m * $x_sum)) / $n;
    
    /*
     * Pearson correlation moment
     * r^2 = 1 - [pcm_error] / [pcm_total]
     */
    $w = $y_sum / $n;
    
    foreach ($y as $n => $z) {
        $pcm_total[] = pow(($z - $w), 2);
        $f = $m * $x[$n] + $b;
        $pcm_error[] = pow(($z - $f), 2);
    }
    
    $r = 1 - array_sum($pcm_error) / array_sum($pcm_total);
    
    $out['slope'] = $m;
    $out['intercept'] = $b;
    $out['pearson'] = $r;
    
    return $out;
}
?>
