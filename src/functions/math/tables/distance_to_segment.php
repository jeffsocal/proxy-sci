<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */

/*
 * returns the distance from a point to a discrete line
 */
function distance_to_segment($x, $y, $p1, $p2)
{
    $p1 = explode(",", $p1);
    $p2 = explode(",", $p2);
    
    $ys = ($p2[1] - $p1[1]);
    $xs = ($p2[0] - $p1[0]);
    
    $d = abs($xs * ($p1[1] - $y) - ($p1[0] - $x) * $ys) / sqrt($xs * $xs + $ys * $ys);
    
    return $d;
}
?>
