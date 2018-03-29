<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */

/*
 * returns the distance from a point to a continous line
 */
function distance_to_line($x, $y, $m, $b)
{
    $a = '-1' * $m;
    $b *= '-1';
    
    $d = abs($a * $x + $y + $b) / sqrt($a * $a + 1);
    
    return $d;
}

?>
