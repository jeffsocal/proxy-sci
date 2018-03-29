<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */

/*
 * returns skewness
 */
function array_skewness($array)
{
    $array = array_rmna($array);
    $i = count($array);
    $average = array_mean($array);
    
    $kTop = 0;
    $kBottom = 0;
    
    for ($x = 0; $x < $i; $x ++) {
        $val = ($array[$x] - $average);
        $kTop += pow($val, 3);
        $kBottom += pow($val, 2);
    }
    $valTop = $kTop / ($i);
    $valBottom = sqrt(pow(($kBottom / ($i)), 3));
    
    if ($valBottom == 0)
        return 0;
    return ($valTop / $valBottom);
}

?>
