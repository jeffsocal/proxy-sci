<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */

/*
 *
 */
function array_normalize($array, $factor = 1)
{
    $array = array_rmna($array);
    if (sizeof($array) == 0) {
        return false;
    }
    
    $max = array_max($array);
    $min = array_min($array);
    if ($max == $min)
        $min = 0;
    
    foreach ($array as $n => $v) {
        $array[$n] = (($v - $min) / ($max - $min)) * $factor;
    }
    return $array;
}

?>
