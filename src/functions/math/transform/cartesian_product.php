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
function cartesian_product($arrays)
{
    $result = array();
    $arrays = array_values($arrays);
    $sizeIn = sizeof($arrays);
    $size = $sizeIn > 0 ? 1 : 0;
    foreach ($arrays as $array) {
        $size = $size * sizeof($array);
    }
    for ($i = 0; $i < $size; $i ++) {
        $result[$i] = array();
        for ($j = 0; $j < $sizeIn; $j ++) {
            array_push($result[$i], current($arrays[$j]));
        }
        for ($j = ($sizeIn - 1); $j >= 0; $j --) {
            if (next($arrays[$j])) {
                break;
            } elseif (isset($arrays[$j])) {
                reset($arrays[$j]);
            }
        }
    }
    return $result;
}

?>
