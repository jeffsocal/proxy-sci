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
function array_align($array1, $array2, $return = 'keys')
{
    $matrix = array();
    foreach ($array1 as $i => $i_v) {
        foreach ($array2 as $j => $j_v) {
            $matrix_d['v'][] = abs($i_v - $j_v);
            $matrix_d['i'][] = $i;
            $matrix_d['j'][] = $j;
        }
    }
    
    $matrix_d = table_sort($matrix_d, 'v', 'ascn', 'number');
    foreach ($matrix_d as $i) {}
    return $matrix_d;
}

?>
