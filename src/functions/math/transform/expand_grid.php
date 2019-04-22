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
use function BenTools\CartesianProduct\cartesian_product;

function expand_grid($arrays)
{
    $array_eg = [];
    for($i = 0; $i < count($arrays); $i++){
        $array_eg[] = $arrays;
    }
    
    $array_eg = cartesian_product($array_eg)->asArray();
    
    $n = count($array_eg);
    
    for ($i = 0; $i < $n; $i ++) {
        $pair = $array_eg[$i];
        unset($array_eg[$i]);
        
        if ($pair[0] == $pair[1]) {
            continue;
        }
//         if (trim($pair[0]) == '' or trim($pair[1]) == '') {
//             continue;
//         }
        
        $pair = array_unique($pair);
        sort($pair);
        
        $array_eg[array_tostring($pair, '-', '')] = $pair;
    }
    
    return array_values($array_eg);
}

?>
