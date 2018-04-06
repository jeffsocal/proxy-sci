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
function expand_grid($arrays)
{
    $arrays = cartesian_product($arrays);
    
    $n = count($arrays);
    
    for ($i = 0; $i < $n; $i ++) {
        $pair = $arrays[$i];
        unset($arrays[$i]);
        
        if ($pair[0] == $pair[1]) {
            continue;
        }
        if (trim($pair[0]) == '' or trim($pair[1]) == '') {
            continue;
        }
        
        sort($pair);
        
        $arrays[array_tostring($pair, '-', '')] = $pair;
    }
    
    return array_values($arrays);
}

?>
