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
function array_range($array, $f = NULL)
{
    if ($f == NULL)
        return $array;
    
    if (! preg_match("/^[\>\<]=*\s+\-*\d+\.*\d*\s*\&*\s*[\>\<]*=*\s*\-*\d*\.*\d*\s*$/", $f))
        return $array;
    
    if (strstr($f, '&')) {
        if (! preg_match("/^[\>\<]=*\s+\-*\d+\.*\d*\s+\&{2}\s+[\>\<]=*\s+\-*\d+\.*\d*\s*$/", $f))
            return $array;
    }
    
    $array = array_rmna($array);
    
    $f = str_replace("<", "\$n <", $f);
    $f = str_replace(">", "\$n >", $f);
    
    $e = "\$out = array_filter(\$array, function (\$n) {return $f ;});";
    
    eval($e);
    
    return $out;
}

?>
