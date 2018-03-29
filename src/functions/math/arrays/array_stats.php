<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */

/*
 * returns an array of all stats
 */
function array_stats($array)
{
    $out['mean'] = array_mean($array);
    $out['median'] = array_median($array);
    $out['sum'] = array_sum($array);
    $out['max'] = array_max($array);
    $out['min'] = array_min($array);
    $out['count'] = sizeof($array);
    
    $out['stdev'] = array_stdev($array);
    $out['absdev'] = array_absdev($array);
    $out['skew'] = array_skewness($array);
    
    return $out;
}

?>
