<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */

/*
 * returns the arthmetic mean
 */
function array_arthmean($array)
{
    $array = array_rmna($array);
    $s = array_sum($array);
    $i = count($array);
    if ($i != 0) {
        $average = $s / $i;
        return ($average);
    } else {
        return ($s);
    }
}

?>
