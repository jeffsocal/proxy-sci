<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */

/*
 * returns a partial string given the start and ending keys of an array
 */
function array_pstr($array, $start, $length)
{
    $new_array = array_slice($array, $start, $length);
    $i = sizeof($new_array);
    $retstr = '';
    for ($x = 0; $x < $i; $x ++) {
        $retstr .= $new_array[$x];
    }
    return ($retstr);
}

?>
