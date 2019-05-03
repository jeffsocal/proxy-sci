<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */

/*
 * remove NAs from 1d array
 */
function array_rmna($array)
{
    foreach ($array as $n => $value) {
        if (preg_match("/^\d+\.?\d*$/", $value))
            $value *= 1;
        if (! is_numeric($value)) {
            unset($array[$n]);
        }
    }
    return array_values($array);
}

?>
