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
function array_quantile($array, $quantile = 0.25)
{
    $quantile = min($quantile, 1);
    $array = array_rmna($array);
    if (sizeof($array) < (1 / $quantile)) {
        return false;
    }
    sort($array, SORT_NUMERIC);
    return $array[floor(sizeof($array) * $quantile) - 1];
}

?>
