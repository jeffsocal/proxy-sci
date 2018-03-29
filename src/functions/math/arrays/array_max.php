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
function array_max($array)
{
    $array = array_rmna($array);
    if (sizeof($array) == 0) {
        return false;
    }
    return max($array);
}

?>
