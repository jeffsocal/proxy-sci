<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */

/*
 * returns the sum of partial array given the start and ending keys
 */
function array_psum($array, $start, $length)
{
    $new_array = array_rmna(array_slice($array, $start, $length));
    return array_sum($new_array);
}

?>
