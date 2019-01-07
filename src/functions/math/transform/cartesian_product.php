<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */

/*
 * provide the all possible iterations value
 */
function cartesian_count(array $cases)
{
    $count = 1;
    foreach ($cases as $c) {
        $count *= count($c);
    }
    return $count;
}

?>
