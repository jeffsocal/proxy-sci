<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */

function getRootSequence($sequence)
{
    return preg_replace("/\[|\]|[0-9]|\./", "", $sequence);
}

?>