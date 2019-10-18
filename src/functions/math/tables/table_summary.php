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
function table_summary($table_array)
{
    $table_cols = table_header($table_array);
    
    $new_table = [];
    foreach ($table_cols as $table_col) {
        $values = $table_array[$table_col];
        
        $new_table['column'][] = $table_col;
        $new_table['n'][] = count($values);
        $new_table['unique'][] = count(array_unique($values));
        $new_table['max'][] = array_max(array_unique($values));
        $new_table['median'][] = array_median(array_unique($values));
        $new_table['mean'][] = array_mean(array_unique($values));
        $new_table['min'][] = array_min(array_unique($values));
        $new_table['q25'][] = array_quantile(array_unique($values), 0.25);
        $new_table['q75'][] = array_quantile(array_unique($values), 0.75);
        $new_table['q005'][] = array_quantile(array_unique($values), 0.005);
        $new_table['q995'][] = array_quantile(array_unique($values), 0.995);
    }
    
    return $new_table;
}


?>