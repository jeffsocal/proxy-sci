<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
namespace ProxySci\Atomic;

use ProxyIO\File\Delim\ReadDelim;

class Mass extends ReadDelim
{

    private $array_atomic_mass;

    function __construct()
    {
        $ini = parse_ini_file(get_include_path() . 'ini/molecular.ini');
        parent::__construct(get_include_path() . $ini['atomic_path']);
        $this->array_atomic_mass = $this->getTableArray();
    }

    public function getMolecularWeight($molfmla)
    {
        preg_match_all("/[A-Z]{1}[a-z]{0,1}\d*/", $molfmla, $molfmla_m);
        
        $molecular_weight = 0;
        foreach ($molfmla_m[0] as $n => $atomic_parts) {
            
            preg_match("/[A-Z]{1}[a-z]{0,1}/", $atomic_parts, $atomic_symbol);
            
            $atomic_quant[0] = 0;
            if (! preg_match("/\d+/", $atomic_parts)) {
                $atomic_quant[0] = 1;
            } else {
                preg_match("/\d+/", $atomic_parts, $atomic_quant);
            }
            
            $atomic_key = array_search($atomic_symbol[0], $this->array_atomic_mass['atomic_sym']);
            if (is_bool($atomic_key) == true) {
                continue;
            }
            
            $molecular_weight += $this->array_atomic_mass['atomic_mass'][$atomic_key] * $atomic_quant[0];
        }
        
        return $molecular_weight;
    }
}

?>