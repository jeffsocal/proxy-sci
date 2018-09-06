<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
namespace ProxySci\Omics\Proteomics;

use ProxyIO\File\Delim\ReadDelim;
use function BenTools\CartesianProduct\cartesian_product;

class PeptideModifications extends Peptide
{

    private $table_unimod;

    public function __construct()
    {
        $ini = parse_ini_file(get_include_path() . 'ini/molecular.ini');
        $unicsv_path = get_include_path() . $ini['unicsv_path'];
        
        $csv = new ReadDelim($unicsv_path);
        
        $this->table_unimod = $csv->getTableArray();
    }

    public function getPeptideMods()
    {
        return $this->table_unimod;
    }

    public function getPTMPeptides($peptide, $n_combs = 2)
    {
        $ps = str_split($peptide);
        $ms = [];
        foreach ($ps as $n => $aa) {
            
            $ms[$n] = [];
            $rows = preg_grep("/" . $aa . "(?!\-term)/", $this->table_unimod['COM_AMINO']);
            
            if ($n == 0)
                $rows = preg_grep("/" . $aa . "|N\-term/", $this->table_unimod['COM_AMINO']);
            
            if ($n == (count($ps) - 1))
                $rows = preg_grep("/" . $aa . "|C\-term/", $this->table_unimod['COM_AMINO']);
            
            if (count($rows) != 0)
                $ms[$n] = array_keys($rows);
            
            $ms[$n][] = '';
        }
        
        $cp = [];
        $cp = cartesian_product($ms);
        
        $pa = [];
        foreach ($cp as $i => $g) {
            $w = array_filter($g);
            if (count($w) > $n_combs)
                continue;
            
            $pst = $ps;
            foreach ($w as $key => $ptm) {
                
                $pst[$key] = "[" . $pst[$key] . round($this->table_unimod['MONO_MASS'][$ptm], 2) . "]";
            }
            
            $pa[] = array_tostring($pst, "", "");
        }
        
        return array_values(array_unique($pa));
    }
}
?>