<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
namespace ProxySci\Omics\Proteomics;

use function BenTools\CartesianProduct\cartesian_product;

class PeptideModifications extends Peptide
{

    private $array_unimod;

    private $max_concur;

    private $max_children;

    public function __construct($ini_array = false)
    {
        $ini = parse_ini_file(get_include_path() . 'ini/molecular.ini');
        $this->array_unimod = parse_json_file($ini['unimod_path']);
        
        $this->max_concur = $ini['max_concurrent'];
        $this->max_children = $ini['max_children'];
        
        if (is_false($ini_array))
            return null;
        
        if (key_exists('max_concurrent', $ini_array))
            $this->max_concur = $ini_array['max_concurrent'];
        
        if (key_exists('max_children', $ini_array))
            $this->max_children = $ini_array['max_children'];
        
        if (key_exists('modification_set', $ini_array)) {
            
            $key = 'applied_sites';
            $mod_set = $ini_array['modification_set'];
            $this->array_unimod = array_intersect_key($this->array_unimod, $mod_set);
            
            foreach ($this->array_unimod as $mod => $array) {
                if (key_exists($key, $mod_set[$mod]))
                    $this->array_unimod[$mod][$key] = $mod_set[$mod][$key];
            }
            
        }
    }

    public function getPTMPeptidesV1($peptide, $n_combs = 2)
    {
        $ps = str_split($peptide);
        $pos = $this->getPeptidePtmPositionArray($peptide);
        
        $as = $pos['seq'];
        $ms = $pos['mod'];
        
        $cp = cartesian_product($ms);
        
        if (count($cp) < 1e7)
            return $this->applyCP($ps, $cp, $n_combs);
        
        return $this->applySingle($ps, $ms);
        
        // print_r($ms);
        // print_r($cp->asArray());
    }

    private function getPeptidePtmPositionArray($peptide)
    {
        $ps = str_split($peptide);
        $ms = [];
        $as = [];
        
        foreach ($ps as $n => $aa) {
            
            $rows = preg_grep("/" . $aa . "(?!\-term)/", $this->table_unimod['COM_AMINO']);
            
            if ($n == 0)
                $rows = preg_grep("/" . $aa . "|N\-term/", $this->table_unimod['COM_AMINO']);
            
            if ($n == (count($ps) - 1))
                $rows = preg_grep("/" . $aa . "|C\-term/", $this->table_unimod['COM_AMINO']);
            
            if (count($rows) != 0) {
                $ms[$n] = array_keys($rows);
                $as[] = $n;
            }
            $ms[$n][] = '';
        }
        
        $as[] = '';
        
        $out['seq'] = $as;
        $out['mod'] = $ms;
        
        return $out;
    }

    public function getPTMPeptides($peptide, $n_combs = 2)
    {
        $ps = str_split($peptide);
        $pos = $this->getPeptidePtmPositionArray($peptide);
        
        $as = $pos['seq'];
        $ms = $pos['mod'];
        
        $n_comb_max = $n_combs;
        
        $num_cp_allpossible = cartesian_count($ms);
        
        // print_message(" - cp all " . $n_combs . " combs", $num_cp_allpossible, ".", 30);
        
        for ($n_combs = $n_comb_max; $n_combs >= 0; $n_combs --) {
            
            $bs = [];
            for ($i = 0; $i < $n_combs; $i ++) {
                $bs[] = $as;
            }
            
            $num_cp_allaminos = cartesian_count($bs);
            // print_message(" - cp aminos " . $n_combs . " combs", $num_cp_allaminos, ".", 30);
            
            if (min($num_cp_allaminos, $num_cp_allpossible) > $this->max_children)
                continue;
            
            if ($num_cp_allpossible <= $num_cp_allaminos) {
                $cp_allpossible = cartesian_product($ms);
                // print_message(" - run ", "all-posible", ".", 30);
                return $this->applyCP($ps, $cp_allpossible, $n_combs);
            }
            // print_message(" - run ", "by-each-amino", ".", 30);
            
            $cp_allaminos = cartesian_product($bs);
            
            $ncp = [];
            foreach ($cp_allaminos->asArray() as $comb) {
                sort($comb);
                $comb = array_unique($comb);
                $ncp['p' . array_tostring($comb, '', '')] = $comb;
            }
            
            $out = [];
            foreach ($ncp as $comb) {
                
                $this_ms = array_intersect_key($ms, array_flip($comb));
                
                $cp = cartesian_product($this_ms);
                
                $this_out = $this->applyCP($ps, $cp, $n_combs);
                $out = array_unique(array_merge($out, $this_out));
            }
            
            return $out;
        }
    }

    private function applyCP($ps, $cp, $n_combs)
    {
        $pa = [];
        foreach ($cp as $i => $g) {
            $w = array_filter($g);
            
            /*
             * check if number of ptms exceeds allowable limit
             */
            if (count($w) > $n_combs)
                continue;
            
            /*
             * check for max allowable concurrent ptms
             */
            $nw = array_count_values($w);
            $over_max_concur = FALSE;
            foreach ($nw as $key => $count) {
                if ($count > $this->table_unimod['MAX_CONCUR'][$key]) {
                    $over_max_concur = TRUE;
                    break;
                }
            }
            if (is_true($over_max_concur))
                continue;
            
            $pst = $ps;
            foreach ($w as $key => $ptm) {
                
                $pst[$key] = "[" . $pst[$key] . round($this->table_unimod['MONO_MASS'][$ptm], 2) . "]";
            }
            
            $pa[] = array_tostring($pst, "", "");
        }
        
        return array_values(array_unique($pa));
    }

    private function applySingle($ps, $ms)
    {
        $pa = [];
        foreach ($ms as $i => $g) {
            $w = array_filter($g);
            
            foreach ($w as $key => $ptm) {
                $pst = $ps;
                
                $pst[$i] = "[" . $pst[$i] . round($this->table_unimod['MONO_MASS'][$ptm], 2) . "]";
                $pa[] = array_tostring($pst, "", "");
            }
        }
        
        return array_values(array_unique($pa));
    }
}
?>