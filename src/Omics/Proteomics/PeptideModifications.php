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

    private $applied_sites;

    private $applied_mods;

    public function __construct($ini_array = false)
    {
        $ini = parse_ini_file(get_include_path() . 'ini/molecular.ini');
        $this->array_unimod = parse_json_file($ini['unimod_path']);
        
        $this->max_concur = $ini['max_concurrent'];
        $this->max_children = $ini['max_children'];
        
        if (! is_false($ini_array)) {
            
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
        /*
         * create an array of where modifications are applied
         */
        
        $this->applied_sites = [];
        $this->applied_mods = [];
        
        $this->applied_mods[] = 'none';
        foreach ($this->array_unimod as $mod => $array) {
            $this->applied_mods[] = $mod;
        }
        
        $applied_mods_n = array_flip($this->applied_mods);
        
        foreach ($this->array_unimod as $mod => $array) {
            foreach ($array['applied_sites'] as $site) {
                $this->applied_sites[$site][] = $applied_mods_n[$mod];
            }
        }
    }

    public function getPTMPeptidesV1($peptide, $max_concur = 2)
    {
        $ps = str_split($peptide);
        $pos = $this->getPeptidePtmPositionArray($peptide);
        
        $as = $pos['seq'];
        $ms = $pos['mod'];
        
        $cp = cartesian_product($ms);
        
        if (count($cp) < 1e7)
            return $this->applyCP($ps, $cp, $max_concur);
        
        return $this->applySingle($ps, $ms);
        
        // print_r($ms);
        // print_r($cp->asArray());
    }

    private function getPeptidePtmPositionArray($peptide)
    {
        $ps = str_split($peptide);
        $es = [];
        
        foreach ($ps as $n => $aa) {
            
            $ps[$n] = [];
            $ps[$n][0] = 0;
            if (key_exists($aa, $this->applied_sites))
                $ps[$n] = array_merge($ps[$n], $this->applied_sites[$aa]);
            
            if ($n == 0 & key_exists('N-TERM', $this->applied_sites))
                $ps[$n] = array_merge($ps[$n], $this->applied_sites['N-TERM']);
            
            if ($n == (count($ps) - 1) & key_exists('C-TERM', $this->applied_sites))
                $ps[$n] = array_merge($ps[$n], $this->applied_sites['C-TERM']);
            
            $ps[$n] = array_unique($ps[$n]);
            
            if (count($ps[$n]) > 1)
                $es[] = $n;
        }
        
        $out['seq'] = $ps;
        $out['loc'] = $es;
        return $out;
    }

    public function getPTMPeptides($peptide)
    {
        $ps = str_split($peptide);
        $out = $this->getPeptidePtmPositionArray($peptide);
        
        $seq = $out['seq'];
        $loc = $out['loc'];
        
        /*
         * calc n peptides from all possible ptm applications
         */
        $num_cp_allpossible = cartesian_count($seq);
        
        for ($max_concur = $this->max_concur; $max_concur >= 0; $max_concur --) {
            
            $bs = [];
            for ($i = 0; $i < $max_concur; $i ++) {
                $bs[] = $loc;
            }
            
            $num_cp_allaminos = cartesian_count($bs);
            
            if (min($num_cp_allaminos, $num_cp_allpossible) > $this->max_children)
                continue;
            
            if ($num_cp_allpossible <= $num_cp_allaminos) {
                $cp_allpossible = cartesian_product($seq)->asArray();
                
                return $this->applyCP($ps, $cp_allpossible, $max_concur);
            }
            
            $cp_allaminos = cartesian_product($bs)->asArray();
            
            $ncp = [];
            foreach ($cp_allaminos as $comb) {
                sort($comb);
                $comb = array_unique($comb);
                $ncp['p' . array_tostring($comb, '', '')] = $comb;
            }
            
            $out = [];
            foreach ($ncp as $comb) {
                
                $this_ms = array_intersect_key($seq, array_flip($comb));
                
                $cp = cartesian_product($this_ms)->asArray();
                
                $this_out = $this->applyCP($ps, $cp, $max_concur);
                $out = array_unique(array_merge($out, $this_out));
            }
            
            return $out;
        }
    }

    private function applyCP($ps, $cp, $max_concur = 3)
    {
        $applied_mods_n = array_flip($this->applied_mods);
        $pa = [];
        foreach ($cp as $i => $g) {
            $w = array_filter($g);
            
            /*
             * check if number of concurrent ptms exceeds allowable global limit
             */
            if (count($w) > $max_concur)
                continue;
            
            $nw = array_count_values($w);
            $over_max_concur = FALSE;
            foreach ($nw as $ptm => $count) {
                
                /*
                 * check if number of concurrent ptms exceeds allowable local limit
                 */
                
                if ($count > $this->array_unimod[$this->applied_mods[$ptm]]['max_concurrent']) {
                    $over_max_concur = TRUE;
                    break;
                }
            }
            if (is_true($over_max_concur))
                continue;
            
            /*
             * apply the ptms
             */
            $pst = $ps;
            foreach ($w as $key => $ptm) {
                
                $pst[$key] = "[" . $pst[$key] . round($this->array_unimod[$this->applied_mods[$ptm]]['mass_monoiso'], 2) . "]";
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