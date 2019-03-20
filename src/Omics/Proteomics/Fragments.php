<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
namespace ProxySci\Omics\Proteomics;

class Fragments extends Peptide
{

    private $frag_series;

    private $frag_charge;

    private $frag_decay;

    private $frag_mass_max;

    private $frag_mass_min;

    private $m_int;

    private $m_slp;

    function __construct($input_series = "abycz", $input_charge = "123", $input_decay = "aw")
    {
        $this->frag_series = $input_series;
        $this->frag_charge = $input_charge;
        $this->frag_decay = $input_decay;
        
        // $this->m_int = 223.1347;
        // $this->m_slp = 1.000512;
        
        // $this->m_int = - 77.01891391;
        // $this->m_slp = 0.01000512049;
        
        $this->m_int = 584.9754;
        $this->m_slp = 0.9995063;
        
        $this->frag_mass_min = 240;
        $this->frag_mass_max = 1800;
        
        parent::__construct();
    }

    public function getSeqMassArray($aa)
    {
        /*
         * new out array
         */
        $array = [];
        
        foreach ($this->getSeqArray('n' . $aa . 'c', $this->seq_regex) as $n => $v) {
            $array['aa'][$n] = $v;
            $array['mass'][$n] = 0;
            foreach ($this->getSeqArray($v, $this->seq_mod_regex) as $nn => $vv) {
                $array['mass'][$n] += $this->getExactMass($vv);
            }
        }
        
        return $array;
    }

    public function getFragmentMassArray($aa)
    {
        /*
         * new out array
         */
        $array_frags = [];
        
        $array_seq = $this->getSeqMassArray($aa);
        $peptide_nmass = array_sum($array_seq['mass']);
        
        $amonia = ($this->mass_H * 3 + $this->mass_N);
        $water = ($this->mass_H * 2 + $this->mass_O);
        
        $j = count($array_seq['mass']) - 1;
        for ($i = 0; $i < $j; $i ++) {
            $b_ion_mass = array_sum(array_slice($array_seq['mass'], 0, $i));
            $y_ion_mass = array_sum(array_slice($array_seq['mass'], $i)) + $this->mass_H * 2;
            
            $b_ion_n = ($i - 1);
            $y_ion_n = ($j - $i);
            
            if ($b_ion_n > 1 && $b_ion_n < ($j - 1)) {
                $array_frags["b$b_ion_n"] = $b_ion_mass;
                if (strstr($this->frag_series, 'a'))
                    $array_frags["a$b_ion_n"] = $b_ion_mass - (12 + $this->mass_O);
                /*
                 * amonia loss
                 */
                if (strstr($this->frag_decay, 'a')) {
                    if (preg_grep('/R|K|Q|N/', array_slice($array_seq['aa'], 1, $i)))
                        $array_frags["b$b_ion_n" . "-a"] = $array_frags["b$b_ion_n"] - $amonia;
                }
                /*
                 * water loss
                 */
                if (strstr($this->frag_decay, 'a')) {
                    if (preg_grep('/S|T|E|D/', array_slice($array_seq['aa'], 1, $i)))
                        $array_frags["b$b_ion_n" . "-w"] = $array_frags["b$b_ion_n"] - $water;
                }
            }
            
            if ($y_ion_n > 0 && $y_ion_n < ($j - 1)) {
                $array_frags["y$y_ion_n"] = ($peptide_nmass - $b_ion_mass) + $this->mass_H * 2;
                /*
                 * amonia loss
                 */
                if (strstr($this->frag_decay, 'a')) {
                    if (preg_grep('/R|K|Q|N/', array_slice($array_seq['aa'], $i)))
                        $array_frags["y$y_ion_n" . "-a"] = $array_frags["y$y_ion_n"] - $amonia;
                }
                /*
                 * water loss
                 */
                if (strstr($this->frag_decay, 'a')) {
                    if (preg_grep('/S|T|E|D/', array_slice($array_seq['aa'], $i)))
                        $array_frags["y$y_ion_n" . "-w"] = $array_frags["y$y_ion_n"] - $water;
                }
            }
        }
        
        /*
         * additional observed charge states
         */
        foreach ($array_frags as $n => $v) {
            
            if (preg_match("/[czaA-Z]\d/", $n))
                continue;
            
            if (preg_match("/a|w/", $n))
                continue;
            
            if (strstr($n, "y") and strstr($this->frag_series, 'z'))
                $array_frags[str_replace('y', 'z', $n)] = $v - ($this->mass_H * 2 + $this->mass_N);
            
            if (strstr($n, "b") and strstr($this->frag_series, 'c'))
                $array_frags[str_replace('b', 'c', $n)] = $v + ($this->mass_H * 3 + $this->mass_N);
            
            if (strstr($this->frag_charge, '2'))
                $array_frags[$n . " 2+"] = ($v + $this->mass_proton) / 2;
            
            if (strstr($this->frag_charge, '3'))
                $array_frags[$n . " 3+"] = ($v + $this->mass_proton) / 3;
        }
        
        asort($array_frags);
        return $array_frags;
    }

    public function getFragmentWordArray($aa)
    {
        $array = $this->getFragmentMassArray($aa);
        
        /*
         * remove values with highly variable mass defect differences
         */
        $array = array_filter($array, function ($x) {
            return ($x > 240 & $x < 1800);
        });
        
        return array_map(array(
            $this,
            'massToHash'
        ), $array);
    }

    public function massToHash($float, $z_detect = false)
    {
        if (is_true($z_detect)) {
            $int = round(($float * $this->m_slp + $this->m_int) * 2) / 2;
            if (! is_int($int))
                $int = round((chargeMass(neutralMass($float, 2), 1) * $this->m_slp + $this->m_int));
        } else {
            $int = round(($float * $this->m_slp + $this->m_int));
        }
        // echo "WTF " . $hash . PHP_EOL;
        return base26_encode($int);
    }

    public function hashToMass($hash)
    {
        $int = base26_decode($hash);
        return ($int - $this->m_int) / $this->m_slp;
    }
}

/*
 * Inclusion of a and immonium ions comes from the observation of a# fragment ions
 * explination is as follows:
 *
 * Observation of y++ ions, may be additionally invlved in CID events to produce immonium and a+ ions
 * i.e. y++ => a+ + immonium
 *
 * "An internal fragment with just a single side chain formed by a combination of a type and y type cleavage is called an immonium ion"
 *
 * AU: Ioannis A. Papayannopoulos
 * TI: The interpretation of collision-induced dissociation tandem mass spectra of peptides
 * SO: Mass Spectrometry Reviews
 * VL: 14
 * NO: 1
 * PG: 49-73
 * YR: 1995
 * CP: Copyright 1995 John Wiley & Sons, Inc.
 * ON: 1098-2787
 * PN: 0277-7037
 * AD: Biogen Inc., 14 Cambridge Center, Cambridge, Massachusetts 02142
 * DOI: 10.1002/mas.1280140104
 * US: http://dx.doi.org/10.1002/mas.1280140104
 */

/*
 * Fragment Immonium Ions
 * Residue 3-letter 1-letter Immonium ion* Related ions*
 * Alanine Ala A 44
 * Arginine Arg R 129 59,70,73,87,100,112
 * Asparagine Asn N 87 70
 * Aspartic acid Asp D 88 70
 * Cysteine Cys C 76
 * Glutamic acid Glu E 102
 * Glutamine Gln Q 101 56,84,129
 * Glycine Gly G 30
 * Histidine His H [110] 82,121,123,138,166
 * Isoleucine Ile I [86] 44,72
 * Leucine Leu L [86] 44,72
 * Lysine Lys K 101 70,84,112,129
 * Methionine Met M 104 61
 * Phenylalanine Phe F [120] 91
 * Proline Pro P [70]
 * Serine Ser S 60
 * Threonine Thr T 74
 * Tryptophan Trp W [159] 77,117,130,132,170,171
 * Tyrosine Tyr Y [136] 91,107
 * Valine Val V 72 41,55,69
 */

/*
 * Fragment Ion Types
 * Ion
 * Type Neutral Mr
 * a [N]+[M]-CHO
 * a* a-NH3
 * a a-H2O
 * b [N]+[M]-H
 * b* b-NH3
 * b b-H2O
 * c [N]+[M]+NH2
 * d a - partial side chain
 * v y - complete side chain
 * w z - partial side chain
 * x [C]+[M]+CO-H
 * y [C]+[M]+H
 * y* y-NH3
 * y y-H2O
 * z [C]+[M]-NH2
 */

?>