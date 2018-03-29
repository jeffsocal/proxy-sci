<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
namespace ProxySci\Omics\Proteomics;

use ProxyIO\File\Delim\ReadDelim;
use ProxySci\Atomic\Mass;

class Peptide extends ReadDelim
{

    protected $mass_isotope;

    protected $mass_proton;

    protected $mass_H;

    protected $mass_O;

    protected $array_aa_data;

    protected $array_aa_metrics;

    protected $array_aa_exactmass;

    private $seq_regex = "/[a-zA-Z]|\[.+?\]/";

    private $seq_mod_regex = "/[a-zA-Z]|[\-\+]*\d+\.*\d*/";

    private $atom;

    function __construct()
    {
        $ini = parse_ini_file(get_include_path() . 'ini/molecular.ini');
        parent::__construct(get_include_path()  . $ini['aminos_path']);
        
        $this->atom = new Mass();
        
        $this->array_aa_data = $this->getTableArray();
        
        $this->mass_proton = 1.00727646688;
        $this->mass_isotope = 1.0025;
        $this->mass_H = $this->atom->getMolecularWeight('H');
        $this->mass_O = $this->atom->getMolecularWeight('O');
        
        $this->setExactMassAminoAcid();
    }

    //
    private function setExactMassAminoAcid()
    {
        
        //
        $this->array_aa_exactmass = array();
        
        //
        $water = $this->mass_O + $this->mass_H * 2;
        foreach ($this->array_aa_data['1_letter'] as $n => $l) {
            $this->array_aa_exactmass[$l] = $this->atom->getMolecularWeight($this->array_aa_data['Res_Formula'][$n]);
        }
        
        // C3H7NO2Se - Selenocysteine
        $this->array_aa_exactmass['U'] = $this->atom->getMolecularWeight('C3H7NO2Se');
        $this->array_aa_exactmass['Z'] = $this->atom->getMolecularWeight('C3H7NO2Se');
        $this->array_aa_exactmass['&'] = $this->array_aa_exactmass['L'];
        $this->array_aa_exactmass['B'] = ($this->array_aa_exactmass['N'] + $this->array_aa_exactmass['D']) / 2;
        $this->array_aa_exactmass['Z'] = ($this->array_aa_exactmass['E'] + $this->array_aa_exactmass['Q']) / 2;
        $this->array_aa_exactmass['J'] = ($this->array_aa_exactmass['I'] + $this->array_aa_exactmass['L']) / 2;
        $this->array_aa_exactmass['X'] = ($this->array_aa_exactmass['Q'] + $this->array_aa_exactmass['K']) / 2;
        $this->array_aa_exactmass['c'] = $this->mass_O + $this->mass_H;
        $this->array_aa_exactmass['n'] = $this->mass_H;
        $this->array_aa_exactmass['p'] = $this->mass_proton;
        
        // common PTMs
        // [C57.021]
        $this->array_aa_exactmass['z'] = $this->array_aa_exactmass['C'] + 57.021;
        // [K42.011]
        $this->array_aa_exactmass['y'] = $this->array_aa_exactmass['K'] + 42.011;
        // [K43.006]
        $this->array_aa_exactmass['x'] = $this->array_aa_exactmass['K'] + 43.006;
        // [M15.995]
        $this->array_aa_exactmass['w'] = $this->array_aa_exactmass['M'] + 15.995;
        
        /*
         * Ambiguous Amino Acids 3-Letter 1-Letter
         * Asparagine or aspartic acid Asx B -> N (CONVERSION FOR OMSSA DB SEARCH)
         * Glutamine or glutamic acid Glx Z -> Q
         * Leucine or Isoleucine Xle J -> L
         * Unspecified or unknown amino acid Xaa X
         */
    }

    public function getExactMass($component)
    {
        if (is_numeric($component))
            return $component;
        
        if (key_exists($component, $this->array_aa_exactmass))
            return $this->array_aa_exactmass[$component];
        
        return 0;
    }

    private function getExactMassArray()
    {
        return $this->array_aa_exactmass;
    }

    protected function getSeqArray($aa, $regex = "[A-Z]")
    {
        preg_match_all($regex, $aa, $aa_m);
        
        if (! is_array($aa_m) or ! is_array($aa_m[0]))
            return false;
        
        return $aa_m[0];
    }

    public function getMolecularWeight($aa)
    {
        $molecular_weight = 0;
        foreach ($this->getSeqArray('n' . $aa . 'c', $this->seq_regex) as $n => $v) {
            
            foreach ($this->getSeqArray($v, $this->seq_mod_regex) as $nn => $vv) {
                $molecular_weight += $this->getExactMass($vv);
            }
        }
        
        return $molecular_weight;
    }

    public function getAfromAaa($Aaa)
    {
        $this_key = array_search($Aaa, $this->array_aa_data['3_letter']);
        if ($this_key == false)
            return false;
        
        return $this->array_aa_data['1_letter'][$this_key];
    }
}
?>