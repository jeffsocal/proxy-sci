<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
namespace ProxySci\Omics\Proteomics;

class Digest extends Protein
{

    private $length_upper_limit;

    private $length_lower_limit;

    private $mass_upper_limit;

    private $mass_lower_limit;

    private $missed_clevage_max;

    private $enzyme_regex;

    private $Peptide;

    public function __construct($missed_clevage_max = 2)
    {
        parent::__construct();
        $this->setLengthLimits();
        $this->setMassLimits();
        $this->setMissedClevageMax();
        $this->setEnzymeRegex();
        
        $this->Peptide = new Peptide();
    }

    public function setEnzymeRegex($regex = "/.*?[KR\#]/")
    {
        $this->enzyme_regex = $regex;
    }

    public function setLengthLimits($low = 6, $high = 53)
    {
        $this->length_lower_limit = $low;
        $this->length_upper_limit = $high;
    }

    public function setMissedClevageMax($missed_clevage_max = 2)
    {
        $this->missed_clevage_max = $missed_clevage_max + 1;
    }

    public function setMassLimits($low = 600, $high = 6000)
    {
        $this->mass_lower_limit = $low;
        $this->mass_upper_limit = $high;
    }

    public function getPeptides($sequence, $sub_il = FALSE)
    {
        $sequence = preg_replace("/[\n\r\t\s\.\-\*]/", '', $sequence);
        
        if (is_true($sub_il))
            $sequence = preg_replace("/L/", "I", $sequence);
        
        /*
         * get all regular peptides
         */
        preg_match_all($this->enzyme_regex, $sequence . "#", $peptides);
        $peptides = preg_replace("/\#/", "", $peptides[0]);
        
        // print_r($peptides);
        // exit;
        
        //
        $i = sizeof($peptides);
        $concat_peptide = "";
        $array_peptides = array();
        
        $s = 0;
        for ($n = 0; $n < $i; $n ++) {
            
            $keep = 'keep';
            $concat_peptide .= $peptides[$n];
            
            $mass = $this->Peptide->getMolecularWeight($concat_peptide);
            $leng = strlen($concat_peptide);
            
            if ($mass < $this->mass_lower_limit || $leng < $this->length_lower_limit)
                $keep = 'null';
            
            if ($mass > $this->mass_upper_limit || $leng > $this->length_upper_limit)
                $keep = 'reject';
            
            if ($s >= $this->missed_clevage_max)
                $keep = 'reject';
            
            switch ($keep) {
                case 'keep':
                    $array_peptides[] = $concat_peptide;
                    $s ++;
                    break;
                
                case 'null':
                    $s ++;
                    break;
                
                case 'reject':
                    $concat_peptide = "";
                    $n -= $s;
                    $s = 0;
                    break;
            }
        }
        for ($n = 0; $n < min(count($array_peptides), $this->missed_clevage_max); $n ++) {
            if (preg_match("/^M/", $array_peptides[$n]))
                $array_peptides[] = preg_replace("/^M/", '', $array_peptides[$n]);
        }
        
        return $array_peptides;
    }
}
?>