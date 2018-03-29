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

    private $lengthUpperLimit;

    private $lengthLowerLimit;

    public function __construct($lengthLowerLimit = 6, $lengthUpperLimit = 36)
    {
        parent::__construct();
        $this->setLengthLimits($lengthLowerLimit, $lengthUpperLimit);
    }

    private function setLengthLimits($low, $high)
    {
        $this->lengthLowerLimit = $low;
        $this->lengthUpperLimit = $high;
    }

    public function getPeptides($sequence)
    {
        $sequence = preg_replace("/[\n\r\t\s\.\-\*]/", '', $sequence);
        
        /*
         * get all regular peptides
         */
        preg_match_all("/.*?[K|R]/", $sequence, $peptides);
        
        /*
         * get the last non-specific tryptic peptide and add it to the end of the list
         */
        preg_match_all("/.*?[K|R]/", strrev($sequence), $peptides_rev);
        
        //
        if (key_exists(0, $peptides_rev[0])) {
            $peptides_rev = strrev(preg_replace("/[K|R]/", "", $peptides_rev[0][0]));
            if ($peptides_rev != "")
                $peptides[0][] = $peptides_rev;
        }
        
        if (sizeof($peptides[0]) == 0) {
            /*
             * systemError("Low Peptide Count: $i");
             */
            $peptides = array(
                $sequence
            );
        } else {
            $peptides = $peptides[0];
        }
        
        //
        $i = sizeof($peptides);
        $concat_peptide = "";
        $array_peptides = array();
        
        $s = 0;
        for ($n = 0; $n < $i; $n ++) {
            $concat_peptide .= $peptides[$n];
            $stl = strlen($concat_peptide);
            if ($stl <= $this->lengthUpperLimit) {
                if ($stl >= $this->lengthLowerLimit) {
                    $array_peptides[] = $concat_peptide;
                }
                $s ++;
            } else {
                $concat_peptide = "";
                $n -= $s;
                $s = 0;
            }
            if ($n == $i - 1 and $s != 0) {
                $concat_peptide = "";
                $n -= ($s - 1);
                $s = 0;
            }
        }
        return $array_peptides;
    }
}
?>