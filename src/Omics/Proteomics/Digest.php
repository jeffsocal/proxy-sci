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
        $this->setMissedClevageMax($missed_clevage_max);
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

        if ($this->enzyme_regex == "." || $this->enzyme_regex == "/./")
            return $this->nonspecificPeptides($sequence);

        /*
         * get all regular peptides
         */
        preg_match_all($this->enzyme_regex, $sequence . "#", $peptides);

        if (sizeof($peptides[0]) == 0)
            return FALSE;

        $peptides = preg_replace("/\#/", "", $peptides[0]);

        // print_r($peptides);
        // exit;

        //
        $i = sizeof($peptides);
        $concat_peptide = "";
        $array_peptides = array();

        $s = 0;
        for ($n = 0; $n < $i; $n ++) {

            for ($mc = 0; $mc <= $this->missed_clevage_max; $mc ++) {

                $concat_peptide = array_tostring(array_slice($peptides, $n, $mc), '', '');

                $leng = strlen($concat_peptide);
                $mass = $this->Peptide->getMolecularWeight($concat_peptide);

                if ($mass < $this->mass_lower_limit || $leng < $this->length_lower_limit)
                    continue;

                if ($mass > $this->mass_upper_limit || $leng > $this->length_upper_limit)
                    break;

                $array_peptides[] = $concat_peptide;
            }
        }
        for ($n = 0; $n < min(count($array_peptides), $this->missed_clevage_max); $n ++) {
            if (preg_match("/^M/", $array_peptides[$n]))
                $array_peptides[] = preg_replace("/^M/", '', $array_peptides[$n]);
        }

        return array_unique($array_peptides);
    }

    public function getPeptidePosition($peptide, $protein)
    {
        preg_match_all('/[a-zA-Z]/', $protein, $protein_aa);
        preg_match_all('/[a-zA-Z]/', $peptide, $peptide_aa);

        return array_align($protein_aa[0], $peptide_aa[0]);
    }

    private function nonspecificPeptides($sequence)
    {
        preg_match_all("/./", $sequence, $aa);
        $aa = $aa[0];
        $l = count($aa);
        $n = 0;
        $array_peptides = array();
        for ($i = 0; $i < $l; $i ++) {
            if ($i > $l - $this->length_upper_limit)
                break;

            for ($j = $this->length_lower_limit; $j <= $this->length_upper_limit; $j ++) {

                $n ++;
                $concat_peptide = array_tostring(array_slice($aa, $i, $j), '', '');
                $leng = strlen($concat_peptide);
                $mass = $this->Peptide->getMolecularWeight($concat_peptide);

                if ($mass < $this->mass_lower_limit || $leng < $this->length_lower_limit)
                    continue;

                if ($mass > $this->mass_upper_limit || $leng > $this->length_upper_limit)
                    break;

                $array_peptides[] = $concat_peptide;
            }
        }

        return array_unique($array_peptides);
    }
}
?>