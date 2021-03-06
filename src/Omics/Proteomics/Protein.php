<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
namespace ProxySci\Omics\Proteomics;

class Protein extends Peptide
{

    private $proteinBitMap;

    private $proteinCoverage;

    private $proteinSeqHTML;

    private $proteinSeqSTR;

    function __construct()
    {
        parent::__construct();
    }

    function calculateCoverage(string $protein, array $peptides = [])
    {
        $this->proteinSeqBitMap = array_fill(0, strlen($protein), 0);
        $peptides = array_values(array_unique($peptides));
        
        $peptides = preg_replace("/I/", "L", $peptides);
        $protein = preg_replace("/I/", "L", $protein);
        
        $i = sizeof($peptides);
        for ($n = 0; $n < $i; $n ++) {
            $peptide = $peptides[$n];
            if (stristr($protein, $peptide)) {
                preg_match("/$peptide/", $protein, $proteinMATCH, PREG_OFFSET_CAPTURE);
                for ($s = 0; $s < strlen($peptide); $s ++) {
                    $this->proteinSeqBitMap[($proteinMATCH[0][1] + $s)] ++;
                }
            }
        }
        $s = 0;
        $i = strlen($protein);
        for ($n = 0; $n < $i; $n ++) {
            if ($this->proteinSeqBitMap[$n] != 0)
                $s ++;
        }
        $this->proteinSeqCoverage = ($s / $i);
        
        return ($this->proteinSeqCoverage);
    }

    private function htmlCoverageMap($protein)
    {
        preg_match_all("/[A-Z]/", $protein, $proteinBit);
        
        $colors = [
            '#F0F0F0',
            '#FFAA00',
            '#FF9900',
            '#FF8800',
            '#FF7700',
            '#FF6600',
            '#FF5500',
            '#FF4400',
            '#FF3300',
            '#FF2200',
            '#FF1100',
            '#FF0000'
        ];
        
        $max_hits = array_max($this->proteinSeqBitMap);
        
        $this->proteinSeqHTML = '';
        foreach ($proteinBit[0] as $n => $aa) {
            
            $v = round(($this->proteinSeqBitMap[$n] / $max_hits) * (count($colors) - 1));
            $aa = '<a style="color:' . $colors[$v] . ';">' . $aa . '</a>';
            
            $this->proteinSeqHTML .= $aa;
            
            if (($n + 1) % 60 === 0)
                $this->proteinSeqHTML .= "<br>\n";
        }
    }

    private function binaryCoverageMap($protein)
    {
        preg_match_all("/[a-z]/", strtolower($protein), $proteinBit);
        
        $max_hits = array_max($this->proteinSeqBitMap);
        
        $this->proteinSeqSTR = '';
        foreach ($proteinBit[0] as $n => $aa) {
            
            if ($this->proteinSeqBitMap[$n] > 0)
                $aa = strtoupper($aa);
            
            $this->proteinSeqSTR .= $aa;
        }
    }

    function getCoverage(string $protein, array $peptides = [], bool $as_html = FALSE)
    {
        $this->calculateCoverage($protein, $peptides);
        
        $map = '';
        if ($as_html == TRUE) {
            $this->htmlCoverageMap($protein);
            $map = '<div class="sequence">' . $this->proteinSeqHTML . "</div>";
        } else {
            $this->binaryCoverageMap($protein);
            $map = $this->proteinSeqSTR;
        }
        
        $out['coverage'] = truncate($this->proteinSeqCoverage * 100, 2);
        $out['map'] = $map;
        
        return $out;
    }
}

?>