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

    function __construct()
    {
        parent::__construct();
    }

    function seqCoverage($protein, $peptides = array())
    {
        $proteinseq = $this->proteinseq;
        $this->proteinseqBitMap = preg_replace("/[A-Z]/", "0", $proteinseq);
        
        $arrPeptideseqs = $this->arrPeptideseqs;
        $i = sizeof($arrPeptideseqs);
        for ($n = 0; $n < $i; $n ++) {
            $peptide = $arrPeptideseqs[$n];
            if (stristr($proteinseq, $peptide)) {
                preg_match("/$peptide/", $proteinseq, $proteinseqMATCH, PREG_OFFSET_CAPTURE);
                for ($s = 0; $s < strlen($peptide); $s ++) {
                    $this->proteinseqBitMap[($proteinseqMATCH[0][1] + $s)] = 1;
                }
            }
        }
        $s = 0;
        $i = strlen($proteinseq);
        for ($n = 0; $n < $i; $n ++) {
            if ($this->proteinseqBitMap[$n] != 0)
                $s ++;
        }
        // $s = array_sum($proteinseqBitMap);
        // $this->seqCoverageBitMap = $proteinseqBitMap;
        return truncate(($s / $i) * 100, 2);
    }
}

?>