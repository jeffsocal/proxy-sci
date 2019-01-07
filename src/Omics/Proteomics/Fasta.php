<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
namespace ProxySci\Omics\Proteomics;

use ProxyTime\ProgressTimer;

class Fasta
{

    private $fasta_file_path;

    private $regex_array;

    public function __construct($source = 'uniprot')
    {
        $ini = parse_ini_file(get_include_path() . 'ini/molecular.ini');
        $this->fasta_file_path = $ini['fasta_path'];
        
        $source_regex = 'uniprot';
        if (preg_match("/uniref/", $source))
            $source_regex = 'uniref';
        
        $this->regex_array = array();
        foreach ($ini as $key => $value) {
            if (strstr($key, $source_regex . '_regex_'))
                $this->regex_array[str_replace($source_regex . '_regex_', '', $key)] = $value;
        }
    }

    function getFastaPath()
    {
        return $this->fasta_file_path;
    }

    function parseFastaEntry($entry)
    {
        $array = array();
        
        foreach ($this->regex_array as $regex_var => $regex_val) {
            
            $array[$regex_var] = '';
            
            preg_match('/' . $regex_val . '/', $entry, $regex);
            
            if (key_exists(0, $regex))
                $array[$regex_var] = trim(preg_replace('/\w+\=|\n/', '', $regex[0]));
        }
        return $array;
    }

    function getProtein($uniprot_name = "ALBU_HUMAN")
    {
        $entry = FALSE;
        $continue = FALSE;
        $handle = fopen($this->fasta_file_path, "r");
        if ($handle) {
            while (($line = fgets($handle)) !== false) {
                // process the line read.
                
                if (preg_match('/^>.+' . $uniprot_name . '/', $line))
                    $continue = TRUE;
                
                if ($continue == TRUE) {
                    
                    if (preg_match('/^>/', $line)) {
                        
                        if ($entry != FALSE)
                            break;
                        
                        $entry = trim($line, '>');
                    } else {
                        $entry .= $line;
                    }
                }
            }
            
            return $this->parseFastaEntry($entry);
            
            fclose($handle);
        } else {
            // error opening the file.
        }
    }

    function getEntry($number = 100000)
    {
        $entry = FALSE;
        $continue = FALSE;
        
        $n = 0;
        
        $handle = fopen($this->fasta_file_path, "r");
        if ($handle) {
            while (($line = fgets($handle)) !== false) {
                // process the line read.
                
                if (preg_match('/^>/', $line)) {
                    $n ++;
                    if ($entry != FALSE & $n == $number)
                        break;
                    
                    $entry = trim($line, '>');
                } else {
                    $entry .= $line;
                }
            }
            return $this->parseFastaEntry($entry);
            
            fclose($handle);
        } else {
            // error opening the file.
        }
    }

    function getHomology($peptide = 'LAQFIKTASHTK')
    {
        $entry = FALSE;
        $continue = FALSE;
        
        $protein_hits = array();
        
        $handle = fopen($this->fasta_file_path, "r");
        if ($handle) {
            while (($line = fgets($handle)) !== false) {
                // process the line read.
                
                if (preg_match('/^>/', $line)) {
                    
                    $protein = $this->parseFastaEntry($entry);
                    if (preg_match('/' . $peptide . '/', $protein['seq']))
                        $protein_hits[] = $protein;
                    
                    $entry = trim($line, '>');
                } else {
                    $entry .= $line;
                }
            }
            $protein = $this->parseFastaEntry($entry);
            if (preg_match('/' . $peptide . '/', $protein['seq']))
                $protein_hits[] = $protein;
            
            fclose($handle);
        } else {
            // error opening the file.
        }
        return $protein_hits;
    }

    function getAllEntries($progress = FALSE)
    {
        $entry = FALSE;
        $parsed = [];
        
        if (! is_false($progress)) {
            $tmr = new ProgressTimer();
            
            $tmr->proTimerSize(556825);
            $tmr->proTimerStart("Parse FASTA");
        }
        
        $handle = fopen($this->fasta_file_path, "r");
        if ($handle) {
            while (($line = fgets($handle)) !== false) {
                /*
                 * process the line read.
                 */
                if (preg_match('/^>/', $line)) {
                    if (! is_false($entry)) {
                        if (! is_false($progress))
                            $tmr->proTimerPrint();
                        
                        $parsed[] = $this->parseFastaEntry($entry);
                    }
                    
                    $entry = trim($line, '>');
                } else {
                    $entry .= $line;
                }
            }
            /*
             * get the last one
             */
            $parsed[] = $this->parseFastaEntry($entry);
            
            fclose($handle);
        } else {
            // error opening the file.
        }
        
        return $parsed;
    }
}
?>