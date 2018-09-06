<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */
namespace ProxySci\MassSpec\Fragmentation;

use ProxyIO\File\Read;

class ReadMGF
{

    private $Read;

    private $data;

    public function __construct($file_path)
    {
        $Read = new Read($file_path);
        
        $this->data = $this->shardFile($Read->getContents());
    }

    private function shardFile($contents)
    {
        $contents = explode("END IONS", $contents);
        
        /*
         * last element is null by nature
         */
        array_pop($contents);
        
        foreach ($contents as $n => $c) {
            
            $e = explode("\n", $c);
            $i = preg_filter("/^[1-9]/", '$0', $e);
            $m = preg_filter("/^[A-Z]/", '$0', $e);
            
            $csv = function ($str) {
                return explode(" ", $str);
            };
            
            $meta = function ($str) {
                $m = explode("=", $str);
                return array(
                    $m[0] => $m[1]
                );
            };
            
            $contents[$n] = [];
            foreach (array_map($meta, $m) as $kvp) {
                $key = key($kvp);
                if ($key != 'BEGIN IONS')
                    $contents[$n][strtolower($key)] = $kvp[$key];
            }
            
            $spec = array_rowtocol(array_map($csv, $i));
            $spec['mz'] = $spec[0];
            $spec['int'] = $spec[1];
            unset($spec[0], $spec[1]);
            $contents[$n]['spec'] = $spec;
        }
        
        return $contents;
    }
    
    function getData(){
        return $this->data;
    }
    
}
?>