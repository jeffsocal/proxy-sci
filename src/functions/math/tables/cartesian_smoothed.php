<?php

/*
 * Written by Jeff Jones (jeff@socalbioinformatics.com)
 * Copyright (2016) SoCal Bioinformatics Inc.
 *
 * See LICENSE.txt for the license.
 */

/*
 * returns the geometric mean
 */
function cartesian_smoothed($x, $y, $coef = 5)
{
    if (! is_numeric($coef))
        systemError("non numeric value for smoothing coef");
    
    /*
     * calculate number points
     */
    $n = count($x);
    
    /*
     * ensure both arrays of points are the same size
     */
    if ($n != count($y)) {
        systemError("smoothed_line(): Number of elements in coordinate arrays do not match.");
    }
    
    $coeh = floor($coef / 2);
    $coef = $coeh * 2 + 1;
    
    $table = array();
    $y = array_values($y);
    $x = array_values($x);
    $s = sizeof($y);
    for ($i = $coeh; $i < $s; $i ++) {
        
        $t_coef = min($coef, $i, $s - $i);
        $t_coeh = max(0, $i - $coeh);
        
        if ($i < $s - $coef)
            $i += $t_coef;
        
        $arr_x = array_slice($x, $t_coeh, $t_coef);
        $arr_y = array_slice($y, $t_coeh, $t_coef);
        
        $table['x'][] = array_mean($arr_x);
        $table['y'][] = array_mean($arr_y);
    }
    
    return $table;
}

// $val = $this->plotSmoothedValue;
// $pSize = sizeof($this->data[$val]);

// if($coef > $pSize/2) $coef = floor($pSize/2);

// $averageSize = $coef;
// $lookBack = ceil($averageSize/2);
// $pYStart = $this->yEnd - $this->plotXAxisMargin;
// $pXStart = $this->xStart + $this->plotYAxisMargin;
// $x = $this->data['x'][0];
// $y = $this->data[$val][0];

// if($this->xDif == 0 or $averageSize == 0) return;

// if ($this->yDif != 0) { $pYLast = $pYStart - (($y - $this->yMin) / $this->yDif) * $this->plotHeightSize; } else { $pYLast = $pYStart; }
// $pXLast = $pXStart + (($x - $this->xMin) / $this->xDif) * $this->plotWidthSize;

// for($get=-$lookBack;$get<($pSize+$lookBack);$get++){
// $population = array();
// for($s=-$lookBack;$s<$lookBack;$s++){
// $key = ($get+$s)."";
// if(!key_exists($key, $this->data['x'])) continue;
// if($key > $pSize-1) continue;
// if(key_exists('x', $population) and sizeof($population['x']) >= $coef) break;
// $population['x'][] = $this->data['x'][$key];
// $population['y'][] = $this->data[$val][$key];
// }
// if(!key_exists('x', $population)) continue;
// if(sizeof($population['x']) == 0) continue;
// $x = array_sum($population['x']) / sizeof($population['x']);
// $y = array_sum($population['y']) / sizeof($population['y']);

// if($x < $this->xMin) {
// $pXLast = $pXStart + ((0) / $this->xDif) * $this->plotWidthSize;
// continue;
// }
// if($x > $this->xMax) {
// $pXLast = $pXStart + ((0) / $this->xDif) * $this->plotWidthSize;
// continue;
// }

// if ($this->yDif != 0) { $pY = $pYStart - (($y - $this->yMin) / $this->yDif) * $this->plotHeightSize; } else { $pY = $pYStart; }
// $pX = $pXStart + (($x - $this->xMin) / $this->xDif) * $this->plotWidthSize;
// $xDistWidth = (((($this->data['x'][1]-$this->data['x'][0]) / $this->xDif) * $this->plotWidthSize)) / 2;
// for($l=0;$l<$this->lineWidth;$l++){
// imageline($this->image, $pX-$l, $pY, $pXLast-$l, $pYLast, $this->imagePointColor);
// imageline($this->image, $pX, $pY-$l, $pXLast, $pYLast-$l, $this->imagePointColor);
// if($bars == true) imagefilledrectangle($this->image, $pX+$xDistWidth, $pYStart, $pX-$xDistWidth, $pY, $this->imagePointColor);
// }
// $pXLast = $pX;
// $pYLast = $pY;
// }
?>
