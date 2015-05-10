<?php
  $ini_array = parse_ini_file("/db/.config.ini");
  $con =mysqli_connect($ini_array['host'], 
                       $ini_array['usrname'], 
                       $ini_array['passwd'], 
                       $ini_array['database']);
  if(!$con) {
    die('Could not connect: '.mysql_error());
  }
?> 
