<?php
  $path = $_GET['path'];
  $url = "http://localhost:18080/delete?path=".$path;
  $file = file_get_contents($url);
  $result = json_decode($file, true);
  $response = "Successfully deleted!\n"
            . "Total number of samples now is ".number_format($result['total']).".\n"
            . "It takes ".$result['time']." seconds\n";
  echo $response;
?>
