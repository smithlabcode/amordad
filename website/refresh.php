<?php
  $url = "http://localhost:18080/refresh";
  $file = file_get_contents($url);
  $result = json_decode($file, true);
  $response = "Successfully refreshed!\n"
            . "It takes ".$result['time']." seconds\n";
  echo $response;
?>
