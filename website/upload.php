<?php
ob_start();
error_reporting(E_ALL);
ini_set('display_errors', TRUE);
$upload_location = '/db/.upload/';
$allowedExts = array("cv");
$uniq_id = uniqid('amordad_upload_');
$sample = $upload_location.$uniq_id.".cv";
move_uploaded_file($_FILES["fileToUpload"]["tmp_name"], $sample);
$url = "http://localhost:18080/query?path=".urlencode($sample);
$file = file_get_contents($url);
$result = json_decode($file, true);
$id = $result["id"];
echo "******  $id  ************<br/>\n";
foreach ($result as $key => $value) {
  if($key != "id")
    echo "$key  ==================>  $value<br/>\n";
}
?>
