<!DOCTYPE html>
<html lang="en">
  <head>
    <!-- Meta, title, CSS, favicons, etc. -->
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <meta name="description" content="Amordad">
    <meta name="keywords" content="Metagenomics, alignment-free, HMP">
    <meta name="author" content="Wenzheng Li">

    <title>

      Amordad

    </title>

    <!-- Bootstrap core CSS -->
    <link href="css/bootstrap.css" rel="stylesheet">
    <link href="assets/css/src/docs.css" rel="stylesheet">


    <!-- Documentation extras -->
    <link href="css/index.css" rel="stylesheet">

    <link rel="icon" href="img/amordad.ico">

    <style>
      a.query_header {
        color: #993333;
        margin-right: 10px;
        margin-bottom:10px;
        font-size: 15px;
        background-color: #F4A460;
      }

      a.query_header:hover {
        background-color: #F4A460;
      }
      div {
        position: relative;
        top: 20px;
        width: 70%;
        margin-left: auto;
        margin-right: auto;
      }

      hr {
        border-top: dashed 1px #8c8b8b;
      }
    </style>
  </head>

<div class="query">
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
$query = "<p><a href=# class=\"btn btn-outline-inverse btn-lg query_header\">$id</a></p>";
echo "$query<hr>";
?>
</div>

<div class="results">
<?php
require 'mysql_login.php';
asort($result);
$count=0;
foreach ($result as $key => $value) {
  if($key != "id") {
    $count += 1;
    echo "No.$count $key<br>\n";
    $statement = "select url from sample where id=\"$key\"";
    $row = mysqli_fetch_array(mysqli_query($con, $statement));
    $url = $row["url"];
    echo "<a href=$url/>$url</a><br>\n<hr>";
  }
}
?>
</div>
</html>
