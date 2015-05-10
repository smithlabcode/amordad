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

      span.query_header {
        margin-right: 10px;
        margin-bottom:10px;
        font-size: 15px;
        background-color: #0066FF;
      }

      div {
        position: relative;
        width: 70%;
        margin-left: auto;
        margin-right: auto;
      }

      hr {
        border-top: dashed 1px #8c8b8b;
      }
    </style>
  </head>

  <body class="bs-docs-home">
    <header class="navbar navbar-static-top bs-docs-nav" id="top" role="banner">
      <div class="container">
        <div class="navbar-header">
          <button class="navbar-toggle" type="button" data-toggle="collapse" data-target=".bs-navbar-collapse">
            <span class="sr-only">Toggle navigation</span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
          </button>
          <a href="index.php" class="navbar-brand">Amordad</a>
        </div>
        <nav class="collapse navbar-collapse bs-navbar-collapse" role="navigation">
          <ul class="nav navbar-nav">
          </ul>
          <ul class="nav navbar-nav navbar-right">
            <li><a href="http://smithlabresearch.org">SmithLab</a></li>
            <li><a href="https://github.com/smithlabcode/amordad">Repo</a></li>
          </ul>
        </nav>
      </div>
    </header>

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
    $query = "<p><span class=\"btn btn-outline-inverse btn-lg query_header\">$id</span></p>";
    echo "$query";
    $num_results = count($result)-1;
    echo "$num_results results found<br>\n";
    ?>
    <p>Results below are sorted by Best Match</p>
    <hr>
    </div>
    <div class="results">
    <?php
    require 'mysql_login.php';
    asort($result);
    foreach ($result as $key => $value) {
      if($key != "id") {
        echo "$key<br>\n";
        echo "$value<br>\n";
        $statement = "select url from sample where id=\"$key\"";
        $row = mysqli_fetch_array(mysqli_query($con, $statement));
        $url = $row["url"];
        echo "<a href=$url/>$url</a><br>\n<hr>";
      }
    }
    ?>
    </div>
  </body>
</html>
