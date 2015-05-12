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

    <!-- HTML5 shim and Respond.js IE8 support of HTML5 elements and media queries -->
    <!--[if lt IE 9]>
    <script src="https://oss.maxcdn.com/html5shiv/3.7.2/html5shiv.min.js"></script>
    <script src="https://oss.maxcdn.com/respond/1.4.2/respond.min.js"></script>
    <![endif]-->

    <link rel="icon" href="img/amordad.ico">

    <style>
      div#refresh {
        width: 50%;
      }

      div#insert, #delete {
        width: 60%;
      }

      div#delete {
        margin-left: 30px;
      }
      button#refresh {
        margin-left: 30px;
      }

      table {
        margin-left: 30px;
      }

      td {
        padding: 30px;
      }

      caption {
        font-size: 20px;
        color: white;
        background-color: black;
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

    <main>
    <hr>
    <!-- refresh section -->
    <div id="refresh">
      <span style="color: blue; margin: 30px">One click to refresh the database</span> 
      <button type="button" id="refresh_btn" class="btn btn-default btn-success">Refresh</button>
    </div>
    <hr>
    <div id="insert">
    <span style="color: blue; margin: 30px;">One click to insert a metagenome</span><br>
    <?php
    $filename = '/db/.admin/new_samples.txt';
    $file = fopen($filename, "r");

    $table = '<table class="table table-striped table-bordered table-condensed">
                <caption> untracked metagenomes</caption>
                <thead>
                  <th style="text-align: center">path</th><th>insert</th>
                </thead>';
    while($line = fgets($file)) {
      $table .= '<tr><td style="text-align: left">'.$line.'</td>';
      $table .= '<td align=center><button>Insert</button></td></tr>';
      }
    $table .= '</table>';
    echo $table;
    ?>
    </div>
    <hr>
    <div id="delete">
    <span style="color: blue;">One click to delete a metagenome</span><br>
    Metagenome ID:<input type="text" style="margin-left: 10px; width: 40%">
    <button type="submit" style="margin-left: 30px">Delete Me</button>
    </div>
    </main>

    <div style="width: 100%; height:200px;"</div>
    <!-- Footer
    ================================================== -->
    <footer class="bs-docs-footer" role="contentinfo">
      <div class="container">
        @copy 2015 SmithLab
      </div>
    </footer>


    <!-- jQuery (necessary for Bootstrap's Javascript plugins) -->
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"></script>
    <!-- Include all compiled plugins (below), or include individual files as needed -->
    <script src="js/bootstrap.min.js"></script>
  </body>
</html>
