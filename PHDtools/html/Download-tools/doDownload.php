<?php
    $filename = $_GET['filename'];
    header('content-disposition:attachment;filename='.basename($filename));
    readfile($filename);
?>
