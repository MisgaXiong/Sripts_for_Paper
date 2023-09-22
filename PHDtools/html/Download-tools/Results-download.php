<?php
    $dtime = date('Y-m-d-h-i-s', time());
    $perl_path = "/home/dell/miniconda3/envs/Perl/bin/perl";
    $delet_date = "/var/www/other_scripts/bio-pipelines/deletdata/deletdata.pl";
    $command = $perl_path." ".$delet_date." -nowdate ".$dtime;
    //shell_exec($command);
    $results_path = "/home/dell/data/Web-Server/Result-store/";
    $_POST['download_code'] = str_replace(array("\r\n", "\r", "\n", "\s"), "", $_POST['download_code']);
    $arr = explode('_', $_POST['download_code']);
	if(strcmp($arr[0], "Genome-assemble")==0){
        if(!file_exists($results_path.$arr[0]."/".$arr[1]."/Final_Assembly.fasta")){
            die('No such file, please check the download code, or the analysis is still under working, or the file was deleted because of the long time gone.<br>');
        }
        $download_path = "doDownload.php?filename=".$results_path.$arr[0]."/".$arr[1]."/Final_Assembly.fasta";
    }else if(strcmp($arr[0], "Microbiome")==0){
        if(!file_exists($results_path.$arr[0]."/".$arr[1]."/Microbiota_Results.zip")){
            die('No such file, please check the download code, or the analysis is still under working, or the file was deleted because of the long time gone.<br>');
        }
        $download_path = "doDownload.php?filename=".$results_path.$arr[0]."/".$arr[1]."/Microbiota_Results.zip";
    }else if(strcmp($arr[0], "Strain-typing")==0){
        if(!file_exists($results_path.$arr[0]."/".$arr[1]."/StrainTyping.tree")){
            die('No such file, please check the download code, or the analysis is still under working, or the file was deleted because of the long time gone.<br>');
        }
        $download_path = "doDownload.php?filename=".$results_path.$arr[0]."/".$arr[1]."/StrainTyping.tree";
    }else if(strcmp($arr[0], "MinorAllele")==0){
        if(!file_exists($results_path.$arr[0]."/".$arr[1]."/iSNV_Results.zip")){
            die('No such file, please check the download code, or the analysis is still under working, or the file was deleted because of the long time gone.<br>');
        }
        $download_path = "doDownload.php?filename=".$results_path.$arr[0]."/".$arr[1]."/iSNV_Results.zip";
    }else if(strcmp($arr[0], "WGS-quality")==0){
        if(!file_exists($results_path.$arr[0]."/".$arr[1]."/Quality_Reports.zip")){
            die('No such file, please check the download code, or the analysis is still under working, or the file was deleted because of the long time gone.<br>');
        }
        $download_path = "doDownload.php?filename=".$results_path.$arr[0]."/".$arr[1]."/Quality_Reports.zip";
    }else if(strcmp($arr[0], "ConserveSeq")==0){
        if(!file_exists($results_path.$arr[0]."/".$arr[1]."/MostConserve_info.zip")){
            die('No such file, please check the download code, or the analysis is still under working, or the file was deleted because of the long time gone.<br>');
        }
        $download_path = "doDownload.php?filename=".$results_path.$arr[0]."/".$arr[1]."/MostConserve_info.zip";
    }else if(strcmp($arr[0], "AlignPrimers")==0){
        if(!file_exists($results_path.$arr[0]."/".$arr[1]."/StrainLevel-PrimerAlignment.zip")){
            die('No such file, please check the download code, or the analysis is still under working, or the file was deleted because of the long time gone.<br>');
        }
        $download_path = "doDownload.php?filename=".$results_path.$arr[0]."/".$arr[1]."/StrainLevel-PrimerAlignment.zip";
	}else if(strcmp($arr[0], "Compare-genome")==0){
        if(!file_exists($results_path.$arr[0]."/".$arr[1]."/comparative-results.tsv")){
            die('No such file, please check the download code, or the analysis is still under working, or the file was deleted because of the long time gone.<br>');
        }
        $download_path = "doDownload.php?filename=".$results_path.$arr[0]."/".$arr[1]."/comparative-results.tsv";
    }else if(strcmp($arr[0], "Uniq-mut")==0){
        if(!file_exists($results_path.$arr[0]."/".$arr[1]."/lineage-mutation-network.tsv")){
            die('No such file, please check the download code, or the analysis is still under working, or the file was deleted because of the long time gone.<br>');
        }
		$download_path = "doDownload.php?filename=".$results_path.$arr[0]."/".$arr[1]."/lineage-mutation-network.tsv";
    }else if(strcmp($arr[0], "Primer-design")==0){
        if(!file_exists($results_path.$arr[0]."/".$arr[1]."/Primers.tsv")){
            die('No such file, please check the download code, or the analysis is still under working, or the file was deleted because of the long time gone.<br>');
        }
        $download_path = "doDownload.php?filename=".$results_path.$arr[0]."/".$arr[1]."/Primers.tsv";
    }
?>

<html lang="en">
    <head>
        <title>Download</title>
        <meta charset="UTF-8"/>
    </head>
    <body>
        <a href="<?php echo $download_path?>">Results download address</a><br>
    </body>
</html>
