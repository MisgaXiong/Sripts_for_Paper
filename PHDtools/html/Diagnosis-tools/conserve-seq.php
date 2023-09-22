<?php
    $visualization = "visibility:visible";
    $dtime = date('Y-m-d-h-i-s', time());
    $perl_path = "/home/dell/miniconda3/envs/Perl/bin/perl";
    $analysis_path = "/home/dell/data/Web-Server/ConserveSeq/".$dtime."/";
    $cons_pipeline = "/var/www/other_scripts/bio-pipelines/Conserve-sequence/conserve-seq.pl";
    $Result_store_path = "/home/dell/data/Web-Server/Result-store/ConserveSeq/".$dtime."/";
    if($_FILES['genomes']['error'] == 0 && $_FILES['reference']['error'] == 0){
        if($_FILES['genome']['size'] > 104857600){
            die('The upload genome sequences file should be smaller than 100Mb!<br>');
        }
        $iso_filename = $_FILES['genomes']['name'];
        $Ref_filename = $_FILES['reference']['name'];
        $iso_tmp = $_FILES['genomes']['tmp_name'];
        $Ref_tmp = $_FILES['reference']['tmp_name'];
	shell_exec("mkdir ".$analysis_path);
	shell_exec("mkdir ".$analysis_path.".../");
        shell_exec("chmod 777 -R ".$analysis_path);
        if(move_uploaded_file($iso_tmp, $analysis_path.$iso_filename)){
            echo "$iso_filename isolate genomes file uploaded!<br>";
            shell_exec("chmod 777 -R ".$analysis_path);
            if($_FILES['genomes']['type'] == "application/octet-stream"){
                shell_exec("mv ".$analysis_path.$iso_filename." ".$analysis_path."isolate_genomes.fasta");
                $iso_filename = "isolate_genomes.fasta";
                shell_exec("chmod -R 777 ".$analysis_path);
            }else if($_FILES['genomes']['type'] == "application/gzip"){
                shell_exec("zcat ".$analysis_path.$iso_filename." > ".$analysis_path."isolate_genomes.fasta");
                $iso_filename = "isolate_genomes.fasta";
                shell_exec("chmod -R 777 ".$analysis_path);
            }else if($_FILES['genomes']['type'] == "application/zip" || $_FILES['genomes']['type'] == "application/x-zip-compressed"){
                shell_exec("unzip -o ".$analysis_path.$iso_filename." -d ".$analysis_path."isolate_genomes.fasta");
                $iso_filename = "isolate_genomes.fasta";
                shell_exec("chmod -R 777 ".$analysis_path);
            }else{
                shell_exec("rm -rf ".$analysis_path);
                die('The pathogen genomes shoud be .fasta/.fasta.gz/.fasta.zip format<br>');
            }
        }else{
            shell_exec("rm -rf ".$analysis_path);
            die('$iso_filename isolate genomes file upload filed!<br>');
        }
        if(move_uploaded_file($Ref_tmp, $analysis_path.$Ref_filename)){
            echo "$Ref_filename Reference sequence uploaded!<br>";
            shell_exec("chmod 777 -R ".$analysis_path);
            if($_FILES['reference']['type'] == "application/octet-stream"){
                shell_exec("mv ".$analysis_path.$Ref_filename." ".$analysis_path."Ref_genome.fasta");
                $Ref_filename = "Ref_genome.fasta";
                shell_exec("chmod -R 777 ".$analysis_path);
            }else if($_FILES['reference']['type'] == "application/gzip"){
                shell_exec("zcat ".$analysis_path.$Ref_filename." > ".$analysis_path."Ref_genome.fasta");
                shell_exec("rm ".$analysis_path.$Ref_filename);
                $Ref_filename = "Ref_genome.fasta";
                shell_exec("chmod -R 777 ".$analysis_path);
            }else if($_FILES['reference']['type'] == "application/zip" || $_FILES['reference']['type'] == "application/x-zip-compressed"){
                shell_exec("unzip -o ".$analysis_path.$Ref_filename." -d ".$analysis_path."Ref_genome.fasta");
                shell_exec("rm ".$analysis_path.$Ref_filename);
                $Ref_filename = "Ref_genome.fasta";
                shell_exec("chmod -R 777 ".$analysis_path);
            }else{
                shell_exec("rm -rf ".$analysis_path);
                die('The reference genome shoud be .fasta/.fasta.gz/.fasta.zip format<br>');
            }
        }else{
            shell_exec("rm -rf ".$analysis_path);
            die('$Ref_filename Reference sequence upload filed!<br>');
        }
        if(is_numeric($_POST['targetLen'])){
            $targetLen = $_POST['targetLen'];
            if($targetLen > 1000){
                shell_exec("rm -rf ".$analysis_path);
                die('The expected target length of conserved sequence must be shorter than 1000<br>');
            }
        }else{
            shell_exec("rm -rf ".$analysis_path);
            die('The input of the expected target length must be a numeric value!<br>');
		}
        $command = $perl_path." ".$cons_pipeline." -WorkPath ".$analysis_path." -reference ".$Ref_filename." -isolates ".$iso_filename." -algo ".$_POST['algo']." -targetLen ".$targetLen." -StorePath ".$Result_store_path;
        shell_exec($command);
        shell_exec("chmod 777 ".$analysis_path.".../shell.sh");
        if($_POST['waiting'] == "simultaneous"){
            $condition = "TRUE";
            shell_exec("bash ".$analysis_path.".../shell.sh");
            echo "The analysis finished, the results can be download from the follows link:<br>";
            $download_path = "../Download-tools/doDownload.php?filename=".$Result_store_path."MostConserve_info.zip";
        }else{
            $visualization = "visibility:hidden";
            echo "Your Download code is: ConserveSeq_".$dtime."<br>";
            echo "Please waiting for the analysis, remember the download code and you can close the window!<br>";
            pclose(popen('bash '.$analysis_path.'.../shell.sh &', 'r'));
        }
    }else{
        die('Files upload error!');
    }
?>

<html lang="en">
    <head>
        <title>The most conserved sequence identification</title>
            <meta charset="UTF-8"/>
    </head>
<body>
    <a href="<?php echo $download_path?>" style="<?php echo $visualization?>">DownloadConserveInfoFile</a><br>
</body>
</html>
