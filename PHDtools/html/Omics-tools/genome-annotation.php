<?php
    $visualization = "visibility:visible";
    $dtime = date('Y-m-d-h-i-s', time());
    $perl_path = "/home/dell/miniconda3/envs/prokka/bin/perl";
    $prokka_path = "/home/dell/miniconda3/envs/prokka/bin/prokka";
    $analysis_path = "/home/dell/data/Web-Server/Genome-annotation/".$dtime."/";
    $tidy_script = "/var/www/other_scripts/bio-pipelines/Genome-annotation/tidy-gff.pl";
    $Result_store_path = "/home/dell/data/Web-Server/Result-store/Genome-annotation/".$dtime."/";
    if(empty($_POST['kingdom']) or empty($_POST['codon'])){
        die('The kingdom and codon table of this species must be chosen!<br>');
    }
    if($_FILES['genome']['error'] == 0){
        shell_exec("mkdir ".$analysis_path);
        shell_exec("chmod -R 777 ".$analysis_path);
        $gme_tmp = $_FILES['genome']['tmp_name'];
        $gme_filename = $_FILES['genome']['name'];
        if(move_uploaded_file($gme_tmp, $analysis_path.$gme_filename)){
            echo "$gme_filename genome sequence uploaded!<br>";
            if($_FILES['genome']['type'] == "application/octet-stream"){
                shell_exec("mv ".$analysis_path.$gme_filename." ".$analysis_path."isolate_genome.fasta");
                $gme_filename = "isolate_genome.fasta";
                shell_exec("chmod -R 777 ".$analysis_path);
            }else if($_FILES['genome']['type'] == "application/gzip"){
                shell_exec("zcat ".$analysis_path.$gme_filename." > ".$analysis_path."isolate_genome.fasta");
                shell_exec("rm -rf ".$analysis_path.$gme_filename);
                $gme_filename = "isolate_genome.fasta";
                shell_exec("chmod -R 777 ".$analysis_path);
            }else if($_FILES['genome']['type'] == "application/zip" || $_FILES['genome']['type'] == "application/x-zip-compressed"){
                shell_exec("unzip -o ".$analysis_path.$gme_filename." -d ".$analysis_path."isolate_genome.fasta");
                shell_exec("rm -rf ".$analysis_path.$gme_filename);
                $gme_filename = "isolate_genome.fasta";
                shell_exec("chmod -R 777 ".$analysis_path);
            }else{
                shell_exec("rm -rf ".$analysis_path);
                die('The genome shoud be .fasta/.fasta.gz/.fasta.zip format<br>');
            }
        }else{
            shell_exec("rm -rf ".$analysis_path);
            die('The genome file upload error!<br>');
        }
        $kingdom = $_POST['kingdom'];
        $codon_tbl = $_POST['codon'];
        $command_1 = $perl_path." ".$prokka_path." ".$analysis_path.$gme_filename." --outdir ".$analysis_path."annotation --addgenes --prefix genome-annotation --kingdom ".$kingdom." --gcode ".$codon_tbl." &> ".$analysis_path."run_prokka.log";
        shell_exec($command_1);
        shell_exec("chmod -R 777 ".$analysis_path);
        $command_2 = $perl_path." ".$tidy_script." -WorkPath ".$analysis_path."annotation/ -gff genome-annotation.gff";
        shell_exec($command_2);
        shell_exec("chmod -R 777 ".$analysis_path);
        $command_3 = "zip -qjr ".$analysis_path."genome-annotations.zip ".$analysis_path."annotation/genome-annotation.ffn ".$analysis_path."annotation/genome-annotation.faa ".$analysis_path."annotation/genome-annotation.gff3";
        shell_exec($command_3);
        shell_exec("chmod -R 777 ".$analysis_path);
        $command_4 = "mkdir ".$Result_store_path;
        shell_exec($command_4);
        shell_exec("chmod -R 777 ".$Result_store_path);
        $command_5 = "mv ".$analysis_path."genome-annotations.zip ".$Result_store_path;
        shell_exec($command_5);
        $command_6 = "chmod -R 777 ".$Result_store_path;
        shell_exec($command_6);
        $command_7 = "rm -rf ".$analysis_path;
        shell_exec($command_7);
        echo "The analysis finished, the annotation results can be download from the follows link:<br>";
        $download_path = "../Download-tools/doDownload.php?filename=".$Result_store_path."genome-annotations.zip";
    }else{
        die('The genome file upload error!<br>');
    }
?>

<html lang="en">
    <head>
        <title>Genome Annotation</title>
        <meta charset="UTF-8"/>
    </head>
<body>
    <a href="<?php echo $download_path?>" style="<?php echo $visualization?>">DownloadAnnotationsFile</a><br>
</body>
</html>
