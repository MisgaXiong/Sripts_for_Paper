<?php
    $visualization = "visibility:visible";
    $dtime = date('Y-m-d-h-i-s', time());
    $perl_path = "/home/dell/miniconda3/envs/Perl/bin/perl";
    $analysis_path = "/home/dell/data/Web-Server/AlignPrimers/".$dtime."/";
    $check_primers = "/var/www/other_scripts/bio-pipelines/Primer-align/Functions/check-primers.pl";
    $align_pipeline = "/var/www/other_scripts/bio-pipelines/Primer-align/primer-align.pl";
    $Result_store_path = "/home/dell/data/Web-Server/Result-store/AlignPrimers/".$dtime."/";
    if($_FILES['genomes']['error'] == 0){
        $iso_filename = $_FILES['genomes']['name'];
        $iso_tmp = $_FILES['genomes']['tmp_name'];
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
        $_POST['primers'] = str_replace(array("\r"), "", $_POST['primers']);
        $save_file = fopen($analysis_path."primers.fasta", 'w+');
        fwrite($save_file, $_POST['primers']);
        fclose($save_file);
        shell_exec("chmod -R 777 ".$analysis_path);
        $command = $perl_path." ".$check_primers." -WorkPath ".$analysis_path." -primers primers.fasta";
        $result = `$command`;
        if(strcmp($result, "FALSE")==0){
            shell_exec("rm -rf ".$analysis_path);
            die('The name of primers or probe must be "forward", "reverse" or "probe"!<br>');
        }
        if($_FILES['amplicon']['error'] == 0){
            $amp_filename = $_FILES['amplicon']['name'];
            $amp_tmp = $_FILES['amplicon']['tmp_name'];
            if(move_uploaded_file($amp_tmp, $analysis_path.$amp_filename)){
                echo "$amp_filename isolate genomes file uploaded!<br>";
                shell_exec("chmod 777 -R ".$analysis_path);
                if($_FILES['amplicon']['type'] == "application/octet-stream"){
                    shell_exec("mv ".$analysis_path.$amp_filename." ".$analysis_path."amplicon.fasta");
                    $amp_filename = "amplicon.fasta";
                    shell_exec("chmod -R 777 ".$analysis_path);
                }else if($_FILES['amplicon']['type'] == "application/gzip"){
                    shell_exec("zcat ".$analysis_path.$amp_filename." > ".$analysis_path."amplicon.fasta");
                    $amp_filename = "amplicon.fasta";
                    shell_exec("chmod -R 777 ".$analysis_path);
                }else if($_FILES['amplicon']['type'] == "application/zip" || $_FILES['amplicon']['type'] == "application/x-zip-compressed"){
                    shell_exec("unzip -o ".$analysis_path.$amp_filename." -d ".$analysis_path."amplicon.fasta");
                    $amp_filename = "amplicon.fasta";
                    shell_exec("chmod -R 777 ".$analysis_path);
                }else{
                    shell_exec("rm -rf ".$analysis_path);
                    die('The amplicon sequence shoud be .fasta/.fasta.gz/.fasta.zip format<br>');
                }
            }
		}
        if($_FILES['amplicon']['error'] == 0){
            $command = $perl_path." ".$align_pipeline." -WorkPath ".$analysis_path." -genomes ".$iso_filename." -primers primers.fasta -amplicon ".$amp_filename." -StorePath ".$Result_store_path;
        }else{
            $command = $perl_path." ".$align_pipeline." -WorkPath ".$analysis_path." -genomes ".$iso_filename." -primers primers.fasta"." -StorePath ".$Result_store_path;
        }
        shell_exec($command);
        shell_exec("chmod 777 ".$analysis_path.".../shell.sh");
        if($_POST['waiting'] == "simultaneous"){
            $condition = "TRUE";
            shell_exec("bash ".$analysis_path.".../shell.sh");
            echo "The analysis finished, the results can be download from the follows link:<br>";
            $download_path = "../Download-tools/doDownload.php?filename=".$Result_store_path."StrainLevel-PrimerAlignment.zip";
        }else{
            $visualization = "visibility:hidden";
            echo "Your Download code is: AlignPrimers_".$dtime."<br>";
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
    <a href="<?php echo $download_path?>" style="<?php echo $visualization?>">DownloadPrimersMappingFile</a><br>
</body>
</html>
