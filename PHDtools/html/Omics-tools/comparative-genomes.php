<?php
    $visualization = "visibility:hidden";
    $dtime = date('Y-m-d-h-i-s', time());
    $perl_path = "/home/dell/miniconda3/envs/Perl/bin/perl";
    $analysis_path = "/home/dell/data/Web-Server/Compare-genome/".$dtime."/";
    $cmp_pipeline = "/var/www/other_scripts/bio-pipelines/Compare-genome/compare-genome.pl";
    $Result_store_path = "/home/dell/data/Web-Server/Result-store/Compare-genome/".$dtime."/";
    if(empty($_POST['datatype'])){
        die('Please choose the Comparative type like Genome vs Genome... <br>');
    }else{
	shell_exec("mkdir ".$analysis_path);
	shell_exec("mkdir ".$analysis_path.".../");
	shell_exec("chmod -R 777 ".$analysis_path);
        if(strcmp($_POST['datatype'], "GmeVSNGS")==0){
            if($_FILES['reference']['error'] == 0 && $_FILES['FQ1']['error'] == 0 && $_FILES['FQ2']['error'] == 0){
                $Ref_filename = $_FILES['reference']['name'];
                $R1_filename = $_FILES['Fq1']['name'];
                $R2_filename = $_FILES['Fq2']['name'];
                $R1_tmp = $_FILES['Fq1']['tmp_name'];
                $R2_tmp = $_FILES['Fq2']['tmp_name'];
                $Ref_tmp = $_FILES['reference']['tmp_name'];
                if(move_uploaded_file($R1_tmp, $analysis_path.$R1_filename)){
                    echo "$R1_filename FASTQ R1 file uploaded!<br>";
                    shell_exec("chmod 777 -R ".$analysis_path);
                }else{
                    shell_exec("rm -rf ".$analysis_path);
                    die('$R1_filename FASTQ R1 file upload filed!<br>');
                }
                if(move_uploaded_file($R2_tmp, $analysis_path.$R2_filename)){
                    echo "$R2_filename FASTQ R2 file uploaded!<br>";
                    shell_exec("chmod 777 -R ".$analysis_path);
                }else{
                    shell_exec("rm -rf ".$analysis_path);
                    die('$R2_filename FASTQ R2 file upload filed!<br>');
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
                }
                if($_FILES['gff']['error'] == 0){
                    $Gff_filename = $_FILES['gff']['name'];
                    $Gff_tmp = $_FILES['gff']['tmp_name'];
                    if(move_uploaded_file($Gff_tmp, $analysis_path.$Gff_filename)){
                        echo "$Gff_filename Reference sequence uploaded!<br>";
                        $codon = $_POST['codon'];
                        $command = $perl_path." ".$cmp_pipeline." -WorkPath ".$analysis_path." -datatype sequencing -reference ".$Ref_filename." -FQ1 ".$R1_filename." -FQ2 ".$R2_filename." -gff ".$Gff_filename." -codon ".$codon." -StorePath ".$Result_store_path;
                        shell_exec($command);
                        shell_exec("chmod 777 ".$analysis_path.".../shell.sh");
                        echo "Your Download code is: Compare-genome_".$dtime."<br>";
                        echo "Please waiting for the analysis, remember the download code and you can close the window!<br>";
                        pclose(popen('bash '.$analysis_path.'.../shell.sh &', 'r'));
                    }else{
                        shell_exec("rm -rf ".$analysis_path);
                        die('The gff file upload error!<br>');
                    }
                }else{
                   $command = $perl_path." ".$cmp_pipeline." -WorkPath ".$analysis_path." -datatype sequencing -reference ".$Ref_filename." -FQ1 ".$R1_filename." -FQ2 ".$R2_filename." -StorePath ".$Result_store_path;
                   shell_exec($command);
                   shell_exec("chmod 777 ".$analysis_path.".../shell.sh");
                   echo "Your Download code is: Compare-genome_".$dtime."<br>";
                   echo "Please waiting for the analysis, remember the download code and you can close the window!<br>";
                   pclose(popen('bash '.$analysis_path.'.../shell.sh &', 'r'));
                }
            }else{
                shell_exec("rm -rf ".$analysis_path);
                die('The sequencing file or the reference file upload error!<br>');
            }
        }else{
            $visualization = "visibility:visible";
            if(strcmp($_POST['datatype'], "GmeVSGme")==0){
                if($_FILES['reference']['error'] == 0 && $_FILES['isolate']['error'] == 0){
                    $Ref_filename = $_FILES['reference']['name'];
                    $iso_filename = $_FILES['isolate']['name'];
                    $Ref_tmp = $_FILES['reference']['tmp_name'];
                    $iso_tmp = $_FILES['isolate']['tmp_name'];
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
						die('The reference file upload error!<br>');
                    }
                    if(move_uploaded_file($iso_tmp, $analysis_path.$iso_filename)){
                        echo "$iso_filename isolate sequence uploaded!<br>";
                        shell_exec("chmod 777 -R ".$analysis_path);
                        if($_FILES['isolate']['type'] == "application/octet-stream"){
                            shell_exec("mv ".$analysis_path.$iso_filename." ".$analysis_path."iso_genome.fasta");
                            $iso_filename = "iso_genome.fasta";
                            shell_exec("chmod -R 777 ".$analysis_path);
                        }else if($_FILES['isolate']['type'] == "application/gzip"){
                            shell_exec("zcat ".$analysis_path.$iso_filename." > ".$analysis_path."iso_genome.fasta");
                            shell_exec("rm ".$analysis_path.$iso_filename);
                            $iso_filename = "iso_genome.fasta";
                            shell_exec("chmod -R 777 ".$analysis_path);
                        }else if($_FILES['isolate']['type'] == "application/zip" || $_FILES['isolate']['type'] == "application/x-zip-compressed"){
                            shell_exec("unzip -o ".$analysis_path.$iso_filename." -d ".$analysis_path."iso_genome.fasta");
                            shell_exec("rm ".$analysis_path.$iso_filename);
                            $iso_filename = "iso_genome.fasta";
                            shell_exec("chmod -R 777 ".$analysis_path);
                        }else{
                            shell_exec("rm -rf ".$analysis_path);
                            die('The isolate genome shoud be .fasta/.fasta.gz/.fasta.zip format<br>');
                        }
                    }else{
                        shell_exec("rm -rf ".$analysis_path);
                        die('The isolate file upload error!<br>');
                    }
                    if($_FILES['gff']['error'] == 0){
                        $Gff_filename = $_FILES['gff']['name'];
                        $Gff_tmp = $_FILES['gff']['tmp_name'];
                        if(move_uploaded_file($Gff_tmp, $analysis_path.$Gff_filename)){
                            echo "$Gff_filename Reference sequence uploaded!<br>";
                            $codon = $_POST['codon'];
                            $command = $perl_path." ".$cmp_pipeline." -WorkPath ".$analysis_path." -datatype fragment --seqtype genome -reference ".$Ref_filename." -isolate ".$iso_filename." -gff ".$Gff_filename." -codon ".$codon." -StorePath ".$Result_store_path;
                            shell_exec($command);
                            shell_exec("chmod 777 ".$analysis_path.".../shell.sh");
                            shell_exec("bash ".$analysis_path.".../shell.sh");
                            echo "The analysis finished, the comparative results can be download from the follows link:<br>";
                            $download_path = "../Download-tools/doDownload.php?filename=".$Result_store_path."comparative-results.tsv";
                       }else{
                            shell_exec("rm -rf ".$analysis_path);
                            die('The gff file upload error!<br>');
                       }
                    }else{
                        $command = $perl_path." ".$cmp_pipeline." -WorkPath ".$analysis_path." -datatype fragment --seqtype genome -reference ".$Ref_filename." -isolate ".$iso_filename." -StorePath ".$Result_store_path;
                        shell_exec($command);
                        shell_exec("chmod 777 ".$analysis_path.".../shell.sh");
                        shell_exec("bash ".$analysis_path.".../shell.sh");
                        echo "The analysis finished, the comparative results can be download from the follows link:<br>";
                        $download_path = "../Download-tools/doDownload.php?filename=".$Result_store_path."comparative-results.tsv";
                    }
                }else{
                    shell_exec("rm -rf ".$analysis_path);
                    die('The reference or isolate file upload error!<br>');
                }
            }else if(strcmp($_POST['datatype'], "GneVSGne")==0){
                if($_FILES['reference']['error'] == 0 && $_FILES['isolate']['error'] == 0){
                    $Ref_filename = $_FILES['reference']['name'];
                    $iso_filename = $_FILES['isolate']['name'];
                    $Ref_tmp = $_FILES['reference']['tmp_name'];
                    $iso_tmp = $_FILES['isolate']['tmp_name'];
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
						die('The reference file upload error!<br>');
                    }
                    if(move_uploaded_file($iso_tmp, $analysis_path.$iso_filename)){
                        echo "$iso_filename isolate sequence uploaded!<br>";
                        shell_exec("chmod 777 -R ".$analysis_path);
                        if($_FILES['isolate']['type'] == "application/octet-stream"){
                            shell_exec("mv ".$analysis_path.$iso_filename." ".$analysis_path."iso_genome.fasta");
                            $iso_filename = "iso_genome.fasta";
                            shell_exec("chmod -R 777 ".$analysis_path);
                        }else if($_FILES['isolate']['type'] == "application/gzip"){
                            shell_exec("zcat ".$analysis_path.$iso_filename." > ".$analysis_path."iso_genome.fasta");
                            shell_exec("rm ".$analysis_path.$iso_filename);
                            $iso_filename = "iso_genome.fasta";
                            shell_exec("chmod -R 777 ".$analysis_path);
                        }else if($_FILES['isolate']['type'] == "application/zip" || $_FILES['isolate']['type'] == "application/x-zip-compressed"){
                            shell_exec("unzip -o ".$analysis_path.$iso_filename." -d ".$analysis_path."iso_genome.fasta");
                            shell_exec("rm ".$analysis_path.$iso_filename);
                            $iso_filename = "iso_genome.fasta";
                            shell_exec("chmod -R 777 ".$analysis_path);
                        }else{
                            shell_exec("rm -rf ".$analysis_path);
                            die('The isolate genome shoud be .fasta/.fasta.gz/.fasta.zip format<br>');
                        }
                    }else{
                        shell_exec("rm -rf ".$analysis_path);
                        die('The isolate file upload error!<br>');
                    }
                    $command = $perl_path." ".$cmp_pipeline." -WorkPath ".$analysis_path." -datatype fragment --seqtype gene -reference ".$Ref_filename." -isolate ".$iso_filename." -StorePath ".$Result_store_path;
                    shell_exec($command);
                    shell_exec("chmod 777 ".$analysis_path.".../shell.sh");
                    shell_exec("bash ".$analysis_path.".../shell.sh");
                    echo "The analysis finished, the comparative results can be download from the follows link:<br>";
                    $download_path = "../Download-tools/doDownload.php?filename=".$Result_store_path."comparative-results.tsv";
                }
            }else if(strcmp($_POST['datatype'], "PrtVSPrt")==0){
                if($_FILES['reference']['error'] == 0 && $_FILES['isolate']['error'] == 0){
                    $Ref_filename = $_FILES['reference']['name'];
                    $iso_filename = $_FILES['isolate']['name'];
                    $Ref_tmp = $_FILES['reference']['tmp_name'];
                    $iso_tmp = $_FILES['isolate']['tmp_name'];
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
						die('The reference file upload error!<br>');
                    }
                    if(move_uploaded_file($iso_tmp, $analysis_path.$iso_filename)){
                        echo "$iso_filename isolate sequence uploaded!<br>";
                        shell_exec("chmod 777 -R ".$analysis_path);
                        if($_FILES['isolate']['type'] == "application/octet-stream"){
                            shell_exec("mv ".$analysis_path.$iso_filename." ".$analysis_path."iso_genome.fasta");
                            $iso_filename = "iso_genome.fasta";
                            shell_exec("chmod -R 777 ".$analysis_path);
                        }else if($_FILES['isolate']['type'] == "application/gzip"){
                            shell_exec("zcat ".$analysis_path.$iso_filename." > ".$analysis_path."iso_genome.fasta");
                            shell_exec("rm ".$analysis_path.$iso_filename);
                            $iso_filename = "iso_genome.fasta";
                            shell_exec("chmod -R 777 ".$analysis_path);
                        }else if($_FILES['isolate']['type'] == "application/zip" || $_FILES['isolate']['type'] == "application/x-zip-compressed"){
                            shell_exec("unzip -o ".$analysis_path.$iso_filename." -d ".$analysis_path."iso_genome.fasta");
                            shell_exec("rm ".$analysis_path.$iso_filename);
                            $iso_filename = "iso_genome.fasta";
                            shell_exec("chmod -R 777 ".$analysis_path);
                        }else{
                            shell_exec("rm -rf ".$analysis_path);
                            die('The isolate genome shoud be .fasta/.fasta.gz/.fasta.zip format<br>');
                        }
                    }else{
                        shell_exec("rm -rf ".$analysis_path);
                        die('The isolate file upload error!<br>');
                    }
                    $command = $perl_path." ".$cmp_pipeline." -WorkPath ".$analysis_path." -datatype fragment --seqtype protein -reference ".$Ref_filename." -isolate ".$iso_filename." -StorePath ".$Result_store_path;
                    shell_exec($command);
                    shell_exec("chmod 777 ".$analysis_path.".../shell.sh");
                    shell_exec("bash ".$analysis_path.".../shell.sh");
                    echo "The analysis finished, the comparative results can be download from the follows link:<br>";
                    $download_path = "../Download-tools/doDownload.php?filename=".$Result_store_path."comparative-results.tsv";
                }
            }
        }
    }
?>

<html lang="en">
    <head>
        <title>Genome assembly</title>
        <meta charset="UTF-8"/>
    </head>
<body>
    <a href="<?php echo $download_path?>" style="<?php echo $visualization?>">DownloadComparativeFile</a><br>
</body>
</html>
