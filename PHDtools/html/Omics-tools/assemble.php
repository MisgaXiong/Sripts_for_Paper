<?php
    $dtime = date('Y-m-d-h-i-s', time());
    $perl_path = "/home/dell/miniconda3/envs/Perl/bin/perl";
    $analysis_path = "/home/dell/data/Web-Server/Genome-assemble/".$dtime."/";
    $assemble_pipeline = "/var/www/other_scripts/bio-pipelines/Genome-assembly/genome-assembly.pl";
    $Result_store_path = "/home/dell/data/Web-Server/Result-store/Genome-assemble/".$dtime."/";
    if($_FILES['Fq1']['error'] == 0 && $_FILES['Fq2']['error'] == 0 && $_FILES['reference']['error'] == 0){
        $R1_filename = $_FILES['Fq1']['name'];
        $R2_filename = $_FILES['Fq2']['name'];
        $Ref_filename = $_FILES['reference']['name'];
        $R1_tmp = $_FILES['Fq1']['tmp_name'];
        $R2_tmp = $_FILES['Fq2']['tmp_name'];
        $Ref_tmp = $_FILES['reference']['tmp_name'];
	shell_exec("mkdir ".$analysis_path);
	shell_exec("mkdir ".$analysis_path.".../");
        shell_exec("chmod 777 -R ".$analysis_path);
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
                shell_exec("rm -rf ".$analysis_path.$Ref_filename);
                $Ref_filename = "Ref_genome.fasta";
                shell_exec("chmod -R 777 ".$analysis_path);
            }else if($_FILES['reference']['type'] == "application/zip" || $_FILES['reference']['type'] == "application/x-zip-compressed"){
                shell_exec("unzip -o ".$analysis_path.$Ref_filename." -d ".$analysis_path."Ref_genome.fasta");
                shell_exec("rm -rf ".$analysis_path.$Ref_filename);
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
        $taxid = str_replace(array("\r\n", "\s", "\s+", "\r", "\n"), "", $_POST['taxid']);
        if(!empty($taxid)){
            $command = $perl_path." ".$assemble_pipeline." -WorkPath ".$analysis_path." -FQ1 ".$R1_filename." -FQ2 ".$R2_filename." -Ref ".$Ref_filename." -taxid ".$taxid." -StorePath ".$Result_store_path;
        }else{
            $command = $perl_path." ".$assemble_pipeline." -WorkPath ".$analysis_path." -FQ1 ".$R1_filename." -FQ2 ".$R2_filename." -Ref ".$Ref_filename." -StorePath ".$Result_store_path;
        }
        shell_exec($command);
        shell_exec("chmod 777 ".$analysis_path.".../shell.sh");
        echo "Your Download code is: Genome-assemble_".$dtime."<br>";
        echo "Please waiting for the analysis, remember the download code and you can close the window!<br>";
        pclose(popen('bash '.$analysis_path.'.../shell.sh &', 'r'));
    }else{
        die('Files upload error!');
    }
?>
