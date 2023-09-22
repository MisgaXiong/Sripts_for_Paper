<?php
    $dtime = date('Y-m-d-h-i-s', time());
    $perl_path = "/home/dell/miniconda3/envs/Perl/bin/perl";
    $iSNV_pipeline = "/var/www/other_scripts/bio-pipelines/Minor-allele/minor-allele.pl";
    $analysis_path = "/home/dell/data/Web-Server/MinorAllele/".$dtime."/";
    $Result_store_path = "/home/dell/data/Web-Server/Result-store/MinorAllele/".$dtime."/";
    if($_FILES['Fq1']['error'] == 0 && $_FILES['Fq2']['error'] == 0 && $_FILES['isolate']['error'] == 0){
        $R1_filename = $_FILES['Fq1']['name'];
        $R2_filename = $_FILES['Fq2']['name'];
        $iso_filename = $_FILES['isolate']['name'];
        $R1_tmp = $_FILES['Fq1']['tmp_name'];
        $R2_tmp = $_FILES['Fq2']['tmp_name'];
        $iso_tmp = $_FILES['isolate']['tmp_name'];
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
        if(move_uploaded_file($iso_tmp, $analysis_path.$iso_filename)){
            echo "$iso_filename isolate genome uploaded!<br>";
            shell_exec("chmod 777 -R ".$analysis_path);
            if($_FILES['isolate']['type'] == "application/octet-stream"){
                shell_exec("mv ".$analysis_path.$iso_filename." ".$analysis_path."Ref_genome.fasta");
                $iso_filename = "Ref_genome.fasta";
                shell_exec("chmod -R 777 ".$analysis_path);
            }else if($_FILES['isolate']['type'] == "application/gzip"){
                shell_exec("zcat ".$analysis_path.$iso_filename." > ".$analysis_path."Ref_genome.fasta");
                shell_exec("rm ".$analysis_path.$iso_filename);
                $iso_filename = "Ref_genome.fasta";
                shell_exec("chmod -R 777 ".$analysis_path);
            }else if($_FILES['isolate']['type'] == "application/zip" || $_FILES['isolate']['type'] == "application/x-zip-compressed"){
                shell_exec("unzip -o ".$analysis_path.$iso_filename." -d ".$analysis_path."Ref_genome.fasta");
                shell_exec("rm ".$analysis_path.$iso_filename);
                $iso_filename = "Ref_genome.fasta";
                shell_exec("chmod -R 777 ".$analysis_path);
            }else{
                shell_exec("rm -rf ".$analysis_path);
                die('The isolate genome shoud be .fasta/.fasta.gz/.fasta.zip format<br>');
            }
        }else{
            shell_exec("rm -rf ".$analysis_path);
            die('$iso_filename Reference sequence upload filed!<br>');
        }
        
    }else{
        die('The pair-end of FASTQ format file and isolate genome must be correctly uploaded!<br>');
    }
    $MuAF = $_POST['MuAF'];
    if($MuAF < 0 || $MuAF > 0.5){
        shell_exec("rm -rf ".$analysis_path);
        die('The minor allele frequency should be larger than 0 and smaller than 0.5!<br>');
    }else if(is_numeric($MuAF)){
        echo "Minor allele frequency chosen: ".$MuAF."<br>";
    }else{
        shell_exec("rm -rf ".$analysis_path);
        die('The minor allele frequency should be a float value!<br>');
    }
    $depth = $_POST['depth'];
    if(is_numeric($depth)){
        echo "Sequencing depth chosen: ".$depth."<br>";
	}else{
        shell_exec("rm -rf ".$analysis_path);
        die('The sequencing depth should be a numeric value!<br>');
    }
    if(empty($_FILES['GFF']['tmp_name']) == "TRUE"){
	    if($_POST['AliQ'] == "default"){
            $command = $perl_path." ".$iSNV_pipeline." -WorkPath ".$analysis_path." -FQ1 ".$R1_filename." -FQ2 ".$R2_filename." -REF ".$iso_filename." -FLAG 4 -MAF ".$MuAF." -MINDEP ".$depth." -StorePath ".$Result_store_path;
        }else if($_POST['AliQ'] == "medium"){
            $command = $perl_path." ".$iSNV_pipeline." -WorkPath ".$analysis_path." -FQ1 ".$R1_filename." -FQ2 ".$R2_filename." -REF ".$iso_filename." -FLAG 3840 -Q 10 -MAF ".$MuAF." -MINDEP ".$depth." -StorePath ".$Result_store_path;
        }else if($_POST['AliQ'] == "high"){
            $command = $perl_path." ".$iSNV_pipeline." -WorkPath ".$analysis_path." -FQ1 ".$R1_filename." -FQ2 ".$R2_filename." -REF ".$iso_filename." -FLAG 3840 -Q 30 -MAF ".$MuAF." -MINDEP ".$depth." -StorePath ".$Result_store_path;
        }else if($_POST['AliQ'] == "strict"){
            $command = $perl_path." ".$iSNV_pipeline." -WorkPath ".$analysis_path." -FQ1 ".$R1_filename." -FQ2 ".$R2_filename." -REF ".$iso_filename." -FLAG 3840 -Q 60 -MAF ".$MuAF." -MINDEP ".$depth." -StorePath ".$Result_store_path;
        }
        shell_exec($command);
        shell_exec("chmod 777 ".$analysis_path.".../shell.sh");
        echo "Your Download code is: MinorAllele_".$dtime."<br>";
        echo "Please waiting for the analysis, remember the download code and you can close the window!<br>";
        pclose(popen('bash '.$analysis_path.'.../shell.sh &', 'r'));
    }else{
        $GFF_filename = $_FILES['GFF']['name'];
        $GFF_tmp = $_FILES['GFF']['tmp_name'];
        if(move_uploaded_file($GFF_tmp, $analysis_path.$GFF_filename)){
            echo "$GFF_filename gene annotation file uploaded!<br>";
            shell_exec("chmod 777 -R ".$analysis_path);
        }else{
            die('$GFF_filename gene annotation file upload filed!<br>');
        }
        $code = $_POST['code'];
        if($_POST['AliQ'] == "default"){
            $command = $perl_path." ".$iSNV_pipeline." -WorkPath ".$analysis_path." -FQ1 ".$R1_filename." -FQ2 ".$R2_filename." -REF ".$iso_filename." -FLAG 4 -MAF ".$MuAF." -MINDEP ".$depth." -GFF ".$GFF_filename." -CODON ".$code." -StorePath ".$Result_store_path;
        }else if($_POST['AliQ'] == "medium"){
            $command = $perl_path." ".$iSNV_pipeline." -WorkPath ".$analysis_path." -FQ1 ".$R1_filename." -FQ2 ".$R2_filename." -REF ".$iso_filename." -FLAG 3840 -Q 10 -MAF ".$MuAF." -MINDEP ".$depth." -GFF ".$GFF_filename." -CODON ".$code." -StorePath ".$Result_store_path;
        }else if($_POST['AliQ'] == "high"){
            $command = $perl_path." ".$iSNV_pipeline." -WorkPath ".$analysis_path." -FQ1 ".$R1_filename." -FQ2 ".$R2_filename." -REF ".$iso_filename." -FLAG 3840 -Q 30 -MAF ".$MuAF." -MINDEP ".$depth." -GFF ".$GFF_filename." -CODON ".$code." -StorePath ".$Result_store_path;
        }else if($_POST['AliQ'] == "strict"){
            $command = $perl_path." ".$iSNV_pipeline." -WorkPath ".$analysis_path." -FQ1 ".$R1_filename." -FQ2 ".$R2_filename." -REF ".$iso_filename." -FLAG 3840 -Q 60 -MAF ".$MuAF." -MINDEP ".$depth." -GFF ".$GFF_filename." -CODON ".$code." -StorePath ".$Result_store_path;
        }
        shell_exec($command);
        shell_exec("chmod 777 ".$analysis_path.".../shell.sh");
        echo "Your Download code is: MinorAllele_".$dtime."<br>";
        echo "Please waiting for the analysis, remember the download code and you can close the window!<br>";
        pclose(popen('bash '.$analysis_path.'.../shell.sh &', 'r'));
    }
?>
