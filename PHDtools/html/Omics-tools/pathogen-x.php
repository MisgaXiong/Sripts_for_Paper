<?php
    $dtime = date('Y-m-d-h-i-s', time());
    $perl_path = "/home/dell/miniconda3/envs/Perl/bin/perl";
    $microbiome_pipeline = "/var/www/other_scripts/bio-pipelines/Microbiome-analysis/next-microbiome.pl";
    $Result_store_path = "/home/dell/data/Web-Server/Result-store/Microbiome/".$dtime."/";
    $host = $_POST['host'];
    $ctglen = $_POST['ctglen'];
    if(empty($host)){
        die('The host information must be chosen!<br>');
    }
    if($_FILES['Fq1']['error'] == 0 && $_FILES['Fq2']['error'] == 0){
        $R1_filename = $_FILES['Fq1']['name'];
        $R2_filename = $_FILES['Fq2']['name'];
        $R1_tmp = $_FILES['Fq1']['tmp_name'];
        $R2_tmp = $_FILES['Fq2']['tmp_name'];
        $analysis_path = "/home/dell/data/Web-Server/Microbiome/".$dtime."/";
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
		if($_FILES['bgm']['error'] == 0){
            $bgm_tmp = $_FILES['bgm']['tmp_name'];
            $bgm_filename = $_FILES['bgm']['name'];
            if(move_uploaded_file($bgm_tmp, $analysis_path.$bgm_filename)){
                if($_FILES['bgm']['type'] == "application/octet-stream"){
                    shell_exec("mv ".$analysis_path.$bgm_filename." ".$analysis_path."bgm_genome.fasta");
                    $bgm_filename = "bgm_genome.fasta";
                    shell_exec("chmod -R 777 ".$analysis_path);
                }else if($_FILES['bgm']['type'] == "application/gzip"){
                    shell_exec("zcat ".$analysis_path.$bgm_filename." > ".$analysis_path."bgm_genome.fasta");
                    shell_exec("rm -rf ".$analysis_path.$bgm_filename);
                    $bgm_filename = "bgm_genome.fasta";
                    shell_exec("chmod -R 777 ".$analysis_path);
                }else if($_FILES['bgm']['type'] == "application/zip" || $_FILES['bgm']['type'] == "application/x-zip-compressed"){
                    shell_exec("unzip -o ".$analysis_path.$bgm_filename." -d ".$analysis_path."bgm_genome.fasta");
                    shell_exec("rm -rf ".$analysis_path.$bgm_filename);
                    $bgm_filename = "bgm_genome.fasta";
                    shell_exec("chmod -R 777 ".$analysis_path);
                }else{
                    shell_exec("rm -rf ".$analysis_path);
                    die('The background microbiota genome shoud be .fasta/.fasta.gz/.fasta.zip format<br>');
                }
            }else{
                shell_exec("rm -rf ".$analysis_path);
                die('$bgm_filename background microbiota file upload filed!<br>');
            }
            $command = $perl_path." ".$microbiome_pipeline." -WorkPath ".$analysis_path." -FQ1 ".$R1_filename." -FQ2 ".$R2_filename." -Host ".$host." -BGM ".$bgm_filename." -contigLen ".$ctglen." -StorePath ".$Result_store_path;
        }else{
            $command = $perl_path." ".$microbiome_pipeline." -WorkPath ".$analysis_path." -FQ1 ".$R1_filename." -FQ2 ".$R2_filename." -Host ".$host." -contigLen ".$ctglen." -StorePath ".$Result_store_path;
        }
        shell_exec($command);
        shell_exec("chmod 777 ".$analysis_path."shell.sh");
		echo "Your Download code is: Microbiome_".$dtime."<br>";
        echo "Please waiting for the analysis, remember the download code and you can close the window!<br>";
        pclose(popen('bash '.$analysis_path.'.../shell.sh &', 'r'));
    }else{
        die('Files upload error!');
    }
?>
