<?php
    $dtime = date('Y-m-d-h-i-s', time());
    $perl_path = "/home/dell/miniconda3/envs/Perl/bin/perl";
    $analysis_path = "/home/dell/data/Web-Server/Uniq-mut/".$dtime."/";
    $uniq_pipeline = "/var/www/other_scripts/bio-pipelines/Uniq-mutation/uniq-mutation.pl";
    $Result_store_path = "/home/dell/data/Web-Server/Result-store/Uniq-mut/".$dtime."/";
    if($_FILES['isolate']['error'] == 0 && $_FILES['reference']['error'] == 0 && $_FILES['gff']['error'] == 0 && $_FILES['lineage']['error'] == 0){
        if($_FILES['isolate']['size'] > 104857600){
            die('The upload isolate genomes file should be smaller than 100Mb!<br>');
        }
        $iso_filename = $_FILES['isolate']['name'];
        $Ref_filename = $_FILES['reference']['name'];
        $Gff_filename = $_FILES['gff']['name'];
        $lin_filename = $_FILES['lineage']['name'];
		$iso_tmp = $_FILES['isolate']['tmp_name'];
        $Ref_tmp = $_FILES['reference']['tmp_name'];
        $Gff_tmp = $_FILES['gff']['tmp_name'];
        $lin_tmp = $_FILES['lineage']['tmp_name'];
	shell_exec("mkdir ".$analysis_path);
	shell_exec("mkdir ".$analysis_path.".../");
        shell_exec("chmod -R 777 ".$analysis_path);
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
                die('The reference genome should be .fasta/.fasta.gz/.fasta.zip format<br>');
            }
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
        if(move_uploaded_file($Gff_tmp, $analysis_path.$Gff_filename)){
            echo "$Gff_filename Reference sequence uploaded!<br>";
            shell_exec("chmod -R 777 ".$analysis_path);
        }else{
            shell_exec("rm -rf ".$analysis_path);
            die('The gff file upload error!<br>');
        }
        if(move_uploaded_file($lin_tmp, $analysis_path.$lin_filename)){
            echo "$lin_filename lineage information uploaded!<br>";
            shell_exec("chmod 777 -R ".$analysis_path);
            if($_FILES['lineage']['type'] == "application/octet-stream"){
                shell_exec("mv ".$analysis_path.$lin_filename." ".$analysis_path."lineage-info.tsv");
                $lin_filename = "lineage-info.tsv";
                shell_exec("chmod -R 777 ".$analysis_path);
            }else if($_FILES['lineage']['type'] == "application/gzip"){
                shell_exec("zcat ".$analysis_path.$lin_filename." > ".$analysis_path."lineage-info.tsv");
                shell_exec("rm ".$analysis_path.$lin_filename);
                $lin_filename = "lineage-info.tsv";
                shell_exec("chmod -R 777 ".$analysis_path);
            }else if($_FILES['lineage']['type'] == "application/zip" || $_FILES['lineage']['type'] == "application/x-zip-compressed"){
                shell_exec("unzip -o ".$analysis_path.$lin_filename." -d ".$analysis_path."lineage-info.tsv");
                shell_exec("rm ".$analysis_path.$lin_filename);
                $lin_filename = "lineage-info.tsv";
                shell_exec("chmod -R 777 ".$analysis_path);
            }else{
                shell_exec("rm -rf ".$analysis_path);
                die('The lineage information file shoud be .tsv/.tsv.gz/.tsv.zip format<br>');
            }
        }else{
            shell_exec("rm -rf ".$analysis_path);
            die('The lineage information file upload error!<br>');
        }
        $codon = $_POST['codon'];
        if(strcmp($_POST['discrim'], "default")==0){
            $percent = 0.8;
            $minNum = 10;
        }else if(strcmp($_POST['discrim'], "low")){
            $percent = 0.5;
            $minNum = 100;
        }else if(strcmp($_POST['discrim'], "high")){
            $percent = 0.7;
            $minNum = 5;
        }
        $command = $perl_path." ".$uniq_pipeline." -WorkPath ".$analysis_path." -reference ".$Ref_filename." -isolate ".$iso_filename." -gff ".$Gff_filename." -codon ".$codon." -p ".$percent." -minNum ".$minNum." -lineage ".$lin_filename." -StorePath ".$Result_store_path;
        shell_exec($command);
        shell_exec("chmod -R 777 ".$analysis_path);
        echo "Your Download code is: Uniq-mut_".$dtime."<br>";
        echo "Please waiting for the analysis, remember the download code and you can close the window!<br>";
        pclose(popen('bash '.$analysis_path.'.../shell.sh &', 'r'));
    }else{
        die('Files upload error!');
    }
?>
