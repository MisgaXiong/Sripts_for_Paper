<?php
    $dtime = date('Y-m-d-h-i-s', time());
    $perl_path = "/home/dell/miniconda3/envs/Perl/bin/perl";
    $analysis_path = "/home/dell/data/Web-Server/Strain-typing/".$dtime."/";
    $typing_pipeline = "/var/www/other_scripts/bio-pipelines/Strain-typing/strain-typing.pl";
    $Result_store_path = "/home/dell/data/Web-Server/Result-store/Strain-typing/".$dtime."/";
    if(!empty($_FILES['large-genomes'])){
        if($_FILES['large-genome']['size'] > 68157440){
            die('Due to the limitation of computing resource, the upload file size should be smaller than 65Mb!<br>');
        }
        if($_FILES['large-genomes']['error'] == 0){
            if(empty($_POST['resolution']) || strcmp($_POST['resolution'], "default")==0){
                $allele_len = 200;
                $allele_p = 0.95;
            }else if(strcmp($_POST['resolution'], "low")==0){
                $allele_len = 250;
                $allele_p = 0.95;
            }else if(strcmp($_POST['resolution'], "high")==0){
                $allele_len = 200;
                $allele_p = 0.9;
            }
            if(empty($_POST['codon'])){
                $codon_tbl = 1;
            }else{
                $codon_tbl = $_POST['codon'];
            }
            $gme_end = end(explode('.', $_FILES['large-genomes']['name']));
            $gme_tmp = $_FILES['large-genomes']['tmp_name'];
            if(strcmp($gme_end, "gz")==0 || strcmp($gme_end, "zip")==0){
		shell_exec("mkdir ".$analysis_path);
                shell_exec("mkdir ".$analysis_path.".../");
                shell_exec("chmod 777 -R ".$analysis_path);
                $gme_filename = "genomes.fa.".$gme_end;
                if(move_uploaded_file($gme_tmp, $analysis_path.$gme_filename)){
			shell_exec("chmod 777 -R ".$analysis_path);
                    if($_FILES['large-genomes']['type'] == "application/zip" || $_FILES['large-genomes']['type'] == "application/x-zip-compressed"){
                        shell_exec("unzip -o ".$analysis_path.$gme_filename." -d ".$analysis_path."genomes.fa");
                        shell_exec("chmod 777 -R ".$analysis_path." && chmod 777 -R ".$analysis_path."genomes.fa");
                    }else if($_FILES['large-genomes']['type'] == "application/gzip"){
                        shell_exec("mkdir ".$analysis_path."genomes.fa");
                        shell_exec("chmod 777 -R ".$analysis_path." && chmod 777 -R ".$analysis_path."genomes.fa");
                        shell_exec("tar -zxvf ".$analysis_path.$gme_filename." --strip-components 1 -C ".$analysis_path."genomes.fa");
                        shell_exec("chmod 777 -R ".$analysis_path." && chmod 777 -R ".$analysis_path."genomes.fa");
                    }else{
                        shell_exec("rm -rf ".$analysis_path);
                        die('Please upload a correct .zip or .gz format for genomes file!<br>');
                    }
                }else{
                    shell_exec("rm -rf ".$analysis_path);
                    die('Genomes file upload filed!<br>');
                }
                $command_1 = "ls ".$analysis_path."genomes.fa/ | wc -l";
                $gme_num = `$command_1`;
                if($gme_num < 10){
                    shell_exec("rm -rf ".$analysis_path);
                    die('Number of genomes for large-genome pathogen typing should at least 10!<br>');
                }
                shell_exec("mkdir ".$analysis_path."refgenome && chmod 777 -R ".$analysis_path);
                $command_2 = "for i in `ls ".$analysis_path."genomes.fa/`; "."do mv ".$analysis_path."genomes.fa/\${i} ".$analysis_path."genomes.fa/\${i%.*}.wlab; done";
                shell_exec($command_2);
                shell_exec("chmod 777 -R ".$analysis_path." && chmod 777 -R ".$analysis_path."genomes.fa/");
                $command_3 = "for i in `ls ".$analysis_path."genomes.fa/`; "."do mv ".$analysis_path."genomes.fa/\${i} ".$analysis_path."genomes.fa/\${i%.*}.fasta; done";
                shell_exec($command_3);
                shell_exec("chmod 777 -R ".$analysis_path." && chmod 777 -R ".$analysis_path."genomes.fa");
                $command_4 = $perl_path." ".$typing_pipeline." -WorkPath ".$analysis_path." -model AlleleCalling -len ".$allele_len." -percent ".$allele_p." -codon ".$codon_tbl." -StorePath ".$Result_store_path;
                shell_exec($command_4);
                shell_exec("chmod 777 ".$analysis_path.".../shell.sh");
                echo "Your Download code is: Strain-typing_".$dtime."<br>";
                echo "Please waiting for the analysis, remember the download code and you can close the window!<br>";
                pclose(popen('bash '.$analysis_path.'.../shell.sh &', 'r'));
            }else{
                die('All uploaded files should be .zip or .gz format!<br>');
            }
        }else{
            die('Reference or genomes files upload filed, Please re-upload them again!<br>');
        }
    }else if(!empty($_FILES['small-genomes'])){
        if($_FILES['small-genomes']['error'] == 0){
            if($_FILES['small-genomes']['size'] > 41943040){
                die('Due to the limitation of computing resource, the upload file size should be smaller than 40Mb!<br>');
            }
            $gme_end = end(explode('.', $_FILES['small-genomes']['name']));
            $gme_tmp = $_FILES['small-genomes']['tmp_name'];
            if(strcmp($gme_end, "gz")==0 || strcmp($gme_end, "zip")==0){
		shell_exec("mkdir ".$analysis_path);
                shell_exec("mkdir ".$analysis_path.".../");
                shell_exec("chmod 777 -R ".$analysis_path);
                $gme_filename = "genomes.fa".$gme_end;
                if(move_uploaded_file($gme_tmp, $analysis_path.$gme_filename)){
                    shell_exec("chmod 777 -R ".$analysis_path);
                    if($_FILES['small-genomes']['type'] == "application/zip" || $_FILES['small-genomes']['type'] == "application/x-zip-compressed"){
                        shell_exec("unzip -o ".$analysis_path.$gme_filename." -d ".$analysis_path."genomes.fa");
                        shell_exec("chmod 777 -R ".$analysis_path." && chmod 777 -R ".$analysis_path."genomes.fa");
                    }else if($_FILES['small-genomes']['type'] == "application/gzip"){
                        shell_exec("mkdir ".$analysis_path."genomes.fa");
                        shell_exec("chmod 777 -R ".$analysis_path." && chmod 777 -R ".$analysis_path."genomes.fa");
                        shell_exec("tar -zxvf ".$analysis_path.$gme_filename." --strip-components 1 -C ".$analysis_path."genomes.fa");
                        shell_exec("chmod 777 -R ".$analysis_path." && chmod 777 -R ".$analysis_path."genomes.fa");
                    }else{
                        shell_exec("rm -rf ".$analysis_path);
                        die('Please upload a correct .zip or .gz format for genomes file!<br>');
                    }
                    $command_1 = "for i in `ls ".$analysis_path."genomes.fa/`; "."do mv ".$analysis_path."genomes.fa/\${i} ".$analysis_path."genomes.fa/\${i%.*}.wlab; done";
                    shell_exec($command_1);
                    shell_exec("chmod 777 -R ".$analysis_path." && chmod 777 -R ".$analysis_path."genomes.fa/");
                    $command_2 = "for i in `ls ".$analysis_path."genomes.fa/`; "."do mv ".$analysis_path."genomes.fa/\${i} ".$analysis_path."genomes.fa/\${i%.*}.fasta; done";
                    shell_exec($command_2);
                    shell_exec("chmod 777 -R ".$analysis_path." && chmod 777 -R ".$analysis_path."genomes.fa");
                    $command_3 = "grep \> ".$analysis_path."genomes.fa/*fa* | wc -l";
                    $gme_num = `$command_3`;
                    if($gme_num < 5){
                        shell_exec("rm -rf ".$analysis_path);
                        die('Number of genomes for small-genomes pathogen typing should be at least 5!<br>');
                    }
                    $command_4 = $perl_path." ".$typing_pipeline." -WorkPath ".$analysis_path." -model WholeGenome -StorePath ".$Result_store_path;
                    shell_exec($command_4);
                    shell_exec("chmod 777 ".$analysis_path.".../shell.sh");
                    echo "Your Download code is: Strain-typing_".$dtime."<br>";
                    echo "Please waiting for the analysis, remember the download code and you can close the window!<br>";
                    pclose(popen('bash '.$analysis_path.'.../shell.sh &', 'r'));
                }else{
                    shell_exec("rm -rf ".$analysis_path);
                    die('Genomes file upload filed!<br>');
                }
            }else{
                die('All uploaded files should be .zip or .gz format!<br>');
            }
        }else{
            die('Genomes file upload filed, Please re-upload them again!<br>');
        }
    }
?>
