<?php
    $dtime = date('Y-m-d-h-i-s', time());
    $perl_path = "/home/dell/miniconda3/envs/Perl/bin/perl";
    $analysis_path = "/home/dell/data/Web-Server/Primer-design/".$dtime."/";
    $checkfa = "/var/www/other_scripts/bio-pipelines/Primer-design/Functions/checkfa.pl";
    $design_pipeline = "/var/www/other_scripts/bio-pipelines/Primer-design/primer-design.pl";
    $Result_store_path = "/home/dell/data/Web-Server/Result-store/Primer-design/".$dtime."/";
    if(empty($_POST['type'])){
        die('Please choose one of the PCR primer design type!<br>');
    }else{
	shell_exec("mkdir ".$analysis_path);
        shell_exec("mkdir ".$analysis_path.".../");
        shell_exec("chmod -R 777 ".$analysis_path);
        $save_file = fopen($analysis_path."target_seq.fasta", 'w+');
        fwrite($save_file, $_POST['target']);
        fclose($save_file);
        $command = $perl_path." ".$checkfa." -WorkPath ".$analysis_path." -target target_seq.fasta";
        $result = `$command`;
        $result = str_replace(array("\n", "\r"), "", $result);
        if(strcmp("TRUE", $result)==0){
            if($_FILES['conserve']['error'] == 0){
                $cons_tmp = $_FILES['conserve']['tmp_name'];
                $cons_rname = $_FILES['conserve']['name'];
                $cons_name = "conserve_info.tsv";
                if(move_uploaded_file($cons_tmp, $analysis_path.$cons_name)){
                    echo "$cons_rname file uploaded!<br>";
                    shell_exec("chmod -R 777 ".$analysis_path);
                    $command = $perl_path." ".$design_pipeline." -WorkPath ".$analysis_path." -target target_seq.fasta -type ".$_POST['type']." -consFILE ".$cons_name." -StorePath ".$Result_store_path;
                    shell_exec($command);
                    shell_exec("chmod 777 ".$analysis_path.".../shell.sh");
                    echo "Your Download code is: Primer-design_".$dtime."<br>";
                    echo "Please waiting for the analysis, remember the download code and you can close the window!<br>";
                    pclose(popen('bash '.$analysis_path.'.../shell.sh &', 'r'));
                }else{
                    shell_exec("rm -rf ".$analysis_path);
                    die('The conserve information file upload error!<br>');
                }
            }else{
                $command = $perl_path." ".$design_pipeline." -WorkPath ".$analysis_path." -target target_seq.fasta -type ".$_POST['type']." -StorePath ".$Result_store_path;
                shell_exec($command);
                shell_exec("chmod 777 ".$analysis_path.".../shell.sh");
                echo "Your Download code is: Primer-design_".$dtime."<br>";
                echo "Please waiting for the analysis, remember the download code and you can close the window!<br>";
                pclose(popen('bash '.$analysis_path.'.../shell.sh &', 'r'));
            }
        }else{
            shell_exec("rm -rf ".$analysis_path);
            die('Plsease input the FASTA format sequence!<br>');
        }
    }
?>
