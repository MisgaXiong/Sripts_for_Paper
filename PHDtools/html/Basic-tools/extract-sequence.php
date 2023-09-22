<?php
    $visualization = "visibility:hidden";
    $analysis_path = "/home/dell/data/Web-Server/Basic-operation/extract-seq/";
    $perl_path = "/home/dell/miniconda3/envs/Perl/bin/perl";
    $perl_script = "/var/www/other_scripts/bio-pipelines/Basic-opteration/basic-operation.pl";
    if(strcmp($_POST['model'], "motif") == 0){
        $subextract = $_POST['subextract'];
        $motiftype = $_POST['motif-type'];
        $inputtype = $_POST['input-type'];
        if($_FILES['input-file']['error'] > 0 or empty($_FILES['input-file'])){
            if($_FILES['motif-seq']['error'] > 0){
                die('The motif file must be in fasta format and please upload again!<br>');
            }else{
                $save_file = fopen($analysis_path."input-seq.fasta", 'w+');
                fwrite($save_file, $_POST['inputseq']);
                fclose($save_file);
                if(move_uploaded_file($_FILES['motif-seq']['tmp_name'], $analysis_path."motif-seq.fasta")){
                    shell_exec("chmod -R 777 ".$analysis_path);
                    $command = $perl_path." ".$perl_script." ".$analysis_path." seqextract input-seq.fasta motif-extract ".$subextract." motif-seq.fasta ".$motiftype." ".$inputtype;
                    $result=`$command`;
                    shell_exec("rm -rf ".$analysis_path."input*");
                    shell_exec("rm -rf ".$analysis_path."motif*");
                    shell_exec("rm -rf ".$analysis_path."makedb*");
                    shell_exec("rm -rf ".$analysis_path."blast*");
                    shell_exec("chmod -R 777 ".$analysis_path);
                }else{
                    shell_exec("rm -rf ".$analysis_path."*");
                    die('The motif file upload failed!<br>');
                }
            }
        }else{
            if(move_uploaded_file($_FILES['input-file']['tmp_name'], $analysis_path."input-seq.fasta")){
                $visualization = "visibility:visiable";
                shell_exec("chmod -R 777 ".$analysis_path);
                if(move_uploaded_file($_FILES['motif-seq']['tmp_name'], $analysis_path."motif-seq.fasta")){
                    $command = $perl_path." ".$perl_script." ".$analysis_path." seqextract input-seq.fasta motif-extract ".$subextract." motif-seq.fasta ".$motiftype." ".$inputtype." > ".$analysis_path."extracted-seqs.fasta";
                    shell_exec($command);
                    shell_exec("chmod -R 777 ".$analysis_path);
                    shell_exec("rm -rf ".$analysis_path."input*");
                    shell_exec("rm -rf ".$analysis_path."motif*");
                    shell_exec("rm -rf ".$analysis_path."makedb*");
                    shell_exec("rm -rf ".$analysis_path."blast*");
                    shell_exec("chmod -R 777 ".$analysis_path);
                    $download_path = "../Download-tools/doDownload.php?filename=".$analysis_path."extracted-seqs.fasta";
                }else{
                    shell_exec("rm -rf ".$analysis_path."*");
                    die('The motif file upload failed!<br>');
                }
            }else{
                shell_exec("rm -rf ".$analysis_path."*");
                die('The input file upload failed!<br>');
            }
        }
    }else if(strcmp($_POST['model'], "sites") == 0){
        $start = $_POST['start'];
        $end = $_POST['end'];
        if(is_numeric($start) && is_numeric($end)){
            if($start > $end){
                shell_exec("rm -rf ".$analysis_path."*");
                die('The start position must be smaller than the end position!<br>');
            }
            $save_file = fopen($analysis_path."input-seq.fasta", 'w+');
            fwrite($save_file, $_POST['inputseq']);
            fclose($save_file);
            shell_exec("chmod -R 777 ".$analysis_path);
            $ids = explode("\n", $_POST['inputseq']);
            $id_arr = explode(" ", $ids[0]);
            $id = str_replace(">", "", $id_arr[0]);
            $id = str_replace("\r", "", $id);
            $site_info = join("\t", array($id, $start, $end))."\n";
            $info_file = fopen($analysis_path."sitesinfo.tsv", 'w+');
            fwrite($info_file, $site_info);
            fclose($info_file);
            shell_exec("chmod -R 777 ".$analysis_path);
            $command = $perl_path." ".$perl_script." ".$analysis_path." seqextract input-seq.fasta site-extract TRUE sitesinfo.tsv";
            $result=`$command`;
            shell_exec("rm -rf ".$analysis_path."input-seq.fasta");
            shell_exec("rm -rf ".$analysis_path."sitesinfo*");
            shell_exec("chmod -R 777 ".$analysis_path);
        }else{
            die('The extraction sites both start and end positions should be numeric!<br>');
        }
    }
?>
<!DOCTYPE HTML>
<html lang="en">
	<head>
		<title>Extract sequence</title>
		<meta charset="UTF-8"/>
		<link rel="stylesheet" href="../css/extract-sequence-style.css">
	</head>
	<body>
		<div class="navbar">
			<a href="#" class="logo">PHDtools: Basic Sequence Operation</a>
			<ul class="nav">
				<li><a href="../PHDtools-home.html">Back to home</a></li>
				<li><a href="../PHDtools-team.html">Team</a></li>
				<li><a href="#publication">Publications</a></li>
			</ul>
		</div>
		<section id="all">
			<div id="center" style="margin-left: 50px;width:700px;height:100%;float:left;">
			<table style="width:100%;">
				<tr>
					<td colspan="2" style="font-size:23px;">Target sequence extraction</td>
				</tr>
				<tr>
					<td colspan="2" style="font-size:15px; text-align:left; padding-left: 40px; padding-top:25px;">Extract the target sequence based on the similarity or sites.</td>
				</tr>
				<tr>
					<td colspan="2">
					<section id="model-select">
						<div class="model-header clearfix">
							<a href="javaScript:;" class="current" style="width:25%;">Motif extraction</a>
							<a href="javaScript:;" style="width:50%;">Motif extraction (Bulk sequences)</a>
							<a href="javaScript:;" style="width:25%;">Site extraction</a>
						</div>
						<div class="model-body">
							<div class="model-menu" style="display:block;">
								<form  method="POST" enctype="multipart/form-data">
									<table style="width:100%;height:650px;">
										<tr>
											<td colspan="2"><div style="padding-left:40px;"><textarea required name="inputseq" placeholder="Input your sequence (fasta format) here... &#10Example:&#10>sequence1&#10ATGGCAATC...&#10>sequence2&#10ATGGCACAT..." style="width:500px;height:250px;resize:vertical;" maxlength="1950000" onchange="this.value=this.value.substring(0, 1950000)" onkeydown="this.value=this.value.substring(0, 1950000)" onkeyup="this.value=this.value.substring(0, 1950000)"></textarea></div></td>
										</tr>
										<tr>
											<td class="cell1">Extract model</td>
											<td>
												<select class="cell2" name="model">
													<option value="motif">motif</option>
												</select>
											</td>
										</tr>
										<tr>
											<td class="cell1">Reference sequence (motif)</td>
											<td>
												<input class="cell2" required type="file" name="motif-seq"/>
											</td>
										</tr>
										<tr>
											<td class="cell1">Input sequence type</td>
											<td>
												<select class="cell2" name="input-type"/>
													<option value="nucleicacid">nucleicacid</option>
													<option value="protein">protein</option>
												</select>
											</td>
										</tr>
										<tr>
											<td class="cell1">Reference sequence type</td>
											<td>
												<select class="cell2" name="motif-type">
													<option value="nucleicacid">nucleicacid</option>
													<option value="protein">protein</option>
												</select>
											</td>
										</tr>
										<tr>
											<td class="cell1">Extract region</td>
											<td>
												<select class="cell2" name="subextract">
													<option value="TRUE">extract aligned region</option>
													<option value="FALSE">extract whole fragment</option>
												</select>
											</td>
										</tr>
										<tr>
											<td class="cell4" border="1" colspan="2">The motif extraction can extract the fragment which contain the target reference (motif) sequence. You can choose the menu to decide extract the whole fragment or a short region that only aligned to the reference sequence. If your only have small number of sequences to extract, you can use current menu, else you can choose the Bulk sequences menu.</td>
										</tr>
										<tr>
											<td style="text-align:center;">
												<input class="cell3" type="submit" value="extract"/>
											</td>
											<td style="text-align:center;">
												<input class="cell3" type="reset" value="cancel"/>
											</td>
										</tr>
									</table>
								</form>
							</div>
							<div class="model-menu">
								<form  method="POST" enctype="multipart/form-data">
									<table style="width:100%;height:450px;">
										<tr>
											<td class="cell1">Bulk sequences file</td>
											<td>
												<input required class="cell2" type="file" name="input-file"/>
											</td>
										</tr>
										<tr>
											<td class="cell1">Extract model</td>
											<td>
												<select class="cell2" name="model">
													<option value="motif">motif</option>
												</select>
											</td>
										</tr>
										<tr>
											<td class="cell1">Reference sequence (motif)</td>
											<td>
												<input class="cell2" required type="file" name="motif-seq"/>
											</td>
										</tr>
										<tr>
											<td class="cell1">Input sequence type</td>
											<td>
												<select class="cell2" name="input-type"/>
													<option value="nucleicacid">nucleicacid</option>
													<option value="protein">protein</option>
												</select>
											</td>
										</tr>
										<tr>
											<td class="cell1">Reference sequence type</td>
											<td>
												<select class="cell2" name="motif-type">
													<option value="nucleicacid">nucleicacid</option>
													<option value="protein">protein</option>
												</select>
											</td>
										</tr>
										<tr>
											<td class="cell1">Extract region</td>
											<td>
												<select class="cell2" name="subextract">
													<option value="TRUE">extract aligned region</option>
													<option value="FALSE">extract whole fragment</option>
												</select>
											</td>
										</tr>
										<tr>
											<td class="cell4" border="1" colspan="2">The motif extraction can extract the fragment which contain the target reference (motif) sequence. You can choose the menu to decide extract the whole fragment or a short region that only aligned to the reference sequence. The current menu is used to extract the large number of sequences.</td>
										</tr>
										<tr>
											<td style="text-align:center;">
												<input class="cell3" type="submit" value="extract"/>
											</td>
											<td style="text-align:center;">
												<input class="cell3" type="reset" value="cancel"/>
											</td>
										</tr>
									</table>
								</form>
							</div>
							<div class="model-menu">
								<form  method="POST" enctype="multipart/form-data">
									<table style="width:100%;height:500px;">
										<tr>
											<td colspan="2"><div style="padding-left:40px;"><textarea required name="inputseq" placeholder="Input one sequence (fasta format) here... &#10Example:&#10>sequence1&#10ATGGCAATCATGGCACATATTCTAGCCATG..." style="width:500px;height:250px;resize:vertical;"></textarea></div></td>
										</tr>
										<tr>
											<td class="cell1">Extract model</td>
											<td>
												<select class="cell2" name="model">
													<option value="sites">site</option>
												</select>
											</td>
										</tr>
										<tr>
											<td class="cell1">Start site</td>
											<td>
												<input required class="cell2" type="text" name="start" placeholder="eg: 1 (smaller than end site)"/>
											</td>
										</tr>
										<tr>
											<td class="cell1">End site</td>
											<td>
												<input required class="cell2" type="text" name="end" placeholder="eg: 100 (greater than start site)"/>
											</td>
										</tr>
										<tr>
											<td class="cell4" colspan="2">The current menu was used to extract the sub-sequence of a long sequence based on the sites.</td>
										</tr>
										<tr>
											<td style="text-align:center;">
												<input class="cell3" type="submit" value="extract"/>
											</td>
											<td style="text-align:center;">
												<input class="cell3" type="reset" value="cancel"/>
											</td>
										</tr>
									</table>
								</form>
							</div>
						</div>
					</section>
					</td>
				</tr>
				<tr>
					<td colspan="2" style="font-size:10px; text-align:justify;border:2px solid #eeeeee;">The sequence extraction function provides two different modes to extract sequences. This function can be used in a variety of scenarios, such as extracting target pathogen contigs from metagenomic data, or extracting specific gene sequences from a draft genome of a pathogen.</td>
				</tr>
			</table>
			</div>
			<div id="output">
				<table style="margin-left:50px;padding-top:6.5rem;width:500px;height:100%;float:left;">
					<tr>
						<td colspan="2">
							<textarea readonly placeholder="Output results..." style="width:550px;height:500px;resize:vertical;"><?php echo $result?></textarea>
						</td>
					</tr>
					<tr>
						<td colspan="2">
							<a href="<?php echo $download_path?>" style="<?php echo $visualization?>">DownloadExtractedFile</a>
						</td>
					</tr>
				</table>
			</div>
		</section>
		<script type="application/javascript" >
			window.onload=function () {
				let as=document.getElementsByClassName('model-header')[0].getElementsByTagName('a');
				let contents=document.getElementsByClassName("model-menu");
				for (let i=0;i<as.length;i++){
					let a=as[i];
					a.id=i;
					a.onclick=function () {
						for(let j=0;j<as.length;j++){
							as[j].className="";
							contents[j].style.display="none";
						}
						this.className='current';
						contents[this.id].style.display='block';
					}
				}
			}
		</script>
	</body>
</html>
