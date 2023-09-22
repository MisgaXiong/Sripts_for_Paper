<?php
    $visualization = "visibility:hidden";
    $analysis_path = "/home/dell/data/Web-Server/Basic-operation/translate-seq/";
    $perl_path = "/home/dell/miniconda3/envs/Perl/bin/perl";
    $perl_script = "/var/www/other_scripts/bio-pipelines/Basic-opteration/basic-operation.pl";
    if(empty($_POST['codon'])){
        $codon = 1;
    }else{
        $codon = $_POST['codon'];
    }
    if($_FILES['input-file']['error'] > 0 or empty($_FILES['input-file'])){
        $save_file = fopen($analysis_path."input-seq.fasta", 'w+');
        fwrite($save_file, $_POST['inputseq']);
        fclose($save_file);
        shell_exec("chmod -R 777 ".$analysis_path);
        $command = $perl_path." ".$perl_script." ".$analysis_path." seqtranslate ".$analysis_path."input-seq.fasta ".$codon;
        $result=`$command`;
        shell_exec("rm -rf ".$analysis_path."input-seq.fasta");
        shell_exec("chmod -R 777 ".$analysis_path);
    }else{
        if(move_uploaded_file($_FILES['input-file']['tmp_name'], $analysis_path."input-seq.fasta")){
            $visualization = "visibility:visiable";
            shell_exec("chmod -R 777 ".$analysis_path);
            $command = $perl_path." ".$perl_script." ".$analysis_path." seqtranslate ".$analysis_path."input-seq.fasta ".$codon." > ".$analysis_path."translate-seq.fasta";
            shell_exec($command);
            shell_exec("chmod -R 777 ".$analysis_path);
            shell_exec("rm -rf ".$analysis_path."input-seq.fasta");
            shell_exec("chmod -R 777 ".$analysis_path);
            $download_path = "../Download-tools/doDownload.php?filename=".$analysis_path."translate-seq.fasta";
        }
    }
?>
<!DOCTYPE HTML>
<html lang="en">
	<head>
		<title>Sequence translation</title>
		<meta charset="UTF-8"/>
		<link rel="stylesheet" href="../css/sequence-translation-style.css">
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
					<td colspan="2" style="font-size:23px;">Translate nucleic acid into protein</td>
				</tr>
				<tr>
					<td colspan="2" style="font-size:15px; text-align:left; padding-left: 40px; padding-top:25px;">Translate functional gene sequences of pathogens into amino acid sequences.</td>
				</tr>
				<tr>
					<td colspan="2">
					<section id="model-select">
						<div class="model-header clearfix">
							<a href="javaScript:;" class="current" style="width:40%;">Sequence translation</a>
							<a href="javaScript:;" style="width:60%;">Sequence translation (Bulk sequences)</a>
						</div>
						<div class="model-body">
							<div class="model-menu" style="display:block;">
								<form  method="POST" enctype="multipart/form-data">
									<table style="width:100%;height:450px;">
										<tr>
											<td colspan="2"><div style="padding-left:40px;"><textarea required name="inputseq" placeholder="Input your sequence (fasta format) here... &#10Example:&#10>sequence1&#10ATGGCAATC...&#10>sequence2&#10ATGGCACAT..." style="width:500px;height:250px;resize:vertical;" maxlength="1950000" onchange="this.value=this.value.substring(0, 1950000)" onkeydown="this.value=this.value.substring(0, 1950000)" onkeyup="this.value=this.value.substring(0, 1950000)"></textarea></div></td>
										</tr>
										<tr>
											<td class="cell1">Translation code</td>
											<td>
												<select class="cell2" name="codon">
													<option value="1">default</option>
													<option value="1">most viruses</option>
													<option value="11">bacteria or plasmid</option>
												</select>
											</td>
										</tr>
										<tr>
											<td class="cell4" border="1" colspan="2">The current page could translate a small number of sequences. If you want translate a large number of sequence, please choose the Bulk sequences menu.</td>
										</tr>
										<tr>
											<td style="text-align:center;">
												<input class="cell3" type="submit" value="translate"/>
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
									<table style="width:100%;height:250px;">
										<tr>
											<td class="cell1">Bulk sequences file</td>
											<td>
												<input required class="cell2" type="file" name="input-file"/>
											</td>
										</tr>
										<tr>
											<td class="cell1">Translation code</td>
											<td>
												<select class="cell2" name="codon">
													<option value="1">default</option>
													<option value="1">most viruses</option>
													<option value="11">bacteria or plasmid</option>
												</select>
											</td>
										</tr>
										<tr>
											<td class="cell4" border="1" colspan="2">The large number of sequences should be in one FASTA format file to be uploaded.</td>
										</tr>
										<tr>
											<td style="text-align:center;">
												<input class="cell3" type="submit" value="translate"/>
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
			</table>
			</div>
			<div id="output">
				<table style="margin-left:50px;padding-top:6.5rem;width:600px;height:100%;float:left;">
					<tr>
						<td colspan="2">
							<textarea readonly placeholder="Output results..." style="width:550px;height:500px;resize:vertical;"><?php echo $result?></textarea>
						</td>
					</tr>
					<tr>
						<td>
							<a href="<?php echo $download_path?>" style="<?php echo $visualization?>">DownloadTranslatedFile</a>
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
