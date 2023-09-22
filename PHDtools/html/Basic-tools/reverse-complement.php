<?php
    $analysis_path = "/home/dell/data/Web-Server/Basic-operation/reverse-complen/";
    $perl_path = "/home/dell/miniconda3/envs/Perl/bin/perl";
    $perl_script = "/var/www/other_scripts/bio-pipelines/Basic-opteration/basic-operation.pl";
    $save_file = fopen($analysis_path."input-seq.fasta", 'w+');
    fwrite($save_file, $_POST['inputseq']);
    fclose($save_file);
    shell_exec("chmod -R 777 ".$analysis_path);
    $command = $perl_path." ".$perl_script." ".$analysis_path." reversecomplent ".$analysis_path."input-seq.fasta ".$_POST['chose'];
    $result=`$command`;
    shell_exec("rm -rf ".$analysis_path."input-seq.fasta");
    shell_exec("chmod -R 777 ".$analysis_path);
?>
<!DOCTYPE HTML>
<html lang="en">
	<head>
		<meta charset="UTF-8">
		<title>Sequence Reverse Complementation</title>
		<link rel="stylesheet" href="../css/reverse-complement-style.css"/>
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
			<div id="center">
				<form method="POST">
					<table>
						<tr>
							<td colspan="2" style="font-size: 24px; padding-left: 20px;">Reverse or complementary processing</td>
						</tr>
						<tr>
							<td colspan="2" style="font-size:15px; text-align:left; padding-left: 40px;">Perform basic operations on nucleotide sequences.</td>
						</tr>
						<tr>
							<td colspan="2">
								<div style="padding-left:40px;"><textarea required type="text" name="inputseq" placeholder="Input your sequence (fasta format) here... &#10Example: &#10>sequence1&#10ATCGGTAA&#10>sequence2&#10CACCTA&#10>sequence3&#10TTCGGACACGTCGA" style="width:550px;height:250px;resize:vertical;"></textarea></div>
							</td>
						</tr>
						<tr>
							<td class="cell1">Operation model</td>
							<td>
								<select class="cell2" name="chose">
									<option value="reverse">reverse</option>
									<option value="complementation">complementation</option>
									<option value="reverse-complement">reverse and complementation</option>
								</select>
							</td>
						</tr>
						<tr>
							<td style="text-align:center;">
								<input class="cell3" type="submit" value="process"/>
							</td>
							<td style="text-align:center;">
								<input class="cell3" type="reset" value="cancel"/>
							</td>
						</tr>
					</table>
				</form>
			</div>
			<div id="output">
				<table style="margin-left:50px;padding-top:6.5rem;width:500px;height:100%;float:left;">
					<tr>
						<td colspan="2">
							<textarea readonly placeholder="Output results..." style="width:500px;height:370px;resize:vertical;"><?php echo $result?></textarea>
						</td>
					</tr>
				</table>
			</div>
		</section>
	</body>
</html>
