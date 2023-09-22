<?php 
    $analysis_path = "/home/dell/data/Web-Server/Basic-operation/align-seq/";
    $perl_path = "/home/dell/miniconda3/envs/Perl/bin/perl";
    $perl_script = "/var/www/other_scripts/bio-pipelines/Basic-opteration/basic-operation.pl";
    $save_file = fopen($analysis_path."input-seq.fasta", 'w+');
    fwrite($save_file, $_POST['inputseq']);
    fclose($save_file);
    shell_exec("chmod -R 777 ".$analysis_path);
    $command = $perl_path." ".$perl_script." ".$analysis_path." seqalignment ".$analysis_path."input-seq.fasta";
    $result=`$command`;
    shell_exec("rm -rf ".$analysis_path."*fasta");
    shell_exec("chmod -R 777 ".$analysis_path);
?>
<!DOCTYPE HTML>
<html lang="en">
	<head>
		<meta charset="UTF-8">
		<title>Sequence Alignment</title>
		<link rel="stylesheet" href="../css/sequence-alignment-style.css"/>
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
							<td colspan="2" style="font-size: 24px; padding-left: 20px;">Sequence alignment</td>
						</tr>
						<tr>
							<td colspan="2" style="font-size:15px; text-align:left; padding-left: 40px;">Perform pair-wise or multiple sequences alignment.</td>
						</tr>
						<tr>
							<td colspan="2">
								<div><textarea name="inputseq" placeholder="Input your sequences (fasta format) ... &#10Example:&#10>sequence1&#10ATGGCAATC&#10>sequence2&#10ATGGCACAT&#10>sequence3&#10ATGGCATGAA" style="width:550px;height:240px;resize:vertical;" maxlength="950000" onchange="this.value=this.value.substring(0,950000)" onkeydown="this.value=this.value.substring(0, 950000)" onkeyup="this.value=this.value.substring(0, 950000)"></textarea></div>
							</td>
						</tr>
						<tr>
							<td style="text-align:center;">
								<input class="cell3" type="submit" value="alignment"/>
							</td>
							<td style="text-align:center;">
								<input class="cell3" type="reset" value="cancel"/>
							</td>
						</tr>
						<tr>
							<td colspan="2" style="font-size:10px; text-align:justify;border:2px solid #eeeeee;">Since most of the computing resources are allocated to other functions of PHDtools, the file size of this function for multiple sequence alignment cannot exceed 950kb.</td>
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
