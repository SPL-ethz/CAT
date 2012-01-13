<?php
// GUI for the administration of solubility data
// (c) Stefan Schorsch 2012 - SPL - IPE - MAVT - ETHZ
error_reporting(E_ALL);
// DB Connection
mysql_connect('129.132.152.27','PBEToolbox','toolbox2000');
mysql_select_db('PBEToolbox');
$action=(isset($_REQUEST['action'])?$_REQUEST['action']:'home');
if($action=='add_new') {
	mysql_query('INSERT INTO solubilityFunctions (function,solvent,solute,tRange,citation,comments) VALUES (
	"'.$_POST['myFunction'].'",
	"'.$_POST['solvents'].'",
	"'.$_POST['solute'].'",
	"'.$_POST['tRange'].'",
	"'.$_POST['citation'].'",
	"'.$_POST['comments'].'"
	)');
	$result = mysql_query("SELECT id FROM solubilityFunctions ORDER BY id desc LIMIT 0,1");
	$array = mysql_fetch_array($result);
	$status = 'New set added to database, the id is: '.$array['id'];
}
if($action=='do_update') {
	mysql_query('UPDATE solubilityFunctions SET 
	function = "'.$_POST['myFunction'].'",
	solvent = "'.$_POST['solvents'].'",
	solute = "'.$_POST['solute'].'",
	tRange = "'.$_POST['tRange'].'",
	citation = "'.$_POST['citation'].'",
	comments = "'.$_POST['comments'].'"
	WHERE id = '.$_POST['myID']);
	$status = 'Set updated in database, the ID is: '.$_POST['myID'];
}
if($action=='do_delte') {
	mysql_query('DELETE FROM solubilityFunctions WHERE id = '.$_POST['myID']);
	$status = 'Set has been deleted';
}

function myFORM($id,$myFunction,$solvent,$solute,$tRange,$citation,$comments,$whatToDo) {
?>
<form action="index.php" method="post">
<input type="hidden" name="action" value="<?php echo ($whatToDo=='new'?'add_new':'do_update'); ?>" />
<input type="hidden" name="myID" value="<?php echo $id; ?>" />
<table border="0">
<tr style="background-color: #fafafa; color: #000000;"><td colspan="2" style="font-size: 14px; font-weight: bold;"><?php echo $whatToDo=='new'?'Add a new set':'Edit set '.$id; ?></td></tr>
<tr style="background-color: #f0f0f0; color: #000000;"><td>Solvent(s)</td><td><input type="text" name="solvents" value="<?php echo $solvent; ?>" /></td></tr>
<tr style="background-color: #fafafa; color: #000000;"><td>Solute</td><td><input type="text" name="solute" value="<?php echo $solute; ?>" /></td></tr>
<tr style="background-color: #f0f0f0; color: #000000;"><td>Temperature Range</td><td><input type="text" name="tRange" value="<?php echo $tRange; ?>" /></td></tr>
<tr style="background-color: #fafafa; color: #000000;"><td colspan="2"><b>Function of solubility [g/kg]</b></td></tr>
<tr style="background-color: #f0f0f0; color: #000000;"><td valign="top">@(T,t,x,L)<br />T is temperature [K]<br />
t is time [s]<br />
x is composition vector [-]<br />
L is length [m]<br /></td><td valign="top"><textarea name="myFunction" rows="5" cols="70"><?php echo $myFunction; ?></textarea></td></tr>
<tr style="background-color: #fafafa; color: #000000;"><td valign="top">Citation</td><td><textarea name="citation" rows="3" cols="70"><?php echo $citation; ?></textarea></td></tr>
<tr style="background-color: #f0f0f0; color: #000000;"><td valign="top">Comments</td><td><textarea name="comments" rows="3" cols="70"><?php echo $comments; ?></textarea></td></tr>
<tr style="background-color: #fafafa; color: #000000;"><td colspan="2"><input type="submit" value="<?php echo $whatToDo=='new'?'add new set':'edit set'; ?>" /></td></tr>
</table>
</form>
<?php } ?>
<html>
<head>
<title>
ETHZ PBEToolbox Crystal Solubility Database
</title>
<style type="text/css" media="all">
  body {
	background-color: #fdfdfd;
	color: #000000;
	font-family: Arial;
  }
  a {
	color: #000000;	
  }
  a:visited {
	color: #000000;	
  }
  a:hover {
	color: #0d0d0d;	
	text-decoration: none;
  }
 </style>
</head>
<body>
<h1>ETHZ PBEToolbox Crystal Solubility Database</h1>
<?php echo isset($status)?'Status: <span style="background-color: #cfcfcf; color: #000000;">'.$status.'<span> | ':''; ?><a href="index.php">HOME & ADD DATA</a> | <a href="index.php?action=browse">BROWSE DATA</a>
<?php
switch ($action) {
	case 'home':
	case 'add_new':
	myFORM(0,'','','','','','','new');
	break;	
	case 'browse':
	case 'do_delte':
	?>
	<table border="0">
	<tr style="background-color: #a0a0a0; color: #ffffff;"><td width="50">ID</td><td width="150">Function</td><td width="150">Solvent</td><td width="150">Solute</td><td width="150">Temperature Range</td><td width="150">Citation</td><td width="150">Comments</td><td width="110">Options</td></tr>
	<?php
	$i = 0;
	$result = mysql_query('SELECT id,function AS myFunction, solvent,solute,tRange,citation,comments FROM solubilityFunctions');
	while($a = mysql_fetch_array($result)) {
		if($i == 0) {
			$i++;
			echo '<tr style="background-color: #fafafa; color: #000000;">';
		} else {
			$i = 0;
			echo '<tr style="background-color: #f0f0f0; color: #000000;">';
		}
		echo '<td>'.$a['id'].'</td>';
		echo '<td>'.((strlen($a['myFunction'])>20)?substr($a['myFunction'],0,17).'...':$a['myFunction']).'</td>';
		echo '<td>'.$a['solvent'].'</td>';
		echo '<td>'.$a['solute'].'</td>';
		echo '<td>'.$a['tRange'].'</td>';
		echo '<td>'.((strlen($a['citation'])>20)?substr($a['citation'],0,17).'...':$a['citation']).'</td>';
		echo '<td>'.((strlen($a['comments'])>20)?substr($a['comments'],0,17).'...':$a['comments']).'</td>';
		echo '<td><a href="index.php?action=update&myID='.$a['id'].'">Edit</a> | <a href="index.php?action=do_delte&myID='.$a['id'].'">Delete</a>';
		echo '</td></tr>';
	}
	?>
	</table>
	<?php
	break;
	case 'update':
	case 'do_update':
		$result = mysql_query('SELECT id,function AS myFunction, solvent,solute,tRange,citation,comments FROM solubilityFunctions WHERE id = '.$_REQUEST['myID']);
		$a = mysql_fetch_array($result);
		myFORM($_REQUEST['myID'],$a['myFunction'],$a['solvent'],$a['solute'],$a['tRange'],$a['citation'],$a['comments'],'update');
	break;	
}
?>
</body>
</html>
<?php
mysql_close();
?>