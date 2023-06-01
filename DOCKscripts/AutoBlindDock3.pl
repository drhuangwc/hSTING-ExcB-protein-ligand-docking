#!/usr/bin/perl   
use strict;

if($#ARGV!=6)
{
      print "**************************************************************************\n";
      print "* Error in input data!                                                   *\n";
      print "* ./AutoDock.pl [Receptor] [Ligand] [PocketNum] [result_folder] [CB-DOCKpath] [numbermodes] [exhaustiveness]               *\n";
      print "* Example: ./AutoBlindDock.pl receptor.pdb ligand.sdf 5 test CB-DOCKpath vine_path *\n";
      print "**************************************************************************\n";
      exit 0;
}

my $pwd_home = `pwd`; chomp $pwd_home; 
print "\nCurrent working DIR = $pwd_home "."\n";

my $protein=$ARGV[0];
my $ligand=$ARGV[1];
my $PocketNum=$ARGV[2];
my $userDirName=$ARGV[3];
my $progPath=$ARGV[4]; 

my $prepare_receptor=$progPath.'/prepare_receptor4.py'; 
my $prepare_ligand='mk_prepare_ligand.py' ;
my $autodockvina=$progPath.'/vina';

my $sec; my $min; my $hour;
my $day; my $mon; my $year;
my $wday; my $yday; my $isdst;
   ($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst)=localtime(time());
   $year+=1900;
   $mon+=1;

my $var = 5;
my $nm=$ARGV[5];
my $ex=$ARGV[6];

my @lig_name=split /[\.\/]/, $ligand;
my @pro_name=split /[\.\/]/, $protein;
my $mol_name=@lig_name[-2]."_".@pro_name[-2];
my $dock_file_home=".\/dock_file";
mkdir $dock_file_home;
my $userDirPath=".\/dock_file\/$userDirName";
mkdir $userDirPath;

my $dock_time=$year.($mon<10?"0$mon":$mon).($day<10?"0$day":$day).($hour<10?"0$hour":$hour).($min<10?"0$min":$min).($sec<10?"0$sec":$sec);

my $dock_file=$mol_name;
my $outf="$userDirPath\/$mol_name"; 

my $logfile=$dock_file."_log.txt";
my $errfile=$dock_file."_err.txt";
my $config="config.txt";
my %conf;
    mkdir "$outf";
    system ("touch $outf\/$errfile");
    system ("touch $userDirPath\/status.txt");
    open (f, ">>$userDirPath\/status.txt");
    print f "$year-$mon-$day $hour:$min:$sec\n";
    close f;
    
    system ("touch $userDirPath\/$dock_time\_run.txt");
    open (r, ">>$userDirPath\/$dock_time\_run.txt");
    print r "$dock_file  $PocketNum  ";
    close r;
=pod    
    #Statistics "run.txt" file number, up to keep the latest three
	my @runs=glob("$userDirPath\/*_run.txt");
	if(@runs > 3)
	{
		my @times;
		foreach(@runs)
		{
			my @tmp=split /\/|\:|\_/,$_;
			push(@times, $tmp[-2]);
			
		}
		sort @times;
		for(my $i=0; $i<=($#times-3); $i++)
		{
			my $del=$times[$i]."_run.txt";
			system ("rm $userDirPath\/$del");
		}
		
	}
=cut
    
    my $pro_ori = $protein;
    $pro_ori =~ s/\.pdb//g;
    my $pro_format=$pro_ori."_format.pdb";
    my $cmd="$progPath/FormatPDB_Simple $protein $pro_format";
    print "\n".$cmd."\n";
	system ($cmd);
	
	if(-e $pro_format)
	{
	system "echo protein format success! >> $outf\/$logfile ";
        print "\n protein format success! \n";
	}
	else
	{
	system "echo protein format error! >> $outf\/$errfile ";
        print "\n protein format error! \n";

	system("perl -pi -e 's/$year-$mon-$day $hour:$min:$sec\n//' $userDirPath\/status.txt");
	exit;
	}
	
	if(@lig_name[-1] ne "sdf")
	{
        my $cmd="babel -i@lig_name[-1] $ligand -osdf $ligand.sdf -p 7";
        print $cmd."\n\n";
	system "$cmd";#ligand format transfer --error info
	$ligand=$ligand.".sdf";
		
	if(-e $ligand)
	{
		system "echo ligand transfer success! >> $outf\/$logfile ";
	}
	else
	{
		system "echo ligand transfer error! >> $outf\/$errfile ";
		system("perl -pi -e 's/$year-$mon-$day $hour:$min:$sec\n//' $userDirPath\/status.txt");
		exit;
	}
	}
	system ("mv $pro_format $outf\/receptor.pdb");
	system ("cp $ligand $outf\/ligand.sdf");	
	
	####################################################################
	#curvatureSurface		
	###################################################################
    	system("$progPath/curvatureSurface/bin/curvatureSurface $outf\/receptor.pdb $outf\/grid.pdb");
	
	my $grid_filepath=$outf."/grid.pdb";
	if(-e $grid_filepath)
	{
		system "echo curvature calculate success! >> $outf\/$logfile ";
	}
	else
	{
		system "echo curvature calculate error! >> $outf\/$errfile ";
		system("perl -pi -e 's/$year-$mon-$day $hour:$min:$sec\n//' $userDirPath\/status.txt");
		exit;
	}

	####################################################################
	#clusters		
	###################################################################
    	my $cmd="$progPath/clusters $outf\/grid.pdb $PocketNum > $outf\/conf.txt";
    	print "\n".$cmd."\n\n";
	system "$cmd"; 

	####################################################################
	#ADT_scripts: prepare_ligand4.py，prepare_receptor4.py, prepare_dpf4.py, eBoxSize.pl;                       
	#prepare_ligand4.py——The imported ligands are converted from sdf format to pdbqt format；          
	#prepare_receptor4.py——Converts the inputted receptor from pdb format to pdbqt format；         
	#prepare_dpf4.py——Find the center of the docking pocket, used in Re-Docking；                                 
	#eBoxSize.pl——Calculate the size of the ligand;		
	####################################################################
	my $cmd="$prepare_ligand -i $outf\/ligand.sdf -o $outf\/ligand.pdbqt --keep_nonpolar_hydrogens";    

	system "$cmd";
    
	my $lig_filepath=$outf."/ligand.pdbqt";
    
	if(-e $lig_filepath)
	{
		system "echo ligand pdbqt transfer success! >> $outf\/$logfile ";
        print "\n ligand pdbqt transfer success! \n";
	}
	else
	{
		system "echo ligand pdbqt transfer error! >> $outf\/$errfile ";
        print "\n ligand pdbqt transfer error! \n";
		system("perl -pi -e 's/$year-$mon-$day $hour:$min:$sec\n//' $userDirPath\/status.txt");
		exit;
	}

        my $cmd="$prepare_receptor -r $outf\/receptor.pdb -o $outf\/receptor.pdbqt";
        print "\n".$cmd."\n\n";
        system "$cmd";

	
	my $pro_filepath=$outf."/receptor.pdbqt";
	if(-e $pro_filepath)
	{
		system "echo receptor pdbqt transfer success! >> $outf\/$logfile ";
        print "\n receptor pdbqt transfer success! \n";
	}
	else
	{
		system "echo receptor pdbqt transfer error! >> $outf\/$errfile ";
        print "\n receptor pdbqt transfer error! \n";
		system("perl -pi -e 's/$year-$mon-$day $hour:$min:$sec\n//' $userDirPath\/status.txt");
		exit;
	}
	
	system "perl $progPath/eBoxSize.pl $outf\/ligand.pdbqt >$outf\/tem_1.txt";
		
#####################ReDocking and GlobalDocking########################
#（2）get the ligand size;
	my $sx;
	open FILE1, "<$outf\/tem_1.txt" or die "Error in opening tem_1.txt";
	while(<FILE1>)
	{
	  chmod;
	  $sx=$_;
	}
	close FILE1;
	system "rm $outf\/tem_1.txt";
	$sx=int($sx)+1;
####################################################################	
#LocalDocking
	open FILE, "<$outf\/conf.txt" or die "opening error!:$!";
	system "echo Cavities  volume  center_x  center_y  center_z  size_x  size_y  size_z  score>> $outf\/$config ";
	
	my @data=<FILE>;
	foreach(@data)
	{
		if($_ =~ /Cavities/)
		{
			next;
		}
		my @array=split /\s+/, $_;
		my $num=$array[0];
		my $_cx=$array[1];
		my $_cy=$array[2];
		my $_cz=$array[3];

		my $_sx=$array[4];	
		if($_sx >= $sx)
		{
			if($_sx - $sx <= 2*$var )
			{
				$_sx = $sx + 2*$var;
			}
			else
			{
				$_sx = $_sx + 5;
			}
		}
		else
		{
			$_sx = $sx + 2*$var;
		}
		
		my $_sy=$array[5];		
		if($_sy >= $sx)
		{
			if($_sy - $sx <= 2*$var )
			{
				$_sy = $sx + 2*$var;
			}
			else
			{
				$_sy = $_sy + 5;
			}	
		}
		else
		{
			$_sy = $sx + 2*$var;
		}	
		
		my $_sz=$array[6];		
		if($_sz >= $sx)
		{
			if($_sz - $sx <= 2*$var )
			{
				$_sz = $sx + 2*$var;
			}
			else
			{
				$_sz = $_sz + 5;
			}
			
		}
		else
		{
			$_sz = $sx + 2*$var;
		}
		

        my $out_ligand=$mol_name."_out_".$num.".pdbqt";


		$conf{$num}="$num  $array[7]  $_cx  $_cy  $_cz  $_sx  $_sy  $_sz  ";
		system "echo Calculate $num --LocalDocking >> $outf\/$logfile ";
        print "\n Calculate $num --LocalDocking >> $outf\/$logfile \n";
        print "\n\n"."$autodockvina --receptor $outf\/receptor.pdbqt --ligand $outf\/ligand.pdbqt --center_x $_cx --center_y $_cy --center_z $_cz --size_x $_sx --size_y $_sy --size_z $_sz --num_modes $nm --exhaustiveness $ex --out $outf\/$out_ligand >> $outf\/$logfile"."\n\n";
		system "$autodockvina --receptor $outf\/receptor.pdbqt --ligand $outf\/ligand.pdbqt --center_x $_cx --center_y $_cy --center_z $_cz --size_x $_sx --size_y $_sy --size_z $_sz --num_modes $nm --exhaustiveness $ex --out $outf\/$out_ligand >> $outf\/$logfile";
		
		my $lig_outfilepath=$outf."/".$out_ligand;
		if(-e $lig_outfilepath)
		{
			system "echo docking in the $num th cavity success! >> $outf\/$logfile ";
            print "docking in the $num th cavity success!  $outf\/$logfile \n\n";
		}
		else
		{
			system "echo docking error! >> $outf\/$errfile ";
            print "docking error \n";
			system("perl -pi -e 's/$year-$mon-$day $hour:$min:$sec\n//' $userDirPath\/status.txt");
			exit;
		}
	}
	
	
#######################################################################
#process the output files
	system ("mv $outf\/receptor.pdbqt $outf\/receptor.pqbk");
	system ("mv $outf\/ligand.pdbqt $outf\/ligand.pqbk");
	system ("mv $outf\/conf.txt $outf\/cavities.txt ");
	
	my @pdbqt = glob("$outf\/*.pdbqt");
	foreach(@pdbqt)
	{
		my $out_ligand = $_;
		$out_ligand =~ s/\.pdbqt/\.sdf/g;
		system "babel -ipdbqt $_ -osdf $out_ligand ---errorlevel 1";

	}
	
###################################################################
#get the vina score of the first pose
my @data;
my @buff;
my $num=0;
open FILE, "<$outf\/$logfile" or die "Error in opening logfile\n";
while(<FILE>)
{	  
	if($_=~/^Calculate/)
	{
		my @temp=split(/\s+/, $_);
		@buff=();
		push @buff, $temp[1];
	}
	if($_=~/^   1 /)
	{     		
		my @temp=split(/\s+/, $_);
		push @buff, $temp[2]; 
		push @data, [@buff];  
		$num++;             
	}
		  
}

for(my $i=0; $i<$num; $i++)
{
	$conf{$i+1}.=$data[$i][1];
	system "echo $conf{$i+1}>> $outf\/$config";
}

my @sdf =  @pdbqt;

for(my $i=0; $i<=$#sdf; $i++)
{	
	$sdf[$i] =~ s/\.pdbqt/\.sdf/g;
	my $out_ligand_score = $sdf[$i];
	my @tmp = split /\.|\_/,$sdf[$i];
	for(my $m=0; $m<$num; $m++)
	{
		if($tmp[-2] == $data[$m][0])
		{
			$out_ligand_score =~ s/\.sdf/\.$data[$m][1]\.sdf/g;
			rename $sdf[$i], $out_ligand_score;
		}
		
	}
	
}

####################################################################
system("perl -pi -e 's/$year-$mon-$day $hour:$min:$sec\n//' $userDirPath\/status.txt");
system ("mv $outf\/receptor.pqbk $outf\/receptor.pdbqt");
system ("mv $outf\/ligand.pqbk $outf\/ligand.pdbqt");

