#!/usr/bin/perl
#
# Author: Adam Zemla 
# Email: zemla1@llnl.gov
# Version: 06/19/2022
#
# runlga.model_sphere_library.pl
# script to evaluate molecule pairs (model.sphere) 
#
# Usage: runlga.model_sphere_library.pl model sphere library
#
# Examples:
#   predicting compound binding sites
#        ./runlga.model_sphere_library.pl  PDBspheres/nCoV_nsp3.6w9c_A.pdb  Sphere.Y97.35.7los_B.404_B.12.pdb          PDBspheres/PDBspheres.HET_compound_library.txt
#   predicting short-peptide binding sites
#        ./runlga.model_sphere_library.pl  PDBspheres/nCoV_nsp3.6w9c_A.pdb  Sphere.6wx4.6wx4_D.6wx4_I.5_3_4_30.12.pdb  PDBspheres/PDBspheres.PEP_peptides_library.txt
#

$|=1;

# subdirectory with executables
#$bin='bin/';
#$bin='./';
#$bin='';
$bin='./';

# executables
$collect_PDB="collect_PDB.pl";
$lgaexec="lga";

# --------------------------------------------------------------------- #

# usage 
if (@ARGV < 3) {
  printf "Usage:\t$0  model  sphere  library\n\n";
  printf "Examples:\n";
  printf "  predicting compound binding sites\n";
  printf "       ./runlga.model_sphere_library.pl  MODELS/nCoV_nsp3.pdb  Sphere.Y97.35.7los_B.404_B.12.pdb  ~/PDB/PDB_SPHERES/Pdb_HET.at_least_1_HETATM_and_5_CA.list\n";
  printf "       ./runlga.model_sphere_library.pl  MODELS/nCoV_nsp3.pdb  Sphere_HET.Y97.35.7los_B.404_B.12.pdb  PDB_SPHERES/Pdb_HET_PEP.spheres_12_5_CA.BindSite.txt\n";
  printf "  predicting short-peptide binding sites\n";
  printf "       ./runlga.model_sphere_library.pl  MODELS/nCoV_nsp3.pdb  Sphere.6wx4.6wx4_D.6wx4_I.5_3_4_30.12.pdb  ~/PDB/PDB_SPHERES/Pdb_PEP.peptide_25_and_sphere_5_CA.list\n";
  printf "       ./runlga.model_sphere_library.pl  MODELS/nCoV_nsp3.pdb  Sphere_PEP.6wx4.6wx4_D.6wx4_I.5_3_4_30.12.pdb  PDB_SPHERES/Pdb_HET_PEP.spheres_12_5_CA.BindSite.txt\n";
  printf "Models with predicted inserted ligands are stored in the folder RESULTS/:\n";
  printf "       RESULTS/nCoV_nsp3.pdb.Sphere_HET.Y97.35.7los_B.404_B.12.pdb.res\n";
  printf "       RESULTS/nCoV_nsp3.pdb.Sphere_PEP.6wx4.6wx4_D.6wx4_I.5_3_4_30.12.pdb.res\n";
  printf "Examples:\n";
  exit;
}

# subdirectories needed for LGA processing:
$pdbspheres = 'PDBspheres/';                                  # PDBspheres folder with lists and structures of template spheres used in calculations
$pdbspheresLib = 'PDBspheres/PDBspheresLibrary/';             # Folder where PDB template spheres used in calculations are stored
$LIBspheresLibraryList = 'PDBspheres/PDBspheresLibrary.list'; # List of PDB template spheres used in calculations
$dirres = 'RESULTS/';
$pdbs = 'MOL2/';
$tmp = 'TMP/';

print "\n";
# subdirectory with outputs from LGA. LGA puts calculation results to TMP directory
if(!-d $tmp) { mkdir($tmp); print "Making $tmp for output files\n"; }
# subdirectory with results from LGA. LGA results are stored in RESULTS directory
if(!-d $dirres) { mkdir($dirres); print "Making $dirres for output (*.res) files\n"; }
# subdirectory with input structures. LGA takes input mol1.mol2 structures from MOL2 directory
if(!-d $pdbs) { mkdir($pdbs); print "Making $pdbs for input files\n"; }
# subdirectory with input spheres structures. Directory with collected spheres structures from PDBspheres library
if(!-d $pdbspheres) { mkdir($pdbspheres); print "Making $pdbspheres for input files\n"; }
if(!-d $pdbspheresLib) { mkdir($pdbspheresLib); print "Making $pdbspheresLib for input files\n"; }

$res='.res';
$pdb='.pdb';
$lga='.lga';
$lig='.lig';
$stral='.stral';
$d='.';

# --------------------------------------------------------------------- #

# PDBspheres parameters:
$distconta="4.5";  # distance to report residues in contact with ligand (4.5)
$distclash="1.0";  # distance test settings for possible steric clashes (1.0)
#

# read input parameters and filenames
$proteinmdl=$ARGV[0];
$tmplsphere=$ARGV[1];
$LIBspheres=$ARGV[2];

@LINE=split(/\//,$proteinmdl);
@R=reverse(@LINE);
$mol1=$R[0];
@LINE=split(/\//,$tmplsphere);
@R=reverse(@LINE);
$mol2=$R[0];

# Check for HET or PEP sphere
chop($HET=`echo "$mol2" | cut -f3 -d'.'`);
$PEP = $HET =~ /[^0-9]/;
if($PEP == 0) {
  chop($ligand=`echo "$mol2" | cut -f5 -d'.'`); 
}
else {
  chop($ligand=`echo "$mol2" | cut -f4 -d'.'`); 
  chop($lchain=`echo "$mol2" | cut -f4 -d'.' | cut -f2 -d'_'`); 
}

# --------------------------------------------------------------------- #

# collect sphere information from the library
if(!-s "$pdbspheresLib$mol2") {
  chop($infsphere=`grep $mol2 $LIBspheres`);
  system "echo '$infsphere' >> $LIBspheresLibraryList";
}
else {
  chop($infsphere=`grep $mol2 $LIBspheresLibraryList`);
}
chop($tmplsphere=`echo "$infsphere" | cut -f2 -d' '`);
chop($sizsphere=`echo "$infsphere" | cut -f13 -d' '`); 
chop($ressphere=`echo "$infsphere" | sed 's/^.*SPHERE formed by positions: //' | cut -f1 -d' '`); 

if(!-s "$pdbspheresLib$mol2") {
  system "echo 'MOLECULE $mol2' > $pdbspheresLib$mol2";
  system "$bin$collect_PDB $tmplsphere -residues:$ressphere | egrep -e '^ATOM' | grep -v 'HOH' >> $pdbspheresLib$mol2";
  system "echo 'TER' >> $pdbspheresLib$mol2";
  if($PEP == 0) {
    system "$bin$collect_PDB $tmplsphere -residues:X$ligand | egrep -e '^ATOM' -e '^HETATM' | grep -v 'HOH' >> $pdbspheresLib$mol2";
  }
  else {
    system "$bin$collect_PDB $ligand | egrep -e '^ATOM' -e '^HETATM' | grep -v 'HOH' >> $pdbspheresLib$mol2";
  }
  system "echo 'END' >> $pdbspheresLib$mol2";
}

# Selection of DIST cutoff for LGA calculations
$dist=0.7;                           
if($sizsphere >= 250)    {$dist=5.0;}    
elsif($sizsphere >= 160) {$dist=4.5;} 
elsif($sizsphere >= 110) {$dist=4.0;} 
elsif($sizsphere >= 80)  {$dist=3.5;}  
elsif($sizsphere >= 70)  {$dist=3.0;}  
elsif($sizsphere >= 60)  {$dist=2.5;}  
elsif($sizsphere >= 50)  {$dist=2.0;}  
elsif($sizsphere >= 40)  {$dist=1.5;}  
elsif($sizsphere >= 30)  {$dist=1.0;}  

# Parameters for LGA calculations
$par="-4 -o1 -gdc -stral -lga_M -d:$dist -ie";

# create input structure file for LGA processing
$model="$mol1$d$mol2";
system "rm -rf $tmp$model$pdb";
system "echo 'MOLECULE $mol2' > $pdbs$model";
system "cat $pdbspheresLib$mol2 | grep -v '^MOLECULE ' >> $pdbs$model";
system "echo 'MOLECULE $mol1' >> $pdbs$model";
system "$bin$collect_PDB $proteinmdl | grep -v '^MOLECULE ' | grep -v '^END' >> $pdbs$model";
system "echo 'END' >> $pdbs$model";
sleep 1;

# run LGA - it knows to look in MOL2 for the input...
print "Running: lga $model $par\n";
system "$bin$lgaexec $model $par > $dirres$model$res";
if (-z "$dirres$model$res") { die "\nError: Problem running LGA program. Check $model\n"; }

# collect SUMMARY results
open(IN_res,"$dirres$model$res");
while(<IN_res>) {
  chomp($line=$_);
  if(/SUMMARY\(RMSD_GDC\): /) {
    @L1=split /[\ ]+/,$line;
    $scoreGDC=$L1[5];
  }
  elsif(/SUMMARY\(LGA\) /) {
    @L2=split /[\ ]+/,$line;
    $scoreNs=$L2[1];
    $scoreNt=$L2[2];
    $scoreN=$L2[4];
    $scoreRMSD=$L2[5];
    $scoreNc=$L2[10];
    $scoreSeqID=$L2[11];
    $scoreLGA=$L2[7];
  }
}
close(IN_res);
if($scoreNs>0) { 
  $scoreGDC=$scoreGDC*$scoreNt/$scoreNs;
  $scoreLGA=$scoreLGA*$scoreNt/$scoreNs;
}
else { 
  print "Warning: LGA didn\'t create SUMMARY results: $model$res (check input data)\n"; 
  exit;
}

# assembling models
system "cat $proteinmdl | grep -v '^END' >> $dirres$model$res";
if($PEP == 0) {
  system "$bin$collect_PDB $tmp$model$pdb -residues:X$ligand | egrep -e '^ATOM' -e '^HETATM' > $tmp$model$lig";
}
else {
  system "$bin$collect_PDB $tmp$model$pdb -chain:$lchain | egrep -e '^ATOM' -e '^HETATM' > $tmp$model$lig";
}
system "cat $tmp$model$lig >> $dirres$model$res";
system "echo 'END' >> $dirres$model$res";

# estimating contacts
$BSresidues=&acheck("$proteinmdl","$tmp$model$lig","$distconta");
$BSresclash=&acheck("$proteinmdl","$tmp$model$lig","$distclash");

@L5 = split /[\ ]+/,$BSresidues; @R5=reverse(@L5); @BSR = split /[\,]+/,$R5[0];
if($R5[0] =~ /,/) {
  $NBSR=$#BSR + 1; 
}
else {
  $NBSR=0;
}
@L6 = split /[\ ]+/,$BSresclash; @R6=reverse(@L6); @CLR = split /[\,]+/,$R6[0];
$clash_res="";
$NB_clash_res=0;
foreach $icl (0..$#CLR) {
  $lCLR=length($CLR[$icl]);
  if($CLR[$icl] !~ /X/ && $lCLR>0 && $CLR[$icl] !~ /Angstroms/) { # clashes with nonstandard aminoacids not counted
    $NB_clash_res=$NB_clash_res+1;
  }
}
$l=length("$model"); $l=$l+2;
$txt="#Structural_modelpdb.Sphere_template";
printf "%-$l\s  Ns   Nt    N   RMSD :  Nc   SeqID    LGA    GDC  N4  cl : Contact_residues\n",$txt;
printf "%-$l\s %3d %4d  %3d %6.2f : %3d  %6.2f %6.2f %6.2f %3d %3d : %s\n",$model,$scoreNs,$scoreNt,$scoreN,$scoreRMSD,$scoreNc,$scoreSeqID,$scoreLGA,$scoreGDC,$NBSR,$NB_clash_res,$R5[0]; 


# removing input and output files after the processing is done ...
system "rm -rf $pdbs$model $tmp$model$lga $tmp$model$pdb $tmp$model$lig $tmp$model$stral";

print "Done! \n";

exit;

sub acheck {
#
$mode=20;
$molecule=$_[0];
$atomfile=$_[1];
$distance=$_[2];

$dist2=$distance*$distance;

@AA_STD=("ALA","VAL","LEU","ILE","PRO","MET","PHE","TRP","GLY","SER","THR","CYS","TYR","ASN","GLN","ASP","GLU","LYS","ARG","HIS");
  
  $AA{"ALA"}='A';
  $AA{"VAL"}='V';
  $AA{"LEU"}='L';
  $AA{"ILE"}='I';
  $AA{"PRO"}='P';
  $AA{"MET"}='M';
  $AA{"PHE"}='F';
  $AA{"TRP"}='W';
  $AA{"GLY"}='G';
  $AA{"SER"}='S';
  $AA{"THR"}='T';
  $AA{"CYS"}='C';
  $AA{"TYR"}='Y';
  $AA{"ASN"}='N';
  $AA{"GLN"}='Q';
  $AA{"ASP"}='D';
  $AA{"GLU"}='E';
  $AA{"LYS"}='K';
  $AA{"ARG"}='R';
  $AA{"HIS"}='H';
  
  $AA{"XXX"}='X';
  
$ok=1;
      open(IN_pdb,$molecule) || {$ok=0};
      if($ok == 0) {
        print("ERROR: PDB molecule: $molecule Check the file name.\n");
        exit;
      }
      $n_atoms=0;
      $n_atom_records=0;
      while(<IN_pdb>) {
        chomp($line=$_);
        $RES_Chain=' ';
        if(/^ATOM / || /^HETATM/ || /^AAMOL1/ || /^AAMOL2/) {
          $line =~ s/[^a-zA-Z0-9\.\-\_\'\ ]/ /g;
          $_=substr($line,17,3);
          ($RES_Name)=/(\S+)/; 
          $RES_Name =~ s/ //g;
          $RES_Name =~ s/[^a-zA-Z0-9]//g;
          if($RES_Name ne 'HOH' && $RES_Name ne '') {
            if("@AA_STD" !~ /$RES_Name/ || length($RES_Name) != 3) {$RES_Name='XXX';}
            $_=substr($line,12,4);
            ($RES_Atom)=/(\S+)/; 
            $_=substr($line,22,5);
            ($RES_Number)=/(\S+)/; 
            $_=substr($line,21,1);
            ($RES_Chain)=/(\S+)/; 
            $RES_Chain =~ s/[^a-zA-Z0-9]//g;
            $n_atom_records1=$n_atom_records+1;
            if($RES_Chain ne ' ' && $RES_Chain ne '') {
              $pdb_resnumber[$n_atom_records1]="$RES_Number\_$RES_Chain";
              $curchain="$RES_Chain";
            }
            else {
              $pdb_resnumber[$n_atom_records1]="$RES_Number";
              $curchain=' ';
            }
            $n_atom_records=$n_atom_records1;
            $pdb_atom_records[$n_atom_records]="$line";

            $n_atoms=$n_atoms+1;
            $pdb_atom_line[$n_atoms]="$line";
            $pdb_atom[$n_atoms]="$RES_Atom";
            $pdb_name1[$n_atoms]="$RES_Name";
            $pdb_number[$n_atoms]="$RES_Number";
            $pdb_name1[$n_atoms]=$AA{$pdb_name1[$n_atoms]};
            if($RES_Chain ne ' ' && $RES_Chain ne '') {
              $pdb_chain[$n_atoms]="$RES_Chain";
              $pdbatom[$n_atoms]="$RES_Atom\.$RES_Name\.$RES_Number\_$RES_Chain";
              $pdbatomres1[$n_atoms]="$RES_Number\_$RES_Chain";
            }
            else {
              $pdb_chain[$n_atoms]=' ';
              $pdbatom[$n_atoms]="$RES_Atom\.$RES_Name\.$RES_Number";
              $pdbatomres1[$n_atoms]="$RES_Number";
            }
            $_=substr($line,30,8);
            ($c_x)=/(\S+)/;
            $pdbcrd_x[$n_atoms]=$c_x;
            $_=substr($line,38,8);
            ($c_y)=/(\S+)/;
            $pdbcrd_y[$n_atoms]=$c_y;
            $_=substr($line,46,8);
            ($c_z)=/(\S+)/;
            $pdbcrd_z[$n_atoms]=$c_z;
          }
        }
      }
      $pdbn_atoms=$n_atoms;
      close(IN_pdb);

$ok=1;
      open(IN_atoms,$atomfile) || {$ok=0};
      if($ok == 0) {
        print("ERROR: ATOM molecule: $atomfile Check the file name.\n");
        exit;
      }
      $n_atoms=0;
      while(<IN_atoms>) {
        chomp($line=$_);
        $RES_Chain=' ';
        if(/^ATOM / || /^HETATM/ || /^AAMOL1/ || /^AAMOL2/) {
          $line =~ s/[^a-zA-Z0-9\.\-\_\'\ ]/ /g;
          $_=substr($line,21,1);
          ($RES_Chain)=/(\S+)/;
          $RES_Chain =~ s/[^a-zA-Z0-9]//g;
          if($RES_Chain ne ' ' && $RES_Chain ne '') {
            $curchain="$RES_Chain";
          }
          else {
            $curchain=' ';
          }
          $n_atoms=$n_atoms+1;
          $_=substr($line,12,4);
          ($RES_Atom)=/(\S+)/; $a_atom[$n_atoms]="$RES_Atom";
          $_=substr($line,17,3);
          ($RES_Name)=/(\S+)/; 
          $RES_Name =~ s/ //g;
          $RES_Name =~ s/[^a-zA-Z0-9]//g;
          if($RES_Name eq 'HOH' || $RES_Name eq '') {
            $n_atoms=$n_atoms-1;
          }
          else {
            if("@AA_STD" !~ /$RES_Name/ || length($RES_Name) != 3) {$RES_Name='XXX';}
            $a_name[$n_atoms]="$RES_Name";
            $a_name1[$n_atoms]=$AA{$a_name[$n_atoms]};
            $_=substr($line,22,5);
            ($RES_Number)=/(\S+)/;
            $a_number[$n_atoms]="$RES_Number";
            if($RES_Chain ne ' ' && $RES_Chain ne '') {
              $a_chain[$n_atoms]="$RES_Chain";
              $aatom[$n_atoms]="$RES_Atom\.$RES_Name\.$RES_Number\_$RES_Chain";
            }
            else {
              $a_chain[$n_atoms]=' ';
              $aatom[$n_atoms]="$RES_Atom\.$RES_Name\.$RES_Number";
            }
            $_=substr($line,30,8);
            ($c_x)=/(\S+)/;
            $acrd_x[$n_atoms]=$c_x;
            $_=substr($line,38,8);
            ($c_y)=/(\S+)/;
            $acrd_y[$n_atoms]=$c_y;
            $_=substr($line,46,8);
            ($c_z)=/(\S+)/;
            $acrd_z[$n_atoms]=$c_z;
          }
        }
      }
      close(IN_atoms);

$k=0;
$all_k=0;
$old=0;
$prev='#';
$all_prev='#';
$res='#';
foreach $p (1 .. $pdbn_atoms) {
  if($res ne $pdb_number[$p]) {
    $old=$old+1;
    $res=$pdb_number[$p];
    $sdist=9999;
  }
  foreach $a (1 .. $n_atoms) {
    $x=$acrd_x[$a]-$pdbcrd_x[$p];
    $y=$acrd_y[$a]-$pdbcrd_y[$p];
    $z=$acrd_z[$a]-$pdbcrd_z[$p];
    $dist=$x*$x+$y*$y+$z*$z;
    if($dist2>=$dist) {
      $sdist_prev=$sdist;
      $sdist=sqrt($dist);
      @tmpcontact=split(/\./,$aatom[$a]);
      $pdb_name2[$p]=$AA{$tmpcontact[1]};
      $pdbatomres2[$p]=$tmpcontact[2];
      $newcontact="$pdb_name2[$p]$pdbatomres2[$p]";
      $sdist=$sdist;
      if($all_prev ne $pdbatom[$p]) {
        $all_k=$all_k+1;
        $all_prev=$pdbatom[$p];
        $all_atom_line[$all_k]=$pdb_atom_line[$p];
      }
      if($prev eq $pdb_number[$p] && $sdist_prev>$sdist) {
        $prev='#';
        $k=$k-1;
      }
      if($prev ne $pdb_number[$p]) {
        $prev=$pdb_number[$p];
        if($old > 1 and $k > 0) {
          $k=$k+1;
          $sequence[$k]='.';
          $number_res[$k]=' ';
          $atom_line[$k]=' ';
        }
        $k=$k+1;
        $sequence[$k]=$pdb_name1[$p];
        $number_res[$k]=$pdbatomres1[$p];
        $atom_line[$k]=$pdb_atom_line[$p];
        $old=0;
      }
    }  
  }  
}

$bsresout="";  
foreach $i (1 .. $k) {
  if($sequence[$i] ne '.') {
    $bsresout="$bsresout$sequence[$i]$number_res[$i],";
  }
}
$bsresout =~ s/ //g;

return $bsresout;
}
           
