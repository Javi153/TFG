#!/usr/bin/perl
#
# Author: Adam Zemla 
# Email: zemla1@llnl.gov
# Version: 06/19/2021
#
# collect_PDB.pl
# script to collect PDB structures.
#
# Usage: ./collect_PDB.pl  molecule [ outfilename ] [ -chain:A | -chain:A_3 ] [ -modelNB:3 ] [ -cleanhetatm ] [ -range:beg:end ] [ -residues:list ] [ -mutant:R.21_A,S.23_A ] [ -serialNB ] [ -allREMARK ]
#
# where: molecule - (input) name of the local molecule file or a code of PDB structure
#     outfilename - (output) name of the file where molecule will be saved (optional)
#           chain - chainID or chainID_modelNB (optional)
#         modelNB - model number (optional)
#     cleanhetatm - remove HETATM coordinates if they have the same  resnumber and chain  as in ATOM coordinates (optional)
#   range:beg:end - range=(beg,end) e.g. (11,121) (optional)
#   residues:list - list of residues e.g. check PDB file pdb1rhu.ent and residue list=F173_A,H174_A,K175_A,S175A_A,T175B_A,G175C_A,M176_A,T177_A,X501_A
#   mutant:R.21_A - modifying PDB by creating mutation e.g.: ARG at the position 21_A (backbone atoms are replaced to a new R amino acid, other atoms are removed) (optional)
#                   Side-chain atoms for a new mutation (R at the position 21_A can be created by running: scwrl -i molecule_mut.pdb -o molecule_mut_sc.pdb > molecule_mut.log
#                   or by running: Create_Fixed_Mutated.InpStr_OutStr.Muts.sh  molecule  mutated.pdb  R.7_A,S.23_A
#                   further "light" refinement can be performed by running: chimera.minimize_structure.inpmol_outmol.sh  mutated.pdb  mutated.Emin.pdb
#                   NOTE: it also accepts mutations in the format:  KA7R,TA23S
#        serialNB - renumber serial number (optional)
#       allREMARK - print all REMARK lines only (optional)
#
#        PDB code:  1abc  or  1abc_B (if chain B is requested)
#           if  1ja1_*  or  1JA1_*  then all chains are collected
#           if  1ja1_A  or  1JA1A   then chain A is collected
#           if  1JA1_   or  1ja1    then chain ' ' is collected
#
#        in the case of NMR models:
#               1bve_B_5   for PDB entry: 1bve, chain: 'B', modelNB: 5
#               1awo___7   for PDB entry: 1awo, chain: ' ', modelNB: 7
#
#        local files can be specified using full path:
#               DIR/molecule
#

$|=1;

$pdblocal='-';
$chain='no';
$mdl=0;
$cleanhetatm=0;
$serialNB=0;
$allREMARK=0;
$mutant=0;
$respos=0;
$resposOK=1;
$beg=-9999;
$end=9999;
$u='_';
if($#ARGV < 0) {
  print"USAGE: collect_PDB.pl  molecule  [ outfilename ] [ -chain:A | -chain:A_3 ] [ -modelNB:3 ] [ -cleanhetatm ] [ -range:beg:end ] [ -residues:list ] [ -mutant:R.21_A,S.23_A ] [ -serialNB ] [ -allREMARK ]\n";
  print"  where: molecule      - (input) name of the local molecule file or a code of PDB structure\n";
  print"         outfilename   - (output) name of the file where molecule will be saved (optional)\n";
  print"         chain         - chainID or chainID_modelNB (optional)\n";
  print"         modelNB       - model number (optional)\n";
  print"         cleanhetatm   - remove HETATM coordinates if they have the same  resnumber and chain  as in ATOM coordinates (optional)\n";
  print"         range:beg:end - range=(beg,end) e.g. (11,121) (optional)\n";
  print"         residues:list - list of residues e.g. check PDB file pdb1rhu.ent and residue list=F173_A,H174_A,K175_A,S175A_A,T175B_A,G175C_A,M176_A,T177_A,X501_A\n";
  print"         mutant:R.21_A - modifying PDB by creating mutation e.g.: LYS at the position 21_A is replaced by ARG (optional)\n";
  print"                         (backbone atoms are replaced to a new ARG amino acid, old LYS side chain atoms are removed)\n";
  print"                         Side-chain atoms for a new mutation (R at the position 21_A can be created by running:\n";
  print"                         scwrl -i molecule_mut.pdb -o molecule_mut_sc.pdb > molecule_mut.log\n";
  print"                         or by running: Create_Fixed_Mutated.InpStr_OutStr.Muts.sh  molecule  mutated.pdb  R.7_A,S.23_A\n";
  print"                         Further 'light' refinement can be performed by running:\n";
  print"                         chimera.minimize_structure.inpmol_outmol.sh  mutated.pdb  mutated.Emin.pdb\n";
  print"                         NOTE: it also accepts mutations in the format:  KA7R,TA23S\n";
  print"         serialNB      - renumber serial number (optional)\n";
  print"         allREMARK     - print all REMARK lines only (optional)\n";
  exit;
}
else {
  $molecule=$ARGV[0];
  if($#ARGV > 0) {
    foreach $K (1 .. $#ARGV) {
      if($ARGV[$K] =~ /\-chain:/) {
        ($tmp,$chain) = split /[\:\0\n]+/,$ARGV[$K];
        if($chain =~ /_/) {
          ($chain,$mdl) = split /[\_\0\n]+/,$chain;
          if($chain eq '' || $chain eq '.' || $chain eq '*') { $chain=' '; }
        }
      }
      elsif($ARGV[$K] =~ /\-modelNB:/) {
        ($tmp,$mdl) = split /[\:\0\n]+/,$ARGV[$K];
        if($mdl !~ /[0-9]/) { $mdl=0; }
      }
      elsif($ARGV[$K] =~ /\-residues:/) {
        ($tmp,$residues) = split /[\:\0\n]+/,$ARGV[$K];
        (@Listres) = split /[\,\0\n]+/,$residues;
        if($#Listres >= 0) {
          foreach $i (0 .. $#Listres) {
            $respos=$respos+1;
            $ListResPos{$Listres[$i]}='ok';
          }
        }
      }
      elsif($ARGV[$K] =~ /\-mutant:/) {
        ($tmp,$newmuts) = split /[\:\0\n]+/,$ARGV[$K];
        (@Listmuts) = split /[\,\0\n]+/,$newmuts;
        if($#Listmuts >= 0) {
          foreach $i (0 .. $#Listmuts) {
            $mutant=$mutant+1;
            (@R) = split /[\.\_\0\n]+/,$Listmuts[$i];
            if($#R>0) {
              $MUT_aa[$mutant]=$R[0];
              $MUT_num[$mutant]=$R[1];
              if($#R > 1) {
                $MUT_chain[$mutant]=$R[2];
              }
              else {
                $MUT_chain[$mutant]=' ';
              }
              $MUT_resnum[$mutant]="$R[1]\_$MUT_chain[$mutant]";
            }
            else {
              $l=length($R[0]);
              $motseq=substr($REFMOT[$i],$beg,$add);
              $MUT_aa_orig_in[$mutant]=substr($R[0],0,1);
              $MUT_chain[$mutant]=substr($R[0],1,1);
              $MUT_num[$mutant]=substr($R[0],2,$l-3);
              $MUT_aa[$mutant]=substr($R[0],$l-1,1);
              $MUT_resnum[$mutant]="$MUT_num[$mutant]\_$MUT_chain[$mutant]";
            }
          }
        }
      }
      elsif($ARGV[$K] =~ /\-cleanhetatm/) {
        $cleanhetatm=1;
      }
      elsif($ARGV[$K] =~ /\-serialNB/) {
        $serialNB=1;
      }
      elsif($ARGV[$K] =~ /\-allREMARK/) {
        $allREMARK=1;
      }
      elsif($ARGV[$K] =~ /\-range:/) {
        ($tmp,$beg,$end) = split /[\:\0\n]+/,$ARGV[$K];
        if($beg < -9999) {
          $beg=-9999;
        }
        if($end > 9999) {
          $end=9999;
        }
      }
      else {
        $pdblocal=$ARGV[$K];
      }
    }
  }
}

# Subdirectories where program searches for structures
$home_dir=$ENV{"HOME"};
$dirpdb="/PDB/structures";                       # default location of standard PDB entries
$dirpdbopt="PDB/structures";                     # optional location of standard PDB entries
$dirlocalpdb='PDB_local';                        # default location of local PDB structures 
$diruserspdb='PDB_users';                        # default location of users PDB structures 
#

if($mutant > 0 || $respos > 0) {
  $resname{'ALA'}='A'; $resname{'A'}='ALA';
  $resname{'VAL'}='V'; $resname{'V'}='VAL';
  $resname{'LEU'}='L'; $resname{'L'}='LEU';
  $resname{'ILE'}='I'; $resname{'I'}='ILE';
  $resname{'PRO'}='P'; $resname{'P'}='PRO';
  $resname{'MET'}='M'; $resname{'M'}='MET';
  $resname{'PHE'}='F'; $resname{'F'}='PHE';
  $resname{'TRP'}='W'; $resname{'W'}='TRP';
  $resname{'GLY'}='G'; $resname{'G'}='GLY';
  $resname{'SER'}='S'; $resname{'S'}='SER';
  $resname{'THR'}='T'; $resname{'T'}='THR';
  $resname{'CYS'}='C'; $resname{'C'}='CYS';
  $resname{'TYR'}='Y'; $resname{'Y'}='TYR';
  $resname{'ASN'}='N'; $resname{'N'}='ASN';
  $resname{'GLN'}='Q'; $resname{'Q'}='GLN';
  $resname{'ASP'}='D'; $resname{'D'}='ASP';
  $resname{'GLU'}='E'; $resname{'E'}='GLU';
  $resname{'LYS'}='K'; $resname{'K'}='LYS';
  $resname{'ARG'}='R'; $resname{'R'}='ARG';
  $resname{'HIS'}='H'; $resname{'H'}='HIS';
  $resname{'UNK'}='X'; $resname{'X'}='UNK';
  $resname{'XXX'}='X';
}

# foreach $i (1 .. $mutant) {
#   printf "MUTANT: $MUT_aa[$i] , $MUT_num[$i] , $MUT_chain[$i] ,\n";
# }
# exit;

$ok=1;
if($pdblocal eq '' || $pdblocal eq ' ' || $pdblocal eq '>'  || $pdblocal eq '>>' || $pdblocal eq '\>' || $pdblocal eq '\>\>') {
  $pdblocal='-';
}
else {
  if(-s $pdblocal > 100) { $ok=0; }
}
if($ok==1) {
  open(OUT,">$pdblocal") || {$ok=0};
  if($ok==0) {
    print "ERROR: (WRITE) bad file name: $pdblocal \n";
    exit;
  }
  else {
    $sp=' ';
    $nmr=$mdl;
    $one='*';
    $hetatm=1;
    $l=length($molecule);
    if($molecule =~ /\// || $l>=10 || $l<4 || ($l==4 && $molecule =~ /\_/)) {
      @LINE=split(/\//,$molecule);
      @R=reverse(@LINE);
      $model=$R[0];
      $sub0="$molecule";
      $ok=1;
      if(-s $sub0 < 100) { $ok=0; } else { open(IN,$sub0) || {$ok=0}; }
    }
    else {             # checking $dirpdb, $dirlocalpdb and ./ for standard PDB files: pdbnccc.ent
      $model="$molecule";
      $pdb='pdb';
      $ent='.ent';
      $un='_';

      $two=substr($model,1,2);
      $two=~tr/A-Z/a-z/;
      $four=substr($model,0,4);
      $four=~tr/A-Z/a-z/;
      $ch5=substr($model,4,1);
      $ch7=substr($model,6,1);

      if($l>6) {
        if($ch5 eq '_' && $ch7 eq '_') {
          $one=substr($model,5,1);
        }
      }
      elsif($l==6) {
        if($ch5 eq '_' || $ch5 eq ':' ||  $ch5 eq '.') {
          $one=substr($model,5,1);
        }
      }
      else {
        if($l==5) {
          $one=substr($model,4,1);
        }
        else {
          $one=' ';
        }
      }
      if($one eq '_' || $one eq '-' || $one eq '.' || $one eq ':') {
        $one=' ';
      }
      if($l>=8) {
        $nmr=substr($model,7,$l-7);
      }
      if($nmr eq '-') {
        $nmr=0;
      }
      if($one =~ /[a-zA-Z0-9]/) { $c=$one; }
      else { $c=' '; }

      if($l<8) {
        $sub0="$dirlocalpdb/$four";            # four letter code (PDB entry) : nccc
        $sub1="$dirlocalpdb/$four$un$c";       # six letter code_chain (PDB entry): nccc_X
        $sub2="$dirlocalpdb/$pdb$four$ent";    # standard PDB files: pdbnccc.ent
        $sub3="$dirpdb/$pdb$four$ent";
        $sub4="$dirpdb/$two/$pdb$four$ent";
        $sub5="$dirpdbopt/$pdb$four$ent";
        $sub6="$dirpdbopt/$two/$pdb$four$ent";
        $sub7="$home_dir$dirpdb/$pdb$four$ent";
        $sub8="$home_dir$dirpdb/$two/$pdb$four$ent";
        $sub9="$pdb$four$ent";                 # standard PDB files: pdbnccc.ent
      }
      else {
        $sub0="$dirlocalpdb/$pdb$four$ent";    # standard PDB files: pdbnccc.ent
        $sub1="$dirpdb/$pdb$four$ent";
        $sub2="$dirpdb/$two/$pdb$four$ent";
        $sub3="$dirpdbopt/$pdb$four$ent";
        $sub4="$dirpdbopt/$two/$pdb$four$ent";
        $sub5="$home_dir$dirpdb/$pdb$four$ent";
        $sub6="$home_dir$dirpdb/$two/$pdb$four$ent";
        $sub7="$dirlocalpdb/$four";            # four letter code (PDB entry) : nccc
        $sub8="$dirlocalpdb/$four$un$c";       # six letter code_chain (PDB entry): nccc_X
        $sub9="$pdb$four$ent";                 # standard PDB files: pdbnccc.ent
      }

      $ok=1;
      if(-s $sub0 < 100) { $ok=0; } else { open(IN,$sub0) || {$ok=0}; }
      if($ok==0) {
        $ok=1;
        if(-s $sub1 < 100) { $ok=0; } else { open(IN,$sub1) || {$ok=0}; }
        if($ok==0) {
          $ok=1;
          if(-s $sub2 < 100) { $ok=0; } else { open(IN,$sub2) || {$ok=0}; }
          if($ok==0) {
            $ok=1;
            if(-s $sub3 < 100) { $ok=0; } else { open(IN,$sub3) || {$ok=0}; }
            if($ok==0) {
              $ok=1;
              if(-s $sub4 < 100) { $ok=0; } else { open(IN,$sub4) || {$ok=0}; }
              if($ok==0) {
                $ok=1;
                if(-s $sub5 < 100) { $ok=0; } else { open(IN,$sub5) || {$ok=0}; }
                if($ok==0) {
                  $ok=1;
                  if(-s $sub6 < 100) { $ok=0; } else { open(IN,$sub6) || {$ok=0}; }
                  if($ok==0) {
                    $ok=1;
                    if(-s $sub7 < 100) { $ok=0; } else { open(IN,$sub7) || {$ok=0}; }
                    if($ok==0) {
                      $ok=1;
                      if(-s $sub8 < 100) { $ok=0; } else { open(IN,$sub8) || {$ok=0}; }
                      if($ok==0) {
                        $ok=1;
                        if(-s $sub9 < 100) { $ok=0; } else { open(IN,$sub9) || {$ok=0}; }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    $model=~s/\*//g;

# checks for models in the current (./) , $dirlocalpdb and $diruserspdb subdirectories

    if($ok==0) {

      $sub1=$model;
      $sub2="$dirpdb/$model";
      $sub3="$home_dir$dirpdb/$model";
      $sub4="$dirlocalpdb/$model";
      $sub5="$diruserspdb/$model";
      $sub6="$diruserspdb/$model\.pdb";
      $sub7="PDB/$dirlocalpdb/$model";
      $sub8="PDB/$diruserspdb/$model";
  
      $ok=1;
      $nmr=$mdl;
      $one='*';

      if(-s $sub1 < 100) { $ok=0; } else { open(IN,$sub1) || {$ok=0}; }
      if($ok==0) {
        $ok=1;
        if(-s $sub2 < 100) { $ok=0; } else { open(IN,$sub2) || {$ok=0}; }
        if($ok==0) {
          $ok=1;
          if(-s $sub3 < 100) { $ok=0; } else { open(IN,$sub3) || {$ok=0}; }
          if($ok==0) {
            $ok=1;
            if(-s $sub4 < 100) { $ok=0; } else { open(IN,$sub4) || {$ok=0}; }
            if($ok==0) {
              $ok=1;
              if(-s $sub5 < 100) { $ok=0; } else { open(IN,$sub5) || {$ok=0}; }
              if($ok==0) {
                $ok=1;
                if(-s $sub6 < 100) { $ok=0; } else { open(IN,$sub6) || {$ok=0}; }
                if($ok==0) {
                  $ok=1;
                  if(-s $sub7 < 100) { $ok=0; } else { open(IN,$sub7) || {$ok=0}; }
                  if($ok==0) {
                    $ok=1;
                    if(-s $sub8 < 100) { $ok=0; } else { open(IN,$sub8) || {$ok=0}; }
                    if($ok==0) {
                      print OUT "ERROR: (READ) bad model file name: $model \n";
                      close(OUT);
                      exit;
                    }
                  }
                }
              }
            }
          }
        }
      }

    }

# reads model

    $NBatm=0;
    if($ok==1) {
      if($chain ne 'no') { $one=$chain; }
      printf OUT "MOLECULE  $model \n";
      $model_NB=0;
      $residue_name_prev="#####";
      while(<IN>) {
        chop($line=$_);
        if(/^MODEL /) {
          $model_NB=$model_NB + 1;
          ($tmp,$model_NB)=/(\S+)\s+(\S+)/;
          if($nmr == 0) {
            $nmr=$model_NB;
          }
          if($model_NB == $nmr) {
            printf OUT "$line\n";
          }
        }
        elsif(/^TER/ && $one eq '*') {
          if($model_NB == $mdl) {
            printf OUT "$line\n";
          }
        }
        elsif(/^ATOM / || /^HETATM/ || /^TER /) {
          $print_ok=0;
          $fixed=0;
          $b_name=substr($line,0,6);
          $s_name=substr($line,6,5);
          $a_name=substr($line,12,4);
          $a_altp=substr($line,16,1);
          $r_name=substr($line,17,3);
          $c_name=substr($line,21,1);
          $r_numb=substr($line,22,5);
          $r_insr=substr($line,26,1);
          $l=length($line);
          if($serialNB == 1) {
            $NBatm=$NBatm+1;
            $left=sprintf("%-6s%5d %4s",$b_name,$NBatm,$a_name);
          }
          else {
            $left=substr($line,0,16);
          }
          $middle=substr($line,17,9);
          $right=substr($line,27,$l-27);
          $atom="$model:$a_name:$r_name:$r_numb:$c_name";
          $atom=~s/ //g;
          $residue_name="$model:$r_name:$c_name";
          $residue_name=~s/ //g;
          $residue_numb_altp="$model:$a_altp:$r_numb:$c_name";
          $residue_numb="$model:$r_numb:$c_name";
          $residue_numb=~s/ //g;
          if($respos>0) {
            $r_npos=$resname{$r_name}; if($r_npos eq '') { $r_npos='X'; }
            if($c_name eq '' || $c_name eq ' ') {
              $residue_npos="$r_npos$r_numb";
            }
            else {
              $residue_npos="$r_npos$r_numb$u$c_name";
            }
            $residue_npos=~s/ //g;
            if($ListResPos{$residue_npos} eq 'ok') { $resposOK=1; } else { $resposOK=0; }
          }
# Check $RESNUMB_OK{$residue_numb} to fix the problem with duplicated numbering (fix by marking "insertion")
          if($ATOMS{$atom} ne $atom || $RESNUMB_OK{$residue_numb_altp} eq 'ok') {
            if($RESNUMB{$residue_numb} eq $residue_numb && $residue_name ne $residue_name_prev) {
              if($r_insr eq ' ' && $a_altp eq ' ') {
                $r_insr='A'; $fixed=1;
              }
            }
            if($RESNUMB{$residue_numb} ne $residue_numb || $residue_name eq $residue_name_prev || $fixed == 1) {
              $line="$left$sp$middle$r_insr$right";
              if(($c_name eq $one || $one eq '*') && $model_NB == $nmr) {
                $print_ok=1;
              }
              elsif(/^HETATM/ && $hetatm == 1 && $c_name eq ' ') {
                $print_ok=1;
              }
            }
            if($cleanhetatm==1) {
              if(/^ATOM/) {
                $RESNUMB_got{$residue_numb}='got';
              }
              if(/^HETATM/ && $RESNUMB_got{$residue_numb} eq 'got') {
                $print_ok=0;
              }
            }
            if($r_name =~ /GLY/ && $a_name =~ /CB/) {
              $print_ok=0;
            }
            if($print_ok==1) {
              $ATOMS{$atom}="$atom";
              if($RESNUMB{$residue_numb} ne $residue_numb) {
                $residue_name_prev=$residue_name;
                $RESNUMB{$residue_numb}=$residue_numb;
              }
              if($mutant>0) {
                (@AT) = split /[\:\0\n]+/,$atom;
                $acheck=0;
                foreach $i (1 .. $mutant) {
                  if($MUT_num[$i] eq $AT[3] && $MUT_chain[$i] eq $AT[4]) {
                    $acheck=$i;
                    $new_r_aa=$MUT_aa[$i];
                  }
                }
                if($acheck>0) {
                  if($AT[1] eq 'N' || $AT[1] eq 'CA' || $AT[1] eq 'C' || $AT[1] eq 'O' || $AT[1] eq 'O1' || $AT[1] eq 'OC1' || $AT[1] eq 'OT1') {
                    $new_r_name=$resname{$new_r_aa};
                    $middle=~s/$r_name/$new_r_name/;
                    $line="$left$sp$middle$r_insr$right";
                    $line=~s/ O1 / O  /;
                    $line=~s/ OC1/ O  /;
                    $line=~s/ OT1/ O  /;
                  }
                  else {
                    $r_numb=-10000;
                  }
                }
              }
              if($r_numb>=$beg && $r_numb<=$end && $resposOK==1) {
                printf OUT "$line\n";
              }
            }
            $RESNUMB_OK{$residue_numb_altp}='ok';
          }
        }
        elsif($allREMARK == 1) {
          printf OUT "$line\n";
        }
        elsif(/^TITLE / ||
              /^LGA / ||
              /^AAMOL/ ||
              /^OBSLTE / ||
              /^DBREF / ||
              /^CONECT / ||
              /^HET/ ||
              /^HETNAM/ ||
              /^HETSYN/ ||
              /^FORMUL / ||
              /^MASTER / ||
              /^EXPDTA / ||
              /^AUTHOR / ||
              /^JRNL / ||
              /^REVDAT / ||
              /^COMPND / ||
              /^SEQRES / ||
              /^SOURCE / ||
              /^KEYWDS / ||
              /^CRYST1 / ||
              /^REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT:/ ||
              /^REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE:/ ||
              /^REMARK   2[ ]+RESOLUTION/ ||
              /^REMARK  Model name: / ||
              /^REMARK  Seq: / ||
              /^REMARK  Templates: / ||
              /^REMARK  AS2TS name: / ||
              /^REMARK  AS2TS score: / ||
              /^HEADER /) {
          printf OUT "$line\n";
        }
      }
      close(IN);
      if($allREMARK == 0) {
        printf OUT "END \n";
      }
    }

  }
}
close(OUT);
exit;

