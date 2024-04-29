#!/bin/bash
#
### Estimating GDT_HA and GDT_TS from comparison of different RNA structures using LGA structural alignment:
#
# To evaluate structure similarity between two RNA structures we use "C4'" atoms instead of CA (calpha) as representatives
#
# Usage: ./run_GDT_for_RNA_structures_with_unknown_residue_residue_correspondences.sh  pdb1  pdb2  [ OUT ] 
#
# Output OUT: 0 - no coordinates, 1 - rotated first structure, 2 - two aligned structures
#

if [ -z "$3" ]; then
  OUT=0
else
  OUT=$3
fi

mol1=`echo $1 | rev | cut -f1 -d'/' | rev`
mol2=`echo $2 | rev | cut -f1 -d'/' | rev`

molsprocess=GDT.$mol1\.$mol2

bin/collect_PDB.pl $1 > MOL2/$molsprocess
bin/collect_PDB.pl $2 >> MOL2/$molsprocess

# Calculate structural alignment on "C4'" atoms
bin/lga -4 -d:4.0 -o0 -atom:C4, -lga_m -stral $molsprocess > RESULTS/$molsprocess.res

# Use LGA alignment records to select residue-residue correspondences for GDT calculations
grep "^LGA " RESULTS/$molsprocess.res >> MOL2/$molsprocess
# Use calculated residue-residue correspondences in GDT evaluation
bin/lga -3 -sia -d:5.0 -atom:C4, -o$OUT -al $molsprocess > RESULTS/$molsprocess.gdt_res

# Collect data from calculations
NB12=`grep "^SUMMARY(LGA)" RESULTS/$molsprocess.res | perl -ane '@L=split(/[\ ]+/,$_); printf "%s:%s",$L[1],$L[2]'`
NB2=`echo $NB12 | cut -f2 -d':'`
NB=`grep "^SUMMARY(GDT)" RESULTS/$molsprocess.gdt_res | perl -ane '@L=split(/[\ ]+/,$_); printf "%s",$L[2]'`
GDT=`cat RESULTS/$molsprocess.gdt_res | grep "GDT PERCENT_AT" | perl -ane '@LINE=(@F);$V1=($LINE[2]+$LINE[3]+$LINE[5]+$LINE[9])/4.0;$V2=($LINE[3]+$LINE[5]+$LINE[9]+$LINE[17])/4.0;printf "%6.2f:%6.2f",$V1,$V2;'`

# Create output
cat RESULTS/$molsprocess.res | grep -v "# END of job"
cat RESULTS/$molsprocess.gdt_res | grep -v "#     " | grep -v "LGA    " | grep -v "ROTATION" | grep -v "_new =" | grep -v "DEG: " | grep -v "# END of job"

echo "#######################################################"
echo "#######################################################"
echo ""
echo "# Values of GDT_HA and GDT_TS from residue-residue correspondences established by LGA sequence independent analysis"
echo "# GDT_HA and GDT_TS values below are normalized according to the size ($NB2) of the reference (target) structure:"
echo "$NB12 $NB $GDT" | perl -ane '@L=split(/[\ \:]+/,$_); printf "GDT_HA = %6.2f  GDT_TS = %6.2f\n\n",$L[3]*$L[2]/$L[1],$L[4]*$L[2]/$L[1];'
echo "# END of job"

if [ $OUT -gt 0 ]; then
  grep -v "^LGA " TMP/$molsprocess.pdb | egrep -v -e "REMARK   GDT and LCS analysis" -e "REMARK   FIXED Atom-Atom correspondence" -e "REMARK   LGA parameters:" -e "REMARK   #CA        N1   N2" -e "REMARK   SUMMARY:" > Superimposed.$molsprocess.pdb
fi

