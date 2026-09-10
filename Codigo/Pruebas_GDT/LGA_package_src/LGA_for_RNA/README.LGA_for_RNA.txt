The LGA structure alignment program can be used to run structure similarity analysis on any selected atom, not just Calpha
(CA) atoms. For example, by using the following option:

"-atom:C4,"

we run structure comparisons not on Calpha atoms but on C4' atoms. 
Please note that in this example for the atom name "C4'" we use the format with a comma (,) instead of apostrophe (')
(the apostrophe is a special character which can cause some problems in scripting so this is why we use comma in such cases). 

Examples of running LGA calculations on "C4'" atoms:

cd LGA_for_RNA

# using a standard script for molecule-molecule comparison (CASP type evaluation when the residue-residue correspondences are fixed):
../runlga.mol_mol.pl RNA_6ymc.pdb RNA_7edt.pdb -3 -d:4 -atom:C4, -stral -o2

# using a standard script for molecule-molecule comparison:
../runlga.mol_mol.pl RNA_6ymc.pdb RNA_7edt.pdb -4 -d:4 -atom:C4, -stral -o2

# using a script prepared for processing RNA structures when residue-residue correspondences will be automatically established from the "optimal" structural alignment:
./run_GDT_for_RNA_structures_with_unknown_residue_residue_correspondences.sh RNA_6ymc.pdb RNA_7edt.pdb 
./run_GDT_for_RNA_structures_with_unknown_residue_residue_correspondences.sh RNA_6ymc.pdb RNA_7edt.pdb 1
./run_GDT_for_RNA_structures_with_unknown_residue_residue_correspondences.sh RNA_7eqj.pdb RNA_7edt.pdb 2 > Superimposed.GDT.RNA_7eqj.pdb.RNA_7edt.pdb.Alignment_GDT_results.txt

###
For more details how to use the "-atom" option please check the README file at:
http://proteinmodel.org/AS2TS/LGA/lga_format.html

It is important to remember, that GDT calculations (option "-3" in the LGA program) have to be applied to already 
known residue-residue correspondences, so when we need to compare structures for which we don't know residue-residue 
correspondences we need to identify these correspondences first. It can be done manually by creating lists of
corresponding residues from both structures, by using parameters: -aa1 -aa2 -er1 -er2 or -sda, or we can try automated 
approach which uses sequence independent option "-4" in the LGA program to establish such correspondences. 
A script listed above "run_GDT_for_RNA_structures_with_unknown_residue_residue_correspondences.sh"
facilitates such "automated" identification of residue-residue correspondences.

However, a user needs to be aware that this automated approach could not be exactly what may be needed (e.g. for CASP
type evaluation; e.g. using option "-3") because it identifies and uses the residue correspondences established by "optimal"
structural alignment done by LGA with option "-4" which searches for the best alignment and superposition between two molecules
regardless of their residue numbering and desired or expected (fixed) residue-residues correspondences.

