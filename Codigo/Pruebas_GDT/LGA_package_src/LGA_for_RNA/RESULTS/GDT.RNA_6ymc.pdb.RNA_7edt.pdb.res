
#######################################################
#                                                     #
#                        LGA                          #
#                  ---------------                    #
#                                                     #
#               Local-Global Alignment                #
#        A Method for Finding 3-D Similarities        #
#               in Protein Structures                 #
#                                                     #
#                  ------------ 09/2019               #
#                                                     #
#      Adam Zemla (zemla1@llnl.gov)                   #
#      Lawrence Livermore National Laboratory, CA     #
#                                                     #
#######################################################

# Molecule1: number of C4' atoms   52 ( 1118),  selected   52 , name RNA_6ymc.pdb
# Molecule2: number of C4' atoms   24 (  516),  selected   24 , name RNA_7edt.pdb
# PARAMETERS: -4  -d:4.0  -o0  -atom:C4,  -lga_m  -stral  GDT.RNA_6ymc.pdb.RNA_7edt.pdb  
# Search for Atom-Atom correspondence
# Structure alignment analysis 

#      Molecule1      Molecule2  DISTANCE  RMSD(lw:1)
LGA    -       -      X     1_A      -         -
LGA    X     1_A      X     2_A      #         -
LGA    X     2_A      X     3_A     3.247      -
LGA    X     3_A      X     4_A     2.529     0.199
LGA    X     4_A      X     5_A     2.835     0.179
LGA    X     5_A      X     6_A     4.287     0.706
LGA    X     6_A      X     7_A     1.029     1.136
LGA    X     7_A      X     8_A     1.794     0.386
LGA    X     8_A      X     9_A     1.019     0.172
LGA    X     9_A      X    10_A     1.168     0.106
LGA    X    10_A      X    11_A     1.307     0.140
LGA    X    11_A      X    12_A     0.908     0.124
LGA    X    12_A      -       -      -         -
LGA    X    13_A      -       -      -         -
LGA    X    14_A      -       -      -         -
LGA    X    15_A      -       -      -         -
LGA    X    16_A      X    13_B     2.680     0.150
LGA    X    17_A      X    14_B     1.976     0.209
LGA    X    18_A      X    15_B     1.774     0.070
LGA    X    19_A      X    16_B     1.947     0.069
LGA    X    20_A      X    17_B     2.514     0.274
LGA    X    21_A      X    18_B     2.186     0.202
LGA    X    22_A      X    19_B     1.221     0.606
LGA    X    23_A      X    20_B     1.415     0.322
LGA    X    24_A      X    21_B     0.916     0.292
LGA    X    25_A      X    22_B     0.980     0.196
LGA    X    26_A      X    23_B     1.287     1.432
LGA    X     1_B      X    24_B     3.032      -
LGA    X     2_B      -       -      -         -
LGA    X     3_B      -       -      -         -
LGA    X     4_B      -       -      -         -
LGA    X     5_B      -       -      -         -
LGA    X     6_B      -       -      -         -
LGA    X     7_B      -       -      -         -
LGA    X     8_B      -       -      -         -
LGA    X     9_B      -       -      -         -
LGA    X    10_B      -       -      -         -
LGA    X    11_B      -       -      -         -
LGA    X    12_B      -       -      -         -
LGA    X    13_B      -       -      -         -
LGA    X    14_B      -       -      -         -
LGA    X    15_B      -       -      -         -
LGA    X    16_B      -       -      -         -
LGA    X    17_B      -       -      -         -
LGA    X    18_B      -       -      -         -
LGA    X    19_B      -       -      -         -
LGA    X    20_B      -       -      -         -
LGA    X    21_B      -       -      -         -
LGA    X    22_B      -       -      -         -
LGA    X    23_B      -       -      -         -
LGA    X    24_B      -       -      -         -
LGA    X    25_B      -       -      -         -
LGA    X    26_B      -       -      -         -

#C4'           N1   N2   DIST      N    RMSD   Seq_Id      LGA_S     LGA_Q  S_nb  S_N   S_Id    LGA_M 
SUMMARY(LGA)   52   24    4.0     22    2.11   100.00     79.999     0.996     3   14 100.00   79.999 

Unitary ROTATION matrix and the SHIFT vector superimpose molecules  (1=>2)
  X_new =   0.959563 * X  +   0.274718 * Y  +  -0.061387 * Z  + -25.825214
  Y_new =  -0.183348 * X  +   0.775434 * Y  +   0.604223 * Z  + -17.899778
  Z_new =   0.213593 * X  +  -0.568535 * Y  +   0.794447 * Z  +  21.908461 

Euler angles from the ROTATION matrix. Conventions XYZ and ZXZ:
           Phi     Theta       Psi   [DEG:       Phi     Theta       Psi ]
XYZ: -0.188799 -0.215251 -0.621143   [DEG:  -10.8174  -12.3330  -35.5888 ]
ZXZ: -3.040343  0.652699  2.782217   [DEG: -174.1988   37.3969  159.4093 ]
 
# END of job
