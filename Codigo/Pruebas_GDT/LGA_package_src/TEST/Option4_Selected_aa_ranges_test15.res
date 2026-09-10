
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

# Molecule1: number of CA atoms   35 ( 1749),  selected   35 , name 1wdn_A
# Molecule2: number of CA atoms   29 ( 1714),  selected   29 , name 1ggg_A
# PARAMETERS: -4  -o2  -gdc  test15  -al  
# Search for Atom-Atom correspondence
# Structure alignment analysis 

# WARNING! The change of the parameter DIST cutoff may give you better result.

# Checking swapping

#      Molecule1      Molecule2  DISTANCE    Mis    MC     All    Dist_max   GDC_mc  GDC_all
LGA    G    52_A      N   141_A     0.703     0    0.150   0.150     1.277   88.214   88.214
LGA    I    53_A      A   142_A     0.686     0    0.034   0.040     0.766   90.476   90.476
LGA    I    54_A      Y   143_A     0.881     0    0.063   0.068     1.086   90.476   88.667
LGA    P    55_A      M   144_A     0.941     0    0.116   0.151     1.404   83.690   85.048
LGA    A    56_A      E   145_A     0.680     0    0.049   0.049     0.714   90.476   90.476
LGA    L    57_A      L   146_A     0.601     0    0.115   0.605     2.356   90.476   82.917
LGA    Q    58_A      G   147_A     1.223     0    0.103   0.103     1.223   90.595   90.595
LGA    T    59_A      T   148_A     0.944     0    0.692   0.930     4.046   82.143   71.905
LGA    K    60_A      N   149_A     1.587     0    0.088   0.088     2.068   75.119   74.667
LGA    N    61_A      R   150_A     1.446     0    0.070   0.071     1.644   79.286   78.000
LGA    V    62_A      A   151_A     0.801     0    0.087   0.128     1.281   95.238   92.476
LGA    D    63_A      D   152_A     0.715     0    0.114   0.254     1.916   90.476   84.881
LGA    L    64_A      A   153_A     0.884     0    0.104   0.138     1.199   85.952   85.048
LGA    A    65_A      V   154_A     1.960     0    0.110   0.129     2.171   75.000   72.952
LGA    L    66_A      L   155_A     1.617     0    0.030   0.987     2.979   66.905   66.964
LGA    A    67_A      H   156_A     3.167     0    0.566   0.638     4.161   50.476   49.048
LGA    G    68_A      D   157_A     3.242     0    0.435   0.435     3.595   46.667   46.667
LGA    I    69_A      -       -      -        -     -       -         -        -        -
LGA    T    70_A      T   158_A     2.814     0    0.608   0.625     6.153   53.571   42.041
LGA    I    71_A      P   159_A     4.560     0    0.639   0.598     6.834   45.476   39.048
LGA    T    72_A      N   160_A     5.738     0    0.631   0.621     6.497   22.738   21.619
LGA    D    73_A      -       -      -        -     -       -         -        -        -
LGA    E    74_A      -       -      -        -     -       -         -        -        -
LGA    R    75_A      -       -      -        -     -       -         -        -        -
LGA    K    76_A      -       -      -        -     -       -         -        -        -
LGA    K    77_A      -       -      -        -     -       -         -        -        -
LGA    A    78_A      Q   172_A     3.590     0    0.056   0.061     4.527   50.119   46.381
LGA    I    79_A      F   173_A     1.322     0    0.023   0.029     2.016   75.119   74.667
LGA    D    80_A      K   174_A     1.128     0    0.037   0.057     1.874   88.214   85.143
LGA    F    81_A      A   175_A     0.490     0    0.043   0.040     0.847   97.619   96.190
LGA    S    82_A      V   176_A     0.851     0    0.120   0.176     2.947   80.119   75.524
LGA    D    83_A      G   177_A     2.685     0    0.697   0.697     2.685   69.048   69.048
LGA    G    84_A      D   178_A     1.417     0    0.336   0.336     2.244   77.381   77.381
LGA    Y    85_A      S   179_A     1.432     0    0.653   0.598     4.266   81.548   72.667
LGA    Y    86_A      L   180_A     3.819     0    0.130   0.135     4.741   44.167   41.619

# RMSD_GDC results:       CA      MC common percent     ALL common percent   GDC_mc  GDC_all
NUMBER_OF_ATOMS_AA:       29     116    116  100.00     220    153   69.55                29
SUMMARY(RMSD_GDC):     2.226          2.182                  2.416           74.372   71.735

#CA            N1   N2   DIST      N    RMSD   Seq_Id      LGA_S     LGA_Q 
SUMMARY(LGA)   35   29    5.0     29    2.23    17.24     79.272     1.247

Unitary ROTATION matrix and the SHIFT vector superimpose molecules  (1=>2)
  X_new =  -0.420000 * X  +  -0.585654 * Y  +  -0.693260 * Z  +  12.659374
  Y_new =  -0.461115 * X  +   0.795663 * Y  +  -0.392803 * Z  + -31.630075
  Z_new =   0.781648 * X  +   0.154695 * Y  +  -0.604232 * Z  +   5.544326 

Euler angles from the ROTATION matrix. Conventions XYZ and ZXZ:
           Phi     Theta       Psi   [DEG:       Phi     Theta       Psi ]
XYZ: -2.309566 -0.897303  2.890956   [DEG: -132.3284  -51.4117  165.6396 ]
ZXZ: -1.055296  2.219598  1.375412   [DEG:  -60.4640  127.1736   78.8053 ]
 
# END of job
