
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
# PARAMETERS: -3  -o2  -gdc  test14  -al  
# FIXED Atom-Atom correspondence
# GDT and LCS analysis 

LCS - RMSD CUTOFF   5.00      length       segment         l_RMS    g_RMS
  LONGEST_CONTINUOUS_SEGMENT:    25     141_A - 176_A       4.51     7.29
  LCS_AVERAGE:     80.14

LCS - RMSD CUTOFF   2.00      length       segment         l_RMS    g_RMS
  LONGEST_CONTINUOUS_SEGMENT:    20     141_A - 160_A       1.93     8.64
  LCS_AVERAGE:     51.49

LCS - RMSD CUTOFF   1.00      length       segment         l_RMS    g_RMS
  LONGEST_CONTINUOUS_SEGMENT:    17     141_A - 157_A       0.77     8.74
  LCS_AVERAGE:     38.76

LCS_GDT    MOLECULE-1    MOLECULE-2     LCS_DETAILS     GDT_DETAILS                                                    TOTAL NUMBER OF RESIDUE PAIRS:   29
LCS_GDT     RESIDUE       RESIDUE       SEGMENT_SIZE    GLOBAL DISTANCE TEST COLUMNS: number of residues under the threshold assigned to each residue pair
LCS_GDT   NAME NUMBER   NAME NUMBER    1.0  2.0  5.0    0.5  1.0  1.5  2.0  2.5  3.0  3.5  4.0  4.5  5.0  5.5  6.0  6.5  7.0  7.5  8.0  8.5  9.0  9.5 10.0
LCS_GDT     G    52_A     N   141_A     17   20   25      5   14   17   17   18   18   18   19   20   21   22   22   23   23   24   24   25   25   25   26 
LCS_GDT     I    53_A     A   142_A     17   20   25      9   14   17   17   18   18   18   19   20   21   22   22   23   23   24   24   25   25   25   26 
LCS_GDT     I    54_A     Y   143_A     17   20   25      9   14   17   17   18   18   18   19   20   21   22   22   23   23   24   24   25   25   25   26 
LCS_GDT     P    55_A     M   144_A     17   20   25      9   14   17   17   18   18   18   19   20   21   22   22   23   23   24   24   25   25   25   26 
LCS_GDT     A    56_A     E   145_A     17   20   25      9   14   17   17   18   18   18   19   20   21   22   22   23   23   24   24   25   25   25   26 
LCS_GDT     L    57_A     L   146_A     17   20   25      9   14   17   17   18   18   18   19   20   21   22   22   23   23   24   24   25   25   25   26 
LCS_GDT     Q    58_A     G   147_A     17   20   25      9   14   17   17   18   18   18   19   20   21   22   22   23   23   24   24   25   25   25   26 
LCS_GDT     T    59_A     T   148_A     17   20   25      6   14   17   17   18   18   18   19   20   21   22   22   23   23   24   24   25   25   25   26 
LCS_GDT     K    60_A     N   149_A     17   20   25      3    3   17   17   18   18   18   19   20   21   22   22   23   23   24   24   25   25   25   26 
LCS_GDT     N    61_A     R   150_A     17   20   25      8   14   17   17   18   18   18   19   20   21   22   22   23   23   24   24   25   25   25   26 
LCS_GDT     V    62_A     A   151_A     17   20   25      9   14   17   17   18   18   18   19   20   21   22   22   23   23   24   24   25   25   25   26 
LCS_GDT     D    63_A     D   152_A     17   20   25      9   14   17   17   18   18   18   19   20   21   22   22   23   23   24   24   25   25   25   26 
LCS_GDT     L    64_A     A   153_A     17   20   25      9   14   17   17   18   18   18   19   20   21   22   22   23   23   24   24   25   25   25   26 
LCS_GDT     A    65_A     V   154_A     17   20   25      3   14   17   17   18   18   18   19   20   21   22   22   23   23   24   24   25   25   25   26 
LCS_GDT     L    66_A     L   155_A     17   20   25      9   14   17   17   18   18   18   19   20   21   22   22   23   23   24   24   25   25   25   26 
LCS_GDT     A    67_A     H   156_A     17   20   25      3   12   17   17   18   18   18   19   20   21   22   22   23   23   24   24   25   25   25   26 
LCS_GDT     G    68_A     D   157_A     17   20   25      3   14   17   17   18   18   18   19   20   21   22   22   23   23   24   24   25   25   25   26 
LCS_GDT     I    69_A     T   158_A      3   20   25      3    3    3    4    7    9   18   19   20   21   22   22   23   23   24   24   25   25   25   26 
LCS_GDT     T    70_A     P   159_A      3   20   25      3    3    3    5   18   18   18   19   20   21   22   22   23   23   24   24   25   25   25   26 
LCS_GDT     I    71_A     N   160_A      3   20   25      3    5    6   12   16   16   17   19   20   20   20   20   23   23   24   24   25   25   25   26 
LCS_GDT     T    72_A     Q   172_A      0    3   25      0    0    3    3    3    3    4    5    6    6    8    8   10   13   18   22   25   25   25   25 
LCS_GDT     D    73_A     F   173_A      3    3   25      1    1    3    3    4    5    5    7    8   12   16   20   21   23   23   24   25   25   25   26 
LCS_GDT     E    74_A     K   174_A      3    3   25      0    3    3    3    4    6    6    8   12   19   22   22   23   23   24   24   25   25   25   26 
LCS_GDT     R    75_A     A   175_A      3    4   25      0    3    3    4    4    6   11   17   20   21   22   22   23   23   24   24   25   25   25   26 
LCS_GDT     K    76_A     V   176_A      3    4   25      0    3    3    4    4    6   11   15   20   21   22   22   23   23   24   24   25   25   25   26 
LCS_GDT     K    77_A     G   177_A      4    4   13      1    4    4    4    4    5    6    7   10   10   11   17   19   20   24   24   25   25   25   26 
LCS_GDT     A    78_A     D   178_A      4    4   13      3    4    4    4    4    5    6    7    8    9   11   12   14   14   14   16   17   17   20   26 
LCS_GDT     I    79_A     S   179_A      4    4   13      3    4    4    4    4    5    6    7    8    9   11   11   12   13   14   16   17   17   20   23 
LCS_GDT     D    80_A     L   180_A      4    4   10      3    4    4    4    4    4    4    7    7    8   11   12   14   16   17   24   25   25   25   26 
LCS_AVERAGE  LCS_A:  56.80  (  38.76   51.49   80.14 )

GLOBAL_DISTANCE_TEST (summary information about detected largest sets of residues (represented by selected AToms) that can fit under specified thresholds)
GDT DIST_CUTOFF  0.50   1.00   1.50   2.00   2.50   3.00   3.50   4.00   4.50   5.00   5.50   6.00   6.50   7.00   7.50   8.00   8.50   9.00   9.50  10.00
GDT NUMBER_AT      9     14     17     17     18     18     18     19     20     21     22     22     23     23     24     24     25     25     25     26 
GDT PERCENT_AT  31.03  48.28  58.62  58.62  62.07  62.07  62.07  65.52  68.97  72.41  75.86  75.86  79.31  79.31  82.76  82.76  86.21  86.21  86.21  89.66
GDT RMS_LOCAL    0.29   0.57   0.77   0.77   1.10   1.10   1.10   1.57   1.93   3.09   3.33   3.33   3.67   3.67   4.40   4.00   4.51   4.51   4.51   5.22
GDT RMS_ALL_AT   8.65   8.65   8.74   8.74   8.58   8.58   8.58   8.73   8.64   7.12   7.02   7.02   6.98   6.98   6.70   7.08   7.29   7.29   7.29   6.63

# Checking swapping

#      Molecule1      Molecule2  DISTANCE    Mis    MC     All    Dist_max   GDC_mc  GDC_all
LGA    G    52_A      N   141_A     2.383     0    0.150   0.150     2.383   66.786   66.786
LGA    I    53_A      A   142_A     2.088     0    0.034   0.040     2.352   66.786   68.000
LGA    I    54_A      Y   143_A     2.550     0    0.063   0.068     3.154   62.857   60.286
LGA    P    55_A      M   144_A     2.710     0    0.116   0.151     2.942   59.048   58.667
LGA    A    56_A      E   145_A     1.541     0    0.049   0.049     1.925   77.143   76.286
LGA    L    57_A      L   146_A     1.475     0    0.115   0.605     3.347   75.119   70.119
LGA    Q    58_A      G   147_A     2.938     0    0.103   0.103     2.938   62.857   62.857
LGA    T    59_A      T   148_A     1.071     0    0.692   0.930     3.867   77.381   68.980
LGA    K    60_A      N   149_A     1.743     0    0.088   0.088     2.180   75.119   74.667
LGA    N    61_A      R   150_A     1.737     0    0.070   0.071     2.612   71.071   69.810
LGA    V    62_A      A   151_A     1.712     0    0.087   0.128     1.984   75.000   74.571
LGA    D    63_A      D   152_A     2.619     0    0.114   0.254     3.569   57.262   53.690
LGA    L    64_A      A   153_A     2.385     0    0.104   0.138     2.633   66.786   64.857
LGA    A    65_A      V   154_A     2.444     0    0.110   0.129     2.838   62.857   63.238
LGA    L    66_A      L   155_A     1.951     0    0.030   0.987     3.341   68.810   69.107
LGA    A    67_A      H   156_A     2.203     0    0.566   0.638     3.572   59.405   60.476
LGA    G    68_A      D   157_A     2.090     0    0.435   0.435     2.418   64.762   64.762
LGA    I    69_A      T   158_A     4.267     0    0.606   0.579     6.961   41.905   36.190
LGA    T    70_A      P   159_A     3.897     0    0.615   0.580     5.217   37.738   38.857
LGA    I    71_A      N   160_A     8.054     0    0.615   0.600    10.268    9.762    7.810
LGA    T    72_A      Q   172_A    12.481     0    0.731   0.663    13.057    0.000    0.000
LGA    D    73_A      F   173_A     9.088     0    0.617   0.596    10.655    9.524    7.619
LGA    E    74_A      K   174_A     7.122     0    0.540   0.558     8.959    9.405    8.095
LGA    R    75_A      A   175_A     5.917     0    0.669   0.624     6.437   19.286   21.714
LGA    K    76_A      V   176_A     7.309     0    0.662   0.607    10.301    5.952    6.190
LGA    K    77_A      G   177_A    14.234     0    0.584   0.584    14.781    0.000    0.000
LGA    A    78_A      D   178_A    16.157     0    0.609   0.610    17.297    0.000    0.000
LGA    I    79_A      S   179_A    15.480     0    0.166   0.225    15.890    0.000    0.000
LGA    D    80_A      L   180_A    14.598     0    0.070   0.065    14.598    0.000    0.000

# RMSD_GDC results:       CA      MC common percent     ALL common percent   GDC_mc  GDC_all
NUMBER_OF_ATOMS_AA:       29     116    116  100.00     220    152   69.09                29
SUMMARY(RMSD_GDC):     6.574          6.454                  6.485           44.228   43.229

#CA            N1   N2   DIST      N    RMSD    GDT_TS    LGA_S3     LGA_Q 
SUMMARY(GDT)   35   29    5.0     21    3.09    63.793    62.267     0.658

LGA_LOCAL      RMSD:   3.090  Number of atoms:   21  under DIST:   5.00
LGA_ASGN_ATOMS RMSD:   7.122  Number of assigned atoms:   29 
Std_ASGN_ATOMS RMSD:   6.574  Standard rmsd on all 29 assigned CA atoms 

Unitary ROTATION matrix and the SHIFT vector superimpose molecules  (1=>2)
  X_new =  -0.449372 * X  +  -0.794446 * Y  +  -0.408560 * Z  +  13.435723
  Y_new =  -0.599394 * X  +   0.607244 * Y  +  -0.521519 * Z  + -26.990023
  Z_new =   0.662414 * X  +   0.010532 * Y  +  -0.749064 * Z  +  10.332670 

Euler angles from the ROTATION matrix. Conventions XYZ and ZXZ:
           Phi     Theta       Psi   [DEG:       Phi     Theta       Psi ]
XYZ: -2.214112 -0.724037  3.127533   [DEG: -126.8593  -41.4843  179.1944 ]
ZXZ: -0.664539  2.417444  1.554898   [DEG:  -38.0753  138.5093   89.0891 ]
 
# END of job
