#!/usr/bin/env python

# quick hack to add dyes without external programs

import sys,os,math,numpy,random

S488_raw = """\
ATOM     20  N   S488   2        4.337  -2.664   4.721  1.00 23.62      CSP 
ATOM     21  H   S488   2        4.065  -3.386   4.066  1.00  0.00      CSP 
ATOM     22  CA  S488   2        3.375  -2.396   5.778  1.00 23.56      CSP 
ATOM     23  HA  S488   2        3.846  -1.858   6.591  1.00  0.00      CSP 
ATOM     24  CB  S488   2        2.207  -1.536   5.225  1.00 22.69      CSP 
ATOM     25  HB1 S488   2        2.648  -0.644   4.734  1.00  0.00      CSP 
ATOM     26  HB2 S488   2        1.707  -2.135   4.434  1.00  0.00      CSP 
ATOM     27  SG  S488   2        0.995  -0.956   6.470  1.00  0.00      CSP 
ATOM     28  C   S488   2        2.878  -3.687   6.402  1.00 23.60      CSP 
ATOM     29  O   S488   2        1.982  -3.690   7.238  1.00 23.11      CSP 
ATOM     30  C1  S488   2       -3.007  -9.561  -6.486  1.00  0.00      CSP 
ATOM     31  C2  S488   2       -2.856 -10.965  -6.349  1.00  0.00      CSP 
ATOM     32  C3  S488   2       -2.815  -8.972  -7.764  1.00  0.00      CSP 
ATOM     33  C4  S488   2       -2.530 -11.768  -7.483  1.00  0.00      CSP 
ATOM     34  C5  S488   2       -2.442  -9.780  -8.879  1.00  0.00      CSP 
ATOM     35  O6  S488   2       -2.313 -11.114  -8.656  1.00  0.00      CSP 
ATOM     36  C7  S488   2       -3.032  -7.601  -7.967  1.00  0.00      CSP 
ATOM     37  C9  S488   2       -2.933  -7.054  -9.239  1.00  0.00      CSP 
ATOM     38  H23 S488   2       -3.341  -6.955  -7.154  1.00  0.00      CSP 
ATOM     39 HH23 S488   2       -3.152  -6.001  -9.350  1.00  0.00      CSP 
ATOM     40  C11 S488   2       -2.579  -7.843 -10.342  1.00  0.00      CSP 
ATOM     41  N15 S488   2       -2.550  -7.055 -11.472  1.00  0.00      CSP 
ATOM     42  H27 S488   2       -2.731  -6.082 -11.440  1.00  0.00      CSP 
ATOM     43 HH27 S488   2       -2.207  -7.524 -12.299  1.00  0.00      CSP 
ATOM     44  C13 S488   2       -2.264  -9.203 -10.173  1.00  0.00      CSP 
ATOM     45  S70 S488   2       -1.823  -9.999 -11.475  1.00  0.00      CSP 
ATOM     46  O72 S488   2       -3.019 -10.459 -12.155  1.00  0.00      CSP 
ATOM     47  O74 S488   2       -1.105  -9.098 -12.372  1.00  0.00      CSP 
ATOM     48  O76 S488   2       -0.927 -11.078 -11.139  1.00  0.00      CSP 
ATOM     49  C8  S488   2       -3.037 -11.592  -5.108  1.00  0.00      CSP 
ATOM     50  C10 S488   2       -2.875 -12.967  -4.985  1.00  0.00      CSP 
ATOM     51  H24 S488   2       -3.311 -11.023  -4.228  1.00  0.00      CSP 
ATOM     52 HH24 S488   2       -3.007 -13.403  -4.005  1.00  0.00      CSP 
ATOM     53  C12 S488   2       -2.544 -13.761  -6.093  1.00  0.00      CSP 
ATOM     54  N16 S488   2       -2.366 -15.064  -5.682  1.00  0.00      CSP 
ATOM     55  H28 S488   2       -2.541 -15.354  -4.751  1.00  0.00      CSP 
ATOM     56 HH28 S488   2       -2.190 -15.715  -6.436  1.00  0.00      CSP 
ATOM     57  C14 S488   2       -2.411 -13.187  -7.369  1.00  0.00      CSP 
ATOM     58  S71 S488   2       -1.970 -14.155  -8.549  1.00  0.00      CSP 
ATOM     59  O73 S488   2       -0.518 -14.190  -8.563  1.00  0.00      CSP 
ATOM     60  O75 S488   2       -2.455 -15.509  -8.289  1.00  0.00      CSP 
ATOM     61  O77 S488   2       -2.518 -13.724  -9.810  1.00  0.00      CSP 
ATOM     62  C45 S488   2       -3.529  -8.706  -5.317  1.00  0.00      CSP 
ATOM     63  C46 S488   2       -2.570  -8.071  -4.503  1.00  0.00      CSP 
ATOM     64  C47 S488   2       -2.952  -7.227  -3.456  1.00  0.00      CSP 
ATOM     65  C48 S488   2       -4.319  -6.988  -3.181  1.00  0.00      CSP 
ATOM     66  C49 S488   2       -5.264  -7.627  -3.999  1.00  0.00      CSP 
ATOM     67  H56 S488   2       -1.518  -8.228  -4.699  1.00  0.00      CSP 
ATOM     68  H57 S488   2       -2.171  -6.764  -2.871  1.00  0.00      CSP 
ATOM     69  H59 S488   2       -6.323  -7.466  -3.824  1.00  0.00      CSP 
ATOM     70  C50 S488   2       -4.903  -8.475  -5.054  1.00  0.00      CSP 
ATOM     71  C51 S488   2       -6.000  -9.096  -5.879  1.00  0.00      CSP 
ATOM     72  O52 S488   2       -5.750  -9.878  -6.834  1.00  0.00      CSP 
ATOM     73  O53 S488   2       -7.203  -8.824  -5.589  1.00  0.00      CSP 
ATOM     74  CI  S488   2       -4.750  -6.153  -2.143  1.00  0.00      CSP 
ATOM     75  OI  S488   2       -5.944  -5.958  -1.932  1.00  0.00      CSP 
ATOM     76  NJ  S488   2       -3.898  -5.522  -1.324  1.00  0.00      CSP 
ATOM     77  HJ  S488   2       -2.916  -5.618  -1.407  1.00  0.00      CSP 
ATOM     78  CK  S488   2       -4.419  -4.702  -0.275  1.00  0.00      CSP 
ATOM     79  HK1 S488   2       -5.067  -5.316   0.396  1.00  0.00      CSP 
ATOM     80  HK2 S488   2       -5.058  -3.895  -0.708  1.00  0.00      CSP 
ATOM     81  CL  S488   2       -3.305  -4.072   0.548  1.00  0.00      CSP 
ATOM     82  HL1 S488   2       -2.661  -3.456  -0.119  1.00  0.00      CSP 
ATOM     83  HL2 S488   2       -2.671  -4.878   0.984  1.00  0.00      CSP 
ATOM     84  CM  S488   2       -3.863  -3.198   1.673  1.00  0.00      CSP 
ATOM     85  HM1 S488   2       -4.508  -3.822   2.333  1.00  0.00      CSP 
ATOM     86  HM2 S488   2       -4.500  -2.400   1.229  1.00  0.00      CSP 
ATOM     87  CN  S488   2       -2.736  -2.568   2.494  1.00  0.00      CSP 
ATOM     88  HN1 S488   2       -2.097  -1.969   1.805  1.00  0.00      CSP 
ATOM     89  HN2 S488   2       -2.101  -3.378   2.920  1.00  0.00      CSP 
ATOM     90  CO  S488   2       -3.254  -1.665   3.624  1.00  0.00      CSP 
ATOM     91  HO1 S488   2       -3.934  -2.286   4.257  1.00  0.00      CSP 
ATOM     92  HO2 S488   2       -3.885  -0.865   3.169  1.00  0.00      CSP 
ATOM     93  NP  S488   2       -2.211  -1.055   4.481  1.00  0.00      CSP 
ATOM     94  CQ  S488   2       -2.611  -0.315   5.489  1.00  0.00      CSP 
ATOM     95  OQ  S488   2       -3.782  -0.009   5.706  1.00  0.00      CSP 
ATOM     96  CR  S488   2       -1.469   0.348   6.178  1.00  0.00      CSP 
ATOM     97  HR1 S488   2       -1.558   1.450   6.109  1.00  0.00      CSP 
ATOM     98  HR2 S488   2       -1.458   0.034   7.243  1.00  0.00      CSP 
ATOM     99  CS  S488   2       -0.248  -0.177   5.388  1.00  0.00      CSP 
ATOM    100  HS  S488   2        0.223   0.655   4.823  1.00  0.00      CSP 
ATOM    101  CF  S488   2       -0.889  -1.091   4.367  1.00  0.00      CSP 
ATOM    102  OF  S488   2       -0.265  -1.403   3.351  1.00  0.00      CSP """

S488_crds = S488_raw.split("\n")

S594_raw = """\
ATOM   1056  N   S594  66       31.644  -5.141  33.101  1.00  0.00      CSP 
ATOM   1057  H   S594  66       31.511  -4.465  33.834  1.00  0.00      CSP 
ATOM   1058  CA  S594  66       32.296  -6.297  33.673  1.00  0.00      CSP 
ATOM   1059  HA  S594  66       31.663  -7.162  33.539  1.00  0.00      CSP 
ATOM   1060  CB  S594  66       33.711  -6.503  33.061  1.00  0.00      CSP 
ATOM   1061  HB1 S594  66       33.593  -6.590  31.962  1.00  0.00      CSP 
ATOM   1062  HB2 S594  66       34.282  -5.571  33.258  1.00  0.00      CSP 
ATOM   1063  SG  S594  66       34.637  -7.959  33.649  1.00  0.00      CSP 
ATOM   1064  C   S594  66       32.386  -6.015  35.192  1.00  0.00      CSP 
ATOM   1065  C1  S594  66       40.813   6.796  35.511  1.00  0.00      CSP 
ATOM   1066  C2  S594  66       40.407   7.447  36.675  1.00  0.00      CSP 
ATOM   1067  C3  S594  66       41.064   7.583  34.386  1.00  0.00      CSP 
ATOM   1068  C4  S594  66       40.213   8.849  36.707  1.00  0.00      CSP 
ATOM   1069  C5  S594  66       40.877   8.985  34.408  1.00  0.00      CSP 
ATOM   1070  O6  S594  66       40.439   9.598  35.566  1.00  0.00      CSP 
ATOM   1071  C7  S594  66       41.522   6.966  33.240  1.00  0.00      CSP 
ATOM   1072  C9  S594  66       41.760   7.644  32.061  1.00  0.00      CSP 
ATOM   1073  C11 S594  66       41.635   9.056  32.074  1.00  0.00      CSP 
ATOM   1074  C13 S594  66       41.176   9.715  33.249  1.00  0.00      CSP 
ATOM   1075  H13 S594  66       41.034  10.779  33.318  1.00  0.00      CSP 
ATOM   1076  H23 S594  66       41.643   5.901  33.286  1.00  0.00      CSP 
ATOM   1077  N15 S594  66       41.926   9.773  30.934  1.00  0.00      CSP 
ATOM   1078  C23 S594  66       41.689  11.222  30.766  1.00  0.00      CSP 
ATOM   1079  H31 S594  66       40.617  11.379  30.533  1.00  0.00      CSP 
ATOM   1080  H32 S594  66       41.959  11.746  31.707  1.00  0.00      CSP 
ATOM   1081  H33 S594  66       42.328  11.622  29.955  1.00  0.00      CSP 
ATOM   1082  C8  S594  66       40.215   6.697  37.818  1.00  0.00      CSP 
ATOM   1083  C10 S594  66       39.747   7.238  39.001  1.00  0.00      CSP 
ATOM   1084  C12 S594  66       39.549   8.640  39.060  1.00  0.00      CSP 
ATOM   1085  C14 S594  66       39.805   9.439  37.911  1.00  0.00      CSP 
ATOM   1086  H14 S594  66       39.689  10.508  37.891  1.00  0.00      CSP 
ATOM   1087  H24 S594  66       40.396   5.642  37.732  1.00  0.00      CSP 
ATOM   1088  N16 S594  66       39.087   9.217  40.224  1.00  0.00      CSP 
ATOM   1089  C24 S594  66       38.964  10.672  40.456  1.00  0.00      CSP 
ATOM   1090  H34 S594  66       39.967  11.074  40.705  1.00  0.00      CSP 
ATOM   1091  H35 S594  66       38.575  11.156  39.537  1.00  0.00      CSP 
ATOM   1092  H36 S594  66       38.251  10.874  41.279  1.00  0.00      CSP 
ATOM   1093  C18 S594  66       38.839   8.387  41.449  1.00  0.00      CSP 
ATOM   1094  C20 S594  66       39.084   6.891  41.259  1.00  0.00      CSP 
ATOM   1095  C22 S594  66       39.443   6.333  40.038  1.00  0.00      CSP 
ATOM   1096  C26 S594  66       39.782   8.874  42.562  1.00  0.00      CSP 
ATOM   1097  C28 S594  66       37.375   8.547  41.890  1.00  0.00      CSP 
ATOM   1098  H48 S594  66       38.742   6.220  42.031  1.00  0.00      CSP 
ATOM   1099  H37 S594  66       39.358   9.714  43.142  1.00  0.00      CSP 
ATOM   1100  H38 S594  66       40.777   9.141  42.144  1.00  0.00      CSP 
ATOM   1101  H39 S594  66       39.964   8.042  43.280  1.00  0.00      CSP 
ATOM   1102  H40 S594  66       36.732   8.816  41.025  1.00  0.00      CSP 
ATOM   1103  H41 S594  66       37.266   9.290  42.701  1.00  0.00      CSP 
ATOM   1104  H42 S594  66       36.976   7.581  42.275  1.00  0.00      CSP 
ATOM   1105  C17 S594  66       42.309   9.076  29.662  1.00  0.00      CSP 
ATOM   1106  C19 S594  66       42.350   7.553  29.763  1.00  0.00      CSP 
ATOM   1107  C21 S594  66       42.053   6.859  30.927  1.00  0.00      CSP 
ATOM   1108  C25 S594  66       41.283   9.420  28.573  1.00  0.00      CSP 
ATOM   1109  C27 S594  66       43.711   9.549  29.242  1.00  0.00      CSP 
ATOM   1110  H47 S594  66       42.554   6.988  28.867  1.00  0.00      CSP 
ATOM   1111  H25 S594  66       41.515  10.356  28.036  1.00  0.00      CSP 
ATOM   1112  H26 S594  66       40.257   9.443  29.004  1.00  0.00      CSP 
ATOM   1113  H27 S594  66       41.263   8.604  27.813  1.00  0.00      CSP 
ATOM   1114  H28 S594  66       44.294   9.898  30.120  1.00  0.00      CSP 
ATOM   1115  H29 S594  66       43.657  10.348  28.478  1.00  0.00      CSP 
ATOM   1116  H30 S594  66       44.287   8.708  28.796  1.00  0.00      CSP 
ATOM   1117  C29 S594  66       42.086   5.310  30.918  1.00  0.00      CSP 
ATOM   1118  S31 S594  66       40.763   4.698  30.721  1.00  0.00      CSP 
ATOM   1119  O33 S594  66       40.421   4.748  29.302  1.00  0.00      CSP 
ATOM   1120  O35 S594  66       39.717   5.369  31.481  1.00  0.00      CSP 
ATOM   1121  O37 S594  66       40.821   3.296  31.132  1.00  0.00      CSP 
ATOM   1122  H45 S594  66       42.792   5.009  30.124  1.00  0.00      CSP 
ATOM   1123  H46 S594  66       42.577   4.963  31.851  1.00  0.00      CSP 
ATOM   1124  C30 S594  66       39.501   4.794  39.897  1.00  0.00      CSP 
ATOM   1125  S32 S594  66       38.291   4.211  39.299  1.00  0.00      CSP 
ATOM   1126  O34 S594  66       37.797   4.987  38.168  1.00  0.00      CSP 
ATOM   1127  O36 S594  66       37.243   4.138  40.315  1.00  0.00      CSP 
ATOM   1128  O38 S594  66       38.592   2.855  38.847  1.00  0.00      CSP 
ATOM   1129  H43 S594  66       39.684   4.386  40.908  1.00  0.00      CSP 
ATOM   1130  H44 S594  66       40.420   4.526  39.335  1.00  0.00      CSP 
ATOM   1131  C45 S594  66       41.092   5.272  35.491  1.00  0.00      CSP 
ATOM   1132  C46 S594  66       40.042   4.420  35.095  1.00  0.00      CSP 
ATOM   1133  C47 S594  66       40.234   3.041  34.980  1.00  0.00      CSP 
ATOM   1134  C48 S594  66       41.487   2.453  35.267  1.00  0.00      CSP 
ATOM   1135  C49 S594  66       42.528   3.324  35.662  1.00  0.00      CSP 
ATOM   1136  H56 S594  66       39.070   4.834  34.854  1.00  0.00      CSP 
ATOM   1137  H57 S594  66       39.387   2.451  34.659  1.00  0.00      CSP 
ATOM   1138  H59 S594  66       43.505   2.918  35.882  1.00  0.00      CSP 
ATOM   1139  C50 S594  66       42.360   4.715  35.777  1.00  0.00      CSP 
ATOM   1140  C51 S594  66       43.543   5.573  36.181  1.00  0.00      CSP 
ATOM   1141  O52 S594  66       43.472   6.837  36.300  1.00  0.00      CSP 
ATOM   1142  O53 S594  66       44.661   5.011  36.412  1.00  0.00      CSP 
ATOM   1143  CI  S594  66       41.712   1.064  35.171  1.00  0.00      CSP 
ATOM   1144  OI  S594  66       42.808   0.561  35.431  1.00  0.00      CSP 
ATOM   1145  NJ  S594  66       40.755   0.186  34.798  1.00  0.00      CSP 
ATOM   1146  HJ  S594  66       39.849   0.519  34.565  1.00  0.00      CSP 
ATOM   1147  CK  S594  66       41.057  -1.219  34.741  1.00  0.00      CSP 
ATOM   1148  HK1 S594  66       41.405  -1.560  35.745  1.00  0.00      CSP 
ATOM   1149  HK2 S594  66       41.896  -1.384  34.024  1.00  0.00      CSP 
ATOM   1150  CL  S594  66       39.868  -2.084  34.316  1.00  0.00      CSP 
ATOM   1151  HL1 S594  66       39.523  -1.750  33.312  1.00  0.00      CSP 
ATOM   1152  HL2 S594  66       39.033  -1.924  35.034  1.00  0.00      CSP 
ATOM   1153  CM  S594  66       40.218  -3.582  34.267  1.00  0.00      CSP 
ATOM   1154  HM1 S594  66       40.563  -3.905  35.274  1.00  0.00      CSP 
ATOM   1155  HM2 S594  66       41.057  -3.732  33.551  1.00  0.00      CSP 
ATOM   1156  CN  S594  66       39.015  -4.434  33.837  1.00  0.00      CSP 
ATOM   1157  HN1 S594  66       38.675  -4.073  32.839  1.00  0.00      CSP 
ATOM   1158  HN2 S594  66       38.179  -4.265  34.553  1.00  0.00      CSP 
ATOM   1159  CO  S594  66       39.322  -5.941  33.751  1.00  0.00      CSP 
ATOM   1160  HO1 S594  66       39.718  -6.251  34.750  1.00  0.00      CSP 
ATOM   1161  HO2 S594  66       40.148  -6.091  33.018  1.00  0.00      CSP 
ATOM   1162  NP  S594  66       38.173  -6.813  33.398  1.00  0.00      CSP 
ATOM   1163  CQ  S594  66       38.369  -8.111  33.381  1.00  0.00      CSP 
ATOM   1164  OQ  S594  66       39.468  -8.649  33.521  1.00  0.00      CSP 
ATOM   1165  CR  S594  66       37.182  -8.846  32.857  1.00  0.00      CSP 
ATOM   1166  HR1 S594  66       37.424  -9.350  31.901  1.00  0.00      CSP 
ATOM   1167  HR2 S594  66       36.848  -9.593  33.606  1.00  0.00      CSP 
ATOM   1168  CS  S594  66       36.150  -7.710  32.657  1.00  0.00      CSP 
ATOM   1169  HS  S594  66       35.884  -7.624  31.581  1.00  0.00      CSP 
ATOM   1170  CF  S594  66       36.945  -6.476  33.023  1.00  0.00      CSP 
ATOM   1171  OF  S594  66       36.625  -5.378  32.564  1.00  0.00      CSP 
ATOM   1172  O   S594  66       32.051  -4.860  35.585  1.00  0.00      CSP """

S594_crds = S594_raw.split("\n")

Usage = """

	build_dyes.py residue_number residue_type pdb_file

	e.g. ./build_dyes.py 33 S488 myprotein.pdb
	will build Alexa488 onto residue 33 of the supplied myprotein.pdb.
	residue 33 should be a cysteine for everything to work correctly.

"""

if len(sys.argv) != 5:
	print usage

resnum = int(sys.argv[1])
dye = sys.argv[2]
pdbinp = sys.argv[3]
pdbout = sys.argv[4]

sys.stdout.write("Adding chromophore %s to cysteine residue %i in pdb %s\n"%(dye,resnum,pdbinp))
sys.stdout.write("...and writing result to pdb %s\n"%(pdbout))

pdblines = filter(lambda x: x[0:6] == "ATOM  ", open(pdbinp).readlines())

def residue(line):
	return int(line[23:26])

def makevec(line):
	return numpy.array(line[30:54].split(),numpy.float64)

firstres = residue(pdblines[0])
lastres = residue(pdblines[-1])

pre = filter(lambda x: residue(x)<resnum, pdblines)
cys = filter(lambda x: residue(x) == resnum, pdblines)
post = filter(lambda x: residue(x)>resnum, pdblines)

CB_cys = filter(lambda x: x[13:15] == "CB", cys)[0]
SG_cys = filter(lambda x: x[13:15] == "SG", cys)[0]
HG_cys = filter(lambda x: x[13:15] == "HG", cys)[0]

if dye=="488":
	dyelines = S488_crds
elif dye== "594":
	dyelines = S594_crds
SG_dye = filter(lambda x: x[13:15] == "SG", dyelines )[0]
CS_dye = filter(lambda x: x[13:15] == "CS", dyelines )[0]

if resnum == firstres:
	dye = "N"+dye
elif resnum == lastres:
	dye = "C"+dye
else:
	dye = "S"+dye

SG_dye_xyz = makevec(SG_dye)
CS_dye_xyz = makevec(CS_dye)
CB_cys_xyz = makevec(CB_cys)
SG_cys_xyz = makevec(SG_cys)
HG_cys_xyz = makevec(HG_cys)

dx = SG_cys_xyz - SG_dye_xyz # translation vector to apply to dye

sb_vec = CB_cys_xyz - SG_cys_xyz
sh_vec = HG_cys_xyz - SG_cys_xyz
cs_vec = CS_dye_xyz - SG_dye_xyz

x = cs_vec/numpy.linalg.norm(cs_vec)
xnew = sh_vec/numpy.linalg.norm(sh_vec)
print cs_vec
print  x
print xnew
zc = numpy.cross(x,xnew)
z = zc/numpy.linalg.norm(zc)
y = numpy.cross(z,x)

M = numpy.array( [x,y,z], numpy.float64 )
M_T = numpy.transpose(M)

costheta = numpy.dot(x,xnew)
sintheta = numpy.linalg.norm(zc)
theta  = math.acos(costheta)
rot1  = numpy.array( [ [ costheta, -sintheta, 0. ], \
			[ sintheta, costheta, 0. ], \
			[ 0., 0., 1. ] ], numpy.float64 )

rot2 = numpy.dot( M_T, rot1 )
rot = numpy.dot( rot2, M )

print rot

xnewtest = numpy.dot(rot,x)
print x
print xnew
print xnewtest


# first we check for clashes

precrd = map(makevec,pre)
postcrd = map(makevec,post)
protcrd = precrd+postcrd
#dyecrd_raw = map(lambda x: x-SG_dye_xyz, map(makevec,dyelines))
#dyecrd = map(lambda x: SG_cys_xyz + numpy.dot(rot,x), dyecrd_raw)

def isclash(protxyz,dyexyz):
	rmin = 3.0
	for k in dyexyz:
		for l in protxyz:
			if numpy.linalg.norm( l-k ) < rmin:
				print l
				print k
				return 1
	return 0

#==================================================
outp = open(pdbout,"w")

atom_idx = 0
res_idx = 0
pres = -1
for line in pre:
	tmp = int(line[21:26])
	if tmp>pres:
		pres = tmp
		res_idx += 1
	atom_idx+=1
	outp.write("ATOM  %5i%10s%5i%s\n"%(atom_idx,line[11:21],res_idx,line[26:66]))

cys_atoms = []


res_idx += 1
for line in cys:
	atom = line[12:16].strip()
	if atom == "HG":
		continue
	cys_atoms.append(atom)
	atom_idx+=1
	outp.write("ATOM  %5i%6s%4s%5i%s\n"%(atom_idx,line[11:17],dye,res_idx,line[26:66]))

dyecrd_nocys  = []
dyelines_nocys = []
for line in dyelines:
	atom = line[12:16].strip()
	if atom in cys_atoms:
		continue
	#if atom in [ "O

	dyelines_nocys.append(line)
	xyz_old = makevec(line)
	dr = xyz_old - SG_dye_xyz 
	xyz_rottrans = SG_cys_xyz + numpy.dot(rot,dr)
	dyecrd_nocys.append(xyz_rottrans)

	#x,y,z = xyz_rottrans
	#outp.write("ATOM  %5i%10s%5i    %8.3f%8.3f%8.3f%6.2f%6.2f\n"%(atom_idx,line[11:21],res_idx,
	#	x,y,z,1.,0.))

ntry = 0
dyecrd_nocys_save = []
for xyz in dyecrd_nocys:
	dyecrd_nocys_save.append(xyz)
n_nocys = len(dyelines_nocys)

maxtry = 20
while ntry<maxtry and isclash(protcrd,dyecrd_nocys):
	print "Trying to remove clashes by rotating chromophore"
	ntry += 1
	print "Try number %i of %i"%(ntry,maxtry)
	zz = xnew
	crossp = numpy.cross(sb_vec,zz)
	xx = crossp / numpy.linalg.norm(crossp)
	yy = numpy.cross(zz,xx)
	print numpy.linalg.norm(xx)
	print numpy.linalg.norm(yy)
	print numpy.linalg.norm(zz)
	rand_theta = 2.* math.pi * random.random()
	costheta = math.cos(rand_theta)
	sintheta = math.sin(rand_theta)
	MM = numpy.array( [xx,yy,zz], numpy.float64 )
	MM_T = numpy.transpose(MM)
	rot1  = numpy.array( [ [ costheta, -sintheta, 0. ], \
				[ sintheta, costheta, 0. ], \
				[ 0., 0., 1. ] ], numpy.float64 )
	
	rot2 = numpy.dot( MM_T, rot1 )
	rot = numpy.dot( rot2, MM )

	dyecrd_nocys = []

	for xyz_old in dyecrd_nocys_save:
		dr = xyz_old - SG_cys_xyz 
		xyz_rottrans = SG_cys_xyz + numpy.dot(rot,dr)
		dyecrd_nocys.append(xyz_rottrans)
	

for k in range(n_nocys):
	line = dyelines_nocys[k]
	x,y,z = dyecrd_nocys[k]
	atom_idx += 1
	outp.write("ATOM  %5i%6s%4s%5i    %8.3f%8.3f%8.3f%6.2f%6.2f\n"%(atom_idx,line[11:17],dye,res_idx,
		x,y,z,1.,0.))

pres = -1
for line in post:
	#outp.write(line)
	tmp = int(line[21:26])
	if tmp>pres:
		pres = tmp
		res_idx += 1
	atom_idx+=1
	outp.write("ATOM  %5i%10s%5i%s\n"%(atom_idx,line[11:21],res_idx,line[26:66]))


outp.close()


