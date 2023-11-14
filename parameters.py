# config.py

workdir = "/eos/user/m/mbuonsan/Tau3mu_2023/Macro_BDT/AMS/"
cat_label = ["A", "B", "C"]
inputfile = "t3mminitree_xgb_2023_120523.root"

# Coefficienti per la normalizzazione del segnale
Dplus_correction_2017 = 1.05
Bs_correction_2017 = 1.12
f_correction_2017 = 1.
Ds_correction_2017 = 0.89
Lumi_data_2017 = [4.79, 9.63, 4.24, 9.3, 9.89]
wNormDs_2017 = 5.344E-04
wNormB0_2017 = 4.991E-04
wNormBp_2017 = 5.096E-04
wNormW_2017 = 1.611E-04

Dplus_correction_2018 = 1.05
Bs_correction_2018 = 1.12
f_correction_2018 = 1.
Ds_correction_2018 = 0.91
Ds_correction_dM4_2018 = 0.81
Lumi_data_2018 = [3.67, 13.98, 7.06, 6.90, 31.75]
wNormDs_2018 = 5.70E-04
wNormB0_2018 = 5.42E-04
wNormBp_2018 = 5.58E-04
wNormW_2018 = 3.05E-04

common_cut_2017 = "(fv_nC>0 && bs_sv_d2Dsig>=3.75 && !(TMath::IsNaN(fv_d3Dsig)))"
common_cut_2018 = "(fv_nC>0 && !(TMath::IsNaN(fv_d3Dsig)) )"

cutB_A = "(tripletMass<=1.753 && tripletMass>=1.62) || (tripletMass<=2.0 && tripletMass>=1.801)"
cutB_B = "(tripletMass<=1.739 && tripletMass>=1.62) || (tripletMass<=2.0 && tripletMass>=1.815)"
cutB_C = "(tripletMass<=1.727 && tripletMass>=1.62) || (tripletMass<=2.0 && tripletMass>=1.827)"

