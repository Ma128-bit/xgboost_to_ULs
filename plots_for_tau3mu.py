import sys, os, subprocess, json
from datetime import datetime
import numpy as np
import pandas as pd
import ROOT
import argparse
from ROOT import *

def bdt_KS_plot(config, fold_index, categories, year):
    
    with open(config, 'r') as file:
        json_file = json.load(file)
    kfold = json_file['number_of_splits']
    output_path = json_file['output_path']
    date = json_file['date']
    inputfile = json_file['Name']
    out_tree_name = json_file['out_tree_name']
    index_branch = json_file['index_branch']
    Y_column = json_file['Y_column']
    pos_dir_xgboost = json_file['data_path']
    weight = json_file['weight_column']
    
    kfold_s = str(kfold)
    file_name = pos_dir_xgboost + output_path +"/" + date+ "/" + inputfile +"_minitree.root"
    f = ROOT.TFile(file_name, "READ")
    t = f.Get(out_tree_name)
    print("Opened input file:", file_name)

    outputdir = pos_dir_xgboost + output_path +"/" + date+ "/" + "KS_plots"
    
    if not os.path.exists(outputdir):
        subprocess.call("mkdir -p %s" % outputdir, shell=True)

    binning = "(25, 0.0, 1.0)"
    varname = "fold_" + str(fold_index) +"_"

    train_sel = "int("+index_branch+")%" + kfold_s + "!=0"
    test_sel = "int("+index_branch+")%" + kfold_s + "==0"

    sig_sel = Y_column+"!=0"
    bkg_sel = Y_column+"==0"

    signal_label_all = "MC #tau#rightarrow3#mu"
    if categories is not None:
        cat_label = categories.split(' ')
    else:
        return -1
    
    for i in range(len(cat_label)):
        h_BDT_Test_S = ROOT.TH1F()
        h_BDT_Train_S = ROOT.TH1F()
        h_BDT_Test_B = ROOT.TH1F()
        h_BDT_Train_B = ROOT.TH1F()

        h2_BDT_Test_S = ROOT.TH1F()
        h2_BDT_Train_S = ROOT.TH1F()
        h2_BDT_Test_B = ROOT.TH1F()
        h2_BDT_Train_B = ROOT.TH1F()

        phiveto = ""
        if i == 0:
            phiveto = "(!(dimu_OS1<1.044 && dimu_OS1>0.994) && !(dimu_OS2<1.044 && dimu_OS2>0.994))"
        elif i == 1:
            phiveto = "(!(dimu_OS1<1.053 && dimu_OS1>0.985) && !(dimu_OS2<1.053 && dimu_OS2>0.985))"
        elif i == 2:
            phiveto = "(!(dimu_OS1<1.064 && dimu_OS1>0.974) && !(dimu_OS2<1.064 && dimu_OS2>0.974))"

        omega_veto = '(!(dimu_OS1<0.79 && dimu_OS1>0.77) && !(dimu_OS2<0.79 & dimu_OS2>0.77))'
        
        if i == 0:
            phiveto = "(!(dimu_OS1<1.044 && dimu_OS1>0.994) && !(dimu_OS2<1.044 && dimu_OS2>0.994))"
        elif i == 1:
            phiveto = "(!(dimu_OS1<1.053 && dimu_OS1>0.985) && !(dimu_OS2<1.053 && dimu_OS2>0.985))"
        elif i == 2:
            phiveto = "(!(dimu_OS1<1.064 && dimu_OS1>0.974) && !(dimu_OS2<1.064 && dimu_OS2>0.974))"

        isSB = ""
        if i == 0:
            isSB = "((tripletMass<=1.753 && tripletMass>=1.62) || (tripletMass<=2.0 && tripletMass>=1.801))"
        elif i == 1:
            isSB = "((tripletMass<=1.739 && tripletMass>=1.62) || (tripletMass<=2.0 && tripletMass>=1.815))"
        elif i == 2:
            isSB = "((tripletMass<=1.727 && tripletMass>=1.62) || (tripletMass<=2.0 && tripletMass>=1.827))"

        categ = "category==" + str(i)

        t.Draw(varname + ">>hTrain_bkg" + binning, bkg_sel + "&&" + isSB + "&&" + train_sel + "&&" + phiveto + "&&" + categ + "&&" + omega_veto)
        h_BDT_Train_B = ROOT.gDirectory.Get("hTrain_bkg")
        t.Draw(varname + ">>hTest_bkg" + binning, bkg_sel + "&&" + isSB + "&&" + test_sel + "&&" + phiveto + "&&" + categ+ "&&" + omega_veto)
        h_BDT_Test_B = ROOT.gDirectory.Get("hTest_bkg")
        t.Draw(varname + ">>hTrain_sig" + binning, weight+"*(" + sig_sel + "&&" + train_sel + "&&" + phiveto + "&&" + categ + "&&" + omega_veto + ")")
        h_BDT_Train_S = ROOT.gDirectory.Get("hTrain_sig")
        t.Draw(varname + ">>hTest_sig" + binning, weight+"*(" + sig_sel + "&&" + test_sel + "&&" + phiveto + "&&" + categ + "&&" + omega_veto + ")")
        h_BDT_Test_S = ROOT.gDirectory.Get("hTest_sig")

        t.Draw(varname + ">>hTrain_bkg2" + binning, varname + ">0.6&&" + bkg_sel + "&&" + isSB + "&&" + train_sel + "&&" + phiveto + "&&" + categ + "&&" + omega_veto)
        h2_BDT_Train_B = ROOT.gDirectory.Get("hTrain_bkg2")
        t.Draw(varname + ">>hTest_bkg2" + binning, varname + ">0.6&&" + bkg_sel + "&&" + isSB + "&&" + test_sel + "&&" + phiveto + "&&" + categ + "&&" + omega_veto)
        h2_BDT_Test_B = ROOT.gDirectory.Get("hTest_bkg2")
        t.Draw(varname + ">>hTrain_sig2" + binning, weight+"*(" + varname + ">0.6&&" + sig_sel + "&&" + train_sel + "&&" + phiveto + "&&" + categ + "&&" + omega_veto + ")")
        h2_BDT_Train_S = ROOT.gDirectory.Get("hTrain_sig2")
        t.Draw(varname + ">>hTest_sig2" + binning, weight+"*(" + varname + ">0.6&&" + sig_sel + "&&" + test_sel + "&&" + phiveto + "&&" + categ + "&&" + omega_veto + ")")
        h2_BDT_Test_S = ROOT.gDirectory.Get("hTest_sig2")

        print(categ)
        print("Background train (test) ", h_BDT_Train_B.GetEntries(), "(", h_BDT_Test_B.GetEntries(), ")")
        print("Signal train (test) ", h_BDT_Train_S.GetEntries(), "(", h_BDT_Test_S.GetEntries(), ")")

        c3 = ROOT.TCanvas("c3", year + " " + cat_label[i], 150, 10, 800, 800)
        ROOT.gStyle.SetOptStat(0)
        leg = ROOT.TLegend(0.30, 0.70, 0.74, 0.82)

        h_BDT_Train_S.SetLineColor(ROOT.kBlue)
        h_BDT_Train_S.SetLineWidth(2)
        h_BDT_Train_S.Scale(1/h_BDT_Train_S.Integral())
        h_BDT_Train_S.GetXaxis().SetTitle(varname)
        h_BDT_Train_S.SetTitle("")
        h_BDT_Test_S.SetLineColor(ROOT.kBlue)
        h_BDT_Test_S.SetLineWidth(2)
        h_BDT_Test_S.Scale(1/h_BDT_Test_S.Integral())

        h_BDT_Train_B.SetLineColor(ROOT.kRed)
        h_BDT_Train_B.SetLineWidth(2)
        h_BDT_Train_B.Scale(1/h_BDT_Train_B.Integral())
        h_BDT_Train_B.GetXaxis().SetTitle(varname)
        h_BDT_Train_B.SetTitle("")
        h_BDT_Test_B.SetLineColor(ROOT.kRed)
        h_BDT_Test_B.SetLineWidth(2)
        h_BDT_Test_B.Scale(1/h_BDT_Test_B.Integral())

        h_BDT_Train_B.Draw("hist")
        h_BDT_Train_B.GetYaxis().SetRangeUser(10e-5, 10e-1)
        h_BDT_Test_B.Draw("lep same")
        h_BDT_Train_S.Draw("hist same")
        h_BDT_Test_S.Draw("lep same")

        leg.AddEntry(h_BDT_Train_S, signal_label_all + " - train", "f")
        leg.AddEntry(h_BDT_Test_S, signal_label_all + " - test", "lep")
        leg.AddEntry(h_BDT_Train_B, year + " data SB - train", "f")
        leg.AddEntry(h_BDT_Test_B, year + " data SB - test", "lep")
        leg.Draw()
        c3.Update()
        c3.SetLogy()

        print("Signal KS test:", h_BDT_Test_S.KolmogorovTest(h_BDT_Train_S))
        KS_value_S = h_BDT_Test_S.KolmogorovTest(h_BDT_Train_S)
        print("BKG KS test:", h_BDT_Test_B.KolmogorovTest(h_BDT_Train_B))
        KS_value_B = h_BDT_Test_B.KolmogorovTest(h_BDT_Train_B)
        KS = "[bdt>0.0] KS test signal(bkg): " + f"{KS_value_S:.2f}" + "(" + f"{KS_value_B:.2f}" + ")"
        text_KS = ROOT.TLatex(0.12, 0.91, KS)
        text_KS.SetTextSize(0.035)
        text_KS.SetNDC(True)
        text_KS.Draw("same")

        print("Signal KS test:", h2_BDT_Test_S.KolmogorovTest(h2_BDT_Train_S))
        KS2_value_S = h2_BDT_Test_S.KolmogorovTest(h2_BDT_Train_S)
        print("BKG KS test:", h2_BDT_Test_B.KolmogorovTest(h2_BDT_Train_B))
        KS2_value_B = h2_BDT_Test_B.KolmogorovTest(h2_BDT_Train_B)
        KS2 = "[bdt>0.6] KS test signal(bkg): " + f"{KS2_value_S:.2f}" + "(" + f"{KS2_value_B:.2f}" + ")"
        text_KS2 = ROOT.TLatex(0.12, 0.86, KS2)
        text_KS2.SetTextSize(0.035)
        text_KS2.SetNDC(True)
        text_KS2.Draw("same")

        c3.SaveAs(outputdir +"/" + cat_label[i] + "_KS_plot.png")

def to_string_with_precision(a_value, n=6):
    return f"{a_value:.{n}f}"


def bdt_taumass_plot(config, categories, year):
    with open(config, 'r') as file:
        json_file = json.load(file)
    kfold = json_file['number_of_splits']
    output_path = json_file['output_path']
    date = json_file['date']
    inputfile = json_file['Name']
    out_tree_name = json_file['out_tree_name']
    index_branch = json_file['index_branch']
    Y_column = json_file['Y_column']
    pos_dir_xgboost = json_file['data_path']
    weight = json_file['weight_column']
    
    kfold_s = str(kfold)
    file_name = pos_dir_xgboost + output_path +"/" + date+ "/" + inputfile +"_minitree.root"
    f = ROOT.TFile(file_name, "READ")
    t = f.Get(out_tree_name)
    print("Opened input file:", file_name)


    outputdir = pos_dir_xgboost + output_path +"/" + date+ "/" + "Mass_plots"
    
    if not os.path.exists(outputdir):
        subprocess.call("mkdir -p %s" % outputdir, shell=True)

    binning = "(40, 0.0, 1.0)"
    varname = "bdt_cv"

    sig_sel = Y_column+"!=0"
    bkg_sel = Y_column+"==0"

    signal_label_all = "MC #tau#rightarrow3#mu"
    if categories is not None:
        cat_label = categories.split(' ')
    else:
        return -1

    ncut = 5
    proba = np.zeros(ncut, dtype=float)
    edges = np.zeros(ncut, dtype=float)
    
    h_mc = ROOT.TH1F("h_mc", "", 200, 0, 0.2)
    t.Draw("abs(tripletMass-1.777)>>h_mc", sig_sel)
    h_mc = ROOT.gDirectory.Get("h_mc")
    for k in range(5):
        proba[k] = k / 5.0
    h_mc.GetQuantiles(ncut, edges, proba)

    bdt_range_signal = [f"abs(tripletMass-1.777)>{to_string_with_precision(edges[i], 3)} && abs(tripletMass-1.777)<={to_string_with_precision(edges[i+1], 3)}" for i in range(ncut-1)]

    h_data = ROOT.TH1F("h_data", "", 200, 0, 0.2)
    t.Draw("abs(tripletMass-1.777)>>h_data", bkg_sel)
    h_data = ROOT.gDirectory.Get("h_data")
    h_data.GetQuantiles(ncut, edges, proba)

    bdt_range_bkg = [f"abs(tripletMass-1.777)>{to_string_with_precision(edges[i], 3)} && abs(tripletMass-1.777)<={to_string_with_precision(edges[i+1], 3)}" for i in range(ncut-1)]

    bdt_range_bkg.append("tripletMass>0")
    bdt_range_signal.append("tripletMass>0")

    n_signal = len(bdt_range_signal)
    n_bkg = len(bdt_range_bkg)

    #for i in range(3):
    for i in range(1):
        hTrain_signal = []
        hTrain_bkg = []

        phiveto = ""
        if i == 0:
            phiveto = "(!(dimu_OS1<1.044 && dimu_OS1>0.994) && !(dimu_OS2<1.044 && dimu_OS2>0.994))"
        if i == 1:
            phiveto = "(!(dimu_OS1<1.053 && dimu_OS1>0.985) && !(dimu_OS2<1.053 && dimu_OS2>0.985))"
        if i == 2:
            phiveto = "(!(dimu_OS1<1.064 && dimu_OS1>0.974) && !(dimu_OS2<1.064 && dimu_OS2>0.974))"

        omega_veto = "(!(dimu_OS1<0.79 && dimu_OS1>0.77) && !(dimu_OS2<0.79 & dimu_OS2>0.77))"

        isSB = ""
        if i == 0:
            isSB = "((tripletMass<=1.753 && tripletMass>=1.62) || (tripletMass<=2.0 && tripletMass>=1.801))"
        if i == 1:
            isSB = "((tripletMass<=1.739 && tripletMass>=1.62) || (tripletMass<=2.0 && tripletMass>=1.815))"
        if i == 2:
            isSB = "((tripletMass<=1.727 && tripletMass>=1.62) || (tripletMass<=2.0 && tripletMass>=1.827))"

        categ = f"category=={i}"
        category = cat_label[i]

        hs_signal = ROOT.THStack("hs_signal", year + " " + signal_label_all + " " + category)
        x1, x2, y1, y2 = 0.20, 0.80, 0.15, 0.40  # bottom center
        if year == "2018" and i == 2:
            y1, y2 = 0.65, 0.90  # top center
        leg_signal = ROOT.TLegend(x1, y1, x2, y2)

        for j in range(n_signal):
            range_ = bdt_range_signal[j]
            s = str(j)
            print("bdt range", range_)
            t.Draw(f"{varname}>>hTrain_signal{s}{binning}", f"{weight}*({sig_sel}&&{bdt_range_signal[j]}&&{phiveto}&&{categ}&&{omega_veto})")
            hTrain_signal.append(ROOT.gDirectory.Get(f"hTrain_signal{s}"))
            hTrain_signal[j].SetTitle(year + " Signal - cat " + category)
            hTrain_signal[j].SetLineColor(1 + j)
            hTrain_signal[j].SetLineWidth(2)
            hTrain_signal[j].Scale(1 / hTrain_signal[j].GetEntries())
            hs_signal.Add(hTrain_signal[j])
            if range_ == "tripletMass>0":
                range_ = "full"
            range_ = range_.replace("tripletMass", "m(3#mu)")
            range_ = range_.replace("-1.777", "-m(#tau)")
            leg_signal.AddEntry(hTrain_signal[j], range_, "l")


        c1 = ROOT.TCanvas("c1", year + " " + signal_label_all + " " + category, 150, 10, 800, 800)
        hs_signal.Draw("hist nostack")
        
        t1 = ROOT.TPaveText(0.1, 0.9, 0.3, 1.0, year + " Signal - cat " + category)
        t1.Draw()
        leg_signal.Draw()
        hs_signal.GetXaxis().SetTitle("BDT score")

        ROOT.gPad.Modified()
        
        c1.SetLogy()
        #c1.Update()
        c1.SaveAs(outputdir + "/" + category + "_" + varname + "correlation_signal.png")

        c1ratio = ROOT.TCanvas("c1ratio", year + " " + signal_label_all + " " + category, 150, 10, 800, 300)
        ROOT.gStyle.SetOptTitle(0)
        ROOT.gStyle.SetOptStat(0)
        hratio = []
        for j in range(n_signal - 1):
            hratio.append(hTrain_signal[j].Clone("hratio"))
            hratio[j].SetLineColor(1 + j)
            hratio[j].Divide(hTrain_signal[n_signal - 1])
            if j == 0:
                hratio[j].GetXaxis().SetTitle("BDT score")
                hratio[j].GetYaxis().SetRangeUser(0.2, 1.8)
                hratio[j].Draw("lep")
            else:
                hratio[j].Draw("lepSame")

        c1ratio.SaveAs(outputdir + "/" + category + "_" + varname + "correlation_signal_ratio.png")


        c2 = ROOT.TCanvas("c2", year + " data SB " + category, 150, 10, 800, 800)
        hs_bkg = ROOT.THStack("hs_bkg", year + " data SB " + category)
        leg_bkg = ROOT.TLegend(0.20, 0.65, 0.80, 0.9)

        for j in range(n_bkg):
            range_ = bdt_range_bkg[j]
            s = str(j)
            print("bdt range", range_)
            t.Draw(f"{varname}>>hTrain_bkg{s}{binning}", f"{bkg_sel}&&{isSB}&&{bdt_range_bkg[j]}&&{phiveto}&&{categ}&&{omega_veto}")
            hTrain_bkg.append(ROOT.gDirectory.Get(f"hTrain_bkg{s}"))
            hTrain_bkg[j].SetTitle(year + " Data SB - cat " + category)
            hTrain_bkg[j].SetLineColor(1 + j)
            hTrain_bkg[j].SetLineWidth(2)
            print("Entries ", hTrain_bkg[j].GetEntries())
            #hTrain_bkg[j].Scale(1 / hTrain_bkg[j].GetEntries())
            hs_bkg.Add(hTrain_bkg[j])
            if range_ == "tripletMass>0":
                range_ = "full"
            range_ = range_.replace("tripletMass", "m(3#mu)")
            range_ = range_.replace("-1.777", "-m(#tau)")
            leg_bkg.AddEntry(hTrain_bkg[j], range_, "l")

        hs_bkg.Draw("hist nostack")
        t2 = ROOT.TPaveText(0.1, 0.9, 0.3, 1.0, year + " Data SB - cat " + category)
        t2.Draw()
        leg_bkg.Draw()
        hs_bkg.GetXaxis().SetTitle("BDT score")

        ROOT.gPad.Modified()
        c2.SetLogy()
        c2.Update()
        c2.SaveAs(outputdir + "/" + category + "_" + varname + "correlation_bkg.png")

        c2ratio = ROOT.TCanvas("c2ratio", year + " " + signal_label_all + " " + category, 150, 10, 800, 300)
        ROOT.gStyle.SetOptTitle(0)
        ROOT.gStyle.SetOptStat(0)
        for j in range(n_bkg - 1):
            hratio = hTrain_bkg[j].Clone("hratio")
            hratio.SetLineColor(1 + j)
            hratio.Divide(hTrain_bkg[n_bkg - 1])
            if (j == 0 and i == 0) or (j == 1 and i > 0):
                hratio.GetXaxis().SetTitle("BDT score")
                hratio.GetYaxis().SetRangeUser(0.2, 1.8)
                hratio.Draw("lep")
            else:
                hratio.Draw("same lep")

        c2ratio.Update()
        c2ratio.SaveAs(outputdir + "/" + category + "_" + varname + "correlation_bkg_ratio.png")
        """

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--config", type=str, help="Path to the copy of the JSON configuration file")
    parser.add_argument("--categories", type=str, help="List of categories separated by ' ' ")
    parser.add_argument("--index", type=int, help="Fold index")

    args = parser.parse_args()
    
    config = args.config
    categories = args.categories
    fold_index = args.index
    
    #bdt_KS_plot(config, fold_index, categories, "")

    bdt_taumass_plot(config, categories, "")

