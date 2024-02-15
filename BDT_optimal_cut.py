import ROOT
import sys, os, subprocess, json
import math
import ctypes
from parameters import *
from ROOT import *
from array import array
import argparse

class BDTcut3d:
    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c

def TH1_integral(h, xmin, xmax):
    axis = h.GetXaxis()
    bmin = axis.FindBin(xmin)
    bmax = axis.FindBin(xmax)
    integral = h.Integral(bmin, bmax)
    integral -= h.GetBinContent(bmin) * (xmin - axis.GetBinLowEdge(bmin)) / axis.GetBinWidth(bmin)
    integral -= h.GetBinContent(bmax) * (axis.GetBinUpEdge(bmax) - xmax) / axis.GetBinWidth(bmax)
    return integral

def log_significance(S, B):
    significance = 0
    significance = math.sqrt(2 * ((S + B) * math.log(1 + S / B) - S))
    return significance

def Get_BDT_cut_3D_v2(categ, year, file_name):
    t = TChain(out_tree_name)
    t.Add(file_name)
    print(f"Opened input file: {file_name}")

    cat = ""
    if categ == "A":
        cat = "0"
    elif categ == "B":
        cat = "1"
    elif categ == "C":
        cat = "2"

    isSB = ""
    if categ == "A":
        isSB = cutB_A
    elif categ == "B":
        isSB = cutB_B
    elif categ == "C":
        isSB = cutB_C

    signal = f"{weight}*(isMC>0 && isMC<4 && category=={cat})"
    bkg = f"{weight}*(isMC==0 && category=={cat} && ({isSB}))"

    N = 400
    N_str = str(N)
    binning = "("+N_str+",0.0,1.0)"

    t.Draw(f"bdt_cv>>h_test_bkg{binning}", bkg)
    h_test_bkg = gDirectory.Get("h_test_bkg")

    t.Draw(f"bdt_cv>>h_test_signal{binning}", signal)
    h_test_signal = gDirectory.Get("h_test_signal")

    h_test_signal.SetDirectory(0)
    h_test_bkg.SetDirectory(0)

    a, b, c = 0, 0, 0
    N_s_1, N_b_1, N_s_2, N_b_2, N_s_3, N_b_3, S1, S2, S3, S = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

    bkg_scale = 1
    if categ.startswith("A"):
        bkg_scale = 4. * 12. / (380.)
    elif categ.startswith("B"):
        bkg_scale = 4. * 19. / (380.)
    elif categ.startswith("C"):
        bkg_scale = 4. * 25. / (380.)
    print(f"bkg_scale = {bkg_scale}")

    X_min = 0.0
    X_max = 1.0
    #X_min = min(h_test_signal.GetXaxis().GetXmin(), h_test_bkg.GetXaxis().GetXmin())
    #X_max = max(h_test_signal.GetXaxis().GetXmax(), h_test_bkg.GetXaxis().GetXmax())
    #h_test_signal.GetXaxis().SetRangeUser(X_min, X_max)

    a_list = []
    b_list = []
    c_list = []
    AMS_list = []
    sig_int = TH1_integral(h_test_signal, X_min, X_max)
    bkg_int = TH1_integral(h_test_bkg, X_min, X_max)
    step = (X_max - X_min) / N
    for i in range(int(N/2), N):
        print("Category", categ, ": ", i-int(N/2), "/", int(N/2), end='\r')
        a = X_min + i * step
        for j in range(int(N/2), N):
            b = X_min + j * step
            if a < b:
                continue
            for k in range(int(N/2), N):
                c = X_min + k * step
                if b < c:
                    continue

                N_s_1 = TH1_integral(h_test_signal, a, X_max)
                N_b_1 = TH1_integral(h_test_bkg, a, X_max)
                if N_s_1 < sig_int * 0.0001:
                    continue
                if N_b_1 < bkg_int * 0.0001:
                    continue

                N_s_2 = TH1_integral(h_test_signal, b, a)
                N_b_2 = TH1_integral(h_test_bkg, b, a)
                if N_s_2 < sig_int * 0.0001:
                    continue
                if N_b_2 < bkg_int * 0.0001:
                    continue

                N_s_3 = TH1_integral(h_test_signal, c, b)
                N_b_3 = TH1_integral(h_test_bkg, c, b)
                if N_s_3 < sig_int * 0.0001:
                    continue
                if N_b_3 < bkg_int * 0.0001:
                    continue

                N_b_1 = N_b_1 * bkg_scale
                N_b_2 = N_b_2 * bkg_scale
                N_b_3 = N_b_3 * bkg_scale

                if N_b_1 > 0 and N_b_2 > 0 and N_b_3 > 0:
                    S = 0
                    S1 = log_significance(N_s_1, N_b_1)
                    S2 = log_significance(N_s_2, N_b_2)
                    S3 = log_significance(N_s_3, N_b_3)
                    S = math.sqrt(S1 * S1 + S2 * S2 + S3 * S3)
                    AMS_list.append(S)
                    a_list.append(a)
                    b_list.append(b)
                    c_list.append(c)

    #Taking absolute maximum of the combined significance
    
    AMS_max = max(AMS_list)
    AMS_maxID = AMS_list.index(AMS_max)

    bcx = a_list[AMS_maxID]
    bcx = round(bcx, 4)
    bcy = b_list[AMS_maxID]
    bcy = round(bcy, 4)
    bcz = c_list[AMS_maxID]
    bcz = round(bcz, 4)
    print(f"bcx={bcx} bcy={bcy} bcz={bcz}")

    print(f"S_max={AMS_max}")

    #Computing cut efficiency on signal
    N_S_12 = TH1_integral(h_test_signal, bcz, X_max)
    N_S_tot = sig_int
    print(f"Signal events kept by BDT {N_S_12} over {N_S_tot} ratio: {N_S_12/N_S_tot}")

    #Computing cut efficiency on backgroup
    N_B_12 = TH1_integral(h_test_bkg, bcz, X_max)
    N_B_tot = bkg_int
    print(f"Background events kept by BDT {N_B_12} over {N_B_tot} ratio: {N_B_12/N_B_tot}")

    return BDTcut3d(bcx, bcy, bcz)

def Get_BDT_cut_3D(categ, year, file_name):
    t = TChain(out_tree_name)
    t.Add(file_name)
    print(f"Opened input file: {file_name}")

    cat = ""
    if categ == "A":
        cat = "0"
    elif categ == "B":
        cat = "1"
    elif categ == "C":
        cat = "2"

    isSB = ""
    if categ == "A":
        isSB = cutB_A
    elif categ == "B":
        isSB = cutB_B
    elif categ == "C":
        isSB = cutB_C

    signal = f"{weight}*(isMC>0 && isMC<4 && category=={cat})"
    bkg = f"{weight}*(isMC==0 && category=={cat} && ({isSB}))"

    N = 1000
    N_str = str(N)
    binning = "("+N_str+",0.0,1.0)"

    t.Draw(f"bdt_cv>>h_test_bkg{binning}", bkg)
    h_test_bkg = gDirectory.Get("h_test_bkg")

    t.Draw(f"bdt_cv>>h_test_signal{binning}", signal)
    h_test_signal = gDirectory.Get("h_test_signal")

    h_test_signal.SetDirectory(0)
    h_test_bkg.SetDirectory(0)

    h3 = TH3F("h3", "test", N, 0.0, 1.0, N, 0.0, 1.0, N, 0.0, 1.0)
    a, b, c = 0, 0, 0
    N_s_1, N_b_1, N_s_2, N_b_2, N_s_3, N_b_3, S1, S2, S3, S = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

    bkg_scale = 1
    if categ.startswith("A"):
        bkg_scale = 4. * 12. / (380.)
    elif categ.startswith("B"):
        bkg_scale = 4. * 19. / (380.)
    elif categ.startswith("C"):
        bkg_scale = 4. * 25. / (380.)
    print(f"bkg_scale = {bkg_scale}")

    #bkg_scale = 1

    X_min = 0.0
    X_max = 1.0
    #X_min = min(h_test_signal.GetXaxis().GetXmin(), h_test_bkg.GetXaxis().GetXmin())
    #X_max = max(h_test_signal.GetXaxis().GetXmax(), h_test_bkg.GetXaxis().GetXmax())
    #h_test_signal.GetXaxis().SetRangeUser(X_min, X_max)

    dim = 0
    step = (X_max - X_min) / N
    for i in range(N):
        print("Category", categ, ": ", i, "/", N, end='\r')
        a = X_min + i * step
        for j in range(N):
            b = X_min + j * step
            if a < b:
                continue
            for k in range(N):
                c = X_min + k * step
                if b < c:
                    continue

                N_s_1 = TH1_integral(h_test_signal, a, X_max)
                N_b_1 = TH1_integral(h_test_bkg, a, X_max)
                if N_s_1 < TH1_integral(h_test_signal, X_min, X_max) * 0.0001:
                    continue
                if N_b_1 < TH1_integral(h_test_bkg, X_min, X_max) * 0.0001:
                    continue

                N_s_2 = TH1_integral(h_test_signal, b, a)
                N_b_2 = TH1_integral(h_test_bkg, b, a)
                if N_s_2 < TH1_integral(h_test_signal, X_min, X_max) * 0.0001:
                    continue
                if N_b_2 < TH1_integral(h_test_bkg, X_min, X_max) * 0.0001:
                    continue

                N_s_3 = TH1_integral(h_test_signal, c, b)
                N_b_3 = TH1_integral(h_test_bkg, c, b)
                if N_s_3 < TH1_integral(h_test_signal, X_min, X_max) * 0.0001:
                    continue
                if N_b_3 < TH1_integral(h_test_bkg, X_min, X_max) * 0.0001:
                    continue

                N_b_1 = N_b_1 * bkg_scale
                N_b_2 = N_b_2 * bkg_scale
                N_b_3 = N_b_3 * bkg_scale

                if N_b_1 > 0 and N_b_2 > 0 and N_b_3 > 0:
                    S = 0
                    S1 = log_significance(N_s_1, N_b_1)
                    S2 = log_significance(N_s_2, N_b_2)
                    S3 = log_significance(N_s_3, N_b_3)
                    S = math.sqrt(S1 * S1 + S2 * S2 + S3 * S3)
                    h3.SetBinContent(i, j, k, S)
                    dim += 1

    #Taking absolute maximum of the combined significance
    nbinx, nbiny, nbinz = ctypes.c_int(), ctypes.c_int(), ctypes.c_int()
    h3.GetMaximumBin(nbinx, nbiny, nbinz)
    nbinx = nbinx.value
    nbiny = nbiny.value
    nbinz = nbinz.value

    bcx = h3.GetXaxis().GetBinCenter(nbinx)
    bcx = round(bcx, 4)
    bcy = h3.GetYaxis().GetBinCenter(nbiny)
    bcy = round(bcy, 4)
    bcz = h3.GetZaxis().GetBinCenter(nbinz)
    bcz = round(bcz, 4)
    print(f"bcx={bcx} bcy={bcy} bcz={bcz}")

    S_max = h3.GetBinContent(h3.GetMaximumBin())
    print(f"S_max={S_max}")

    #Computing cut efficiency on signal
    N_S_12 = TH1_integral(h_test_signal, bcz, X_max)
    N_S_tot = TH1_integral(h_test_signal, X_min, X_max)
    print(f"Signal events kept by BDT {N_S_12} over {N_S_tot} ratio: {N_S_12/N_S_tot}")

    #Computing cut efficiency on backgroup
    N_B_12 = TH1_integral(h_test_bkg, bcz, X_max)
    N_B_tot = TH1_integral(h_test_bkg, X_min, X_max)
    print(f"Background events kept by BDT {N_B_12} over {N_B_tot} ratio: {N_B_12/N_B_tot}")

    return BDTcut3d(bcx, bcy, bcz)


def BDT_optimal_cut_v3(inputfile, year):
    ncat = len(cat_label)
    outputfile = "_" + year + "_BDT.txt"

    outputfile2 = "config_" + label_out + "_" + year +".txt"

    log = open(workdir + outputfile, "w")

    cuts = []
    for k in range(ncat):
        log.write("category {}\n".format(cat_label[k]))
        print("category {}".format(cat_label[k]))

        file_name = workdir + inputfile + ".root"

        # Open input files
        t = TChain(out_tree_name)
        t.Add(file_name)
        print("Opened input file: {}".format(file_name))

        isSB = ""

        if k == 0:
            isSB = cutB_A
        elif k == 1:
            isSB = cutB_B
        elif k == 2:
            isSB = cutB_C

        c = str(k)
        signal = "*(isMC>0 && isMC<4 && category=={})".format(c)
        bkg = "*(isMC==0 && category=={} && ({}))".format(c, isSB)

        signal = weight+signal
        bkg = weight+bkg
        
        # bdt score distribution
        binning = "(40, 0.0, 1.0)"
        t.Draw(f"bdt_cv>>h_test_bkg"+binning, bkg)
        h_test_bkg = gDirectory.Get("h_test_bkg").Clone("h_test_bkg")
        t.Draw(f"bdt_cv>>h_test_signal"+binning, signal)
        h_test_signal = gDirectory.Get("h_test_signal").Clone("h_test_signal")

        h_test_signal.SetDirectory(0)
        h_test_bkg.SetDirectory(0)
        
        X_min = min(h_test_signal.GetXaxis().GetXmin(), h_test_bkg.GetXaxis().GetXmin())
        X_max = max(h_test_signal.GetXaxis().GetXmax(), h_test_bkg.GetXaxis().GetXmax())

        #X_min = 0.8
        #X_max = 1.0

        # Compute cut and make colz plot
        cut_value = Get_BDT_cut_3D_v2(cat_label[k], year, file_name)
        cuts.append(cut_value)
        log.write("{},{},{}\n".format(cut_value.a, cut_value.b, cut_value.c))

        l = TLine()
        l.SetLineStyle(2)
        l.SetLineColor(1)
        
        # Drawing BDT score from scratch without signal normalization
        c2 = TCanvas("c2", "c2", 150, 10, 840, 540)
        gStyle.SetOptStat(0)
        gStyle.SetOptTitle(0)
        h_test_signal.GetXaxis().SetRangeUser(X_min, X_max)
        h_test_bkg.GetXaxis().SetRangeUser(X_min, X_max)
        h_test_bkg.SetLineColor(207);
        h_test_bkg.SetFillColor(207);
        h_test_bkg.SetFillStyle(3345);
        
        h_test_signal.SetLineColor(59);
        h_test_signal.SetFillColor(59);
        h_test_signal.SetFillStyle(3354);
        
        h_test_signal.SetLineWidth(2)
        h_test_bkg.SetLineWidth(2)
        h_test_bkg.GetXaxis().SetTitle("BDT score")
        h_test_signal.Scale(1 / h_test_signal.Integral())
        h_test_bkg.Scale(1 / h_test_bkg.Integral())
        #h_test_signal.Rebin(2)
        #h_test_bkg.Rebin(2)

        Y_max = 3 * h_test_signal.GetMaximum()

        h_test_bkg.Draw("HISTE")
        h_test_bkg.GetYaxis().SetRangeUser(1E-3, Y_max)
        h_test_bkg.GetXaxis().SetRangeUser(X_min, X_max)
        h_test_signal.Draw("same HISTE")
        c2.Update()
        l.DrawLine(cut_value.a, 1E-3, cut_value.a, Y_max)
        l.DrawLine(cut_value.b, 1E-3, cut_value.b, Y_max)
        l.DrawLine(cut_value.c, 1E-3, cut_value.c, Y_max)
        ta = TLatex(cut_value.a, 1E-3, "a")
        ta.Draw()
        tb = TLatex(cut_value.b, 1E-3, "b")
        tb.Draw()
        tc = TLatex(cut_value.c, 1E-3, "c")
        tc.Draw()

        Ltext = ROOT.TLatex(0.5, 0.83, "#bf{#it{CMS Preliminary}}")
        Ltext.SetTextFont(42);
        Ltext2 = ROOT.TLatex(0.5, 0.78, year+" category "+ cat_label[k])
        Ltext2.SetTextFont(42);
        Ltext2.SetTextSize(0.04);
        Ltext.SetNDC()
        Ltext2.SetNDC()
        Ltext.Draw("same")
        Ltext2.Draw("same")
        
        leg2 = TLegend(0.1, 0.75, 0.4, 0.9)
        #leg2.AddEntry(h_test_signal, "{} {} - signal".format(year, cat_label[k]), "f")
        #leg2.AddEntry(h_test_bkg, "{} {} - bkg".format(year, cat_label[k]), "f")
        leg2.AddEntry(h_test_signal, "MC -- signal", "f")
        leg2.AddEntry(h_test_bkg, "Data sidebands -- bkg", "f")
        leg2.Draw()
        c2.Update()
        c2.SaveAs(workdir + "_Cat_" + cat_label[k] + "_" + year + "_BDT_newnorm.png")

        # Same plot in log scale
        c2.SetLogy()
        c2.Update()
        c2.SaveAs(workdir +  "_Cat_" + cat_label[k] + "_" + year + "_BDT_log_newnorm.png")

    log.close()

    log2 = open(workdir + outputfile2, "w")
    log2.write("{},tripletMass,bdt_cv,category,isMC,combine_weight,dimu_OS1,dimu_OS2\n".format(out_tree_name))
    log2.write("A1,B1,C1,A2,B2,C2,A3,B3,C3\n")
    log2.write("{},{},{},".format(cuts[0].a, cuts[1].a, cuts[2].a))
    log2.write("{},{},{},".format(cuts[0].b, cuts[1].b, cuts[2].b))
    log2.write("{},{},{}\n".format(cuts[0].c, cuts[1].c, cuts[2].c))
    log2.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--config", type=str, help="Path to the copy of the JSON configuration file")
    args = parser.parse_args()
    config = args.config
    
    with open(config, 'r') as file:
        json_file = json.load(file)
    output_path = json_file['output_folder']
    date = json_file['date']
    date2 = json_file['date']
    label_out = json_file['label']
    out_tree_name = json_file['out_tree_name']
    pos_dir_xgboost = config.split(output_path)[0]
    weight = "combine_weight"

    if not output_path.endswith("/"):
        output_path += "/"

    if not pos_dir_xgboost.endswith("/"):
        pos_dir_xgboost += "/"

    if not date.endswith("/"):
        date += "/"

    if date.startswith("/"):
        date = date[1:]

    if output_path.startswith("/"):
        output_path = output_path[1:]
    
    inputfile = "t3mminitree_"+date
    workdir = pos_dir_xgboost + output_path + label_out + "_" + date
    print(workdir)

    year="2022"
    BDT_optimal_cut_v3(inputfile, year)
