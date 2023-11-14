import ROOT
import math
import ctypes
from parameters import *
from ROOT import *
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

    reversedphiveto = ""
    if categ == "A":
        reversedphiveto = "((dimu_OS1<1.044 && dimu_OS1>0.994) || (dimu_OS2<1.044 && dimu_OS2>0.994))"
    elif categ == "B":
        reversedphiveto = "((dimu_OS1<1.053 && dimu_OS1>0.985) || (dimu_OS2<1.053 && dimu_OS2>0.985))"
    elif categ == "C":
        reversedphiveto = "((dimu_OS1<1.064 && dimu_OS1>0.974) || (dimu_OS2<1.064 && dimu_OS2>0.974))"

    phiveto = ""
    if categ == "A":
        phiveto = "(!(dimu_OS1<1.044 && dimu_OS1>0.994) && !(dimu_OS2<1.044 && dimu_OS2>0.994))"
    elif categ == "B":
        phiveto = "(!(dimu_OS1<1.053 && dimu_OS1>0.985) && !(dimu_OS2<1.053 && dimu_OS2>0.985))"
    elif categ == "C":
        phiveto = "(!(dimu_OS1<1.064 && dimu_OS1>0.974) && !(dimu_OS2<1.064 && dimu_OS2>0.974))"

    signal = f"{weight}*(isMC>0 && isMC<4 && category=={cat} && {phiveto})"
    bkg = f"{weight}*(isMC==0 && category=={cat} && ({isSB}) && {phiveto})"

    N = 300
    N_str = str(N)
    binning = "("+N_str+",0.0,1.0)"

    t.Draw(f"bdt_cv>>h_test_bkg{binning}", bkg)
    h_test_bkg = gDirectory.Get("h_test_bkg")
    h_test_bkg2 = gDirectory.Get("h_test_bkg")
    t.Draw(f"bdt_cv>>h_test_signal{binning}", signal)
    h_test_signal = gDirectory.Get("h_test_signal")
    h_test_signal2 = gDirectory.Get("h_test_signal")

    h_test_signal.SetDirectory(0)
    h_test_bkg.SetDirectory(0)
    h_test_signal2.SetDirectory(0)
    h_test_bkg2.SetDirectory(0)

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

    X_min = min(h_test_signal.GetXaxis().GetXmin(), h_test_signal.GetXaxis().GetXmin())
    X_max = max(h_test_signal.GetXaxis().GetXmax(), h_test_signal.GetXaxis().GetXmax())
    h_test_signal.GetXaxis().SetRangeUser(X_min, X_max)

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
    print(f"nbinx={nbinx}, nbiny={nbiny}, nbinz={nbinz}")
    print("h3.GetMaximumBin(nbinx, nbiny, nbinz)")
    h3.GetMaximumBin(nbinx, nbiny, nbinz)
    nbinx = nbinx.value
    nbiny = nbiny.value
    nbinz = nbinz.value
    print(f"nbinx={nbinx}, nbiny={nbiny}, nbinz={nbinz}")

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

from ROOT import *
from array import array

def BDT_optimal_cut_v3(inputfile, year):
    ncat = len(cat_label)
    outputfile = inputfile.replace(".root", "") + "_" + year + "_BDT.txt"

    log = open(workdir + "/" + outputfile, "w")

    for k in range(ncat):
        log.write("category {}\n".format(cat_label[k]))
        print("category {}".format(cat_label[k]))

        file_name = workdir + "/" + inputfile

        # Open input files
        t = TChain(out_tree_name)
        t.Add(file_name)
        print("Opened input file: {}".format(file_name))

        isSB = ""
        reversedphiveto = ""
        phiveto = ""

        if k == 0:
            isSB = cutB_A
            reversedphiveto = "((dimu_OS1<1.044 && dimu_OS1>0.994) || (dimu_OS2<1.044 && dimu_OS2>0.994))"
            phiveto = "(!(dimu_OS1<1.044 && dimu_OS1>0.994) && !(dimu_OS2<1.044 && dimu_OS2>0.994))"
        elif k == 1:
            isSB = cutB_B
            reversedphiveto = "((dimu_OS1<1.053 && dimu_OS1>0.985) || (dimu_OS2<1.053 && dimu_OS2>0.985))"
            phiveto = "(!(dimu_OS1<1.053 && dimu_OS1>0.985) && !(dimu_OS2<1.053 && dimu_OS2>0.985))"
        elif k == 2:
            isSB = cutB_C
            reversedphiveto = "((dimu_OS1<1.064 && dimu_OS1>0.974) || (dimu_OS2<1.064 && dimu_OS2>0.974))"
            phiveto = "(!(dimu_OS1<1.064 && dimu_OS1>0.974) && !(dimu_OS2<1.064 && dimu_OS2>0.974))"

        c = str(k)
        signal = "*(isMC>0 && isMC<4 && category=={} && {})".format(c, phiveto)
        bkg = "*(isMC==0 && category=={} && ({}) && {})".format(c, isSB, phiveto)

        signal = weight+signal
        bkg = weight+bkg
        
        # bdt score distribution
        binning = (500, 0.0, 1.0)
        t.Draw(f"bdt_cv>>h_test_bkg{binning}", bkg)
        h_test_bkg = gDirectory.Get("h_test_bkg").Clone("h_test_bkg")
        h_test_bkg2 = gDirectory.Get("h_test_bkg").Clone("h_test_bkg2")
        t.Draw(f"bdt_cv>>h_test_signal{binning}", signal)
        h_test_signal = gDirectory.Get("h_test_signal").Clone("h_test_signal")
        h_test_signal2 = gDirectory.Get("h_test_signal").Clone("h_test_signal2")

        h_test_signal.SetDirectory(0)
        h_test_bkg.SetDirectory(0)
        h_test_signal2.SetDirectory(0)
        h_test_bkg2.SetDirectory(0)

        h_test_bkg.SetLineColor(kBlack);
        h_test_signal.SetLineColor(kRed);

        X_min = min(h_test_signal.GetXaxis().GetXmin(), h_test_signal.GetXaxis().GetXmin())
        X_max = max(h_test_signal.GetXaxis().GetXmax(), h_test_signal.GetXaxis().GetXmax())

        X_min = 0.8
        X_max = 1.0

        # Compute cut and make colz plot
        cut_value = Get_BDT_cut_3D(cat_label[k], year, file_name)
        log.write("{},{},{}\n".format(cut_value.a, cut_value.b, cut_value.c))

        c1 = TCanvas("c1", "c1", 150, 10, 800, 800)
        h_test_bkg.Draw("HISTE")
        h_test_signal.Draw("same HISTE")
        c1.Update()
        l = TLine()
        l.DrawLine(cut_value.a, 0, cut_value.a, 0.1)
        l.DrawLine(cut_value.b, 0, cut_value.b, 0.1)
        l.DrawLine(cut_value.c, 0, cut_value.c, 0.1)

        leg = TLegend(0.35, 0.75, 0.65, 0.9)
        leg.AddEntry(h_test_signal, "{}_signal".format(cat_label[k]), "f")
        leg.AddEntry(h_test_bkg, "{}_bkg".format(cat_label[k]), "f")
        leg.Draw()
        c1.Update()
        c1.SaveAs(workdir + "/" + outputfile + "_" + cat_label[k] + "_" + year + "_normBDT_newnorm.png")

        # Drawing BDT score from scratch without signal normalization
        c2 = TCanvas("c2", "c2", 150, 10, 800, 800)
        gStyle.SetOptStat(0)
        gStyle.SetOptTitle(0)
        h_test_signal2.GetXaxis().SetRangeUser(X_min, X_max)
        h_test_bkg2.GetXaxis().SetRangeUser(X_min, X_max)
        h_test_bkg2.SetLineColor(kBlack)
        h_test_signal2.SetLineColor(kRed)
        h_test_signal2.SetLineWidth(2)
        h_test_bkg2.SetLineWidth(2)
        h_test_bkg2.GetXaxis().SetTitle("BDT score")
        h_test_signal2.Scale(1 / h_test_signal2.Integral())
        h_test_bkg2.Scale(1 / h_test_bkg2.Integral())
        h_test_signal2.Rebin(4)
        h_test_bkg2.Rebin(4)

        Y_max = 0.3

        h_test_bkg2.Draw("HISTE")
        h_test_bkg2.GetYaxis().SetRangeUser(1E-5, Y_max)
        h_test_bkg2.GetXaxis().SetRangeUser(X_min, X_max)
        h_test_signal2.Draw("same HISTE")
        c2.Update()
        l.DrawLine(cut_value.a, 1E-5, cut_value.a, Y_max)
        l.DrawLine(cut_value.b, 1E-5, cut_value.b, Y_max)
        l.DrawLine(cut_value.c, 1E-5, cut_value.c, Y_max)
        ta = TLatex(cut_value.a, 1E-5, "a")
        ta.Draw()
        tb = TLatex(cut_value.b, 1E-5, "b")
        tb.Draw()
        tc = TLatex(cut_value.c, 1E-5, "c")
        tc.Draw()

        leg2 = TLegend(0.35, 0.75, 0.65, 0.9)
        leg2.AddEntry(h_test_signal2, "{} {} - signal".format(year, cat_label[k]), "f")
        leg2.AddEntry(h_test_bkg2, "{} {} - bkg".format(year, cat_label[k]), "f")
        leg2.Draw()
        c2.Update()
        c2.SaveAs(workdir + "/" + outputfile + "_" + cat_label[k] + "_" + year + "_BDT_newnorm.png")

        # Same plot in log scale
        c2.SetLogy()
        c2.Update()
        c2.SaveAs(workdir + "/" + outputfile + "_" + cat_label[k] + "_" + year + "_BDT_log_newnorm.png")

    log.close()
    print("Exiting ROOT")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--config", type=str, help="Path to the copy of the JSON configuration file")
    with open(config, 'r') as file:
        json_file = json.load(file)
    output_path = json_file['output_path']
    date = json_file['date']
    inputfile = json_file['Name']
    out_tree_name = json_file['out_tree_name']
    pos_dir_xgboost = config.split(output_path)[0]
    weight = json_file['weight_column']

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
    
    inputfile = inputfile +"_minitree.root"
    workdir = pos_dir_xgboost + output_path + date

    year="2022"
    BDT_optimal_cut_v3(inputfile, year)
