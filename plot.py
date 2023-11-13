import ROOT

def bdt_KS_plot(inputfile, year):
    kfold = 10
    kfold_s = str(kfold)

    file_name = " /analysis_xgboost_Run3/" + inputfile
    f = ROOT.TFile(file_name, "READ")
    t = f.Get("OutputTree")
    print("Opened input file:", file_name)

    outputfile = inputfile.replace(".root", "")

    binning = "(25, 0.0, 1.0)"
    varname = "bdt_fold1"

    train_sel = "int(evt)%" + kfold_s + "!=0"
    test_sel = "int(evt)%" + kfold_s + "==0"

    sig_sel = "isMC>0 && isMC<4"
    bkg_sel = "isMC==0"

    signal_label_all = "MC #tau#rightarrow3#mu"
    cat_label = ["categoryA", "categoryB", "categoryC"]

    for i in range(3):
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

        isSB = ""
        if i == 0:
            isSB = "((tripletMass<=1.753 && tripletMass>=1.62) || (tripletMass<=2.0 && tripletMass>=1.801))"
        elif i == 1:
            isSB = "((tripletMass<=1.739 && tripletMass>=1.62) || (tripletMass<=2.0 && tripletMass>=1.815))"
        elif i == 2:
            isSB = "((tripletMass<=1.727 && tripletMass>=1.62) || (tripletMass<=2.0 && tripletMass>=1.827))"

        categ = "category==" + str(i)

        t.Draw(varname + ">>hTrain_bkg" + binning, bkg_sel + "&&" + isSB + "&&" + train_sel + "&&" + phiveto + "&&" + categ)
        h_BDT_Train_B = ROOT.gDirectory.Get("hTrain_bkg")
        t.Draw(varname + ">>hTest_bkg" + binning, bkg_sel + "&&" + isSB + "&&" + test_sel + "&&" + phiveto + "&&" + categ)
        h_BDT_Test_B = ROOT.gDirectory.Get("hTest_bkg")
        t.Draw(varname + ">>hTrain_sig" + binning, "weight*(" + sig_sel + "&&" + train_sel + "&&" + phiveto + "&&" + categ + ")")
        h_BDT_Train_S = ROOT.gDirectory.Get("hTrain_sig")
        t.Draw(varname + ">>hTest_sig" + binning, "weight*(" + sig_sel + "&&" + test_sel + "&&" + phiveto + "&&" + categ + ")")
        h_BDT_Test_S = ROOT.gDirectory.Get("hTest_sig")

        t.Draw(varname + ">>hTrain_bkg2" + binning, varname + ">0.6&&" + bkg_sel + "&&" + isSB + "&&" + train_sel + "&&" + phiveto + "&&" + categ)
        h2_BDT_Train_B = ROOT.gDirectory.Get("hTrain_bkg2")
        t.Draw(varname + ">>hTest_bkg2" + binning, varname + ">0.6&&" + bkg_sel + "&&" + isSB + "&&" + test_sel + "&&" + phiveto + "&&" + categ)
        h2_BDT_Test_B = ROOT.gDirectory.Get("hTest_bkg2")
        t.Draw(varname + ">>hTrain_sig2" + binning, "weight*(" + varname + ">0.6&&" + sig_sel + "&&" + train_sel + "&&" + phiveto + "&&" + categ + ")")
        h2_BDT_Train_S = ROOT.gDirectory.Get("hTrain_sig2")
        t.Draw(varname + ">>hTest_sig2" + binning, "weight*(" + varname + ">0.6&&" + sig_sel + "&&" + test_sel + "&&" + phiveto + "&&" + categ + ")")
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

        c3.SaveAs(outputfile + "_" + cat_label[i] + "_KS.png")

bdt_KS_plot("your_input_file.root", "your_year")

