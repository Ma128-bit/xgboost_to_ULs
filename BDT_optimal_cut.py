import ROOT
import math

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
    t = ROOT.TChain("OutputTree")
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

    signal = f"weight_MC*(isMC>0 && isMC<4 && category=={cat} && {phiveto})"
    bkg = f"weight*(isMC==0 && category=={cat} && ({isSB}) && {phiveto})"

    N = 500
    binning = "(500,0.0,1.0)"

    t.Draw(f"bdt_cv>>h_test_bkg{binning}", bkg)
    h_test_bkg = ROOT.gDirectory.Get("h_test_bkg")
    h_test_bkg2 = ROOT.gDirectory.Get("h_test_bkg")
    t.Draw(f"bdt_cv>>h_test_signal{binning}", signal)
    h_test_signal = ROOT.gDirectory.Get("h_test_signal")
    h_test_signal2 = ROOT.gDirectory.Get("h_test_signal")

    h_test_signal.SetDirectory(0)
    h_test_bkg.SetDirectory(0)
    h_test_signal2.SetDirectory(0)
    h_test_bkg2.SetDirectory(0)

    h3 = ROOT.TH3F("h3", "test", N, 0.0, 1.0, N, 0.0, 1.0, N, 0.0, 1.0)
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
                    S = math.sqrt(S1*S1 + S2*S2 + S3*S3)
                    h3.SetBinContent(i, j, k, S)
                    dim += 1

    nbinx, nbiny, nbinz = ROOT.Long(), ROOT.Long(), ROOT.Long()
    print("h3.GetMaximumBin(nbinx, nbiny, nbinz)")
    h3.GetMaximumBin(nbinx, nbiny, nbinz)
    print(f"nbinx={nbinx}, nbiny={nbiny}, nbinz={nbinz}")

    bcx = h3.GetXaxis().GetBinCenter(nbinx)
    bcy = h3.GetYaxis().GetBinCenter(nbiny)
    bcz = h3.GetZaxis().GetBinCenter(nbinz)
    print(f"bcx={bcx} bcy={bcy} bcz={bcz}")

    S_max = h3.GetBinContent(h3.GetMaximumBin())
    print(f"S_max={S_max}")

    N_S_12 = TH1_integral(h_test_signal, bcz, X_max)
    N_S_tot = TH1_integral(h_test_signal, X_min, X_max)
    print(f"Signal events kept by BDT {N_S_12} over {N_S_tot} ratio: {N_S_12/N_S_tot}")

    N_B_12 = TH1_integral(h_test_bkg, bcz, X_max)
    N_B_tot = TH1_integral(h_test_bkg, X_min, X_max)
    print(f"Background events kept by BDT {N_B_12} over {N_B_tot} ratio: {N_B_12/N_B_tot}")

    return BDTcut3d(bcx, bcy, bcz)


# Example usage:
file_name = "your_input_file.root"
categ = "A"
year = "2023"
cut_value = Get_BDT_cut_3D(categ, year, file_name)
print(f"BDT cut values: {cut_value.a}, {cut_value.b}, {cut_value.c}")
