
from ROOT import ROOT, gROOT, TFile, TTree, TH1F, TPad, TCanvas, TLine, TLegend, THStack, TGraph, gPad, TGraphErrors, TColor, TLatex, TArrow, TH2F, TMath,    TF1
from ROOT import kBlue, kRed, kDashed, kGreen, kGray, kBlack
from ROOT import gRandom
import sys
import math
import os
from array import array

gROOT.LoadMacro("AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

gROOT.SetBatch(True)

def plot(fres = 'result.txt', ex='ex0', ul_max = 20, tag=''):
    lines = []
    with open(fres, 'r') as f:
        lines = f.readlines()
        pass
    ul_exp_toy = 0.
    ul_exp_classic = 0.
    a_exp_new = array('f', [])

    ul_obs_toy = 0.
    ul_obs_classic = 0.
    a_obs_new = array('f', [])

    for line in lines:
        line = line.strip().split()
        if line[0] != ex:
            continue
        if line[1] == 'exp':
            ul_exp_toy = float(line[2])
            ul_exp_classic = abs(float(line[3])-ul_exp_toy)/ul_exp_toy*100
            for i in range(4, len(line)):
                a_exp_new.append(abs(float(line[i])-ul_exp_toy)/ul_exp_toy*100)
        if line[1] == 'obs':
            ul_obs_toy = float(line[2])
            ul_obs_classic = abs(float(line[3])-ul_obs_toy)/ul_obs_toy*100
            for i in range(4, len(line)):
                a_obs_new.append(abs(float(line[i])-ul_obs_toy)/ul_obs_toy*100)


    n = len(a_exp_new)
    #a_nsmall = array('f', [0.5, 1.5, 3.5, 5.5])
    a_nsmall = array('f', [-1, 0, 1, 3, 5])
    g_exp = TGraph(n, a_nsmall, a_exp_new)
    g_obs = TGraph(n, a_nsmall, a_obs_new)
    h2 = TH2F('h2', '', int(a_nsmall[-1]-a_nsmall[0])+1, int(a_nsmall[0])-0.5, int(a_nsmall[-1])+0.5, 100, 0, ul_max)
    Cs_nsmall = TCanvas('Cs_nsmall', '', 800, 600)
    h2.Draw()
    g_obs.Draw('PLsame')
    g_obs.SetLineWidth(2)
    g_exp.Draw('PLsame')
    g_exp.SetMarkerStyle(24)
    g_exp.SetLineStyle(2)
    g_exp.SetLineWidth(2)
    #g_obs.GetYaxis().SetRangeUser(0, ul_max)
    for i in range(h2.GetXaxis().GetNbins()):
        dn = -1+i
        #if dn < 0:
        #    h2.GetXaxis().SetBinLabel(i+1, 'n_{small}%i' % (dn))
        #elif dn == 0:
        #    h2.GetXaxis().SetBinLabel(i+1, 'n_{small}')
        #else:
        #    h2.GetXaxis().SetBinLabel(i+1, 'n_{small}+%i' % (dn))
        h2.GetXaxis().SetBinLabel(i+1, '%i' % (dn))
    h2.GetXaxis().SetLabelSize(0.08)
    h2.GetXaxis().SetTitle('n_{small}-n_{small}^{0}')
    h2.GetYaxis().SetTitle('|UL(A.F.)/UL(toy)-1|#times100 [%]')
    Cs_nsmall.Update()
    pxmin = Cs_nsmall.GetUxmin()
    pxmax = Cs_nsmall.GetUxmax()
    pymin = Cs_nsmall.GetUymin()
    pymax = Cs_nsmall.GetUymax()
    line_exp = TLine(pxmin, ul_exp_classic, pxmax, ul_exp_classic)
    line_exp.SetLineStyle(2)
    line_exp.SetLineWidth(2)
    line_exp.SetLineColor(kRed)
    line_exp.Draw()
    line_obs = TLine(pxmin, ul_obs_classic, pxmax, ul_obs_classic)
    #line_obs.SetLineStyle(2)
    line_obs.SetLineColor(kRed)
    line_obs.SetLineWidth(2)
    line_obs.Draw()
    ltx = TLatex()
    nsig = 1
    nbkg = 0
    if ex=='ex0':
        nbkg = 3.
    elif ex == 'ex1':
        nbkg = 8.11
    elif ex == 'sbottom':
        nsig = (3.758+3.561+2.010)*4
        nbkg = 2.817 + 5.325 + 9.616
    ltx.DrawLatex(pxmin+0.05*(pxmax-pxmin), pymin+0.9*(pymax-pymin), 's=%.1f, b=%.1f, n_{small}^{0}=b+#mus-1' % (nsig, nbkg))
    leg = TLegend(0.65, 0.7, 0.95, 0.92)
    leg.AddEntry(line_exp, 'Exp. (classic)', 'L')
    leg.AddEntry(line_obs, 'Obs. (classic)', 'L')
    leg.AddEntry(g_exp, 'Exp. (new)', 'LP')
    leg.AddEntry(g_obs, 'Obs. (new)', 'LP')
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.Draw()
    Cs_nsmall.SaveAs('Cs_nsmall_'+ex+'_'+tag+'.png')
    Cs_nsmall.SaveAs('Cs_nsmall_'+ex+'_'+tag+'.pdf')
    return

if __name__ == '__main__':
    #plot('result.txt', 'ex0', 15, 'data2')
    #plot('result.txt', 'ex1', 15, 'data2')
    plot('result.txt', 'sbottom', 30, 'sbottom')
