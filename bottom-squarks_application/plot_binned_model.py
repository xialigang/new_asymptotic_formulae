
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

from aux_funcs import get_sigmas



def plot_s0s1b0b1(fname_sigmas='output_sigmas_900_v1.txt', list_mu=[0.065], tag='0.065', mass='900'):
    list_hist_s_mu = []
    list_hist_b_mu = []
    list_hist_s_0 = []
    list_hist_b_0 = []
    for mu in list_mu:
        tag += '_'+str(mu)
    ymax = 0
    Nmax = 3
    for mu in list_mu:
        dict_mu_qmu, dict_mu_calib, dict_mumumu_sb = get_sigmas(mu, fname_sigmas, 1)
        print('dict_mumumu_sb =', dict_mumumu_sb)
        hist_s_mu = TH1F('hist_s_'+str(mu)+'_mu', '', Nmax, 0, Nmax)
        hist_b_mu = TH1F('hist_b_'+str(mu)+'_mu', '', Nmax, 0, Nmax)
        hist_s_0 = TH1F('hist_s_'+str(mu)+'_0', '', Nmax, 0, Nmax)
        hist_b_0 = TH1F('hist_b_'+str(mu)+'_0', '', Nmax, 0, Nmax)
        for i in range(3):
            hist_s_mu.SetBinContent(i+1, dict_mumumu_sb[(mu,mu,0)][0][i]*mu)
            hist_b_mu.SetBinContent(i+1, dict_mumumu_sb[(mu,mu,0)][1][i])
            hist_s_0.SetBinContent(i+1, dict_mumumu_sb[(mu,0,0)][0][i]*mu)
            hist_b_0.SetBinContent(i+1, dict_mumumu_sb[(mu,0,0)][1][i])
            _ymax = max(dict_mumumu_sb[(mu,mu,0)][0][i]*mu, dict_mumumu_sb[(mu,mu,0)][1][i], dict_mumumu_sb[(mu,0,0)][0][i]*mu,dict_mumumu_sb[(mu,0,0)][1][i])
            ymax = max(ymax, _ymax)
        list_hist_s_mu.append(hist_s_mu)
        list_hist_b_mu.append(hist_b_mu)
        list_hist_s_0.append(hist_s_0)
        list_hist_b_0.append(hist_b_0)
    Cs_sb = TCanvas('Cs_sb_'+tag, '', 800, 600)
    list_hist_s_mu[0].Draw()
    for i in range(len(list_mu)):
        list_hist_s_mu[i].Draw('same')
        list_hist_b_mu[i].Draw('same')
        list_hist_s_0[i].Draw('same')
        list_hist_b_0[i].Draw('same')
        list_hist_s_mu[i].SetLineColor(kRed)
        list_hist_s_0[i].SetLineColor(kRed)
        list_hist_s_0[i].SetLineStyle(kDashed)
        list_hist_b_0[i].SetLineStyle(kDashed)
    list_hist_s_mu[0].GetXaxis().SetTitle('bin index (m_{X}=%s)' % (mass))
    list_hist_s_mu[0].GetYaxis().SetTitle('Events')
    list_hist_s_mu[0].GetYaxis().SetRangeUser(0, 1.5*ymax)
    for i in range(Nmax):
        list_hist_s_mu[0].GetXaxis().SetBinLabel(i+1, '%i' % (i))
    leg = TLegend(0.2, 0.5, 0.6, 0.92)
    leg.AddEntry(list_hist_s_mu[0], 'sig.(#mu=%s,#mu_{H}=%s)' % (list_mu[0], list_mu[0]), 'L')
    leg.AddEntry(list_hist_b_mu[0], 'bkg.(#mu=%s,#mu_{H}=%s)' % (list_mu[0], list_mu[0]), 'L')
    leg.AddEntry(list_hist_s_0[0], 'sig.(#mu=%s,#mu_{H}=0)' % (list_mu[0]), 'L')
    leg.AddEntry(list_hist_b_0[0], 'bkg.(#mu=%s,#mu_{H}=0)' % (list_mu[0]), 'L')
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.Draw()
    Cs_sb.SaveAs('Cs_sb_'+tag+'.png')
    Cs_sb.SaveAs('Cs_sb_'+tag+'.pdf')

    return

option = 0
if len(sys.argv)>1:
    option = int(sys.argv[1])


if __name__ == '__main__':
    plot_s0s1b0b1('my_6bin_model.csv', [0.20], 'sbottom', 'sbottom')

