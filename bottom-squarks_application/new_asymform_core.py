from ROOT import ROOT, gROOT, TFile, TTree, TH1F, TPad, TCanvas, TLine, TLegend, THStack, TGraph, gPad, TGraphErrors, TColor, TLatex, TArrow, TH2F, TMath, TF1
from ROOT import kBlue, kRed, kDashed, kGreen, kGray, kBlack
from ROOT import gRandom
import sys
import math
import os
from array import array
#from scipy import optimize
#from scipy import integrate
#from fit import dofit
#from CDFexact import get_C, get_nll2, get_muhat_alpha_beta, get_qmu_from_mu_alpha_beta
#from CDFhist import read_obs_hist, plot_llscan
#from solvefunc import newton_method

from aux_funcs import *

gROOT.LoadMacro("AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

gROOT.SetBatch(True)

yaxis_low = 1e-6

def do_one_mu_for_one_mass(mass='ex0', mu=0.045, nbins = 100, xmin = 0, xmax = 20, logy=0, fast=0, do_q0=False, fname_toy_mu='', fname_toy_0='', fname_sigmas='', fname_interval='', fname_calib='', list_special=[], doTilde=1):

    #print('***************len(dict_toy_muhat) = ', len(dict_toy_muhat))
    qmu_exp = 0.
    qmu_obs = 0.

    if is_missing([fname_toy_mu, fname_toy_0]):
        return 0

    mupre = -1
    qmu_mupre = -1
    #qmu_exp, qmu_obs, qmu_0 = get_sigmas(mu, fname_sigmas, 0)
    dict_mu_qmu, dict_mu_calib, dict_mumumu_sb = get_sigmas(mu, fname_sigmas, 1)

    print('dict_mu_qmu =', dict_mu_qmu)
    print('dict_mu_calib =', dict_mu_calib)
    print('dict_mumumu_sb =', dict_mumumu_sb)

    qmu_exp = dict_mu_qmu[mu][0]
    qmu_obs = dict_mu_qmu[mu][1]
    # = {}
    #qmu_exp, qmu_obs, sigma_muH, sigma_muH0, sigma_muH_in_muH0, mu_key,  = get_sigmas3(mu, fname_sigmas, 0)
    #qmu_exp, qmu_obs, sigma_muH, sigma_muH0, mu_key, qmu_key = get_sigmas3(mu, fname_sigmas, 0)
    #qmu_exp, qmu_obs, sigma_muH, sigma_muH0, mu_key, qmu_key, mupre, qmu_mupre = get_sigmas4(mu,'output_sigmas_'+mass+'_v17.txt', 0)
    #qmu_exp, qmu_obs, sigma_muH, sigma_muH0, mu_key, qmu_key, mupre, qmu_mupre = get_sigmas4(mu,'output_sigmas_'+mass+'_v17.txt', 0)

    if qmu_obs < 0:
        qmu_obs = 0
    sigma = 0
    list_q0_calib = []
    if do_q0:
        list_q0_calib = get_sigmas3(mu, fname_sigmas, 1)
        print('list_q0_calib =', list_q0_calib)
    if do_q0 == 0:
        sigma = mu/math.sqrt(qmu_exp)


    #print('dict_mu_sb =', dict_mu_sb)
    #print('dict_mumu_calib =', dict_mumu_calib)

    print('qmu_exp, qmu_obs =', qmu_exp, qmu_obs)

    #sigma_muH0 = dict_mumu_calib[(0,mu)][2]
    #sigma_muH = dict_mumu_calib[(mu,mu)][1]

    sigma_muH0 = 0
    sigma_muH = 0
    if len(dict_mu_calib)>0:
        sigma_muH0 = dict_mu_calib[mu][0]
        sigma_muH = dict_mu_calib[mu][1]
    
    sigma = mu/math.sqrt(qmu_exp)
    #sigma_muH_in_muH0 = dict_mumu_calib[(0,mu)][1]
    if do_q0:
        if mass == 'ex0':
            sigma = sigma_muH0
        else:
            sigma = sigma_muH0
    print('sigma(wald) =', sigma)

    h_toy_mu = TH1F('h_toy_mu', '', nbins, xmin, xmax)
    h_toy_0 = TH1F('h_toy_0', '', nbins, xmin, xmax)

    list_sigma = []
    if len(list_special) > 0:
        list_sigma += list_special
    list_sigma += [sigma, sigma_muH0, sigma_muH]
    #list_sigma.append(qmu_0)

    if fname_toy_0.endswith('.root'):
        p0_toy_exp, dp0_toy_exp, p0_toy_obs, dp0_toy_obs = get_toys(h_toy_0, fname_toy_0, 'tree_stat', qmu_exp, qmu_obs, -1, mu, 0, list_sigma, -0.5*mu, mu*1.5, -1, 3*pow(mu/sigma,2), tag=mass, fast = fast, do_q0 = do_q0, doTilde = doTilde)
    else:
        p0_toy_exp, dp0_toy_exp, p0_toy_obs, dp0_toy_obs = get_toys_csv(h_toy_0, fname_toy_0, 'tree_stat', qmu_exp, qmu_obs, -1, mu, 0, list_sigma, -0.5*mu, mu*1.5, -1, 3*pow(mu/sigma,2), tag=mass, fast = fast, do_q0 = do_q0, doTilde = doTilde)
    if not do_q0:
        if fname_toy_mu.endswith('.root'):
            pmu_toy_exp, dpmu_toy_exp, pmu_toy_obs, dpmu_toy_obs = get_toys(h_toy_mu, fname_toy_mu, 'tree_stat', qmu_exp, qmu_obs, -1, mu, mu, list_sigma, -0.5*mu, mu*1.5, -1, 3*pow(mu/sigma,2), tag=mass, fast = fast, doTilde = doTilde)
        else:
            pmu_toy_exp, dpmu_toy_exp, pmu_toy_obs, dpmu_toy_obs = get_toys_csv(h_toy_mu, fname_toy_mu, 'tree_stat', qmu_exp, qmu_obs, -1, mu, mu, list_sigma, -0.5*mu, mu*1.5, -1, 3*pow(mu/sigma,2), tag=mass, fast = fast, doTilde = doTilde)
    if do_q0:
        p0_new_obs = 1.-Fnewq0(qmu_obs, mu, 0, [sigma], dict_mu_sb, list_q0_calib=list_q0_calib)
        p0_classic_obs = 1. - CDFasymq0(qmu_obs, 0, sigma, sigma)
        print('debug: classic(1), new(1) =', CDFasymq0(1000,0, sigma, sigma), Fnewq0(1000, mu, 0, [sigma], dict_mu_sb, list_q0_calib=list_q0_calib))
        h_classic_0 = TH1F('h_classic_0', '', nbins, xmin, xmax)
        h_new_0 = TH1F('h_new_0', '', nbins, xmin, xmax)
        p_0_pre = 0.
        for i in range(nbins):
            print('--------------------------------------------working on bin', i)
            x1 = h_classic_0.GetBinLowEdge(i+1)
            x2 = h_classic_0.GetBinLowEdge(i+2)
            # classic formulae
            p_0 = CDFasymq0(x2, 0, sigma, sigma)
            if i>0:
                p_0 -= CDFasymq0(x1, 0, sigma, sigma)
            h_classic_0.SetBinContent(i+1, p_0)

            p_0 = Fnewq0(x2, mu, 0, [sigma], dict_mu_sb, list_q0_calib=list_q0_calib)
            h_new_0.SetBinContent(i+1, p_0 - p_0_pre)
            p_0_pre = p_0
        print('debug: h_classic_0, h_new_0 =', h_classic_0.Integral(), h_new_0.Integral())
        Cs_q0 = TCanvas('Cs_q0', '', 800, 600)
        if logy:
            Cs_q0.SetLogy()
        h_toy_0.Draw('PE')
        h_toy_0.SetMarkerStyle(20)
        if logy:
            h_toy_0.GetYaxis().SetRangeUser(yaxis_low, 1000)
        else:
            h_toy_0.GetYaxis().SetRangeUser(0, h_toy_mu.GetMaximum()*1.2)
        h_toy_0.GetYaxis().SetTitle('Probability')
        #h_toy_0.GetXaxis().SetTitle('#tilde{q}_{#mu} (%s, #mu=%.3f)' % (mass, mu))
        h_toy_0.GetXaxis().SetTitle('q_{0} (%s)' % (mass))
        #h_toy_0.GetXaxis().SetTitle('q_{#mu} (%s, #mu=%.3f)' % (mass, mu))
        h_toy_0.Draw('PE,same')
        #classic
        h_classic_0.Draw('hist,same')
        h_classic_0.SetLineColor(kRed-7)
        h_classic_0.SetLineStyle(2)
        #new
        h_new_0.Draw('hist,same')
        h_new_0.SetLineColor(kRed)

        Cs_q0.Update()
        Uymin = Cs_q0.GetUymin()
        Uymax = Cs_q0.GetUymax()
        print('Uymin, Uymax =', Uymin, Uymax)
        #a_exp = TArrow(qmu_exp, Uymin+0.8*(Uymax-Uymin), qmu_exp, Uymin, 0.017, '|>')
        #a_exp.Draw()
        #a_exp.SetFillColor(kGray)
        #a_exp.SetLineColor(kGray)
        #a_exp.SetLineWidth(2)

        a_obs = TArrow(qmu_obs, Uymin+0.8*(Uymax-Uymin), qmu_obs, Uymin, 0.017, '|>')
        #a_obs.Draw()
        a_obs.SetLineWidth(2)
        ltx = TLatex()
        #ltx.DrawLatex(qmu_obs, turn2log(Uymin+0.1*(Uymax-Uymin), logy), 'q_{0}^{obs}')
        #ltx.DrawLatex(xmin + 0.05*(xmax-xmin),  turn2log(Uymin+0.9*(Uymax-Uymin), logy), 'p_{0} obs.(toy)=%.3f#pm%.3f%%' % (p0_toy_obs*100, dp0_toy_obs*100))
        #ltx.DrawLatex(xmin + 0.05*(xmax-xmin),  turn2log(Uymin+0.85*(Uymax-Uymin), logy), 'p_{0} obs.(classic)=%.3f%%' % (p0_classic_obs*100))
        #ltx.DrawLatex(xmin + 0.05*(xmax-xmin),  turn2log(Uymin+0.80*(Uymax-Uymin), logy), 'p_{0} obs.(new)=%.3f%%' % (p0_new_obs*100))

        leg = TLegend(0.65, 0.5, 0.95, 0.92)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.Draw()
        leg.AddEntry(h_toy_0, 'toy(#mu_{H}=0)', 'PE')
        #leg.AddEntry(h_toy_mu, 'toy(#mu_{H}=#mu)', 'PE')
        leg.AddEntry(h_classic_0, 'classic(#mu_{H}=0)', 'L')
        #leg.AddEntry(h_classic_mu, 'classic(#mu_{H}=#mu)', 'L')
        leg.AddEntry(h_new_0, 'new(#mu_{H}=0)', 'L')
        #leg.AddEntry(h_new_mu, 'new(#mu_{H}=#mu)', 'L')
        plotname = 'Cs_q0_%s_mu%.3f_logy%i' % (mass, mu, logy)
        Cs_q0.SaveAs(plotname+'.png')
        Cs_q0.SaveAs(plotname+'.pdf')
        Cs_q0.SaveAs(plotname+'.root')


        sigma_muH0 = sigma
        a_q0 = array('f', [])
        a_dq0 = array('f', [])
        a_Z_classic = array('f', [])
        a_Z_new = array('f', [])
        a_Z_toys = array('f', [])
        a_dZ_toys = array('f', [])
        list_q0 = [0.1, 0.5, 1, 2, 4, 6, 8, 10, 12, 14, 16]
        if mass == 'ex1' and 0:
            list_q0 = [0.1, 0.5, 1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
        for q0_obs in list_q0:
            p0_toy_exp, dp0_toy_exp, p0_toy_obs, dp0_toy_obs = get_toys(h_toy_0, fname_toy_0, 'tree_stat', qmu_exp, q0_obs, -1., mu, 0, list_sigma, -0.5*mu, mu*1.5, -1, 3*pow(mu/sigma,2), tag=mass, fast = 1, do_q0 = do_q0, doTilde = doTilde)
            Z_classic = Phi_quantile(CDFasymq0(q0_obs,0, sigma_muH0, sigma_muH0))
            Z_new = Phi_quantile(Fnewq0(q0_obs, mu, 0, [sigma_muH0], dict_mu_sb, list_q0_calib=list_q0_calib))
            Z_toys = Phi_quantile(1.-p0_toy_obs)
            Z_toys_up = Phi_quantile(1.-p0_toy_obs-dp0_toy_obs)
            Z_toys_dn = Phi_quantile(1.-p0_toy_obs+dp0_toy_obs)
            a_q0.append(q0_obs)
            a_dq0.append(0)
            a_Z_classic.append(Z_classic)
            a_Z_new.append(Z_new)
            a_Z_toys.append(Z_toys)
            a_dZ_toys.append(max(abs(Z_toys_up-Z_toys), abs(Z_toys_dn-Z_toys)))
            print('q0_obs, Z_classic, Z_new, Z_toys =', p0_toy_obs, dp0_toy_obs, q0_obs, Z_classic, Z_new, Z_toys, Z_toys_up, Z_toys_dn)
        nZ = len(a_q0)
        g_Z_classic = TGraph(nZ, a_q0, a_Z_classic)
        g_Z_new = TGraph(nZ, a_q0, a_Z_new)
        g_Z_toys = TGraphErrors(nZ, a_q0, a_Z_toys, a_dq0, a_dZ_toys)
        Cs_Z = TCanvas('Cs_Z', '', 800, 600)
        g_Z_toys.Draw('APE')
        g_Z_toys.SetMarkerStyle(20)
        g_Z_toys.SetMarkerSize(1)
        g_Z_classic.Draw('Lsame')
        g_Z_classic.SetLineStyle(2)
        g_Z_classic.SetLineWidth(2)
        g_Z_classic.SetLineColor(kRed)
        g_Z_new.Draw('Lsame')
        g_Z_new.SetLineWidth(2)
        g_Z_new.SetLineColor(kGreen)
        #g_Z_toys.GetXaxis().SetTitle('q_{0}^{obs}')
        g_Z_toys.GetXaxis().SetTitle('q_{0}')
        g_Z_toys.GetYaxis().SetTitle('Z')
        g_Z_toys.GetYaxis().SetRangeUser(0, 5)
        leg_Z = TLegend(0.5,0.2, 0.95, 0.5)
        leg_Z.AddEntry(g_Z_toys, 'toys', 'PE')
        leg_Z.AddEntry(g_Z_classic, 'classic', 'L')
        leg_Z.AddEntry(g_Z_new, 'new', 'L')
        leg_Z.SetFillStyle(0)
        leg_Z.SetBorderSize(0)
        leg_Z.Draw()
        plotname = 'Cs_Z_%s_mu%.3f_logy%i' % (mass, mu, logy)
        Cs_Z.SaveAs(plotname+'.png')
        Cs_Z.SaveAs(plotname+'.pdf')
        return 1


    CLs_toy_exp = pmu_toy_exp / p0_toy_exp
    CLs_toy_obs = pmu_toy_obs / p0_toy_obs


    dCLs_toy_exp_pmu = CLs_toy_exp * dpmu_toy_exp/pmu_toy_exp
    dCLs_toy_exp_p0 = CLs_toy_exp * dp0_toy_exp/p0_toy_exp
    dCLs_toy_exp = math.sqrt(pow( dCLs_toy_exp_pmu,2) + pow( dCLs_toy_exp_p0,2))

    dCLs_toy_obs_pmu = CLs_toy_obs * dpmu_toy_obs/pmu_toy_obs
    dCLs_toy_obs_p0 = CLs_toy_obs * dp0_toy_obs/p0_toy_obs
    dCLs_toy_obs = math.sqrt(pow( dCLs_toy_obs_pmu,2) + pow( dCLs_toy_obs_p0,2))

    pmu_classic_exp = 1.-Fclassic(qmu_exp, mu, mu, sigma, doTilde = doTilde)
    p0_classic_exp = 1. - Fclassic(qmu_exp, mu, 0, sigma, doTilde = doTilde)
    CLs_classic_exp = pmu_classic_exp / p0_classic_exp

    pmu_classic_obs = 1.-Fclassic(qmu_obs, mu, mu, sigma, doTilde = doTilde)
    p0_classic_obs = 1. - Fclassic(qmu_obs, mu, 0, sigma, doTilde = doTilde)
    CLs_classic_obs = pmu_classic_obs / p0_classic_obs

    ############### determine the optimal nsmall below
    #list_s = dict_mu_sb[(mu,mu)][0]
    #list_b = dict_mu_sb[(mu,mu)][1]
    #stot = 0
    #btot = 0
    #b5 = list_b[5]
    #if b5>10:
    #    b5 = 10
    #for i in range(6):
    #    si = list_s[i]
    #    bi = list_b[i]
    #    if i == 5:
    #        si = 10/bi*si
    #        bi = 10
    #    stot += si
    #    btot += bi
    #n0 = 0
    #prob0 = 0.
    #k = 0
    #while 1:
    #    prob0 += TMath.Poisson(k, btot + mu*stot)
    #    if prob0 > 0.001:
    #        n0 = k
    #        break
    #    k += 1
    #    pass
    #print('n0, prob0 =', n0, prob0)
    #pmu0_new_exp = 1.-Fnew(qmu_exp, mu, mu, list_sigma, dict_mu_sb,   = dict_mu_qmu, dict_mumu_calib = dict_mumu_calib, doTilde = doTilde, nsmall=n0)
    #direction = 1
    #if pmu0_new_exp < pmu_classic_exp:
    #    direction = -1
    #best_nsmall_i_think = int(btot+mu*stot)-1
    #if best_nsmall_i_think>max(6,b5) and 0:
    #    best_nsmall_i_think = max(6,b5)
    nsmall = -1
    #_pmu_new_exp = pmu_classic_exp
    #for n in range(best_nsmall_i_think+1):
    #    pmu_new_exp = 1.-Fnew(qmu_exp, mu, mu, list_sigma, dict_mu_sb,   = dict_mu_qmu, dict_mumu_calib = dict_mumu_calib, doTilde = doTilde, nsmall = n)
    #    print("*****************", n, pmu_toy_exp, pmu_classic_exp, pmu_new_exp, pmu_classic_exp-pmu_toy_exp, pmu_new_exp-pmu_toy_exp)
    #    if (direction == 1 and pmu_new_exp >= _pmu_new_exp) or (direction==-1 and pmu_new_exp <= _pmu_new_exp):
    #        _pmu_nex_exp = pmu_new_exp
    #    if (direction == 1 and pmu_new_exp < _pmu_new_exp) or (direction==-1 and pmu_new_exp > _pmu_new_exp):
    #        nsmall = n-1
    #        break
    print('>>>>>>>>>>>>>>>>>>>best nsmall =',nsmall)
    nsmall0 = nsmall
    #for n in range(best_nsmall_i_think+1):
    #for n in range(nsmall+1):
    #    p0_new_exp  = 1.-Fnew(qmu_exp, mu, 0, list_sigma, dict_mu_sb,   = dict_mu_qmu, dict_mumu_calib = dict_mumu_calib, doTilde = doTilde, nsmall = n)
    #    CLs_new_exp = _pmu_new_exp / p0_new_exp
    #    #print("*****************", n, pmu_toy_exp, pmu_classic_exp, pmu_new_exp, pmu_classic_exp-pmu_toy_exp, pmu_new_exp-pmu_toy_exp)
    #    print("*****************", n, p0_toy_exp, p0_classic_exp, p0_new_exp, p0_classic_exp-p0_toy_exp, p0_new_exp-p0_toy_exp)
    #    if direction == 1 and CLs_new_exp >= CLs_classic_exp:
    #        nsmall0 = n-1
    #        #break
    #    elif direction == -1 and CLs_new_exp <= CLs_classic_exp:
    #        nsmall0 = n-1
    #        #break
    print('>>>>>>>>>>>>>>>>>>>best nsmall0 =',nsmall0)

    ########################################

    nsmall = -9999
    nsmall0= -9999
    print('qmu_exp, qmu_obs =', qmu_exp, qmu_obs)
    #pmu_new_exp = 1.-Fnew(qmu_exp, mu, mu, list_sigma, dict_mu_sb,   = dict_mu_qmu, dict_mumu_calib = dict_mumu_calib, doTilde = doTilde)
    #print('pmu_new_exp =', pmu_new_exp)
    pmu_new_exp = 1.-Fnew(qmu_exp, mu, mu, list_sigma, dict_mumumu_sb,  dict_mu_qmu = dict_mu_qmu, dict_mu_calib = dict_mu_calib, doTilde = doTilde, nsmall = nsmall)
    p0_new_exp  = 1.-Fnew(qmu_exp, mu, 0, list_sigma, dict_mumumu_sb,  dict_mu_qmu = dict_mu_qmu, dict_mu_calib = dict_mu_calib, doTilde = doTilde, nsmall = nsmall0)
    CLs_new_exp = 1
    if p0_new_exp > 0:
        CLs_new_exp = pmu_new_exp / p0_new_exp

    pmu_new_obs = 1.-Fnew(qmu_obs, mu, mu, list_sigma, dict_mumumu_sb,  dict_mu_qmu = dict_mu_qmu, dict_mu_calib = dict_mu_calib, doTilde = doTilde, nsmall = nsmall)
    p0_new_obs  = 1.-Fnew(qmu_obs, mu, 0, list_sigma, dict_mumumu_sb,  dict_mu_qmu = dict_mu_qmu, dict_mu_calib = dict_mu_calib, doTilde = doTilde, nsmall = nsmall0)
    CLs_new_obs = 1
    if p0_new_obs == 0:
        print('!!!p0_new_obs =', p0_new_obs)
    else:
        CLs_new_obs = pmu_new_obs / p0_new_obs


    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>')
    print('mu =', mu)
    print('toy exp:', CLs_toy_exp, pmu_toy_exp, p0_toy_exp)
    print('toy obs:', CLs_toy_obs, pmu_toy_obs, p0_toy_obs)
    print('classic exp:', CLs_classic_exp, pmu_classic_exp, p0_classic_exp)
    print('classic obs:', CLs_classic_obs, pmu_classic_obs, p0_classic_obs)
    print('new exp:', CLs_new_exp, pmu_new_exp, p0_new_exp)
    print('new obs:', CLs_new_obs, pmu_new_obs, p0_new_obs)

    print('CLs_toy_exp =', CLs_toy_exp)
    print('CLs_toy_obs =', CLs_toy_obs)
    print('CLs_classic_exp =', CLs_classic_exp)
    print('CLs_classic_obs =', CLs_classic_obs)
    print('CLs_new_exp =', CLs_new_exp)
    print('CLs_new_obs =', CLs_new_obs)

    list_CLs = [
            CLs_toy_exp,
            CLs_toy_obs,
            CLs_classic_exp,
            CLs_classic_obs,
            CLs_new_exp,
            CLs_new_obs,
            ]
    if fast:
        return list_CLs

    h_classic_mu = TH1F('h_classic_mu', '', nbins, xmin, xmax)
    h_classic_0 = TH1F('h_classic_0', '', nbins, xmin, xmax)
    h_new_mu = TH1F('h_new_mu', '', nbins, xmin, xmax)
    h_new_0 = TH1F('h_new_0', '', nbins, xmin, xmax)
    hr_classic_mu = TH1F('hr_classic_mu', '', nbins, xmin, xmax)
    hr_classic_0 = TH1F('hr_classic_0', '', nbins, xmin, xmax)
    hr_new_mu = TH1F('hr_new_mu', '', nbins, xmin, xmax)
    hr_new_0 = TH1F('hr_new_0', '', nbins, xmin, xmax)
    p_mu_pre = 0
    p_0_pre = 0

    for i in range(nbins):
        print('--------------------------------------------working on bin', i)
        x1 = h_classic_mu.GetBinLowEdge(i+1)
        x2 = h_classic_mu.GetBinLowEdge(i+2)
        # classic formulae
        p_mu = Fclassic(x2, mu, mu, sigma, doTilde = doTilde) - Fclassic(x1, mu, mu, sigma, doTilde = doTilde)
        p_0 = Fclassic(x2, mu, 0, sigma, doTilde = doTilde) - Fclassic(x1, mu, 0, sigma, doTilde = doTilde)
        if i==0 or x1==0:
            p_mu = Fclassic(x2, mu, mu, sigma, doTilde = doTilde)
            p_0 = Fclassic(x2, mu, 0, sigma, doTilde = doTilde)
        h_classic_mu.SetBinContent(i+1, p_mu)
        h_classic_0.SetBinContent(i+1, p_0)
        # debug try 2bin
        p_mu = Fnew(x2, mu, mu, list_sigma, dict_mumumu_sb,  dict_mu_qmu = dict_mu_qmu, dict_mu_calib = dict_mu_calib, doTilde = doTilde, nsmall = nsmall)
        p_0 = Fnew(x2, mu, 0, list_sigma, dict_mumumu_sb,  dict_mu_qmu = dict_mu_qmu, dict_mu_calib = dict_mu_calib, doTilde = doTilde, nsmall = nsmall0)
        h_new_mu.SetBinContent(i+1, p_mu - p_mu_pre)
        h_new_0.SetBinContent(i+1, p_0 - p_0_pre)
        p_mu_pre = p_mu
        p_0_pre = p_0

    fsc_pmu = 5.
    ymax_hr = 0.
    for i in range(nbins):
        p_toy_mu = h_toy_mu.Integral(i+1, nbins+1)
        hr_classic_mu.SetBinContent(i+1, (-p_toy_mu + h_classic_mu.Integral(i+1, nbins+1))*fsc_pmu)
        hr_new_mu.SetBinContent(i+1, (-p_toy_mu + h_new_mu.Integral(i+1, nbins+1))*fsc_pmu)
        p_toy_0 = h_toy_0.Integral(i+1, nbins+1)
        hr_classic_0.SetBinContent(i+1, -p_toy_0 + h_classic_0.Integral(i+1, nbins+1))
        hr_new_0.SetBinContent(i+1, -p_toy_0 + h_new_0.Integral(i+1, nbins+1))

        ymax_hr = max(ymax_hr, abs(hr_classic_mu.GetBinContent(i+1)), abs(hr_new_mu.GetBinContent(i+1)), abs(hr_classic_0.GetBinContent(i+1)), abs(hr_new_0.  GetBinContent(i+1)))

    print('h_classic_mu integral =', h_classic_mu.Integral())
    print('h_classic_0 integral =', h_classic_0.Integral())
    print('h_new_mu integral =', h_new_mu.Integral() ) #,
    print('h_new_0 integral =', h_new_0.Integral() ) #

    Cs_rqmu = TCanvas('Cs_rqmu', '', 800, 600)
    hr_classic_0.Draw('hist')
    hr_classic_0.SetLineColor(kRed-7)
    hr_classic_0.SetLineStyle(2)
    hr_classic_0.GetYaxis().SetRangeUser(-(int(ymax_hr*10)+1)*0.1, (int(ymax_hr*10)+1)*0.1)
    hr_classic_0.GetYaxis().SetTitle('Tail Probability Difference')
    if doTilde:
        hr_classic_0.GetXaxis().SetTitle('#tilde{q}_{#mu} (%s, #mu=%.3f)' % (mass, mu))
    else:
        hr_classic_0.GetXaxis().SetTitle('q_{#mu} (%s, #mu=%.3f)' % (mass, mu))
    #new
    hr_classic_mu.Draw('hist,same')
    hr_classic_mu.SetLineColor(kBlue-7)
    hr_classic_mu.SetLineStyle(2)
    hr_new_0.Draw('hist,same')
    hr_new_0.SetLineColor(kRed)
    hr_new_mu.Draw('hist,same')
    hr_new_mu.SetLineColor(kBlue)
    Cs_rqmu.Update()
    l1 = TLine(Cs_rqmu.GetUxmin(), 0, Cs_rqmu.GetUxmax(), 0)
    l1.Draw()
    plotname = 'Cs_rqmu_%s_mu%.3f_logy%i' % (mass, mu, logy)
    Cs_rqmu.SaveAs(plotname+'.png')
    Cs_rqmu.SaveAs(plotname+'.pdf')


    Cs_qmu = TCanvas('Cs_qmu', '', 800, 600)
    if logy:
        Cs_qmu.SetLogy()
    h_toy_mu.Draw('PE')
    h_toy_mu.SetMarkerStyle(20)
    #h_toy_mu.SetLineColor(kBlue)
    if logy:
        h_toy_mu.GetYaxis().SetRangeUser(yaxis_low, 1000)
    else:
        h_toy_mu.GetYaxis().SetRangeUser(0, h_toy_mu.GetMaximum()*1.2)
    h_toy_mu.GetYaxis().SetTitle('Probability')
    if doTilde:
        h_toy_mu.GetXaxis().SetTitle('#tilde{q}_{#mu} (%s, #mu=%.3f)' % (mass, mu))
    else:
        h_toy_mu.GetXaxis().SetTitle('q_{#mu} (%s, #mu=%.3f)' % (mass, mu))
    h_toy_0.Draw('PE, same')
    h_toy_mu.Draw('PE,same')
    h_toy_0.SetMarkerStyle(24)
    #h_toy_0.SetLineColor(kRed)
    #classic
    h_classic_0.Draw('hist,same')
    h_classic_0.SetLineColor(kRed-7)
    h_classic_0.SetLineStyle(2)
    h_classic_mu.Draw('hist,same')
    h_classic_mu.SetLineColor(kBlue-7)
    h_classic_mu.SetLineStyle(2)
    #new
    h_new_0.Draw('hist,same')
    h_new_0.SetLineColor(kRed)
    h_new_mu.Draw('hist,same')
    h_new_mu.SetLineColor(kBlue)

    Cs_qmu.Update()
    Uymin = Cs_qmu.GetUymin()
    Uymax = Cs_qmu.GetUymax()
    print('Uymin, Uymax =', Uymin, Uymax)
    a_exp = TArrow(qmu_exp, Uymin+0.8*(Uymax-Uymin), qmu_exp, Uymin, 0.017, '|>')
    a_exp.Draw()
    a_exp.SetFillColor(kGray)
    a_exp.SetLineColor(kGray)
    a_exp.SetLineWidth(2)

    a_obs = TArrow(qmu_obs, Uymin+0.8*(Uymax-Uymin), qmu_obs, Uymin, 0.017, '|>')
    a_obs.Draw()
    a_obs.SetLineWidth(2)
    ltx = TLatex()
    if doTilde:
        ltx.DrawLatex(qmu_exp, turn2log(Uymin+0.1*(Uymax-Uymin), logy), '#tilde{q}_{#mu}^{exp}')
        ltx.DrawLatex(qmu_obs, turn2log(Uymin+0.1*(Uymax-Uymin), logy), '#tilde{q}_{#mu}^{obs}')
    else:
        ltx.DrawLatex(qmu_exp, turn2log(Uymin+0.1*(Uymax-Uymin), logy), 'q_{#mu}^{exp}')
        ltx.DrawLatex(qmu_obs, turn2log(Uymin+0.1*(Uymax-Uymin), logy), 'q_{#mu}^{obs}')
    ltx.DrawLatex(xmin + 0.05*(xmax-xmin),  turn2log(Uymin+0.9*(Uymax-Uymin), logy), 'CLs exp.(toy)=%.3f#pm%.3f' % (CLs_toy_exp, dCLs_toy_exp))
    ltx.DrawLatex(xmin + 0.05*(xmax-xmin),  turn2log(Uymin+0.85*(Uymax-Uymin), logy), 'CLs exp.(classic)=%.3f' % CLs_classic_exp)
    ltx.DrawLatex(xmin + 0.05*(xmax-xmin),  turn2log(Uymin+0.80*(Uymax-Uymin), logy), 'CLs exp.(new)=%.3f' % CLs_new_exp)
    #ltx.DrawLatex(xmin + 0.05*(xmax-xmin),  turn2log(Uymin+0.75*(Uymax-Uymin), logy), 'CLs obs.(toy)=%.3f' % CLs_toy_obs)
    ltx.DrawLatex(xmin + 0.05*(xmax-xmin),  turn2log(Uymin+0.75*(Uymax-Uymin), logy), 'CLs obs.(toy)=%.3f#pm%.3f' % (CLs_toy_obs, dCLs_toy_obs))
    ltx.DrawLatex(xmin + 0.05*(xmax-xmin),  turn2log(Uymin+0.7*(Uymax-Uymin), logy), 'CLs obs.(classic)=%.3f' % CLs_classic_obs)
    ltx.DrawLatex(xmin + 0.05*(xmax-xmin),  turn2log(Uymin+0.65*(Uymax-Uymin), logy), 'CLs obs.(new)=%.3f' % CLs_new_obs)

    leg = TLegend(0.65, 0.5, 0.95, 0.92)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.Draw()
    leg.AddEntry(h_toy_0, 'toy(#mu_{H}=0)', 'PE')
    leg.AddEntry(h_toy_mu, 'toy(#mu_{H}=#mu)', 'PE')
    leg.AddEntry(h_classic_0, 'classic(#mu_{H}=0)', 'L')
    leg.AddEntry(h_classic_mu, 'classic(#mu_{H}=#mu)', 'L')
    leg.AddEntry(h_new_0, 'new(#mu_{H}=0)', 'L')
    leg.AddEntry(h_new_mu, 'new(#mu_{H}=#mu)', 'L')
    plotname = 'Cs_qmu_%s_mu%.3f_logy%i' % (mass, mu, logy)
    Cs_qmu.SaveAs(plotname+'.png')
    Cs_qmu.SaveAs(plotname+'.pdf')
    #Cs_qmu.SaveAs('~/public/'+plotname+'.png')
    #Cs_qmu.SaveAs('~/public/'+plotname+'.pdf')
    #Cs_qmu.SaveAs(plotname+'.root')
    return list_CLs


def get_ul(list_mu = [], list_cls = [], target=0.05):
    n = len(list_mu)
    mu = 0.
    for i in range(n-1):
        mu0 = list_mu[i]
        mu1 = list_mu[i+1]
        cls0 = list_cls[i]
        cls1 = list_cls[i+1]
        if (cls0-target)*(cls1-target)<=0:
            mu =  mu0 + (mu1-mu0)/(cls1-cls0)*(target-cls0)
            pass
        pass
    return mu

def plot_CLs(list_mu_CLs, mass=900):
    if len(list_mu_CLs) == 0:
        return
    n = len(list_mu_CLs)
    list_g = []
    a_mu = array('f', [])
    a_CLs_toy_exp = array('f', [])
    a_CLs_toy_obs = array('f', [])
    a_CLs_classic_exp = array('f', [])
    a_CLs_classic_obs = array('f', [])
    a_CLs_new_exp = array('f', [])
    a_CLs_new_obs = array('f', [])

    list_a_Cls = []
    for i in range(n):
        mu = list_mu_CLs[i][0]
        a_mu.append(mu)
        list_CLs = list_mu_CLs[i][1]
        a_CLs_toy_exp.append(list_CLs[0])
        a_CLs_toy_obs.append(list_CLs[1])
        a_CLs_classic_exp.append(list_CLs[2])
        a_CLs_classic_obs.append(list_CLs[3])
        a_CLs_new_exp.append(list_CLs[4])
        a_CLs_new_obs.append(list_CLs[5])
        pass

    ul_toy_exp = get_ul(a_mu, a_CLs_toy_exp)
    ul_classic_exp = get_ul(a_mu, a_CLs_classic_exp)
    ul_new_exp = get_ul(a_mu, a_CLs_new_exp)
    ul_toy_obs = get_ul(a_mu, a_CLs_toy_obs)
    ul_classic_obs = get_ul(a_mu, a_CLs_classic_obs)
    ul_new_obs = get_ul(a_mu, a_CLs_new_obs)

    print('ul_toy_exp =', ul_toy_exp)
    print('ul_toy_obs =', ul_toy_obs)


    g_CLs_toy_exp = TGraph(n, a_mu, a_CLs_toy_exp)
    g_CLs_toy_obs = TGraph(n, a_mu, a_CLs_toy_obs)
    g_CLs_classic_exp = TGraph(n, a_mu, a_CLs_classic_exp)
    g_CLs_classic_obs = TGraph(n, a_mu, a_CLs_classic_obs)
    g_CLs_new_exp = TGraph(n, a_mu, a_CLs_new_exp)
    g_CLs_new_obs = TGraph(n, a_mu, a_CLs_new_obs)

    mu_min = a_mu[0] - 0.1*(a_mu[-1] - a_mu[0])
    mu_max = a_mu[-1] + 0.1*(a_mu[-1] - a_mu[0])

    CLs_max = max(
            a_CLs_toy_exp[0],
            a_CLs_classic_exp[0],
            a_CLs_new_exp[0],
            a_CLs_toy_obs[0],
            a_CLs_classic_obs[0],
            a_CLs_new_obs[0],
            )*1.5
    Cs_CLs = TCanvas('Cs_CLs', '', 800, 600)
    h2 = TH2F('h2', '', 100, mu_min, mu_max, 100, 0, CLs_max)
    h2.Draw()
    h2.GetXaxis().SetTitle("#mu (%s)" % mass )
    h2.GetYaxis().SetTitle("CLs")

    g_CLs_toy_exp.Draw('LPsame')
    g_CLs_toy_exp.SetLineStyle(2)
    g_CLs_toy_exp.SetLineWidth(2)
    g_CLs_toy_exp.SetMarkerStyle(24)
    g_CLs_classic_exp.Draw('Lsame')
    g_CLs_classic_exp.SetLineStyle(2)
    g_CLs_classic_exp.SetLineWidth(2)
    g_CLs_classic_exp.SetMarkerStyle(24)
    g_CLs_classic_exp.SetLineColor(kRed)
    g_CLs_classic_exp.SetMarkerColor(kRed)

    g_CLs_new_exp.Draw('Lsame')
    g_CLs_new_exp.SetLineStyle(2)
    g_CLs_new_exp.SetLineWidth(2)
    g_CLs_new_exp.SetMarkerStyle(24)
    g_CLs_new_exp.SetLineColor(kGreen)
    g_CLs_new_exp.SetMarkerColor(kGreen)


    g_CLs_toy_obs.Draw('LPsame')
    g_CLs_toy_obs.SetLineStyle(1)
    g_CLs_toy_obs.SetLineWidth(2)
    g_CLs_classic_obs.Draw('Lsame')
    g_CLs_classic_obs.SetLineStyle(1)
    g_CLs_classic_obs.SetLineWidth(2)
    g_CLs_classic_obs.SetLineColor(kRed)
    g_CLs_classic_obs.SetMarkerColor(kRed)

    g_CLs_new_obs.Draw('Lsame')
    g_CLs_new_obs.SetLineStyle(1)
    g_CLs_new_obs.SetLineWidth(2)
    g_CLs_new_obs.SetLineColor(kGreen)
    g_CLs_new_obs.SetMarkerColor(kGreen)


    leg = TLegend(0.6, 0.55, 0.95, 0.92)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.Draw()
    leg.AddEntry(g_CLs_toy_exp, 'toy exp.', 'LP')
    leg.AddEntry(g_CLs_classic_exp, 'classic exp.', 'L')
    leg.AddEntry(g_CLs_new_exp, 'new exp.', 'L')
    leg.AddEntry(g_CLs_toy_obs, 'toy obs.', 'LP')
    leg.AddEntry(g_CLs_classic_obs, 'classic obs.', 'L')
    leg.AddEntry(g_CLs_new_obs, 'new obs.', 'L')

    Cs_CLs.Update()
    xmin = Cs_CLs.GetUxmin()
    xmax = Cs_CLs.GetUxmax()
    ymin = Cs_CLs.GetUymin()
    ymax = Cs_CLs.GetUymax()

    ltx = TLatex()
    x_ul = 0.55
    y_ul = 0.45
    ltx.DrawLatex(xmin+x_ul*(xmax-xmin), ymin+(y_ul)*(ymax-ymin), 'UL exp. (toy) = %.3f' % ul_toy_exp )
    ltx.DrawLatex(xmin+x_ul*(xmax-xmin), ymin+(y_ul-0.05)*(ymax-ymin), 'UL exp. (classic) = %.3f' % ul_classic_exp )
    ltx.DrawLatex(xmin+x_ul*(xmax-xmin), ymin+(y_ul-0.10)*(ymax-ymin), 'UL exp. (new) = %.3f' % ul_new_exp )
    ltx.DrawLatex(xmin+x_ul*(xmax-xmin), ymin+(y_ul-0.15)*(ymax-ymin), 'UL obs. (toy) = %.3f' % ul_toy_obs )
    ltx.DrawLatex(xmin+x_ul*(xmax-xmin), ymin+(y_ul-0.20)*(ymax-ymin), 'UL obs. (classic) = %.3f' % ul_classic_obs )
    ltx.DrawLatex(xmin+x_ul*(xmax-xmin), ymin+(y_ul-0.25)*(ymax-ymin), 'UL obs. (new) = %.3f' % ul_new_obs )

    line0 = TLine(xmin, 0.05, xmax, 0.05)
    line0.Draw()
    line0.SetLineColor(kBlue)

    plotname = 'Cs_CLs_%s' % (mass)
    Cs_CLs.SaveAs(plotname + '.png')
    Cs_CLs.SaveAs('~/public/' + plotname + '.png')
    Cs_CLs.SaveAs(plotname + '.pdf')
    Cs_CLs.SaveAs('~/public/' + plotname + '.pdf')

    with open('ul_'+str(mass)+'.txt', 'w') as f:
        f.write('%s exp %.3f %.3f %.3f\n' % (mass, ul_toy_exp, ul_classic_exp, ul_new_exp))
        f.write('%s obs %.3f %.3f %.3f\n' % (mass, ul_toy_obs, ul_classic_obs, ul_new_obs))



