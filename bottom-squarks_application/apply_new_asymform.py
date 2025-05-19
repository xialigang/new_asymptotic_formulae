

from ROOT import ROOT, gROOT, TFile, TTree, TH1F, TPad, TCanvas, TLine, TLegend, THStack, TGraph, gPad, TGraphErrors, TColor, TLatex, TArrow, TH2F, TMath
from ROOT import kBlue, kRed, kDashed, kGreen, kGray
import sys
import math
import os
from array import array


#gROOT.LoadMacro("AtlasStyle.C")
#from ROOT import SetAtlasStyle
#SetAtlasStyle()

gROOT.SetBatch(True)


from new_asymform_core import do_one_mu_for_one_mass, plot_CLs

dict_mass_mu = {
        #'ex0': [ 1, 2, 2.5, 3, 3.5, 4, 5],
        #'ex1': [ 2, 3, 3.5, 4, 4.5, 5, 6],
        'ex0': [2,  3, 3.5, 4, 5, 6, 7],
        'ex1': [ 2, 3, 3.5, 4, 5, 6, 7],
        #'ex2': [ 2, 2.5, 3, 3.5, 4, 5, 6 ],
        'ex2': [ 1, 2, 2.5, 3, 3.5, 4, 5],
        #'3500': [0.3, 0.5, 0.7]
        #'3500': [0.3, 0.35, 0.4, 0.5, 0.6],
        #'4000': [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        #'4500': [0.3, 0.4, 0.5, 0.6, 0.7],
        #'3400': [0.3, 0.4, 0.5, 0.6, 0.7],
        #'3000': [0.3, 0.4, 0.5, 0.6, 0.7],
        #'2500': [0.3, 0.4, 0.5, 0.6, 0.7],
        #'1900': [0.6, 0.7, 0.8, 0.9],
        #'3500': [0.3, 0.4, 0.5],
        #'4000': [0.4, 0.6, 0.8, 1.0],
        #'4500': [0.4, 0.5, 0.6, 0.7],
        '251': [0.200, 0.225, 0.250, 0.275, 0.300, 0.325, 0.350, 0.375, 0.400, 0.425, 0.450],
        '550': [0.080, 0.085, 0.090, 0.095, 0.100, 0.105, 0.110 ],
        '500': [0.125, 0.131, 0.137, 0.143, 0.149, 0.155, 0.161, 0.167, 0.173, 0.179, 0.185],
        '700': [0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080, 0.085],
        '900': [0.040, 0.045, 0.050, 0.055, 0.060, 0.065, 0.070, 0.075, 0.080, 0.085, 0.090],
        'sbottom': [0.15, 0.18, 0.20, 0.22, 0.25, 0.30, 0.35],
        }

def do_one_mass(m='ex0', fast=0, do_q0=False):
    mass = m
    nbins, xmax = 48, 12
    if do_q0:
        do_one_mu_for_one_mass(m, 2, nbins, 0, xmax, 1, do_q0 = do_q0)
        return 1
    list_mu_CLs = []
    for mu in dict_mass_mu[m]:

        ###############zy
        nbins, xmax = 40, 10
        if fast == 2:
            fast0 = 0
        else:
            fast0 = fast
        toydir = os.path.join('toys_from_yan', m)
        fname_toy_mu = os.path.join(toydir, 'mu_%i-muprime_%i.root' % (mu*10, mu*10))
        fname_toy_0 = os.path.join(toydir, 'mu_%i-muprime_%i.root' % (mu*10, 0))
        fname_sigmas = 'output_sigmas_zy_%s_v1.txt' % (m)
        fname_interval = 'signal_interval_zy_%s_v4.txt' % (m)
        list_special = []
        if mass == '3400':
            list_special = [0.000893931 , 0.045441]
        elif mass == '3000':
            list_special = [0.00288978 , 0.0486466]
        elif mass == '2500':
            list_special = [0.0102432 , 0.0617949]
        elif mass == '1900':
            list_special = [0.046229 , 0.0575648]
        ################bbbbb
        toydir = os.path.join('toys', 'X'+mass)
        if do_q0:
            fname_toy_mu = os.path.join(toydir, 'HH_bb_gamgam_resonant_X%s_result_toys_mu%s_muprime_0_combined.root' % (mass, 'x'))
            fname_toy_0 = fname_toy_mu
        else:
            fname_toy_mu = os.path.join(toydir, 'HH_bbbb_resonant_X%s_result_toys_mu_%.6f_muprime_%.6f_combined.root' % (mass, mu, mu))
            fname_toy_0 = os.path.join(toydir, 'HH_bbbb_resonant_X%s_result_toys_mu_%.6f_muprime_0_combined.root' % (mass, mu))
        fname_sigmas = 'output_sigmas_bbbb_%s_v14.txt' % (m)
        fname_interval = 'signal_interval_bbbb_%s.txt' % (m)
        list_special = []
        if mass == '3500':
            list_special += [0.0114309 , 0.22346]
        elif mass == '4500':
            list_special += [0.00695206 , 0.273886]
        elif mass == '4000':
            list_special += [0.00814194 , 0.244827]
        nbins, xmax = 40, 20
        if m == '3500':
            nbins, xmax = 60, 15
        ###############bbyy
        nbins, xmax = 40, 20
        if mass == '700':
            nbins, xmax = 30, 15
        elif mass == '251':
            nbins, xmax = 25, 25
        toydir = 'toys'
        fname_toy_mu = os.path.join(toydir, 'HH_bb_gamgam_resonant_X%s_result_toys_mu_%.6f_muprime_%.6f_combined.root' % (mass, mu, mu))
        fname_toy_0 = os.path.join(toydir, 'HH_bb_gamgam_resonant_X%s_result_toys_mu_%.6f_muprime_0_combined.root' % (mass, mu))
        fname_sigmas = 'output_sigmas_%s_v1.txt' % (m)
        fname_interval = 'signal_interval_%s.txt' % (m)
        list_special = []
        if mass == '251':
            list_special = [ 0.0111898, 0.174794 ]
        elif mass == '500':
            list_special = [ 0.00703924, 0.223689 ]
        elif mass == '700':
            list_special = [0.00980635 , 0.147181]
        elif mass == '900':
            list_special = [0.00474105 , 0.205993]
        #######my own toy example
        nbins, xmax = 60, 15
        toydir = os.path.join('toys', 'X'+mass)
        if do_q0:
            fname_toy_mu = os.path.join(toydir, 'HH_bb_gamgam_resonant_X%s_result_toys_mu%s_muprime_0_combined.root' % (mass, 'x'))
            fname_toy_0 = fname_toy_mu
        else:
            fname_toy_mu = os.path.join(toydir, 'HH_bb_gamgam_resonant_X%s_result_toys_mu_%.6f_muprime_%.6f_combined.root' % (mass, mu, mu))
            fname_toy_0 = os.path.join(toydir, 'HH_bb_gamgam_resonant_X%s_result_toys_mu_%.6f_muprime_0_combined.root' % (mass, mu))
        fname_sigmas = 'output_sigmas_%s_v2.txt' %(mass)
        fname_interval= ''
        #fname_interval = 'output_s0s1b0b1_%s_v2.txt' %(mass)
        list_special = []
        if mass == 'ex1':
            list_special = [0.149239 , 0.017662]
        elif mass == 'ex0':
            list_special = [0.14925 , 0.0178638]
    
        #######ATLAS searching for bottom-squarks
        nbins, xmax = 48, 12 
        toydir = os.path.join('toys', 'X'+mass)
        fname_toy_mu = os.path.join(toydir, 'toys_RegionA_mu%.2f_muH%.2f_combined.csv' % (mu, mu))
        fname_toy_0 = os.path.join(toydir, 'toys_RegionA_mu%.2f_muH%.2f_combined.csv' % (mu, 0))
        fname_sigmas = 'my_6bin_model.csv'
        fname_interval= ''
        #fname_interval = 'output_s0s1b0b1_%s_v2.txt' %(mass)
        list_special = []
        if mass == 'sbottom':
            list_special = [0.05228640755407068,0.20992515772040818]
    

        ###################
        list_CLs = do_one_mu_for_one_mass(str(m), mu, nbins, 0, xmax, 1, fast=fast0, fname_toy_mu=fname_toy_mu, fname_toy_0=fname_toy_0, fname_sigmas=fname_sigmas, fname_interval=fname_interval, list_special = list_special, doTilde=1)

        print(mu, list_CLs)
        if list_CLs == 0:
            continue
        list_mu_CLs.append([mu, list_CLs])
        pass
    if len(list_mu_CLs) == len(dict_mass_mu[m]) and len(list_mu_CLs)>1:
        plot_CLs(list_mu_CLs, m)
    pass

fast=0
mass='sbottom'

if len(sys.argv)>1:
    mass = str(sys.argv[1])
if len(sys.argv)>2:
    fast = int(sys.argv[2])
#do_one_mass('ex2', fast)
#do_one_mass('ex1', fast)
#do_one_mass('ex0', fast)

do_one_mass(mass, fast)


