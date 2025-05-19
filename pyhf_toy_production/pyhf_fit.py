import time
import os
import ROOT
import sys
import json
import cabinetry
import pyhf
from cabinetry.model_utils import prediction
from pyhf.contrib.utils import download
import numpy as np
import matplotlib.pyplot as plt
import math

#ROOT.gROOT.LoadMacro("AtlasStyle.C")
#from ROOT import SetAtlasStyle
#SetAtlasStyle()

ROOT.gROOT.SetBatch(True)

#plt.rcParams.update()
if not os.path.isdir('bottom-squarks') and 0:
    download('https://www.hepdata.net/record/resource/1935437?view=true','bottom-squarks')

time0 = time.ctime()

action = ''
test_poi = 0.19
toy_mode = 1
init_mode = 0
n_toys = 100
rndseed = 0

if len(sys.argv)>1:
    action = sys.argv[1]
if len(sys.argv)>2:
    test_poi = float(sys.argv[2])
if action == 'toys':
    if len(sys.argv)>3:
        toy_mode = int(sys.argv[3])
    if len(sys.argv)>4:
        n_toys = int(sys.argv[4])
    if len(sys.argv)>5:
        rndseed = int(sys.argv[5])
if action == 'mymodel':
    if len(sys.argv)>3:
        init_mode = int(sys.argv[3])


np.random.seed(rndseed)
SR = 'RegionA'
nbins = 6
mumax = 1
#SR = 'RegionB'
#nbins = 2

bkg_only_workspace = pyhf.Workspace(json.load(open('bottom-squarks/'+SR+'/BkgOnly.json')))
patchset = pyhf.PatchSet(json.load(open('bottom-squarks/'+SR+'/patchset.json')))
workspace = patchset.apply(bkg_only_workspace, 'sbottom_600_280_150')

model, data = cabinetry.model_utils.model_and_data(workspace)

par_names = model.config.par_names
init_params = model.config.suggested_init()
fixed_params = model.config.suggested_fixed()
poi_index = model.config.poi_index
poi_name = model.config.poi_name
unbounded_bounds = model.config.suggested_bounds()
unbounded_bounds[poi_index] = (-mumax, mumax)
print('par_names =', par_names, len(par_names))
print('poi_index, poi_name, bound =', poi_index, poi_name, unbounded_bounds[poi_index])
print('init_params =', init_params)
print('fixed_params =', fixed_params)
print('unbounded_bounds =', unbounded_bounds)

#prefit_model = prediction(model)
#cabinetry.visualize.data_mc(prefit_model, data)
#pre_params = model.config.suggested_init()
#print('pre_params =', pre_params)
# bkg-only fit
init_poi = 0
if init_mode == 1:
    init_poi = test_poi
init_params[poi_index] = init_poi
fixed_params[poi_index] = True
fit_results_mu0 = cabinetry.fit.fit(model, data, init_pars=init_params, fix_pars=fixed_params)
fixed_params[poi_index] = False
print('fit_results_mu0 =', fit_results_mu0)
postfit_model_mu0 = prediction(model, fit_results=fit_results_mu0)
print('postfit_model_mu0 =', postfit_model_mu0)
cabinetry.visualize.data_mc(postfit_model_mu0, data)
#model_yield_mu0 = postfit_model_mu0.model_yields
#print('model_yield_mu0 =', model_yield_mu0)
#n_reg = len(model_yield_mu0)
#n_sam = len(model_yield_mu0[0])
#n_bins = len(model_yield_mu0[0][0])
#print('n_reg, n_sam, n_bins =', n_reg, n_sam, n_bins)

def getZ(b=1,s=1):
    Z = math.sqrt(2*((s+b)*math.log(1.+s/b)-s))
    return Z



# unconditional fit
bestfit_pars_obs, twice_nll_obs = pyhf.infer.mle.fit(data, model, return_fitted_val=True)
print('twice_nll_obs =', twice_nll_obs)
npars = len(par_names)
for i in range(npars):
    print(i,par_names[i], bestfit_pars_obs[i])
poi_obs = bestfit_pars_obs[poi_index]
print('poi_obs =', poi_obs)
# conditional fit with mu=test_poi
bestfit_pars_mu1, twice_nll_mu1 = pyhf.infer.mle.fixed_poi_fit(test_poi, data, model, return_fitted_val=True)
print('twice_nll_mu1 =', twice_nll_mu1)
# conditional (bkg-only) fit with mu=0
bestfit_pars_mu0, twice_nll_mu0 = pyhf.infer.mle.fixed_poi_fit(0., data, model, return_fitted_val=True)
for i in range(npars):
    print(i,par_names[i], bestfit_pars_mu0[i])
print('twice_nll_mu0 =', twice_nll_mu0)
qmu_tilde_obs = twice_nll_mu1 - twice_nll_obs
if poi_obs<0:
    qmu_tilde_obs = twice_nll_mu1 - twice_nll_mu0

if action == 'mymodel':
    asimov_data_mu0 = cabinetry.model_utils.asimov_data(model, fit_results=fit_results_mu0, poi_name=poi_name, poi_value=0)
    asimov_data_mu1 = cabinetry.model_utils.asimov_data(model, fit_results=fit_results_mu0, poi_name=poi_name, poi_value=test_poi)
    print('data =', data, len(data))
    print('asimov_data_mu0 =', asimov_data_mu0, len(asimov_data_mu0))
    print('asimov_data_mu1 =', asimov_data_mu1, len(asimov_data_mu1))

    # prepare the 6-bin model
    my_6bin_model = []
    for ibin in range(nbins):
        nobs = data[ibin]
        b = asimov_data_mu0[ibin]
        s = asimov_data_mu1[ibin]-asimov_data_mu0[ibin]
        Z = getZ(b,s)
        print('bin, nobs, b, s, Z =', ibin, nobs, b, s, Z)
        my_6bin_model.append({'bin':[ibin, nobs, b, s], 'Z':Z})
        pass
    # re-order it with decreasing significance
    my_6bin_model_sorted = sorted(my_6bin_model, key=lambda x: x['Z'], reverse=True)
    print('my_6bin_model_sorted =', my_6bin_model_sorted)
    # to obtain qmu_tilde_exp, we will fit to asimov_data_0
    bestfit_pars_mu1_asimov0, twice_nll_mu1 = pyhf.infer.mle.fixed_poi_fit(test_poi, asimov_data_mu0, model, return_fitted_val=True)
    bestfit_pars_mu0_asimov0, twice_nll_mu0 = pyhf.infer.mle.fixed_poi_fit(0, asimov_data_mu0, model, return_fitted_val=True)
    print('cross check bestfit_pars_mu0_asimov0 =', bestfit_pars_mu0_asimov0)
    qmu_tilde_exp = twice_nll_mu1 - twice_nll_mu0
    sigma_wald = abs(test_poi-0)/math.sqrt(qmu_tilde_exp)
    # to obtain sigma(mu) coefficients
    # floating parameters: poi, mu_ttbar
    init_params = model.config.suggested_init()
    print('check1, init_params =', init_params, len(init_params))
    #for i in range(len(par_names)):
    #    if par_names[i] in [poi_name, 'mu_ttbar']:
    #        init_params[i] = 0.1 #make it far from best-fit to get a reasonable fit
    #init_params = fit_results_mu0.bestfit 
    #print('check2, init_params =', init_params, len(init_params))
    #print('check3, init_params =', init_params, len(init_params))
    init_params[poi_index] = 0.001
    fixed_params = model.config.suggested_fixed()
    _fit_results = cabinetry.fit.fit(model, asimov_data_mu0, init_pars=init_params, fix_pars=fixed_params, par_bounds=unbounded_bounds)
    mu_unc_asimov0 = _fit_results.uncertainty[poi_index]
    for i in range(len(par_names)): 
        if par_names[i] not in [poi_name, 'mu_ttbar']:
            fixed_params[i] = True
        else:
            fixed_params[i] = False
    init_params = []
    for a in _fit_results.bestfit:
        init_params.append(float(a))
    init_params[poi_index] = 0.001
    _fit_results = cabinetry.fit.fit(model, asimov_data_mu0, init_pars=init_params, fix_pars=fixed_params, par_bounds=unbounded_bounds)
    mu_unc_asimov0_nosyst = _fit_results.uncertainty[poi_index]
    init_params = model.config.suggested_init()
    fixed_params = model.config.suggested_fixed()
    _fit_results = cabinetry.fit.fit(model, asimov_data_mu1, init_pars=init_params, fix_pars=fixed_params)
    mu_unc_asimov1 = _fit_results.uncertainty[poi_index]
    for i in range(len(par_names)): 
        if par_names[i] not in [poi_name, 'mu_ttbar']:
            fixed_params[i] = True
        else:
            fixed_params[i] = False
    init_params = []
    for a in _fit_results.bestfit:
        init_params.append(float(a))
    _fit_results = cabinetry.fit.fit(model, asimov_data_mu1, init_pars=init_params, fix_pars=fixed_params)
    mu_unc_asimov1_nosyst = _fit_results.uncertainty[poi_index]
    sigma0 = 0
    if mu_unc_asimov0 > mu_unc_asimov0_nosyst:
        sigma0 = math.sqrt(pow(mu_unc_asimov0,2)-pow(mu_unc_asimov0_nosyst,2))
    sigma1 = 0
    if mu_unc_asimov1 > mu_unc_asimov1_nosyst:
        sigma1 = math.sqrt(pow(mu_unc_asimov1,2)-pow(mu_unc_asimov1_nosyst,2)) 
    kappa = 0
    if sigma1 > sigma0:
        kappa = math.sqrt(pow(sigma1,2)-pow(sigma0,2))/test_poi 
    print('mu_unc_asimov0, mu_unc_asimov0_nosyst =', mu_unc_asimov0, mu_unc_asimov0_nosyst)
    print('mu_unc_asimov1, mu_unc_asimov1_nosyst =', mu_unc_asimov1, mu_unc_asimov1_nosyst)
    print('sigma0, sigma1, kappa =', sigma0, sigma1, kappa)
    filename = 'my_%s_6bin_model_mu%.2f_muH%.2f_seed%i.csv' % (SR, test_poi, init_poi, rndseed)
    with open(filename, 'w') as f:
        for abin in my_6bin_model_sorted:
            f.write('nominal %.2f %.2f %i %.6f %.6f\n' % (test_poi, init_poi, abin['bin'][1], abin['bin'][2], abin['bin'][3]))
        f.write('%.2f %e %.6f %.6f %.6f\n' % (test_poi, poi_obs, sigma_wald, qmu_tilde_exp, qmu_tilde_obs))


def tmu_asymp(mutest, muhat, sigma):
    return (mutest - muhat) ** 2 / sigma**2


def tmu_tilde_asymp(mutest, muhat, sigma):
    a = tmu_asymp(mutest, muhat, sigma)
    b = tmu_asymp(mutest, muhat, sigma) - tmu_asymp(0.0, muhat, sigma)
    return np.where(muhat > 0, a, b)

def qmu_asymp(mutest, muhat, sigma):
    return np.where(
        muhat < mutest, tmu_asymp(mutest, muhat, sigma), np.zeros_like(muhat)
    )


def qmu_tilde_asymp(mutest, muhat, sigma):
    return np.where(
        muhat < mutest, tmu_tilde_asymp(mutest, muhat, sigma), np.zeros_like(muhat)
    )

if action == 'toys':
    unbounded_bounds = model.config.suggested_bounds()
    unbounded_bounds[model.config.poi_index] = (-mumax, mumax)
    bounded_bounds = model.config.suggested_bounds()
    fixed_params = model.config.suggested_fixed()
    true_poi = test_poi
    if toy_mode == 0:
        true_poi = 0

    # according to the ATLAS convention, we use bestfit_pars_mu0 to produce toys
    pars_init = bestfit_pars_mu0
    print('pars_init =', pars_init)
    pars_init[poi_index] = true_poi
    print('pars_init =', pars_init)
    print('bestfit_pars_mu0 =', bestfit_pars_mu0)
    toys = model.make_pdf(pyhf.tensorlib.astensor(pars_init)).sample((n_toys,))
    pars = np.asarray([pyhf.infer.mle.fit(toy, model, par_bounds=unbounded_bounds) for toy in toys])
    muhat = pars[:, model.config.poi_index]
    print('muhat =', muhat)
    qmu_tilde = np.asarray(
            [
                pyhf.infer.test_statistics.qmu_tilde(
                    test_poi,
                    toy,
                    model,
                    init_pars=model.config.suggested_init(),
                    par_bounds=bounded_bounds,
                    fixed_params=fixed_params,
                    )
                for toy in toys
                ]
            )

    Cs_muhat = ROOT.TCanvas('Cs_muhat', '', 800, 600)
    h_muhat = ROOT.TH1F('h_muhat', '', 100, -mumax, mumax)
    for a in muhat:
        h_muhat.Fill(a)
        pass
    h_muhat.Scale(1./h_muhat.Integral())
    h_muhat.Draw('hist')
    h_muhat.GetXaxis().SetTitle('#hat{#mu}')
    h_muhat.GetYaxis().SetTitle('Probability')
    Cs_muhat.SaveAs('Cs_sbottom_muhat_mu%.2f_muH%.2f_seed%i.png' % (test_poi, true_poi, rndseed))
    Cs_muhat.SaveAs('Cs_sbottom_muhat_mu%.2f_muH%.2f_seed%i.pdf' % (test_poi, true_poi, rndseed))

    Cs_qmu_tilde = ROOT.TCanvas('Cs_qmu_tilde', '', 800, 600)
    h_qmu_tilde = ROOT.TH1F('h_qmu_tilde', '', 100, 0, 10)
    for a in qmu_tilde:
        h_qmu_tilde.Fill(a)
        pass
    h_qmu_tilde.Scale(1./h_qmu_tilde.Integral())
    h_qmu_tilde.Draw('hist')
    h_qmu_tilde.GetXaxis().SetTitle('#tilde{q}_{#mu}')
    h_qmu_tilde.GetYaxis().SetTitle('Probability')
    Cs_qmu_tilde.SaveAs('Cs_sbottom_qmu_tilde_mu%.2f_muH%.2f_seed%i.png' % (test_poi, true_poi, rndseed))
    Cs_qmu_tilde.SaveAs('Cs_sbottom_qmu_tilde_mu%.2f_muH%.2f_seed%i.pdf' % (test_poi, true_poi, rndseed))

    Cs_check_wald = ROOT.TCanvas('Cs_check_wald', '', 800, 600)
    h_check_wald = ROOT.TH2F('h_check_wald', '', 100, -mumax, mumax, 100, 0, 10)
    for i in range(len(muhat)):
        h_check_wald.Fill(muhat[i], qmu_tilde[i])
    h_check_wald.Scale(1./h_check_wald.Integral())
    h_check_wald.Draw('colz')
    h_check_wald.GetXaxis().SetTitle('#hat{#mu}')
    h_check_wald.GetYaxis().SetTitle('#tilde{q}_{#mu}')
    Cs_check_wald.SaveAs('Cs_sbottom_check_wald_mu%.2f_muH%.2f_seed%i.png' % (test_poi, true_poi, rndseed))
    Cs_check_wald.SaveAs('Cs_sbottom_check_wald_mu%.2f_muH%.2f_seed%i.pdf' % (test_poi, true_poi, rndseed))

    toys_filename = 'toys_%s_mu%.2f_muH%.2f_seed%i.csv' % (SR, test_poi, true_poi, rndseed)
    with open(toys_filename, 'w') as f:
        for i in range(len(muhat)):
            f.write('%.6f %.6f\n' % (muhat[i], qmu_tilde[i]))

    #print('muhat =', muhat, len(muhat))
    #weights = []
    #for i in range(len(muhat)):
    #    weights.append(1./len(muhat))
    #muhat_sigma = np.std(muhat)
    #fig, ax = plt.subplots()
    #fig.set_size_inches(7, 5)

    #ax.set_xlabel(r'$\hat{\mu}$')
    #ax.set_ylabel('Density')
    #ax.set_ylim(top=0.5)

    #ax.hist(muhat, bins=np.linspace(-1, 1, 50), weights = weights)
    #ax.axvline(true_poi, label='true poi', color='black', linestyle='dashed')
    ##ax.axvline(np.mean(muhat), label='empirical mean', color='red', linestyle='dashed')
    #ax.legend()

    #fig.savefig('Cs_muhat.png')
    #fig.savefig('Cs_muhat.pdf')

    #muhat_asymp = np.linspace(-1, 1)
    #fig, ax = plt.subplots()

    #ax.set_xlabel(r'$\hat{\mu}$')
    #ax.set_ylabel(r'$\tilde{q}_{\mu}$')
    #ax.scatter(muhat, qmu_tilde, alpha=0.2, label='toys')
    #ax.plot(
    #        muhat_asymp,
    #        qmu_tilde_asymp(test_poi, muhat_asymp, muhat_sigma),
    #        label='asymptotic',
    #        #label=r'$\\tilde{q}_{\\mu}$ asymptotic',
    #        c='r',
    #        )
    #ax.legend()
    #fig.savefig('Cs_qmutilde_muhat.png')
    #fig.savefig('Cs_qmutilde_muhat.pdf')





if action=='scan':
    for _test_poi in [0.1, 0.15, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.35]:
        CLs_obs, CLs_exp_band = pyhf.infer.hypotest(_test_poi, data, model, return_expected_set=True)
        #print('_test_poi =', _test_poi)
        print('CLs: test_poi, obs, exp =', _test_poi, CLs_obs, CLs_exp_band[2])
        #print('CLs_obs =', CLs_obs)
        #print('CLs_exp_band =', CLs_exp_band)
if action=='limit':
    limitmax = 0.5
    if SR == 'RegionB':
        limitmax = 6
    scan = np.linspace(0., limitmax, 31)
    obs_limit, exp_limits, (scan, results) = pyhf.infer.intervals.upper_limits.upper_limit(data, model, scan, return_results=True, level=0.05)
    print('obs_limit =', obs_limit)
    print('exp_limits =', exp_limits)

time1 = time.ctime()
print('time spent:', time0, '--->', time1)
