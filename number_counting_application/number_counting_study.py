

from ROOT import ROOT, gROOT, TFile, TTree, TH1F, TPad, TCanvas, TLine, TLegend, THStack, TGraph, gPad, TGraphErrors, TColor, TLatex, TArrow, TH2F, TMath, TF1
from ROOT import kBlue, kRed, kDashed, kGreen, kGray, kBlack
from ROOT import gRandom
import sys
import math
import os
from array import array


#gROOT.LoadMacro("AtlasStyle.C")
#from ROOT import SetAtlasStyle
#SetAtlasStyle()

gROOT.SetBatch(True)


def Phi(x):
    return ROOT.Math.normal_cdf(x, 1)

def get_nll2(n,b,s,mu,Delta=0,theta=0, theta_aux=0):
    nll2 = 9999
    if n==0:
        nll2 = -2*(-(b+Delta*theta+mu*s))+pow(theta-theta_aux,2)
    else:
        nll2 = -2*(n*math.log(b+Delta*theta+mu*s)-(b+Delta*theta+mu*s))+pow(theta-theta_aux,2)
    return nll2

def get_mu_theta(n,b,s,mu,Delta=0,theta_aux=0):
    if n==0:
        theta_h = -Delta+theta_aux
        theta_hh = -Delta+theta_aux
        muhat = -(b+Delta*theta_h)/s
        return muhat, theta_h, theta_hh
    muhat = (n-b)/s
    theta_h = theta_aux
    theta_hh = theta_aux
    if Delta > 0:
        theta_hh = 0.5*(-(Delta-theta_aux+(b+mu*s)/Delta)+math.sqrt(pow(Delta-theta_aux+(b+mu*s)/Delta,2)-4*((b+mu*s)/Delta*(Delta-theta_aux)-n)))
    return muhat,theta_h,theta_hh
def get_qmu(n,b,s,mu,Delta=0, doTilde=0, option=0, theta_aux=0):
    muhat, theta_h, theta_hh = get_mu_theta(n,b,s,mu,Delta,theta_aux)
    if muhat>mu and option==0:
        return 0
    elif muhat>=0 or doTilde==0:
        nll2 = get_nll2(n,b,s,muhat,Delta,theta_h, theta_aux)
        nll2_mu = get_nll2(n,b,s,mu,Delta,theta_hh, theta_aux)
        return nll2_mu - nll2
    else:
        nll2_mu = get_nll2(n,b,s,mu,Delta,theta_hh, theta_aux)
        _mu,_th,theta_hh0 = get_mu_theta(n,b,s,0,Delta, theta_aux)
        nll2_0 = get_nll2(n,b,s,0,Delta,theta_hh0, theta_aux)
        if n==0 and 0:
            print('nll2_mu, nll2_0, _mu, _th, theta_hh0 =', nll2_mu, nll2_0, _mu, _th, theta_hh0, nll2_mu - nll2_0)
        return nll2_mu - nll2_0

def CDFasym(qmu, mu, muH, sigma, sigmap=0, doTilde=0):
    if sigmap<=0:
        sigmap = sigma
    if qmu <= pow(mu/sigma,2) or doTilde==False:
        return Phi(math.sqrt(qmu)*sigma/sigmap-(mu-muH)/sigmap)
    else:
        return Phi(((qmu-mu*mu/sigma/sigma)/(2*mu/sigma/sigma)+muH)/sigmap)
        #return Phi((qmu-pow(mu/sigma,2)+2*mu*muH/pow(sigma,2))/(2.*mu*sigma/pow(sigma,2)))

def newCDFasym(qmu, mu, muH, b, s, Delta=0, doTilde=0, big_nsmall=0):
    # expected qmu
    qmu_exp = get_qmu(b,b,s,mu,Delta)
    # sigma from Wald's theorem
    sigma = mu/math.sqrt(qmu_exp)
    # R factor
    _qmu_exp  = get_qmu(b,b,s,mu,Delta=0)
    rqmu = qmu_exp / _qmu_exp
    #print('rqmu =', rqmu)
    #print('b,s,mu,muH,qmu_exp,sigma =', b,s,mu, muH, qmu_exp, sigma)
    nsmall = int(b+mu*s)-1
    if nsmall < 3:
        nsmall = 3
    if b>10 and nsmall>10:
        nsmall = 10
    if big_nsmall == 1:
        nsmall = int(b+(mu+5*sigma)*s)
    #print('nsmall in newCDFasym =', nsmall)
    Prob_small = 0.
    Prob_big = 0.
    Poisson_small = 0.
    nu = b+muH*s
    for k in range(0,nsmall+1):
        Poisson_small += TMath.Poisson(k,b+muH*s)
    Poisson_big = 1- Poisson_small
    muH_corr = 0.
    for k in range(0,nsmall+1):
        Poisson_k = TMath.Poisson(k,b+muH*s)
        _qmu = get_qmu(k, b, s, mu, option=1, doTilde=doTilde)
        _qmu *= rqmu
        mup = float((k-b)/s)
        muH_corr += mup * Poisson_k
        if _qmu == 0 or mup==mu:
            #print('mu, mup, _qmu =', mu, mup, _qmu)
            sigma1 = sigma
        elif mup>=0:
            sigma1 = abs(mup-mu)/math.sqrt(_qmu)
        else:
            sigma1 = math.sqrt((-2*mu*mup+mu*mu)/_qmu)
        sigmap = float(Delta/s)
        #print('k, mup, _qmu, sigma1, sigmap, rqmu =', k, mup, _qmu, sigma1, sigmap, rqmu)
        Prob_small += Poisson_k*CDFasym(qmu, mu, mup, sigma1, sigmap, doTilde=doTilde)
        pass
    muH_corr = (muH - muH_corr)/Poisson_big
    Prob_big = Poisson_big * CDFasym(qmu, mu, muH_corr, sigma, doTilde=1)
    Prob = Prob_big + Prob_small
    #Prob = Prob_small
    return Prob
def cdf_mu(muobs, mu, muH, b, s, Delta=0,doTilde=0, big_nsmall=0):
    qmu_exp = get_qmu(b,b,s,mu, Delta)
    sigma = mu/math.sqrt(qmu_exp)
    #print('b,s,mu,muH,qmu_exp,sigma =', b,s,mu, muH, qmu_exp, sigma)
    nsmall = int(b+muH*s)-1
    if nsmall < 3:
        nsmall = 3
    if b>10 and nsmall > 10:
        nsmall = 10
    if big_nsmall == 1:
        nsmall = int(b+(mu+5*sigma)*s)
    #print('in cdf_mu, nsmall =', nsmall)
    Prob_small = 0.
    Prob_big = 0.
    Poisson_small = 0.
    for k in range(0,nsmall+1):
        Poisson_small += TMath.Poisson(k,b+muH*s)
    Poisson_big = 1- Poisson_small
    muH_corr = 0.
    for k in range(0,nsmall+1):
        Poisson_k = TMath.Poisson(k,b+muH*s)
        mup = (k-b)/s
        muH_corr += mup * Poisson_k
        sigmap = float(Delta/s)
        Prob_small += Poisson_k*Phi((muobs -mup)/sigmap)
        pass
    muH_corr = (muH - muH_corr)/Poisson_big
    Prob_big = Poisson_big * Phi((muobs-muH_corr)/sigma)
    Prob = Prob_big + Prob_small
    #Prob = Prob_small
    return Prob


def CDFwald(tmu, mu, muH, sigma0, sigma1):
    mup = mu+math.sqrt(tmu)*sigma0
    mum = mu-math.sqrt(tmu)*sigma0
    return Phi((mup-muH)/sigma1) - Phi((mum-muH)/sigma1)

def get_sigma(n,b,s,Delta=0):
    sigma = math.sqrt(n+Delta*Delta)/s
    #sigma = math.sqrt(b+mu*s)/s
    return sigma
def asymp_wald(X,P):
    muhat = X[0]
    mu = P[0]
    sigma = P[1]
    if muhat>mu:
        return 0
    elif muhat>0 and muhat<=mu:
        return pow((mu-muhat)/sigma,2)
    else:
        return (mu*mu-2*mu*muhat)/sigma/sigma

def get_toys(b=10,s=1,mu=1,muH=1,ntoys=1000, nbins=100, xmin=0, xmax=10, big_nsmall=0):
    # suppose a background uncertainty of 10%
    Delta = 0.1*b
    if b>10:
        Delta = 0.05*b
    # sigma from fit for muH
    sigma_muH = get_sigma(b+muH*s,b,s,Delta)
    # sigma from Wald's theorem
    qmu_exp = get_qmu(b,b,s,mu,Delta)
    sigma = mu/math.sqrt(qmu_exp)
    # define histograms
    h = TH1F('h_qmu', '', nbins, xmin, xmax)
    h.Sumw2()
    h_muhat = TH1F('h_muhat', '', nbins, muH-5*sigma_muH,muH+5*sigma_muH)
    h_muhat_asymp = TH1F('h_muhat_asymp', '', nbins, muH-5*sigma_muH,muH+5*sigma_muH)
    h_muhat_new = TH1F('h_muhat_new', '', nbins, muH-5*sigma_muH,muH+5*sigma_muH)
    h_qmu_muhat = TH2F('h_qmu_muhat', '', 100, min(-0.5*mu,muH-3*sigma_muH), max(1.1*mu,muH+3*sigma_muH), 100, 0, xmax)

    mu_mean = 0.
    mu2_mean = 0.
    thetaH = 0. #nominal theta values (usually 0 before unblinding and fitted to observed data after unblinding)
    # now produce toys
    for itoy in range(ntoys):
        # produce auxiliary measurements for systematic variations
        theta_aux = gRandom.Gaus(0,1)+thetaH
        # produce toys for observables, here it is just the number of events
        n = gRandom.Poisson(b+thetaH*Delta+muH*s)
        # get best-fit results
        muhat, theta_h, theta_hh = get_mu_theta(n,b,s,mu,Delta, theta_aux)
        qmu = get_qmu(n,b,s,mu,Delta,doTilde=1,theta_aux=theta_aux)
        # fill histograms
        h_muhat.Fill(muhat)
        h_qmu_muhat.Fill(muhat, qmu)
        h.Fill(qmu)
        # n=0 is special and printed for monitoring
        if n==0:
            print('n,theta_aux, muhat, qmu =', n,theta_aux, muhat, qmu)
        # just want to check if the mean is still 0 or not
        mu_mean += muhat
        mu2_mean += muhat*muhat
        pass
    mu_mean /= ntoys
    mu2_mean /= ntoys
    sigma_mean = math.sqrt(mu2_mean - mu_mean*mu_mean)
    print('mu_mean =', mu_mean, '+/-', sigma_mean)
    h_muhat.Scale(1./h_muhat.Integral())
    h.Scale(1./h.Integral())
    h_af0 = TH1F('h_af0', '', nbins, xmin, xmax)
    h_af1 = TH1F('h_af1', '', nbins, xmin, xmax)
    h_af = TH1F('h_af', '', nbins, xmin, xmax)
    h_afnew = TH1F('h_afnew', '', nbins, xmin, xmax)
    h_ptail= TH1F('h_ptail', '', nbins, xmin, xmax)
    h_ptailnew= TH1F('h_ptailnew', '', nbins, xmin, xmax)
    for ib in range(nbins):
        qmuh = h_af0.GetBinLowEdge(ib+2)
        qmul = h_af0.GetBinLowEdge(ib+1)
        Ph0 = CDFasym(qmuh, mu, muH, sigma, doTilde=1)
        Pl0 = CDFasym(qmul, mu, muH, sigma, doTilde=1)
        #xxxxxxx
        Phnew = newCDFasym(qmuh, mu, muH, b, s, Delta, 1, big_nsmall=big_nsmall)
        Plnew = newCDFasym(qmul, mu, muH, b, s, Delta, 1, big_nsmall=big_nsmall)
        if ib==0:
            h_af0.SetBinContent(ib+1, Ph0)
            h_afnew.SetBinContent(ib+1, Phnew)
        else:
            h_af0.SetBinContent(ib+1, Ph0-Pl0)
            h_afnew.SetBinContent(ib+1, Phnew-Plnew)
        pass
    for ib in range(nbins):
        CDF_toy = h.Integral(ib+1,nbins+1)
        CDF_af = h_af0.Integral(ib+1,nbins+1)
        CDF_afnew = h_afnew.Integral(ib+1,nbins+1)
        h_ptail.SetBinContent(ib+1, CDF_af - CDF_toy)
        h_ptailnew.SetBinContent(ib+1, CDF_afnew - CDF_toy)
    print('PDF intergral =', h_af0.Integral(), h_afnew.Integral())
    Cs_af = TCanvas('Cs_af', '', 800, 600)
    Cs_af.SetLogy()
    h.Draw('PE')
    h.SetMarkerStyle(20)
    h.SetLineWidth(2)
    h_af0.Draw('hist,same')
    h_af0.SetLineColor(kRed)
    h_af0.SetLineStyle(kDashed)
    h_afnew.Draw('hist,same')
    h_afnew.SetLineColor(kBlue)
    h.Draw('PEsame')
    h.GetYaxis().SetRangeUser(1e-5,1)
    #h.GetYaxis().SetRangeUser(1e-3,1)
    h.GetYaxis().SetTitle('Probability')
    h.GetXaxis().SetTitle('#tilde{q}_{#mu}')
    leg_af = TLegend(0.6, 0.7, 0.95, 0.92)
    leg_af.AddEntry(h, 'Toy', 'PE')
    leg_af.AddEntry(h_af0, 'Classic', 'L')
    leg_af.AddEntry(h_afnew, 'New', 'L')
    leg_af.SetFillStyle(0)
    leg_af.SetBorderSize(0)
    leg_af.Draw()
    plotname = 'Cs_af_b%i_s%i_mu%.2f_muH%i' % (b, s, mu, muH)
    Cs_af.SaveAs(plotname+'.png')
    Cs_af.SaveAs(plotname+'.pdf')

    Cs_ptail = TCanvas('Cs_ptail', '', 800, 600)
    h_ptail.Draw('hist')
    h_ptail.SetLineColor(kRed)
    h_ptail.SetLineStyle(kDashed)
    h_ptailnew.Draw('hist,same')
    h_ptailnew.SetLineColor(kBlue)
    L0 = TLine(xmin, 0, xmax, 0)
    L0.SetLineWidth(1)
    L0.SetLineStyle(2)
    L0.Draw()
    ymin = min(h_ptail.GetMinimum(), h_ptailnew.GetMinimum())
    ymax = max(h_ptail.GetMaximum(), h_ptailnew.GetMaximum())
    dy = ymax-ymin
    ymin = ymin - 0.5*dy
    ymax = ymax + 0.5*dy
    h_ptail.GetYaxis().SetRangeUser(ymin, ymax)
    h_ptail.GetYaxis().SetTitle('Tail Probability Difference')
    h_ptail.GetXaxis().SetTitle('#tilde{q}_{#mu}')
    leg_ptail = TLegend(0.6, 0.7, 0.95, 0.92)
    leg_ptail.AddEntry(h_ptail, 'Classic', 'L')
    leg_ptail.AddEntry(h_ptailnew, 'New', 'L')
    leg_ptail.SetFillStyle(0)
    leg_ptail.SetBorderSize(0)
    leg_ptail.Draw()
    plotname = 'Cs_ptail_b%i_s%i_mu%.2f_muH%i' % (b, s, mu, muH)
    Cs_ptail.SaveAs(plotname+'.png')
    Cs_ptail.SaveAs(plotname+'.pdf')

    for i in range(nbins):
        muobs0 = h_muhat.GetBinLowEdge(i+1)
        muobs1 = h_muhat.GetBinLowEdge(i+2)
        p0 = cdf_mu(muobs0, mu, muH, b, s, Delta=Delta, doTilde=1, big_nsmall=big_nsmall)
        p1 = cdf_mu(muobs1, mu, muH, b, s, Delta=Delta, doTilde=1, big_nsmall=big_nsmall)
        #print('i,muobs0, muobs1,p0,p1 =', i,muobs0, muobs1,p0,p1)
        _sigma = sigma
        if mu==muH:
            _sigma = sigma_muH
        if i==0:
            h_muhat_new.SetBinContent(i+1,p1)
            h_muhat_asymp.SetBinContent(i+1, Phi((muobs1-muH)/_sigma))
        else:
            h_muhat_new.SetBinContent(i+1,p1-p0)
            h_muhat_asymp.SetBinContent(i+1, Phi((muobs1-muH)/_sigma)-Phi((muobs0-muH)/_sigma))

    print('check integral: h_muhat_asymp, h_muhat_new =', h_muhat_asymp.Integral(), h_muhat_new.Integral())
    Cs_muhat = TCanvas('Cs_muhat', '', 800, 600)
    Cs_muhat.SetLogy()
    h_muhat.Draw('PE')
    h_muhat.SetMarkerStyle(20)
    h_muhat.SetLineWidth(2)
    h_muhat_asymp.Draw('hist,same')
    h_muhat_asymp.SetLineColor(kRed)
    h_muhat_asymp.SetLineStyle(kDashed)
    h_muhat_new.Draw('hist,same')
    h_muhat_new.SetLineColor(kBlue)
    h_muhat.Draw('PEsame')
    h_muhat.GetYaxis().SetRangeUser(1e-5,1)
    h_muhat.GetYaxis().SetTitle('PDF')
    h_muhat.GetXaxis().SetTitle('#hat{#mu}')
    plotname = 'Cs_muhat_b%i_s%i_mu%.2f_muH%i' % (b, s, mu, muH)
    Cs_muhat.SaveAs(plotname+'.png')

    Cs_qmu_muhat = TCanvas('Cs_qmu_muhat', '', 800, 600)
    h_qmu_muhat.Draw('colz')
    h_qmu_muhat.GetYaxis().SetTitle('q_{#mu}')
    h_qmu_muhat.GetXaxis().SetTitle('#hat{#mu}')
    fasym = TF1('fasym', asymp_wald, min(-0.5*mu,muH-5*sigma_muH), max(1.5*mu,muH+5*sigma_muH),2)
    if mu == muH:
        fasym.SetParameters(mu, sigma_muH)
    else:
        fasym.SetParameters(mu, sigma)
    fasym.Draw('same')
    fasym.SetLineStyle(2)
    fasym.SetLineWidth(2)
    plotname = 'Cs_qmu_muhat_b%i_s%i_mu%.2f_muH%i' % (b, s, mu, muH)
    Cs_qmu_muhat.SaveAs(plotname+'.png')
    Cs_qmu_muhat.SaveAs(plotname+'.pdf')

#get_toys(100, 10, 3, 3, ntoys=100000, nbins=50, xmin=0, xmax=15, big_nsmall=0)
#get_toys(100, 10, 3, 0, ntoys=100000, nbins=50, xmin=0, xmax=30, big_nsmall=0)
#get_toys(10, 3, 3, 3, ntoys=100000, nbins=50, xmin=0, xmax=20, big_nsmall=0)
#get_toys(10, 3, 3, 0, ntoys=100000, nbins=50, xmin=0, xmax=20, big_nsmall=0)
get_toys(3, 1, 3, 3, ntoys=100000, nbins=50, xmin=0, xmax=10, big_nsmall=0)
#get_toys(3, 1, 3, 0, ntoys=100000, nbins=50, xmin=0, xmax=10, big_nsmall=0)



