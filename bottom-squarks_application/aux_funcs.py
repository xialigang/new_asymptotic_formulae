from ROOT import ROOT, gROOT, TFile, TTree, TH1F, TPad, TCanvas, TLine, TLegend, THStack, TGraph, gPad, TGraphErrors, TColor, TLatex, TArrow, TH2F, TMath, TF1
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

yaxis_low = 1e-6

dict_toy_muhat = {}
dict_toy_thetahat_fixmu = {}
dict_toy_thetahat_floatmu = {}
def get_sigmas(mu, filepath, option = 0):
    lines = None
    lines0 = []
    lines1 = []
    with open(filepath,'r') as f:
        while 1:
            line = f.readline().strip()
            if line == '':
                break
            if line.startswith('nominal') or line.startswith('calib'):
                lines1.append(line)
            else:
                lines0.append(line)
                pass
            pass
        pass
    dict_mu_calib = {}
    dict_mu_qmu = {}
    dict_mumumu_sb = {}

    #0.25 -0.000297236 2.05234 0.927242 0 1.70645 1.7607 1.59571 0.0398229 0.233417
    for line in lines0:
        line = line.strip().split()
        mu = float(line[0])
        if option == 0:
            dict_mu_qmu[mu] = [float(line[2]), float(line[3])] # qmu_exp, qmu_obs
            dict_mu_calib[mu] = [float(line[5]), float(line[6])]# sigmad2L_muH0, sigmad2L_muH
        else:
            dict_mu_qmu[mu] = [float(line[3]), float(line[4])] # qmu_exp, qmu_obs

    a = None
    list_s = []
    list_b = []
    for line in lines1:
        line = line.strip().split()
        if option == 0:
            mu  = float(line[1])
            muH = float(line[2])
            muH1= float(line[3])
        else:
            mu  = float(line[1])
            muH = float(line[2])
            muH1= 0
        if a == None:
            a = (mu, muH, muH1)
        a = (mu, muH, muH1)
        _b = float(line[4])
        _s = float(line[5])/mu
        if len(list_s)>=3 and 1: # we have only 3 bins in SR, remove the CR 
            _s = 0
            _b = 0
        list_s.append(_s)
        list_b.append(_b)
        
        if len(list_s)==6:
            dict_mumumu_sb[a] = [list_s, list_b]
            list_s = []
            list_b = []
            pass
        pass
    return dict_mu_qmu, dict_mu_calib, dict_mumumu_sb


def Phi_quantile(p):
    return ROOT.Math.normal_quantile(p, 1)
def Phi(x):
    return ROOT.Math.normal_cdf(x, 1)

def fasym(X,P):
    x = X[0]
    mu = P[0]
    sigma = P[1]
    doTilde = P[2]
    if x>= mu:
        return 0
    elif x>0 or doTilde == 0:
        return pow((x-mu)/sigma,2)
    else:
        return (-2*x*mu+mu*mu)/(sigma*sigma)
def fasymq0(X,P):
    x = X[0]
    mu = P[0]
    sigma = P[1]
    if x<0:
        return 0
    else:
        return pow(x/sigma,2)
def get_toys(h, fname, treename='tree_stat', qmu_exp=1., qmu_obs=1., qmu_cut = -1., mu=1, muH=1, list_sigma=[], xmin=-2, xmax=10, ymin=-1, ymax=10, tag='', fast=0, do_q0=False, doTilde=0):
    h_qmu_muhat = TH2F('h_qmu_muhat', '', 100, xmin, xmax, 100, ymin, ymax)
    f = TFile(fname, 'read')
    treename = 'tree_stat'
    t = f.Get(treename)
    nentries = t.GetEntries()
    neff = 0.
    nobs = 0.
    nexp = 0.
    psum = 0.
    for i in range(nentries):
        t.GetEntry(i)
        psum += 1.
        muhat = getattr(t, 'POI_muhat', muH)
        if do_q0 and 0:
            if t.raw_status_fit_muhat != 0 or t.raw_status_fit_mu_0 != 0:
                continue
        #nb_events_toy = t.nb_events_toy
        if do_q0:
            qmu = t.q0
            #qmu = t.t0
            #if nb_events_toy != 1 and 0:
             #   continue
        else:
            if doTilde:
                qmu = t.qmutilde
            else:
                qmu = t.qmu
        if math.isnan(qmu):
            #print('qmu =', qmu)
            continue
        h_qmu_muhat.Fill(muhat, qmu)
        if do_q0:
            if qmu<0:
                continue
        if qmu < qmu_cut:
            continue
        if qmu < 0:
            #continue
            qmu = 0
        neff += 1.
        if qmu > qmu_exp:
            nexp += 1.
        if qmu > qmu_obs:
            nobs += 1.
        h.Fill(qmu)
        pass
    print('psum =', psum)
    psum = psum/t.GetEntries()
    print('psum =', psum)
    if do_q0:
        fmu = TF1('fmu', fasymq0, xmin, xmax, 2)
    else:
        fmu = TF1('fmu', fasym, xmin, xmax, 2)
    fmu.SetParameters(mu, list_sigma[2], doTilde)
    print('debug qmu_obs, nobs, neff =', qmu_obs, nobs, neff)
    pobs = nobs/neff
    pexp = nexp/neff
    dpobs = math.sqrt(pobs*(1.-pobs)/neff)
    dpexp = math.sqrt(pexp*(1.-pexp)/neff)
    if fast:
        return pexp, dpexp, pobs, dpobs

    Cs = TCanvas('Cs', '', 600, 600)
    h_qmu_muhat.Draw('colz')
    #h_qmu_muhat.Draw('scatt')
    h_qmu_muhat.SetMarkerStyle(6)
    fmu.Draw('same')
    fmu.SetLineWidth(3)
    fmu.SetLineColor(kBlack)
    fmu.SetLineStyle(2)

    #fnewmu.Draw('same')
    #fnewmu.SetLineWidth(3)
    #fnewmu.SetLineColor(kBlack)
    #fnewmu.SetLineStyle(2)

    h_qmu_muhat.GetXaxis().SetTitle("#hat{#mu}")
    if do_q0:
        h_qmu_muhat.GetYaxis().SetTitle("q_{0}")
    else:
        if doTilde:
            h_qmu_muhat.GetYaxis().SetTitle("#tilde{q}_{#mu}")
        else:
            h_qmu_muhat.GetYaxis().SetTitle("q_{#mu}")
    L0 = TLine(0, ymin, 0, ymax)
    L0.Draw()
    L0.SetLineWidth(2)
    #L0.SetLineColor(kGreen)
    L0.SetLineStyle(kDashed)

    Ly0 = TLine(xmin, 0, xmax, 0)
    Ly0.Draw()
    Ly0.SetLineWidth(2)
    #Ly0.SetLineColor(kGreen)
    Ly0.SetLineStyle(kDashed)

    Lmu = TLine(mu, ymin, mu, ymax)
    if not do_q0:
        Lmu.Draw()
    Lmu.SetLineWidth(2)
    #Lmu.SetLineColor(kGreen)
    Lmu.SetLineStyle(kDashed)

    Cs.SaveAs('Cs_h_qmu_muhat_%s_mu%.6f_muH%.6f.png' % (tag, mu, muH))
    Cs.SaveAs('Cs_h_qmu_muhat_%s_mu%.6f_muH%.6f.pdf' % (tag, mu, muH) )
    h.Scale(1./h.Integral(0,h.GetNbinsX()+1) * psum)
    pobs = nobs/neff
    pexp = nexp/neff
    dpobs = math.sqrt(pobs*(1.-pobs)/neff)
    dpexp = math.sqrt(pexp*(1.-pexp)/neff)
    return pexp, dpexp, pobs, dpobs

def get_toys_csv(h, fname, treename='tree_stat', qmu_exp=1., qmu_obs=1., qmu_cut = -1., mu=1, muH=1, list_sigma=[], xmin=-2, xmax=10, ymin=-1, ymax=10, tag='', fast=0, do_q0=False, doTilde=0):
    h_qmu_muhat = TH2F('h_qmu_muhat', '', 100, xmin, xmax, 100, ymin, ymax)
    lines = []
    with open(fname, 'r') as f:
        lines = f.readlines()
    nentries = len(lines)
    neff = 0.
    nobs = 0.
    nexp = 0.
    psum = 0.
    for i in range(nentries):
        line = lines[i].strip().split()
        psum += 1.
        muhat = float(line[0])
        qmu = float(line[1])
        if math.isnan(qmu):
            #print('qmu =', qmu)
            continue
        h_qmu_muhat.Fill(muhat, qmu)
        if qmu < qmu_cut:
            continue
        if qmu < 0:
            qmu = 0
        neff += 1.
        if qmu > qmu_exp:
            nexp += 1.
        if qmu > qmu_obs:
            nobs += 1.
        h.Fill(qmu)
        pass
    print('psum =', psum)
    psum = psum/nentries
    print('psum =', psum)
    if do_q0:
        fmu = TF1('fmu', fasymq0, xmin, xmax, 2)
    else:
        fmu = TF1('fmu', fasym, xmin, xmax, 2)
    fmu.SetParameters(mu, list_sigma[2], doTilde)
    print('debug qmu_obs, nobs, neff =', qmu_obs, nobs, neff)
    pobs = nobs/neff
    pexp = nexp/neff
    dpobs = math.sqrt(pobs*(1.-pobs)/neff)
    dpexp = math.sqrt(pexp*(1.-pexp)/neff)
    if fast:
        return pexp, dpexp, pobs, dpobs

    Cs = TCanvas('Cs', '', 600, 600)
    h_qmu_muhat.Draw('colz')
    #h_qmu_muhat.Draw('scatt')
    h_qmu_muhat.SetMarkerStyle(6)
    fmu.Draw('same')
    fmu.SetLineWidth(3)
    fmu.SetLineColor(kBlack)
    fmu.SetLineStyle(2)

    #fnewmu.Draw('same')
    #fnewmu.SetLineWidth(3)
    #fnewmu.SetLineColor(kBlack)
    #fnewmu.SetLineStyle(2)

    h_qmu_muhat.GetXaxis().SetTitle("#hat{#mu}")
    if do_q0:
        h_qmu_muhat.GetYaxis().SetTitle("q_{0}")
    else:
        if doTilde:
            h_qmu_muhat.GetYaxis().SetTitle("#tilde{q}_{#mu}")
        else:
            h_qmu_muhat.GetYaxis().SetTitle("q_{#mu}")
    L0 = TLine(0, ymin, 0, ymax)
    L0.Draw()
    L0.SetLineWidth(2)
    #L0.SetLineColor(kGreen)
    L0.SetLineStyle(kDashed)

    Ly0 = TLine(xmin, 0, xmax, 0)
    Ly0.Draw()
    Ly0.SetLineWidth(2)
    #Ly0.SetLineColor(kGreen)
    Ly0.SetLineStyle(kDashed)

    Lmu = TLine(mu, ymin, mu, ymax)
    if not do_q0:
        Lmu.Draw()
    Lmu.SetLineWidth(2)
    #Lmu.SetLineColor(kGreen)
    Lmu.SetLineStyle(kDashed)

    Cs.SaveAs('Cs_h_qmu_muhat_%s_mu%.6f_muH%.6f.png' % (tag, mu, muH))
    Cs.SaveAs('Cs_h_qmu_muhat_%s_mu%.6f_muH%.6f.pdf' % (tag, mu, muH) )
    h.Scale(1./h.Integral(0,h.GetNbinsX()+1) * psum)
    pobs = nobs/neff
    pexp = nexp/neff
    dpobs = math.sqrt(pobs*(1.-pobs)/neff)
    dpexp = math.sqrt(pexp*(1.-pexp)/neff)
    return pexp, dpexp, pobs, dpobs

def Fclassic(qmu=1., mu=1., mup=1., sigma=1., doTilde=True):
    if qmu <= pow(mu/sigma,2) or doTilde == False:
        return ROOT.Math.normal_cdf(math.sqrt(qmu)-(mu-mup)/sigma)
    else:
        return ROOT.Math.normal_cdf((qmu-(mu*mu-2*mu*mup)/sigma/sigma)/(2*mu/sigma))
    pass




def get_qmu_from_musb(muhat, mu, list_s, list_b, doTilde = 0, mu1=-9999):
    nbins = len(list_s)
    qmu = 0
    if mu1 == -9999:
        mu1 = muhat
    for i in range(nbins):
        si = list_s[i]
        bi = list_b[i]
        if i==nbins-1:
            if bi>10:
                si = 10/bi*si
                bi = 10
        ni = bi + mu1*si
        if doTilde and mu1<0:
            if ni<=0:
                qmu += 2*mu*si
            else:
                qmu += -2*(ni*math.log((bi+mu*si)/(bi+0*si))-(mu-0)*si)
        else:
            if ni<=0:
                qmu += 2*(mu-muhat)*si
            else:
                qmu += -2*(ni*math.log((bi+mu*si)/(bi+muhat*si))-(mu-muhat)*si)
    return qmu

def get_sigma_from_musb(muH, list_s, list_b):
    nbins = len(list_s)
    sigma = 0.
    for i in range(nbins):
        si = list_s[i]
        bi = list_b[i]
        ni = bi + muH*si
        if si == 0:
            continue
        sigma += si*si/ni
    sigma = 1./math.sqrt(sigma)
    return sigma


def CDFasymq0(q0, muH, sigma, sigmap):
    return Phi((sigma*math.sqrt(q0)-muH)/sigmap)
    #return Phi((sigma*math.sqrt(q0)-muH)/sigmap)-Phi((-sigma*math.sqrt(q0)-muH)/sigmap)


def CDFasym(qmu, mu, muH, sigma, sigmastar, sigmap, doTilde=1, doBound=0, d1L_sigma = 0):
    if doBound == 1:
        if doTilde == 1:
            _qmu = (-2*mu*muH+mu*mu)/pow(sigma,2)
            delta_qmu = _qmu*2*sigmap/(mu-2*muH)
            return Phi((qmu-_qmu)/delta_qmu)
        else:
            muhat = mu-sigma*math.sqrt(qmu)
            Bhalf = 0.5+Phi(-d1L_sigma)*math.exp(0.5*pow(d1L_sigma,2))
            muhat_min = muH
            if muhat >= muhat_min:
                return Bhalf/(0.5)*Phi(sigma/sigmap*math.sqrt(qmu)-(mu-muH)/sigmap)
            else:
                return Bhalf+(1-Bhalf)/(0.5)*(0.5-Phi(-sigma/sigmap*math.sqrt(qmu)+(mu-muH)/sigmap))

    if qmu <= pow(mu/sigma,2) or doTilde==False:
        return Phi(sigma/sigmap*math.sqrt(qmu)-(mu-muH)/sigmap)
    else:
        #return Phi((qmu-(mu*mu-2*mu*muH)/sigma/sigma)/(2.*mu*sigmap/sigma/sigma))
        return Phi((qmu-pow(mu/sigma,2)+2*mu*muH/pow(sigmastar,2))/(2.*mu*sigmap/pow(sigmastar,2)))




def get_muhat(x=[], test_mu=1., init_mu=0.1, do_q0=False, flag_bigb5=False, doTilde=True, fmup=0):
    N = len(x)
    s = 0.
    n = 0.
    istar = -1
    for i in range(N):
        #print('debug :', i, x[i])
        s += x[i][0]
        n += x[i][2]
        if x[i][2]>0 and istar <0:
            istar = i
        pass
    #print('debug istar =', istar)
    A = 0.
    B = 0.
    if n==0:
        mu= -x[0][1]/x[0][0]
        if mu < -x[0][1]/x[0][0]:
            mu= -x[0][1]/x[0][0]
        qmu = 0.
        q0 = 0.
        for i in range(N):
            if mu>0 or doTilde==0:
                qmu += 2.*(-(mu-test_mu)*x[i][0])
                q0 += 2.*(-(mu-0.)*x[i][0])
            else:
                qmu += 2.*(-(0-test_mu)*x[i][0])
                q0 += 2.*(-(mu-0.)*x[i][0])
            pass
        if do_q0:
            return mu, q0
        else:
            return mu, qmu
    qmu = 0.
    q0 = 0.
    itr = 0
    mu = init_mu
    mu_pre = init_mu
    if 1:
        #print('debug qmu<0 before: init_mu, mu_pre, mu, qmu =', init_mu, mu_pre, mu, qmu)
        mumin = -x[0][1]/x[0][0]*0.999999
        #mumin = -x[0][1]/x[0][0]*0.5
        #mumin = -0.01335
        #mumin = test_mu*1.
        #mumax = test_mu*1.0
        mumax = test_mu*5.0
        if do_q0:
            mumax = 6
        Nmu = 1000
        dmu = (mumax-mumin)/Nmu
        logLmax = -9999
        init_mu = 0.
        for imu in range(-1, Nmu+1):
            logL = 0.
            _mu = mumin + imu*dmu
            if imu == -1:
                _mu = 0
            for i in range(N):
                if x[i][1] + _mu*x[i][0] > 0:
                    logL += x[i][2]*math.log(x[i][1] + _mu*x[i][0]) - (x[i][1] + _mu*x[i][0])
                elif x[i][1] + _mu*x[i][0] == 0 and x[i][2] == 0:
                    logL +=  - (x[i][1] + _mu*x[i][0])
                else:
                    logL += -9999
                if x[i][2]>0:
                    logL -= x[i][2]*math.log(x[i][2]) - x[i][2]
                pass
            #print('debug : _mu, logL =', _mu, logL)
            if logL > logLmax:
                logLmax = logL
                init_mu = _mu
        mu = init_mu
        qmu = 0.
        q0 = 0.
        for i in range(N):
            q0 += 2.*(x[i][2]*math.log((x[i][1]+mu*x[i][0])/(x[i][1]+0*x[i][0]))-(mu-0)*x[i][0])
            if mu>0 or doTilde==0:
                if x[i][2] == 0 and x[i][1]+mu*x[i][0]==0:
                    qmu += 2.0*(-(mu-test_mu)*x[i][0])
                elif x[i][2]>0 and x[i][1]+mu*x[i][0]<=0:
                    for j in range(10):
                        print('ERROR! negative expected n!!!')
                else:
                    qmu += 2.*(x[i][2]*math.log((x[i][1]+mu*x[i][0])/(x[i][1]+test_mu*x[i][0]))-(mu-test_mu)*x[i][0])
            else:
                qmu += 2.*(x[i][2]*math.log((x[i][1]+0*x[i][0])/(x[i][1]+test_mu*x[i][0]))-(0-test_mu)*x[i][0])
            pass
    if do_q0:
        return mu, q0
    else:
        return mu, qmu




def get_nsmall(mu, muH, dict_mumumu_sb):
    _b = 0
    _s = 0
    _n = 0
    _b5= 0
    for i in range(6):
        #si = dict_mumumu_sb[(mu,muH,muH)][0][i]
        #bi = dict_mumumu_sb[(mu,muH,muH)][1][i]
        si = dict_mumumu_sb[(mu,muH,0)][0][i]
        bi = dict_mumumu_sb[(mu,muH,0)][1][i]
        if i==5 and bi>10:
            si = si*10/bi
            bi = 10
        if i==5:
            _b5 = bi
        _b += bi
        _s += si
    _n = int(_b+mu*_s)-1
    if _n > 10:
        _n = 10
    if _n < 3:
        _n = 3
    if muH > 0 and 0:
        _n += 5
    return _b,_s,_n

def get_n_for_0(mu, dict_mumumu_sb, sigma_muH, sigma_muH0):
    prob_improve_mu = 0.
    _b, _s, _n = get_nsmall(mu, mu, dict_mumumu_sb)
    for k in range(0,_n+1):
        prob_improve_mu += TMath.Poisson(k, _b + mu*_s)
    prob_improve_0 = 0.
    b, s, n = get_nsmall(mu, 0, dict_mumumu_sb)
    n_for_0 = 0
    for k in range(0,max(n,_n)+1):
        prob_improve_0 += TMath.Poisson(k, b + 0*s)
        if prob_improve_0 >= prob_improve_mu:
            n_for_0 = k
            break
        pass
    #print('prob(mu), prob(0) =', Phi(((_n-_b)/_s-mu)/sigma_muH), Phi((n_for_0-b)/s/sigma_muH0))
    #print('prob_improve_mu, prob_improve_0, mu, muH, n, n_for_0 =', prob_improve_mu, prob_improve_0, mu, 0, n, n_for_0)
    if n_for_0 < 3:
        n_for_0 = 3
    return n_for_0

def Fnew(qmu=1, mu=1, muH=1, list_sigma=[], dict_mumumu_sb={}, dict_mu_qmu={}, dict_mu_calib={}, doTilde=1, nsmall = -9999):

    debug = 0


    _muH = (mu, muH, muH)
    _muH = (mu, muH, 0) # for sbottom search
    if nsmall>-9999:
        print('mu, muH, dict_mumumu_sb =', mu, muH, dict_mumumu_sb[_muH])
    
    s0 = dict_mumumu_sb[_muH][0][0]
    s1 = dict_mumumu_sb[_muH][0][1]
    s2 = dict_mumumu_sb[_muH][0][2]
    s3 = dict_mumumu_sb[_muH][0][3]
    s4 = dict_mumumu_sb[_muH][0][4]
    s5 = dict_mumumu_sb[_muH][0][5]

    b0 = dict_mumumu_sb[_muH][1][0]
    b1 = dict_mumumu_sb[_muH][1][1]
    b2 = dict_mumumu_sb[_muH][1][2]
    b3 = dict_mumumu_sb[_muH][1][3]
    b4 = dict_mumumu_sb[_muH][1][4]
    b5 = dict_mumumu_sb[_muH][1][5]

    bmax = 10.
    if b5>bmax and bmax>0:
        s5 = s5 * bmax/b5
        b5 = bmax

    s = s0 + s1 + s2 + s3 + s4 + s5
    b = b0 + b1 + b2 + b3 + b4 + b5

    stot = s

    muhat_min = -b0/s0
    list_s = [s0, s1, s2, s3, s4, s5]
    list_b = [b0, b1, b2, b3, b4, b5]
    if debug:
        print('s =', s)
        print('b =', b)
        print('list_s =', list_s)
        print('list_b =', list_b)

    #stot_muH0 = 0.
    #stot_muH = 0.
    #for i in range(6):
    #    stot_muH0 += dict_mumumu_sb[(mu, 0, 0)][0][i]
    #    stot_muH += dict_mumumu_sb[(mu, mu, mu)][0][i]
    #print('stot_muH0, stot_muH =', stot_muH0, stot_muH)

    special_sigma0 = list_sigma[0]
    special_k = list_sigma[1]
    sigma = list_sigma[2]
    sigma_muH0 = list_sigma[3]
    sigma_muH = list_sigma[4]
    sigma_mu = sigma_muH

    #special_sigma0 = math.sqrt(pow(special_sigma0,2) + pow(extra_dmu,2))
    #special_k = math.sqrt(pow(special_k,2) + pow(extra_dmu/mu,2))

    rqmu_negative = 1
    rqmu_0 = 1
    rqmu_mu = 1

    #_qmu_0 = get_qmu_from_musb(0, mu, dict_mumumu_sb[(mu, muH, 0)][0], dict_mumumu_sb[(mu, muH, 0)][1], doTilde=doTilde)
    #qmu_0 = dict_mu_qmu[mu][0] if muH == 0 else dict_mu_calib[mu][1]
    #rqmu_0 = qmu_0 / _qmu_0

    _qmu_0 = get_qmu_from_musb(0, mu, dict_mumumu_sb[(mu, 0, 0)][0], dict_mumumu_sb[(mu, 0, 0)][1], doTilde=doTilde)
    qmu_0 = dict_mu_qmu[mu][0]
    rqmu_0 = qmu_0 / _qmu_0


    #list_qmu_calib = [[_muhat_min, rqmu_negative], [0, rqmu_0]]
    list_qmu_calib = [[0, rqmu_0]]
    list_qmu_calib.sort()

    if nsmall >-9999 or 1:
        print('mu, muH, calib, =', mu, muH, list_qmu_calib)

    Poisson_big = 1.
    #n = int(b+mu*s)-1
    #if n > max(6, int(b5)) and 1:
    #    n = max(6, int(b5))

    _b, _s, n = get_nsmall(mu, muH, dict_mumumu_sb)
    #n -= 1
    if muH == 0:
        n = get_n_for_0(mu, dict_mumumu_sb, sigma_muH, sigma_muH0)
        if n < 3 and 0:
            n = 3

    if nsmall > -9999:
        n = nsmall
        print('Manually setting nsmall =', nsmall)
    print('nsmall =', n)
    muH_corr = 0.
    sigma_corr = 0.

    for k in range(0,n+1):
        Poisson_big -= TMath.Poisson(k, b + muH*s)
    Poisson_small = 1. - Poisson_big
    Poisson_small = 0.
    Prob_small = 0.

    muH_corr = 0.
    sigma_corr = 0.

    qmu_exp = dict_mu_qmu[mu][0]

    for k in range(0,n+1):
        if qmu == qmu_exp and nsmall > -9999:
            print('we are working on k =',k)
        Poisson_k = TMath.Poisson(k, b + muH*s)
        for i0 in range(0, k+1):
            for i1 in range(0, k-i0+1):
                for i2 in range(0, k-i0-i1+1):
                    for i3 in range(0, k-i0-i1-i2+1):
                        if b3+s3==0 and i3>0:
                            continue
                        for i4 in range(0, k-i0-i1-i2-i3+1):
                            if b4+s4==0 and i4>0:
                                continue
                            i5 = k-i0-i1-i2-i3-i4
                            if b5+s5 == 0 and i5 != 0:
                                continue
                            x = [[s0, b0, i0]]
                            if s1+b1>0:
                                x.append([s1, b1, i1])
                            if s2+b2>0:
                                x.append([s2, b2, i2])
                            if s3+b3>0:
                                x.append([s3, b3, i3])
                            if s4+b4>0:
                                x.append([s4, b4, i4])
                            if s5>0:
                                #x.append([s5, b5 + db5, i5 + db5])
                                x.append([s5, b5, i5])
                            list_n = [i0, i1, i2, i3, i4, i5]
                            toy = '%.2f,%.2f,%i,%i,%i,%i,%i,%i' % (muH, mu, i0, i1, i2, i3, i4, i5)
                            mup = -9999
                            _qmu = -9999
                            flag_n0 = 0
                            if len(dict_toy_muhat)==0 or toy not in dict_toy_muhat:
                                mup, _qmu = get_muhat(x, mu, muH, flag_bigb5=False, doTilde=doTilde)
                                dict_toy_muhat[toy] = [mup]
                            else:
                                res = dict_toy_muhat[toy]
                                mup = res[0]
                                _qmu = get_qmu_from_muhat(x, mu, mup, 0, 0, 0, doTilde=doTilde)
                            sigmap = math.sqrt(pow(special_sigma0,2)+pow(special_k*mup,2)) #251
                            Poisson_small += Poisson_k * math.factorial(k)/math.factorial(i5)/math.factorial(i0)/math.factorial(i1)/math.         factorial(i2)/math.factorial(i3)/math.factorial(i4) * pow((muH*s5+b5)/(muH*s+b), i5) * pow((muH*s0+b0)/(muH*s+b),i0) * pow((muH*s1+b1)/(muH*s+ b),i1) * pow((muH*s2+b2)/(muH*s+b),i2) * pow((muH*s3+b3)/(muH*s+b),i3) * pow((muH*s4+b4)/(muH*s+b),i4)
                            sigma_corr += pow(mup,2) * Poisson_k * math.factorial(k)/math.factorial(i5)/math.factorial(i0)/math.factorial(i1)/math.         factorial(i2)/math.factorial(i3)/math.factorial(i4) * pow((muH*s5+b5)/(muH*s+b), i5) * pow((muH*s0+b0)/(muH*s+b),i0) * pow((muH*s1+b1)/(muH*s+ b),i1) * pow((muH*s2+b2)/(muH*s+b),i2) * pow((muH*s3+b3)/(muH*s+b),i3) * pow((muH*s4+b4)/(muH*s+b),i4)
                            sigma1 = sigma
                            sigma2 = sigma
                            _mu = mu
                            doBound = False
                            if _qmu>0 and mup != mu:
                                _rqmu = 1
                                n_mu_calib = len(list_qmu_calib)
                                if mup < list_qmu_calib[0][0]:
                                    _rqmu = list_qmu_calib[0][1]
                                elif mup > list_qmu_calib[n_mu_calib-1][0]:
                                    _rqmu = list_qmu_calib[n_mu_calib-1][1]
                                else:
                                    for j in range(n_mu_calib-1):
                                        aj = list_qmu_calib[j]
                                        ajp1 = list_qmu_calib[j+1]
                                        if aj[0]<= mup and ajp1[0]>mup:
                                            _rqmu = aj[1]+(ajp1[1]-aj[1])/(ajp1[0]-aj[0])*(mup-aj[0])
                                            break

                                if i0+i1+i2+i3+i4==0 and 0:
                                    if doTilde==0:
                                        #_qmu *= 1-(mu+mup)/2*stot*pow(special_k,2)
                                        _qmu -= (mu*mu-mup*mup)*pow(stot*special_k,2)
                                    else:
                                        #_qmu *= 1-mu/2*stot*pow(special_k,2)
                                        _qmu -= (mu*mu)*pow(stot*special_k,2)
                                else:
                                    _qmu *= _rqmu
                                if mup>=0 or doTilde == 0:
                                    sigma1 = abs(mup-mu)/math.sqrt(_qmu)
                                    sigma2 = sigma1

                                    delta_qmu = 2*stot*math.sqrt(pow(special_sigma0,2)+pow((mu-muhat_min)*special_k,2))
                                    rel_delta_qmu = delta_qmu / (2*(mu-muhat_min)*stot)
                                    qmu_bound = 2*(mu-muhat_min)*stot
                                    if doTilde==0 and i0+i1+i2+i3+i4 == 0: #doTilde = 1
                                        sigmap = (mu-mup)/2*rel_delta_qmu
                                        #sigmap = math.sqrt(pow((mu-mup)/2*rel_delta_qmu,2)+pow(sigmap,2))
                                    elif doTilde==0 and abs((mup-muhat_min)/muhat_min)<0.001:
                                        sigmap = math.sqrt(pow((mu-mup)/2*rel_delta_qmu,2)+pow(sigmap,2))
                                else:
                                    sigma2 = math.sqrt((-2.*mu*mup+mu*mu)/_qmu)
                                    sigma1 = sigma2
                                    if i0+i1+i2+i3+i4 == 0: #doTilde = 1
                                        delta_qmu = 2*mu*stot*special_k
                                        rel_delta_qmu = delta_qmu / (2*mu*stot)
                                        sigmap = (mu-2*mup)/2.*delta_qmu/_qmu
                                    elif abs((mup-muhat_min)/muhat_min)<0.001:
                                        delta_qmu = 2*mu*stot*special_k
                                        rel_delta_qmu = delta_qmu / (2*mu*stot)
                                        sigmap = math.sqrt(pow(sigmap,2)+pow((mu-2*mup)/2.*delta_qmu/_qmu,2))
                            else:
                                _qmu = 0
                                mup = muH
                                sigma1 = sigma
                                sigma2 = sigma
                            mnp = pow((muH*s0+b0)/(muH*s+b),i0)
                            if s1+b1>0:
                                mnp *= pow((muH*s1+b1)/(muH*s+b),i1)
                            if s2+b2>0:
                                mnp *= pow((muH*s2+b2)/(muH*s+b),i2)
                            if s3+b3>0:
                                mnp *= pow((muH*s3+b3)/(muH*s+b),i3)
                            if s4+b4>0:
                                mnp *= pow((muH*s4+b4)/(muH*s+b),i4)
                            if s5+b5>0:
                                mnp *= pow((muH*s5+b5)/(muH*s+b),i5)
                            Prob_small += Poisson_k * math.factorial(k)/math.factorial(i5)/math.factorial(i0)/math.factorial(i1)/math.factorial(i2)/math.factorial(i3)/math.factorial(i4)*mnp*CDFasym(qmu, _mu, mup, sigma1, sigma2, sigmap, doTilde=doTilde) 
                            #Prob_small += Poisson_k * math.factorial(k)/math.factorial(i5)/math.factorial(i0)/math.factorial(i1)/math.factorial(i2)/math.factorial(i3)/math.factorial(i4) * pow((muH*s5+b5)/(muH*s+b), i5) * pow((muH*s0+b0)/(muH*s+b),i0) * pow((muH*s1+b1)/(muH*s+b),i1) * pow((muH*s2+b2)/(muH*s+b),i2) * pow((muH*s3+b3)/(muH*s+b),i3) * pow((muH*s4+b4)/(muH*s+b),i4)* CDFasym(qmu, _mu, mup, sigma1, sigma2, sigmap, doTilde=doTilde)
                            muH_corr += mup*Poisson_k * math.factorial(k)/math.factorial(i5)/math.factorial(i0)/math.factorial(i1)/math.factorial(i2)/math.factorial(i3)/math.factorial(i4)*mnp
                            #muH_corr += mup * Poisson_k * math.factorial(k)/math.factorial(i5)/math.factorial(i0)/math.factorial(i1)/math.         factorial(i2)/math.factorial(i3)/math.factorial(i4) * pow((muH*s5+b5)/(muH*s+b), i5) * pow((muH*s0+b0)/(muH*s+b),i0) * pow((muH*s1+b1)/(muH*s+ b),i1) * pow((muH*s2+b2)/(muH*s+b),i2) * pow((muH*s3+b3)/(muH*s+b),i3) * pow((muH*s4+b4)/(muH*s+b),i4)
            pass
        pass
    #print('size of dict_toy_muhat =', len(dict_toy_muhat))
    #print('size of dict_toy_thetahat_fixmu =', len(dict_toy_thetahat_fixmu))
    #print('size of dict_toy_thetahat_floatmu =', len(dict_toy_thetahat_floatmu))

    if Poisson_small == 0 or Prob_small == 0:
        muH_corr = muH
    else:
        muH_corr = muH_corr/Poisson_small # average
        muH_corr = (muH - muH_corr*Poisson_small)/Poisson_big

    #sigma_corr = math.sqrt((sigma*sigma+muH*muH-sigma_corr)/Poisson_big-muH_corr*muH_corr)
    if debug:
        print('nsmall, b+muH*s, mu, muH =',n, b+muH*s, mu, muH)
        print('second: muH_corr =', muH_corr, muH)
        print('second: sigma_corr =', sigma_corr, sigma)

    sigma_wald = sigma
    sigma_corr = sigma
    print('mu, muH, muH_corr, sigma, sigma_corr =', mu, muH, muH_corr, sigma, sigma_corr)
    Prob_big = Poisson_big * CDFasym(qmu, mu, muH_corr, sigma_corr, sigma_corr, sigma_corr,  doTilde=doTilde)
    if debug == 1:
        print('Prob_big, Prob_small =', Prob_big, Prob_small)
    Prob = Prob_big + Prob_small
    #Prob = Prob_small
    #Prob = Prob_big
    return Prob

def Fnewq0(q0, mu, muH, list_sigma, list_sb, list_q0_calib=[]):
    #mu = 0.
    debug = 0

    s0 = list_sb[0]
    s1 = list_sb[1]
    s2 = list_sb[2]
    s3 = list_sb[3]
    s4 = list_sb[4]
    s5 = list_sb[5]

    b0 = list_sb[6]
    b1 = list_sb[7]
    b2 = list_sb[8]
    b3 = list_sb[9]
    b4 = list_sb[10]
    b5 = list_sb[11]

    special_sigma0 = list_sb[12]
    special_k = list_sb[13]

    db5 = 0.
    if b5>10:
        db5 = b5 - 10.
        b5 = 10.
        s5 = 0.


    s = s0 + s1 + s2 + s3 + s4 + s5
    b = b0 + b1 + b2 + b3 + b4 + b5
    #print('s0, s1, b0, b1, b+muH*s =', s0, s1, b0, b1, b + muH*s)

    stot = s
    muhat_min = -b0/s0
    sigma = list_sigma[0]
    #sigma0= list_sigma[1]
    Poisson_big = 1.

    #n = int(b+mu*s)-1
    n = int(b+muH*s)-1

    #n = int(b+mu*s)
    #n = int(b+muH*s)-1
    #n = 6
    if n > max(6, int(b5)):
        n = max(6, int(b5))

    #n -= 3
    #n -= 2
    #n -= 1
    #n += 1
    #n += 2
    #n = int(b)
    #n = int(b5)
    #n = 10
    print('debug nsmall =', n)
    #n = 9
    #n = 3
    muH_corr = 0.
    sigma_corr = 0.
    #for k in range(0, n+1):
    #    muH_corr += (n-b)/s*TMath.Poisson(k, b+muH*s)
    #    sigma_corr += pow((n-b)/s,2)*TMath.Poisson(k, b+muH*s)
    #    pass
    for k in range(0,n+1):
        Poisson_big -= TMath.Poisson(k, b + muH*s)
    Poisson_small = 1. - Poisson_big
    Prob_small = 0.
    sigma_muH0 = sigma
    sigma = sigma_muH0

    list_s = [s0, s1, s2, s3, s4, s5]
    list_b = [b0, b1, b2, b3, b4, b5]

    _qmu_exp = get_qmu_from_musb(0, mu, list_s, list_b)
    _sigma = get_sigma_from_musb(0, list_s, list_b)
    rq0 = pow(_sigma/sigma_muH0,2)
    list_mu_q0_exp = [[0, rq0]]
    for a in list_q0_calib:
        mu_key = a[0]
        q0_key = a[1]
        _q0_key = get_qmu_from_musb(mu_key, 0, list_s, list_b)
        rq0_key = 1
        if q0_key > 0 and _q0_key>0:
            rq0_key = q0_key / _q0_key
        list_mu_q0_exp.append([mu_key, rq0_key])
    list_mu_q0_exp.sort()

    #rq0 = 1.
    print('debug rq0, rq0_key =', rq0, list_mu_q0_exp)



    muH_corr = 0.
    sigma_corr = 0.
    flag_check = 0
    for k in range(0,n+1):
        Poisson_k = TMath.Poisson(k, b + muH*s)
        if flag_check:
            if k != 1:
                continue
        for i0 in range(0, k+1):
            for i1 in range(0, k-i0+1):
                for i2 in range(0, k-i0-i1+1):
                    for i3 in range(0, k-i0-i1-i2+1):
                        for i4 in range(0, k-i0-i1-i2-i3+1):
                            i5 = k-i0-i1-i2-i3-i4
                            x = [
                                    [s0, b0, i0],
                                    [s1, b1, i1],
                                    [s2, b2, i2],
                                    [s3, b3, i3],
                                    [s4, b4, i4],
                                    ]
                            if s5>0:
                                #x.append([s5, b5 + db5, i5 + db5])
                                x.append([s5, b5, i5])
                            list_n = [i0, i1, i2, i3, i4, i5]
                            toy = '%i,%i,%i,%i,%i,%i' % (i0, i1, i2, i3, i4, i5)
                            mup = -9999
                            _q0 = -9999
                            #mup, _q0 = get_muhat(x, 0, muH, do_q0=True)
                            if len(dict_toy_muhat)==0 or toy not in dict_toy_muhat:
                                mup, _q0 = get_muhat(x, 0, muH, do_q0=True)
                                dict_toy_muhat[toy] = [mup]
                            else:
                                res = dict_toy_muhat[toy]
                                mup = res[0]
                                _q0 = get_qmu_from_muhat(x, 0, mup, 0, 0, 0, doTilde=0)
                            #mup, _q0 = get_muhat(x, 0, muH, do_q0=True)
                            if (abs((mup-muhat_min)/muhat_min)<0.001 or i0+i1+i2+i3+i4==0) and 0:
                                _q0 = 0-2*mup*stot + pow(mup*special_k*stot,2)
                            else:
                                n_mu_calib = len(list_mu_q0_exp)
                                _rq0 = 1
                                if mup < list_mu_q0_exp[0][0]:
                                    _rq0 = list_mu_q0_exp[0][1]
                                elif mup > list_mu_q0_exp[n_mu_calib-1][0]:
                                    _rq0 = list_mu_q0_exp[n_mu_calib-1][1]
                                else:
                                    for j in range(n_mu_calib-1):
                                        aj = list_mu_q0_exp[j]
                                        ajp1 = list_mu_q0_exp[j+1]
                                        if aj[0]<= mup and ajp1[0]>mup:
                                            _rq0 = aj[1]+(ajp1[1]-aj[1])/(ajp1[0]-aj[0])*(mup-aj[0])
                                            break
                                _q0 *= _rq0
                            sigmap = math.sqrt(pow(special_sigma0,2)+pow(special_k*mup,2)) #251
                            muH_corr += mup * Poisson_k * math.factorial(k)/math.factorial(i5)/math.factorial(i0)/math.factorial(i1)/math.         factorial(i2)/math.factorial(i3)/math.factorial(i4) * pow((muH*s5+b5)/(muH*s+b), i5) * pow((muH*s0+b0)/(muH*s+b),i0) * pow((muH*s1+b1)/(muH*s+ b),i1) * pow((muH*s2+b2)/(muH*s+b),i2) * pow((muH*s3+b3)/(muH*s+b),i3) * pow((muH*s4+b4)/(muH*s+b),i4)
                            sigma_corr += (pow(mup,2)+pow(sigmap,2))* Poisson_k * math.factorial(k)/math.factorial(i5)/math.factorial(i0)/math.factorial(i1)/math.         factorial(i2)/math.factorial(i3)/math.factorial(i4) * pow((muH*s5+b5)/(muH*s+b), i5) * pow((muH*s0+b0)/(muH*s+b),i0) * pow((muH*s1+b1)/(muH*s+ b),i1) * pow((muH*s2+b2)/(muH*s+b),i2) * pow((muH*s3+b3)/(muH*s+b),i3) * pow((muH*s4+b4)/(muH*s+b),i4)
                            if _q0>0:
                                sigma1 = abs(mup)/math.sqrt(_q0)
                                #if abs((mup-muhat_min)/muhat_min)<0.001 and i0+i1+i2+i3+i4==0:
                                if abs((mup-muhat_min)/muhat_min)<0.001 or i0+i1+i2+i3+i4==0:
                                    #sigmap = special_k*(0-2.*mup)/2.
                                    #sigmap /= 1.-(0+mup)*pow(special_k,2)*stot/2.
                                    delta_qmu = 2*stot*math.sqrt(pow(special_sigma0,2)+pow((0-mup)*special_k,2))
                                    rel_delta_qmu = delta_qmu / (2*abs(0-mup)*stot)
                                    sigmap = abs(0-mup)/2*rel_delta_qmu
                            else:
                                sigma1 = sigma
                            Prob_small += Poisson_k * math.factorial(k)/math.factorial(i5)/math.factorial(i0)/math.factorial(i1)/math.factorial(i2)/math.factorial(i3)/math.factorial(i4) * pow((muH*s5+b5)/(muH*s+b), i5) * pow((muH*s0+b0)/(muH*s+b),i0) * pow((muH*s1+b1)/(muH*s+b),i1) * pow((muH*s2+b2)/(muH*s+b),i2) * pow((muH*s3+b3)/(muH*s+b),i3) * pow((muH*s4+b4)/(muH*s+b),i4)* CDFasymq0(q0, mup, sigma1, sigmap)
                            #Prob_small += CDFasymq0(q0, mup, sigma1, sigmap)
            pass
        pass
    muH_corr_small = muH_corr/Poisson_small # average
    sigma_corr_small = math.sqrt(sigma_corr/Poisson_small-muH_corr_small*muH_corr_small) # standard deviation, sqrt(<X^2>-<X>^2)

    muH_corr = (muH - muH_corr_small*Poisson_small)/Poisson_big
    sigma_corr_big = math.sqrt((pow(muH,2)+pow(sigma,2)-Poisson_small*(pow(muH_corr_small,2)+pow(sigma_corr_small,2)))-pow(muH_corr,2))
    muH_corr_big = muH + Poisson_small/Poisson_big*(muH-muH_corr_small)*math.exp(-0.5*(pow(muH-muH_corr_small,2)-pow(muH-muH_corr,2))/pow(sigma,2))
    if debug or 1:
        print('nsmall, b+muH*s, mu, muH =',n, b+muH*s, mu, muH)
        print('second: muH_corr_small, muH_corr, muH_corr_big, muH =', muH_corr_small, muH_corr)
        print('second: sigma_corr_small, sigma_corr_big, sigma =', sigma_corr_small, sigma_corr_big, sigma)

    Prob_big = Poisson_big * CDFasymq0(q0, muH_corr, sigma, sigma)

    if debug == 1:
        print('Prob_big, Prob_small =', Prob_big, Prob_small)
    Prob = Prob_big + Prob_small
    #Prob = Prob_small
    return Prob



def turn2log(x, logy=0):
    if logy==0:
        return x
    else:
        return pow(10,x)
def is_missing(list_f=[]):
    flag = 0
    for f in list_f:
        if not os.path.exists(f):
            print(f,'is missing !!')
            flag = 1
            break
    return flag







def get_qmu_from_muhat(x=[], test_mu=1., mu=1., delta=0., theta=0., theta_cond=0., doTilde=True):
    qmu = 0.
    N = len(x)
    s = 0.
    n = 0.
    for i in range(N):
        s += x[i][0]
        n += x[i][2]
    if n==0 or 0:
        mu= -x[0][1]/x[0][0]
        qmu = 0.
        for i in range(N):
            if mu>0 or doTilde==0:
                qmu += 2.*(-(mu-test_mu)*x[i][0])
            else:
                qmu += 2.*(-(0-test_mu)*x[i][0])
        return qmu
    for i in range(N):
        n = x[i][2]
        s = x[i][0]
        b = x[i][1]
        v = b*(1+delta*theta) + mu*s
        if mu<0 and doTilde==1:
            v = b*(1+delta*theta) + 0*s
        v_cond = b*(1+delta*theta_cond) + test_mu*s
        if n==0:
            qmu += 2.*(-(v-v_cond))
        elif n>0 and v_cond == 0:
            for j in range(10):
                print('ERROR! v_cond =', v_cond)
        else:
            qmu += 2.*(n*math.log(v/v_cond)-(v-v_cond))
    if delta>0:
        qmu += pow(theta_cond,2) - pow(theta,2)
    return qmu

