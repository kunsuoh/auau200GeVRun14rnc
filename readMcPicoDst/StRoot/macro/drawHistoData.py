import ROOT
ROOT.gROOT.SetBatch(True)

import rootpy.ROOT as ROOT

from rootpy.tree import Cut
from rootpy.tree.categories import Categories
from rootpy.io import root_open
from rootpy.plotting import Hist, Hist2D, Hist3D, histogram, Canvas, Legend

import math

from ROOT import TTree, TFile, AddressOf, gROOT


# load input files
input = 'Gamma'
if input=='Pi0':
    infileSim = root_open('root/out_pi0_15.root')
    #infileSim = root_open('test.pi0.root')
    infileWeight = root_open('/Users/kunsu/auau200GeVRun10/Efficiency/PartnerFinding/Decay/Out/Pion_STAR.root')
    wt = infileWeight.fitfun_pt_MB
    gid = 10007


if input=='Gamma':
    #infileSim = root_open('root/out_gamma_15.root')
    infileSim = root_open('root/out_gamma_SL15k.root')
    infileWeight = root_open('/Users/kunsu/auau200GeVRun10/Efficiency/PartnerFinding/Decay/Out/pi02gamma_dir_3.root')
    wt = infileWeight.funDirGamma_MB
    gid = 1

tree = infileSim.eT


canvas = Canvas()
#canvas.SetLogy()
def drawRatio(xval, xmin, xmax, xvalname,
              xbin=100,ybin=100,
              #  cut='rchfthit.pxl1 && rchfthit.pxl2 && rchfthit.ist ',
              #  cut2=' && truth.pxl1 && truth.pxl2 && truth.ist',
              cut='rcHftHit_pxl1 && rcHftHit_pxl2 && rcHftHit_ist && sqrt(mc_x**2+mc_y**2) > 1.9 && sqrt(mc_x**2+mc_y**2) < 2.1 && parentGid ==1',
              cut2='',
              histoname='histoRatio', drawOption='E0',inTree=tree):
    histo1 = inTree.Draw(xval, hist=Hist(xbin,xmin,xmax), selection=cut)
    histo1.color = 'red'
    histo1.SetMarkerStyle(24)
    histo2 = inTree.Draw(xval, hist=Hist(xbin,xmin,xmax), selection=cut + cut2)
    histo2.color = 'blue'
    histo2.Divide(histo1)
    histo2.SetMaximum(1)
    histo2.Draw(drawOption)
    canvas.SaveAs('Eff/'+histoname+'_'+xvalname+'.pdf')

drawRatio('rcPt',0,10,'rcPt',cut2='&& truth_pxl1 && truth_pxl2 && truth_ist',histoname='histoRatio'+input+'BeamPipe');
drawRatio('rcPt',0,10,'rcPt',cut2='&& truth_pxl1 && truth_pxl2 && truth_ist',
          cut='rcHftHit_pxl1 && rcHftHit_pxl2 && rcHftHit_ist && sqrt(mc_x**2+mc_y**2) < 1.9 && parentGid ==1',
          histoname='histoRatio'+input+'Inner');
drawRatio('rcPt',0,10,'rcPt',cut2='&& truth_pxl1 && truth_pxl2 && truth_ist',
          cut='rcHftHit_pxl1 && rcHftHit_pxl2 && rcHftHit_ist && sqrt(mc_x**2+mc_y**2) >2.1 && parentGid ==1',
          histoname='histoRatio'+input+'Outer');



tree = infileSim.T


# define histograms
hRcPt0 = ROOT.TH1F("hRcPt0","hRcPt0",100, 0, 10) 
hRcPt1 = ROOT.TH1F("hRcPt1","hRcPt1",100, 0, 10) 
hRcPt2 = ROOT.TH1F("hRcPt2","hRcPt2",100, 0, 10) 
hRcPt3 = ROOT.TH1F("hRcPt3","hRcPt3",100, 0, 10) 
hRcPt4 = ROOT.TH1F("hRcPt4","hRcPt4",100, 0, 10) 
hRcPt5 = ROOT.TH1F("hRcPt5","hRcPt5",100, 0, 10)
hConvR0 = ROOT.TH1F("hConvR0","hConvR0",100, 0, 10)
hConvR1 = ROOT.TH1F("hConvR1","hConvR1",100, 0, 10)
hConvR2 = ROOT.TH1F("hConvR2","hConvR2",100, 0, 10)
hConvR3 = ROOT.TH1F("hConvR3","hConvR3",100, 0, 10)
hConvR4 = ROOT.TH1F("hConvR4","hConvR4",100, 0, 10)
hConvR5 = ROOT.TH1F("hConvR5","hConvR5",100, 0, 10)

hRcPt0Wt = ROOT.TH1F("hRcPt0Wt","hRcPt0Wt",100, 0, 10) # pt > 0.6
hRcPt1Wt = ROOT.TH1F("hRcPt1Wt","hRcPt1Wt",100, 0, 10) # 0.6 < pt < 0.8
hRcPt2Wt = ROOT.TH1F("hRcPt2Wt","hRcPt2Wt",100, 0, 10) # 0.8 < pt < 1.0
hRcPt3Wt = ROOT.TH1F("hRcPt3Wt","hRcPt3Wt",100, 0, 10) # 1.0 < pt < 1.5
hRcPt4Wt = ROOT.TH1F("hRcPt4Wt","hRcPt4Wt",100, 0, 10) # 1.5 < pt < 2.0
hRcPt5Wt = ROOT.TH1F("hRcPt5Wt","hRcPt5Wt",100, 0, 10) # pt > 2.0
hConvR0Wt = ROOT.TH1F("hConvR0Wt","hConvR0Wt",100, 0, 10)
hConvR1Wt = ROOT.TH1F("hConvR1Wt","hConvR1Wt",100, 0, 10)
hConvR2Wt = ROOT.TH1F("hConvR2Wt","hConvR2Wt",100, 0, 10)
hConvR3Wt = ROOT.TH1F("hConvR3Wt","hConvR3Wt",100, 0, 10)
hConvR4Wt = ROOT.TH1F("hConvR4Wt","hConvR4Wt",100, 0, 10)
hConvR5Wt = ROOT.TH1F("hConvR5Wt","hConvR5Wt",100, 0, 10)


hMomPt = ROOT.TH1F("hMomPt","hMomPt",100, 0, 10)
hMomPtWt = ROOT.TH1F("hMomPtWt","hMomPtWt",100, 0, 10)

# pair loop
for pair in tree:
    if pair.parentGid == gid: # and pair.mcPairPt > 0.2 and pair.mcPairPt < 10:
        momPt = pair.mcPairPt
        pt1 =  abs(pair.pt1)
        pt2 =  abs(pair.pt2)
        rcConvR = math.sqrt(pair.rc_x**2 + pair.rc_y**2)
        hft1 = pair.rcHftHit1_pxl1 and pair.rcHftHit1_pxl2 and pair.rcHftHit1_ist
        hft2 = pair.rcHftHit2_pxl1 and pair.rcHftHit2_pxl2 and pair.rcHftHit2_ist
        #cut = pair.truth1_pxl1 and pair.truth2_pxl1 and pair.truth1_pxl2 and pair.truth2_pxl2
        #cut = (pair.truth1_pxl1 or pair.truth2_pxl1) and pair.truth1_pxl2 and pair.truth2_pxl2
        cut = pair.mass < 0.05 and pair.pairDca < 0.01
        if pt1 > 0.6 and pt2 > 0.6 and hft1 and hft2 and cut:
            hRcPt0.Fill(pt1)
            hRcPt0.Fill(pt2)
            hRcPt0Wt.Fill(pt1, wt.Eval(momPt))
            hRcPt0Wt.Fill(pt2, wt.Eval(momPt))
            hConvR0.Fill(rcConvR)
            hConvR0Wt.Fill(rcConvR, wt.Eval(momPt))
        
        if pt1 > 0.6 and pt1 < 0.8 and pt2 > 0.6 and pt2 < 0.8 and hft1 and hft2 and cut:
            hRcPt1.Fill(pt1)
            hRcPt1.Fill(pt2)
            hRcPt1Wt.Fill(pt1, wt.Eval(momPt))
            hRcPt1Wt.Fill(pt2, wt.Eval(momPt))
            hConvR1.Fill(rcConvR)
            hConvR1Wt.Fill(rcConvR, wt.Eval(momPt))
        
        if pt1 > 0.8 and pt1 < 1 and pt2 > 0.8 and pt2 < 1 and hft1 and hft2 and cut:
            hRcPt2.Fill(pt1)
            hRcPt2.Fill(pt2)
            hRcPt2Wt.Fill(pt1, wt.Eval(momPt))
            hRcPt2Wt.Fill(pt2, wt.Eval(momPt))
            hConvR2.Fill(rcConvR)
            hConvR2Wt.Fill(rcConvR, wt.Eval(momPt))
        
        if pt1 > 1 and pt1 < 1.5 and pt2 > 1 and pt2 < 1.5 and hft1 and hft2 and cut:
            hRcPt3.Fill(pt1)
            hRcPt3.Fill(pt2)
            hRcPt3Wt.Fill(pt1, wt.Eval(momPt))
            hRcPt3Wt.Fill(pt2, wt.Eval(momPt))
            hConvR3.Fill(rcConvR)
            hConvR3Wt.Fill(rcConvR, wt.Eval(momPt))
        
        if pt1 > 1.5 and pt1 < 2. and pt2 > 1.5 and pt2 < 2. and hft1 and hft2 and cut:
            hRcPt4.Fill(pt1)
            hRcPt4.Fill(pt2)
            hRcPt4Wt.Fill(pt1, wt.Eval(momPt))
            hRcPt4Wt.Fill(pt2, wt.Eval(momPt))
            hConvR4.Fill(rcConvR)
            hConvR4Wt.Fill(rcConvR, wt.Eval(momPt))
        
        if pt1 > 2 and pt2 > 2 and hft1 and hft2 and cut:
            hRcPt5.Fill(pt1)
            hRcPt5.Fill(pt2)
            hRcPt5Wt.Fill(pt1, wt.Eval(momPt))
            hRcPt5Wt.Fill(pt2, wt.Eval(momPt))
            hConvR5.Fill(rcConvR)
            hConvR5Wt.Fill(rcConvR, wt.Eval(momPt))
        
        hMomPt.Fill(momPt)
        hMomPtWt.Fill(momPt, wt.Eval(momPt))

fout = ROOT.TFile('histo_All'+input+'Cut.root','RECREATE')
fout.cd()
hRcPt0.Write()
hRcPt1.Write()
hRcPt2.Write()
hRcPt3.Write()
hRcPt4.Write()
hRcPt5.Write()
hConvR0.Write()
hConvR1.Write()
hConvR2.Write()
hConvR3.Write()
hConvR4.Write()
hConvR5.Write()
hRcPt0Wt.Write()
hRcPt1Wt.Write()
hRcPt2Wt.Write()
hRcPt3Wt.Write()
hRcPt4Wt.Write()
hRcPt5Wt.Write()
hConvR0Wt.Write()
hConvR1Wt.Write()
hConvR2Wt.Write()
hConvR3Wt.Write()
hConvR4Wt.Write()
hConvR5Wt.Write()

hMomPt.Write()
hMomPtWt.Write()
wt.Write()




'''


canvas = Canvas()
#canvas.SetLogy()


def drawRatio(xval, xmin, xmax, xvalname,
              xbin=100,ybin=100,
              #  cut='rchfthit.pxl1 && rchfthit.pxl2 && rchfthit.ist ',
              #  cut2=' && truth.pxl1 && truth.pxl2 && truth.ist',
              cut='',
              cut2='',
              histoname='histoRatio', drawOption='E0',inTree=tree):
    histo1 = inTree.Draw(xval, hist=Hist(xbin,xmin,xmax), selection=cut)
    histo1.color = 'red'
    histo1.SetMarkerStyle(24)
    histo2 = inTree.Draw(xval, hist=Hist(xbin,xmin,xmax), selection=cut + cut2)
    histo2.color = 'blue'
    histo2.Divide(histo1)
    histo2.SetMaximum(1)
    histo2.SetMinimum(0)
    histo2.Draw(drawOption)
    canvas.SaveAs('Eff/'+histoname+'_'+xvalname+'.pdf')

drawRatio('rcPt',0,5,'rcPt',
          cut='rchfthit.pxl1 && rchfthit.pxl2 && rchfthit.ist && parentGid == 1 && sqrt(mc.x**2+mc.y**2) < 2.1 && sqrt(mc.x**2+mc.y**2) > 1.9',
          cut2='&& truth.pxl1 && truth.pxl2 && truth.ist',
          histoname='histoRatio'+input+'BeamPipe');

drawRatio('rcPt',0,5,'rcPt',
          cut='rchfthit.pxl1 && rchfthit.pxl2 && rchfthit.ist && parentGid == 1 && sqrt(mc.x**2+mc.y**2) < 1.9',
          cut2='&& truth.pxl1 && truth.pxl2 && truth.ist',
          histoname='histoRatio'+input+'InnerBeamPipe');

drawRatio('rcPt',0,5,'rcPt',
          cut='rchfthit.pxl1 && rchfthit.pxl2 && rchfthit.ist && parentGid == 1 && sqrt(mc.x**2+mc.y**2) > 2.1',
          cut2='&& truth.pxl1 && truth.pxl2 && truth.ist',
          histoname='histoRatio'+input+'OuterBeamPipe');

drawRatio('rcPt',0,5,'rcPt',
          cut='rchfthit.pxl1 && rchfthit.pxl2 && rchfthit.ist && parentGid == 10007 && sqrt(mc.x**2+mc.y**2) < .1',
          cut2='&& truth.pxl1 && truth.pxl2 && truth.ist',
          histoname='histoRatioPi0');



def draw1D(xval='sqrt(rc.x**2+rc.y**2)', xmin=0, xmax=10, xvalname='Radius',
           xbin=100,ybin=100,
           #  cut='rchfthit.pxl1 && rchfthit.pxl2 && rchfthit.ist ',
           #  cut2=' && truth.pxl1 && truth.pxl2 && truth.ist',
           cut='rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit1.ist && rchfthit2.pxl1 && rchfthit2.pxl2 && rchfthit2.ist && parentGid ==1',
           cut2='',
           histoname='histoRatio', drawOption='E0',inTree=tree,ndata=1):
    inData = root_open('outHist.root')
    histoData = inData.Get('h'+str(ndata))
    histoData.Sumw2()
    histoDataTemp = inData.Get('hLS'+str(ndata))
    histoDataTemp.Sumw2()
    histoData.Add(histoDataTemp,-1)
    histoData.scale(1/histoData.GetMaximum())
    histo2 = inTree.Draw(xval, hist=Hist(xbin,xmin,xmax), selection=cut + cut2)
    histo2.color = 'blue'
    histo2.scale(1/histo2.GetMaximum())
    histoData.SetMaximum(10)
    histoData.SetMinimum(0.001)
    histoData.Draw(drawOption)
    histo2.Draw(drawOption+'same')
    canvas.SaveAs('Eff/'+histoname+'_'+xvalname+str(ndata)+'.pdf')


def draw1DBoth(xval='sqrt(rc.x**2+rc.y**2)', xmin=0, xmax=10, xvalname='RadiusBoth',
           xbin=100,ybin=100,
           #  cut='rchfthit.pxl1 && rchfthit.pxl2 && rchfthit.ist ',
           #  cut2=' && truth.pxl1 && truth.pxl2 && truth.ist',
           cut='rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit1.ist && rchfthit2.pxl1 && rchfthit2.pxl2 && rchfthit2.ist',
           cut2='',
           histoname='histoRatio', drawOption='E0',inTree=tree,inTree2=tree2,ndata=1):
    inData = root_open('outHist.root')
    histoData = inData.Get('h'+str(ndata))
    histoData.Sumw2()
    histoDataTemp = inData.Get('hLS'+str(ndata))
    histoDataTemp.Sumw2()
    histoData.Add(histoDataTemp,-1)
    histoData.scale(1/histoData.GetMaximum())
    histo2 = inTree.Draw(xval, hist=Hist(xbin,xmin,xmax), selection=cut + cut2)
    histo2.color = 'blue'
    histo2.scale(1/histo2.GetMaximum())
    histo3 = inTree2.Draw(xval, hist=Hist(xbin,xmin,xmax), selection=cut + cut2)
    histo3.color = 'red'
    histo3.scale(1/histo3.GetMaximum())
    histoData.SetMaximum(10)
    histoData.SetMinimum(0.001)
    histoData.Draw(drawOption)
    histo2.Draw(drawOption+'same')
    histo3.Draw(drawOption+'same')
    canvas.SaveAs('Eff/'+histoname+'_'+xvalname+str(ndata)+'.pdf')


draw1DBoth(cut2='&& abs(pt1) > 0.6 && abs(pt2) > 0.6',ndata = 0)

draw1D(cut2='&& abs(pt1) > 0.6 && abs(pt2) > 0.6',ndata = 0)
draw1D(cut2='&& abs(pt1) > 0.6 && abs(pt2) > 0.6 && abs(pt1) < 0.8 && abs(pt2) < 0.8',ndata = 1)
draw1D(cut2='&& abs(pt1) > 0.8 && abs(pt2) > 0.6 && abs(pt1) < 1. && abs(pt2) < 1.',ndata = 2)
draw1D(cut2='&& abs(pt1) > 1. && abs(pt2) > 1. && abs(pt1) < 1.5 && abs(pt2) < 1.5',ndata = 3)
draw1D(cut2='&& abs(pt1) > 1.5 && abs(pt2) > 1.5 && abs(pt1) < 2. && abs(pt2) < 2.',ndata = 4)
draw1D(cut2='&& abs(pt1) > 2. && abs(pt2) > 2.',ndata = 5)

draw1D(cut2='&& abs(pt1) > 0.6 && abs(pt2) > 0.6 && rc.z < -3',ndata = 6)
draw1D(cut2='&& abs(pt1) > 0.6 && abs(pt2) > 0.6 && rc.z > -3 && rc.z < 0',ndata = 7)
draw1D(cut2='&& abs(pt1) > 0.6 && abs(pt2) > 0.6 && rc.z > 0 && rc.z < 3',ndata = 8)
draw1D(cut2='&& abs(pt1) > 0.6 && abs(pt2) > 0.6 && rc.z > 3',ndata = 9)


    # Data
    infile = root_open('../out_gamma10_5.root')
    ntuple = infile.nt
    xval = 'sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)'
    xbin = 100
    xmin = 0
    xmax = 10
    cut='mass < 0.05 && pairDca < 0.01 && abs(pt1) > 0.6 && abs(pt2) > 0.6'
    histoUS = ntuple.Draw(xval, hist=Hist(xbin,xmin,xmax), selection=cut+'&& pt1*pt2<0')
    histoLS = ntuple.Draw(xval, hist=Hist(xbin,xmin,xmax), selection=cut+'&& pt1*pt2>0')
    histoData = histoUS - histoLS                                                               # data UL-LS
    histoData.scale(1/histoData.GetMaximum())
    histoData.color = 'black'
    histoData.SetMarkerStyle(20)
    
    
    # Simulation
def drawRatio(xval, xmin, xmax, xvalname,
              xbin=100,ybin=100,
              #  cut='rchfthit.pxl1 && rchfthit.pxl2 && rchfthit.ist ',
              #  cut2=' && truth.pxl1 && truth.pxl2 && truth.ist',
              cut='rchfthit.pxl1 && rchfthit.pxl2 && rchfthit.ist && sqrt(mc.x**2+mc.y**2) > 1.8 && parentGid ==1',
              cut2='',
              histoname='histoRatio', drawOption='E0',inTree=tree):
    histo1 = inTree.Draw(xval, hist=Hist(xbin,xmin,xmax), selection=cut)
    histo1.color = 'red'
    histo1.SetMarkerStyle(24)
    histo2 = inTree.Draw(xval, hist=Hist(xbin,xmin,xmax), selection=cut + cut2)
    histo2.color = 'blue'
    histo2.Divide(histo1)
    histo2.SetMaximum(1)
    histo2.Draw(drawOption)
    canvas.SaveAs('Eff/'+histoname+'_'+xvalname+'.pdf')

drawRatio('rcPt',0,10,'rcPt',cut2='&& truth.pxl1 && truth.pxl2 && truth.ist',histoname='histoRatio'+input+'3');

tree = infile2.T
drawRatio('pt1',0,10,'pt1',cut='rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit1.ist && sqrt(mc1.x**2+mc1.y**2) > 1.8',cut2='&& truth1.pxl1 && truth1.pxl2 && truth1.ist',histoname='histoRatioPair'+input+'3',inTree=tree);



ntuple2 = infile2.T
def draw2D(
xval, xmin, xmax, xvalname,
yval='sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)', ymin=0, ymax=10, yvalname='ConvR',
xbin=100,ybin=300,
cut='mass < 0.1 && pairDca < 0.01 && pt1*pt2<0 && abs(pt1) > 0.8 && abs(pt2) > 0.8',
histoname='histoData'
):
histo1 = ntuple.Draw(xval+':'+yval, hist=Hist2D(xbin,xmin,xmax,ybin,ymin,ymax), selection=cut)
histo1.Draw('col2');
histo1.xaxis.SetTitle(xvalname)
histo1.yaxis.SetTitle(yvalname)
canvas.SaveAs(histoname+'_'+yvalname+'_'+xvalname+'.pdf')

def draw1D(
xval, xmin, xmax, xvalname, xbin=100,
cut='mass < 0.05 && pairDca < 0.01 && nsige1>-1 && nsige2>-1 && abs(pt1) > 0.6 && abs(pt2) > 0.6',
histoname='histoPi0Dalitz',
cutSim='abs(pt1) > 0.6 && abs(pt2) > 0.6', cutSim2 = '',
xvalSim='',xvalSim2=''
):
histo1 = ntuple.Draw(xval, hist=Hist(xbin,xmin,xmax), selection=cut+' && pt1*pt2<0') # data UL
histo1.color = 'red'
histo1.SetMarkerStyle(24)
histo2 = ntuple.Draw(xval, hist=Hist(xbin,xmin,xmax), selection=cut+' && pt1*pt2>0') # data LS
#histo2.scale(histo1.integral(90,99)/histo2.integral(90,99))
histo2.color = 'blue'
histo2.SetMarkerStyle(24)
histo3 = histo1-histo2                                                               # data UL-LS
histo3.scale(1/histo3.GetMaximum())
histo3.color = 'black'
histo3.SetMarkerStyle(20)


if xvalSim == '':
xvalSim = xval
histo4 = ntuple2.Draw(xvalSim, hist=Hist(xbin,xmin,xmax), selection=cutSim + cutSim2)          # sim #1
histo4.scale(1/histo4.GetMaximum())
histo4.color = 'red'
histo4.scale(histo3.GetMaximum()/histo4.GetMaximum())
if xvalSim2 == '':
xvalSim2 = xval
histo5 = ntuple2.Draw(xvalSim2, hist=Hist(xbin,xmin,xmax), selection=cutSim + cutSim2)         # sim #2
histo5.scale(1/histo5.GetMaximum())
histo5.color = 'blue'
histo5.scale(histo3.GetMaximum()/histo4.GetMaximum())

histo3.SetMinimum(1e-4)
#histo1.Draw('E0')
#histo2.Draw('sameE0')
histo3.Draw('E0')
canvas.SaveAs('data.pdf')
canvas.SaveAs('data.root')

histo4.Draw('sameE0')
histo5.Draw('sameE0')
canvas.SaveAs('4bitsLog1GeV'+histoname+'_'+xvalname+'.pdf')



draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR0000',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2'+
       ' && truth1.pxl1==0 && truth1.pxl2==0 && truth2.pxl1==0 && truth2.pxl2==0')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR0001',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2'+
       ' && truth1.pxl1==0 && truth1.pxl2==0 && truth2.pxl1==0 && truth2.pxl2')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR0010',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2'+
       ' && truth1.pxl1==0 && truth1.pxl2==0 && truth2.pxl1 && truth2.pxl2==0')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR0011',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2'+
       ' && truth1.pxl1==0 && truth1.pxl2==0 && truth2.pxl1 && truth2.pxl2')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR0100',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2'+
       ' && truth1.pxl1==0 && truth1.pxl2 && truth2.pxl1==0 && truth2.pxl2==0')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR0101',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2'+
       ' && truth1.pxl1==0 && truth1.pxl2 && truth2.pxl1==0 && truth2.pxl2')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR0110',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2'+
       ' && truth1.pxl1==0 && truth1.pxl2 && truth2.pxl1 && truth2.pxl2==0')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR0111',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2'+
       ' && truth1.pxl1==0 && truth1.pxl2 && truth2.pxl1 && truth2.pxl2')

draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR1000',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2'+
       ' && truth1.pxl1 && truth1.pxl2==0 && truth2.pxl1==0 && truth2.pxl2==0')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR1001',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2'+
       ' && truth1.pxl1 && truth1.pxl2==0 && truth2.pxl1==0 && truth2.pxl2')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR1010',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2'+
       ' && truth1.pxl1 && truth1.pxl2==0 && truth2.pxl1 && truth2.pxl2==0')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR1011',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2'+
       ' && truth1.pxl1 && truth1.pxl2==0 && truth2.pxl1 && truth2.pxl2')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR1100',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2'+
       ' && truth1.pxl1 && truth1.pxl2 && truth2.pxl1==0 && truth2.pxl2==0')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR1101',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2'+
       ' && truth1.pxl1 && truth1.pxl2 && truth2.pxl1==0 && truth2.pxl2')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR1110',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2'+
       ' && truth1.pxl1 && truth1.pxl2 && truth2.pxl1 && truth2.pxl2==0')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR1111',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2'+
       ' && truth1.pxl1 && truth1.pxl2 && truth2.pxl1 && truth2.pxl2')






draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR3_1',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
    cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2 && rchfthit1.ist && rchfthit2.ist && (truth1.pxl2==0 || truth2.pxl2==0) && (truth1.pxl1 && truth2.pxl1)')


draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR12',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2 && (truth1.pxl2 || truth2.pxl2) && !(truth1.pxl2==0 && truth2.pxl2==0) && truth1.pxl1 && truth2.pxl1')


draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2 && truth1.pxl1 && truth1.pxl2 && truth2.pxl1 && truth2.pxl2')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR1',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2 && (truth1.pxl1==0 || truth1.pxl2==0 || truth2.pxl1==0 || truth2.pxl2==0)')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR2',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2 && (truth1.pxl1==0 || truth2.pxl1==0) && (truth1.pxl2 && truth2.pxl2)')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR3',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2 && (truth1.pxl2==0 || truth2.pxl2==0) && (truth1.pxl1 && truth2.pxl1)')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR4',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR9',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2 && truth1.pxl1 && truth2.pxl1 && truth1.pxl2==0 && truth2.pxl2==0')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR10',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2 && truth1.pxl1==0 && truth2.pxl1==0 && truth1.pxl2 && truth2.pxl2')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR11',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && rchfthit1.pxl1 && rchfthit1.pxl2 && rchfthit2.pxl1 && rchfthit2.pxl2 && (truth1.pxl1 || truth2.pxl1) && !(truth1.pxl1==0 && truth2.pxl1==0) && truth1.pxl2 && truth2.pxl2')


# nHits
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR5',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && nHits1.pxl1 && nHits1.pxl2 && nHits2.pxl1 && nHits2.pxl2 && (truth1.pxl1==0 || truth1.pxl2==0 || truth2.pxl1==0 || truth2.pxl2==0)')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR6',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && nHits1.pxl1 && nHits1.pxl2 && nHits2.pxl1 && nHits2.pxl2 && (truth1.pxl1==0 || truth2.pxl1==0)')
draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR7',100,xvalSim='sqrt((mc.x)**2+(mc.y)**2)',xvalSim2='sqrt((rc.x)**2+(rc.y)**2)',
       cutSim2 =' && nHits1.pxl1 && nHits1.pxl2 && nHits2.pxl1 && nHits2.pxl2 && truth1.pxl1 && truth1.pxl2 && truth2.pxl1 && truth2.pxl2')



draw1D('mass',0, 0.1,'Mass',100)
'''
'''
draw2D('phi',-4,4,'Phi')
draw2D('eta',-1.1,1.1,'Eta')
draw2D('theta',0,0.05,'Theta')
draw2D('mass',-0.01,0.09,'Mass')
draw2D('abs(pt1)',0,20,'Pt1')
draw2D('abs(pt2)',0,20,'Pt2')
draw2D('v0z',-10,10,'V0z')
draw2D('pairDca',0,0.01,'pairDca')
draw2D('mass',-0.01,0.19,'Mass',
       'theta',0,0.4,'Theta',
       xbin=100)
draw2D('pairDca',0,0.01,'pairDca',
       'theta',0,0.05,'Theta',
       xbin=100)
draw2D('phi',-4,4,'Phi',
       'eta',-1.1,1.1,'Eta',
       xbin=100,ybin=100)
draw2D('abs(pt1)',0,20,'Pt1',
       'abs(pt2)',0,20,'Pt2',
       xbin=100,ybin=100)
draw2D('v0x',-10,10,'V0x',
       'v0y',-10,10,'V0y',
       xbin=400,ybin=400)
'''

