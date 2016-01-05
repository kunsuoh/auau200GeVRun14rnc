import ROOT
ROOT.gROOT.SetBatch(True)

import rootpy.ROOT as ROOT

from rootpy.tree import Cut
from rootpy.tree.categories import Categories
from rootpy.io import root_open
from rootpy.plotting import Hist, Hist2D, Hist3D, histogram, Canvas, Legend

# Simulation
infile = root_open('out_gamma_5.root')
#infile = root_open('out_k0s_2.root')
T = infile.T

canvas = Canvas()
#canvas.SetLogz()

def draw2D(
           xval, xmin, xmax, xvalname,
           yval='sqrt(rc.x**2+rc.y**2)', ymin=0, ymax=10, yvalname='ConvR',
           xbin=100,ybin=300,
           cut='',
           histoname='histo',xTitle='',yTitle=''
           ):
    histo1 = T.Draw(xval+':'+yval, hist=Hist2D(xbin,xmin,xmax,ybin,ymin,ymax), selection=cut)
    histo1.Draw('col2');
    if xTitle=='' :
        histo1.xaxis.SetTitle(xvalname)
    else :
        histo1.xaxis.SetTitle(xTitle)
    if yTitle=='' :
        histo1.yaxis.SetTitle(yvalname)
    else :
        histo1.yaxis.SetTitle(yTitle)
    canvas.SaveAs(histoname+'_'+yvalname+'_'+xvalname+'.pdf')

def draw1D(
           xval, xmin, xmax, xvalname, xbin=100,
           cut='abs(pt1) > 0.6 && abs(pt2) > 0.6 &&sqrt(rc.x**2+rc.y**2) > 2 && sqrt(rc.x**2+rc.y**2) < 4',
           histoname='histo1D'
           ):
    histo1 = T.Draw(xval, hist=Hist(xbin,xmin,xmax), selection=cut)
    histo1.color = 'red'
    histo1.Draw('E0')
    histo1.xaxis.SetTitle(xvalname)
    canvas.SaveAs(histoname+'_'+xvalname+'.pdf')

def draw1DCompare(
                  xval, xmin, xmax, xvalname, xbin=100,
                  cut='abs(pt1) > 0.2 && abs(pt2) > 0.2',
                  histoname='histo1DCompare_gamma'
                  ):
    canvas1 = Canvas()
    histo1 = T.Draw(xval, hist=Hist(xbin,xmin,xmax), selection=cut + '&& sqrt(mc.x**2+mc.y**2) < 2.1')
    histo1.scale(1/histo1.integral())
    histo1.color = 'red'
    histo2 = T.Draw(xval, hist=Hist(xbin,xmin,xmax), selection=cut + '&& sqrt(mc.x**2+mc.y**2) > 2.1 && sqrt(mc.x**2+mc.y**2) < 4')
    histo2.scale(1/histo2.integral())
    histo2.color = 'blue'
    histo3 = T.Draw(xval, hist=Hist(xbin,xmin,xmax), selection=cut + '&& sqrt(mc.x**2+mc.y**2) > 4')
    histo3.scale(1/histo3.integral())
    histo3.color = 'black'
    
    histo1.SetMaximum(histo1.GetMaximum()*1.5)
    histo1.xaxis.SetTitle(xvalname)
    histo1.Draw('E0')
    histo2.Draw('E0same')
    histo3.Draw('E0same')
    
    
    
    canvas1.SaveAs(histoname+'_'+xvalname+'.pdf')


draw1DCompare('mass',0., 0.1,'Mass',100)
draw1DCompare('sqrt(rc.x**2+rc.y**2)',0, 10,'ConvR',100)
draw1DCompare('sqrt(mc.x**2+mc.y**2)',0, 10,'mcConvR',100)
draw1DCompare('pairDca',0, 0.01,'pairDca',100)

draw1DCompare('length',0, .03,'length',100)
draw1DCompare('angle',0, .01,'angle',100)
draw1DCompare('openangle',0, 0.02,'OpenAngle',100)
draw1DCompare('phiV',0, 0.8,'phiV',100)


'''
    draw2D('sqrt(v0x**2+v0y**2)',0,10,'RcConvR',
    'sqrt(mcv0x**2+mcv0y**2)',0,10,'McTruthConvR',
    xTitle='Rc Conversion Radius [cm]',
    yTitle='Mc Truth Conversion Radius [cm]',
    xbin=100,ybin=100,
    histoname='histoPi0',
    cut='')
    draw2D('sqrt(v0x**2+v0y**2)',0,10,'RcConvR',
    'mass',0,1,'Mas',
    xTitle='Rc Conversion Radius [cm]',
    yTitle='',
    xbin=100,ybin=100,
    histoname='histoPi0',
    cut='')
    draw1D('mass',0, 0.05,'Mass',100,cut='')
    draw1D('mass',0, 0.1,'Mass',100)
    draw1D('sqrt(v0x**2+v0y**2)-sqrt(mcv0x**2+mcv0y**2)',-10,10,'ConvR_{Rc-Mc}')
    draw1D('sqrt((v0x-mcv0x)**2+(v0y-mcv0y)**2)',0,10,'positionXY_{Rc-Mc}',cut='sqrt(v0x**2+v0y**2)>2.1')
    
    draw2D('phi',-4,4,'Phi')
    draw2D('eta',-1.1,1.1,'Eta')
    draw2D('phiV',0,0.8,'phiV')
    draw2D('mass',0.0,0.2,'Mass')
    draw2D('abs(pt1)',0,20,'Pt1')
    draw2D('abs(pt2)',0,20,'Pt2')
    draw2D('v0z',-10,10,'V0z')
    draw2D('mcPairPt',0,25,'PairPt')
    draw2D('pairDca',0,0.01,'pairDca')
    draw2D('mass',0,0.2,'Mass',
    'phiV',0,0.8,'phiV',
    xbin=100)
    draw2D('mass',0,0.2,'Mass',
    'mcPairPt',0,20,'mcPairPt',
    xbin=100)
    draw2D('pairDca',0,0.01,'pairDca',
    'phiV',0,0.8,'phiV',
    xbin=100)
    draw2D('phi',-4,4,'Phi',
    'eta',-1.1,1.1,'Eta',
    xbin=100,ybin=100)
    draw2D('abs(pt1)',0,20,'Pt1',
    'abs(pt2)',0,20,'Pt2',
    xbin=100,ybin=100)
    draw2D('abs(mcPairPt)',0,20,'PairPt',
    'abs(pt1)',0,20,'Pt1',
    xbin=100,ybin=100)
    draw2D('abs(mcPairPt)',0,20,'PairPt',
    'abs(pt2)',0,20,'Pt2',
    xbin=100,ybin=100)
    draw2D('v0x',-10,10,'V0x',
    'v0y',-10,10,'V0y',
    xbin=400,ybin=400,
    histoname='histo',
    cut='mass<0.1 && pairDca<0.01')
    draw2D('sqrt(mcv0x**2+mcv0y**2)',0,10,'McConvR',
    'mass', 0, 0.1, 'mass',
    xbin=500,ybin=100)
    draw2D('sqrt(v0x**2+v0y**2)',0,10,'ConvR',
    'mass', 0, 0.1, 'mass',
    xbin=200,ybin=100)
    
    draw2D('v0x',-10,10,'V0x',
    'v0y',-10,10,'V0y',
    xbin=400,ybin=400,
    cut='sqrt(mcv0x**2+mcv0y**2) < 2.1 && pairDca < 0.001',
    histoname='histo_PairDca0.001')
    draw2D('v0x',-10,10,'V0x',
    'v0y',-10,10,'V0y',
    xbin=400,ybin=400,
    cut='sqrt(mcv0x**2+mcv0y**2) < 2.1 && pairDca < 0.001 && mass > 0.01',
    histoname='histo_PairDca0.001Mass0.01')
    
    draw2D('v0x',-10,10,'V0x',
    'v0y',-10,10,'V0y',
    xbin=400,ybin=400,
    cut='sqrt(mcv0x**2+mcv0y**2) < 2.1 && theta>0.01',
    histoname='histo_Thteta0.01')
    draw2D('mcv0x',-10,10,'MCV0x',
    'mcv0y',-10,10,'MCV0y',
    xbin=400,ybin=400,
    histoname='histo')
    
    nhisto = 0
    categories = Categories.from_string('{abs(pt1)|1.5,1.6,1.8,2.0,2.2,2.5,3.0,4.0,5.0,7.0}')
    #categories = Categories.from_string('{sqrt(mcv0x**2+mcv0y**2)|2.1}x{v0z|-6,0,6}x{abs(pt1)|1.5,2.5,3.5,4.5,5.5,6.5}')
    for cut in categories:
    cut2 = cut & "mass < 0.01 && pairDca < 0.01 && pt1*pt2<0 "
    cutSim = cut2 & "sqrt(mcv0x**2+mcv0y**2) >.1"
    
    canvas = Canvas()
    canvas.SetLogy()
    h2d = ntuple.Draw('sqrt(rc.x**2+rc.y**2)', hist=Hist(100, 0, 10), selection=cutSim)
    h2d1 = ntuple.Draw('sqrt(rc.x**2+rc.y**2)', hist=Hist(100, 0, 10), selection=cutSim&"sqrt(mcv0x**2+mcv0y**2) < 2.1")
    h2d2 = ntuple.Draw('sqrt(rc.x**2+rc.y**2)', hist=Hist(100, 0, 10), selection=cutSim&"sqrt(mcv0x**2+mcv0y**2) > 2.1")
    
    h2d.SetTitle(str(cut))
    h2d.xaxis.SetTitle("Conversion radius [cm]")
    h2d1.SetMarkerColor(2)
    h2d2.SetMarkerColor(4)
    h2d1.SetMarkerStyle(24)
    h2d2.SetMarkerStyle(24)
    
    h2d.SetMaximum(h2d.GetMaximum()*20)
    
    h2d.Draw('')
    h2d1.Draw('same')
    h2d2.Draw('same')
    
    legend = Legend(4, leftmargin=0.35, margin=0.3)
    legend.AddEntry(h2d, style='LEP', label='Simulation')
    legend.AddEntry(h2d1, style='LEP', label='  only beam pipe')
    legend.AddEntry(h2d2, style='LEP', label='  others')
    legend.Draw()
    
    label = ROOT.TText(0.01, 0.01, cut.str)
    label.SetTextFont(43)
    label.SetTextSize(6)
    label.SetNDC()
    label.Draw()
    canvas.SaveAs('histSimConvR_' + str(nhisto) + '.pdf')
    
    
    canvas2 = Canvas()
    canvas2.SetLogy()
    
    hTheta = ntuple.Draw('theta', hist=Hist(100, 0, 0.1), selection=cutSim)
    hTheta1 = ntuple.Draw('theta', hist=Hist(100, 0, 0.1), selection=cutSim&"sqrt(mcv0x**2+mcv0y**2) < 2.1")
    hTheta2 = ntuple.Draw('theta', hist=Hist(100, 0, 0.1), selection=cutSim&"sqrt(mcv0x**2+mcv0y**2) > 2.1")
    
    hTheta1.SetMarkerColor(2)
    hTheta2.SetMarkerColor(4)
    hTheta1.SetMarkerStyle(24)
    hTheta2.SetMarkerStyle(24)
    
    hTheta.Draw();
    hTheta1.Draw('same');
    hTheta2.Draw('same');
    
    canvas2.SaveAs('histSimTheta_' + str(nhisto) + '.pdf')
    
    
    nhisto += 1
    '''