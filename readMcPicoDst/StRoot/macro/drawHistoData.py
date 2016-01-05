import ROOT
ROOT.gROOT.SetBatch(True)

import rootpy.ROOT as ROOT

from rootpy.tree import Cut
from rootpy.tree.categories import Categories
from rootpy.io import root_open
from rootpy.plotting import Hist, Hist2D, Hist3D, histogram, Canvas, Legend

# Data
infile = root_open('../out_gamma9.root')
ntuple = infile.nt
# Simulation
infile2 = root_open('out.root')
ntuple2 = infile2.nt2

canvas = Canvas()
canvas.SetLogz()

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
           cut='mass < 0.02 && pairDca < 0.01 && abs(pt1) > 0.6 && abs(pt2) > 0.6 && abs(pt1) < 0.8 && abs(pt2) < 0.8 && nsige1>-1 && nsige2>-1',
           histoname='histo1DData',
           cutSim='abs(pt1) > 0.6 && abs(pt2) > 0.6 && abs(pt1) < 0.8 && abs(pt2) < 0.8 ',
           xvalSim='',xvalSim2=''
           ):
    histo1 = ntuple.Draw(xval, hist=Hist(xbin,xmin,xmax), selection=cut+' && pt1*pt2<0')
    histo1.color = 'red'
    histo1.SetMarkerStyle(24)
    histo2 = ntuple.Draw(xval, hist=Hist(xbin,xmin,xmax), selection=cut+' && pt1*pt2>0')
    #histo2.scale(histo1.integral(90,99)/histo2.integral(90,99))
    histo2.color = 'blue'
    histo2.SetMarkerStyle(24)
    histo3 = histo1-histo2
    histo3.scale(1/histo3.GetMaximum())
    histo3.color = 'black'
    histo3.SetMarkerStyle(20)
    if xvalSim == '':
        xvalSim = xval
    histo4 = ntuple2.Draw(xvalSim, hist=Hist(xbin,xmin,xmax), selection=cutSim)
    histo4.scale(1/histo4.GetMaximum())
    histo4.color = 'red'
    histo4.scale(histo3.GetMaximum()/histo4.GetMaximum())
    if xvalSim2 == '':
        xvalSim2 = xval
    histo5 = ntuple2.Draw(xvalSim2, hist=Hist(xbin,xmin,xmax), selection=cutSim)
    histo5.scale(1/histo5.GetMaximum())
    histo5.color = 'blue'
    histo5.scale(histo3.GetMaximum()/histo4.GetMaximum())

    histo3.SetMinimum(0)
    #histo1.Draw('E0')
    #histo2.Draw('sameE0')
    histo3.Draw('E0')
    histo4.Draw('sameE0')
    histo5.Draw('sameE0')
    histo1.xaxis.SetTitle(xvalname)
    canvas.SaveAs(histoname+'_'+xvalname+'.pdf')

draw1D('sqrt((v0x+0.2383)**2+(v0y+0.1734)**2)',0, 10,'ConvR',100,xvalSim='sqrt((v0x)**2+(v0y)**2)',xvalSim2='sqrt((mcv0x)**2+(mcv0y)**2)')
draw1D('mass',0, 0.1,'Mass',100)

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

