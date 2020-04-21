import pickle, os 
import math 
import ROOT 

ss = ["shared25","shared50","shared75","shared100"]
samples = {"shared25":ROOT.kBlack, "shared50":ROOT.kRed, "shared75":ROOT.kBlue, "shared100":ROOT.kGreen+1}
tags = {"shared25": "shared 25% hits", "shared50": "shared 50% hits", "shared75":"shared 75% hits" , "shared100":"shared 100% hits"}

plots=[]
titles={}  
for var in ["Phi","PhiB","Chi2","Bx","Time"]:
    for st in ["MB1","MB2","MB3","MB4"]:
        for q in ["q1","q3","q5","q8"]:
            plots.append("h%sRes_%s_%s" %(var,st,q))
            titles["h%sRes_%s_%s" %(var,st,q)] = '; %s_{bay}-%s_{std};' %(var,var)


outpath = "/afs/cern.ch/user/f/folguera/www/private/L1TPhase2/DTTP/200421_Groupings/"
if not os.path.exists(outpath):
    os.mkdir(outpath)
    print "cp ~folguera/public/utils/index.php %s/" %outpath
    os.system("cp ~folguera/public/utils/index.php %s/" %outpath)
os.system("cp EventDumpList.log %s/" %outpath)

outpath = outpath + "BayesToStd/"
if not os.path.exists(outpath):
    os.mkdir(outpath)
    print "cp ~folguera/public/utils/index.php %s/" %outpath
    os.system("cp ~folguera/public/utils/index.php %s/" %outpath)

outFile = ROOT.TFile("GroupingComparison_BayesToStd_Apr21.root","RECREATE")
outFile.cd()

ROOT.gROOT.ProcessLine('.L PlotTemplate.C+')
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

with open('GroupingComparison_BayesToStd_Apr21.pickle', 'rb') as handle:
    b = pickle.load(handle)

leg = ROOT.TLegend(0.6,0.6,0.88,0.4);
leg.SetTextSize(0.03);

for s in ss:        
    for plot in plots:
        b[s][plot].SetTitle(titles[plot])
        b[s][plot].Write()
        b[s][plot].SetLineColor(samples[s])
        b[s][plot].SetMarkerStyle(8)
        b[s][plot].SetMarkerColor(samples[s])
        b[s][plot].SetMarkerSize(0.8)
        b[s][plot].SetLineWidth(2)

    leg.AddEntry(b[s][plots[0]], tags[s],'l');

canvas = ROOT.CreateCanvas('name',False,True)


for plot in plots:
    drawn = False
    counts=0
    disp = ROOT.TLatex()
    disp.SetTextSize(0.025)
    for s in ss: 
        if not (drawn): 
            b[s][plot].Draw()
            drawn=True
        else:
            b[s][plot].Draw('same')            
        
        disp.DrawLatexNDC(0.13,0.87-counts*0.03,"#color[%d]{%s: %2.2f (%2.2f)}" %(samples[s],tags[s],b[s][plot].GetMean(),b[s][plot].GetRMS()))
#        disp.DrawLatexNDC(0.13,0.87-counts*0.03,"%s: %2.2f (%2.2f)" %(tags[s],b[s][plot].GetMean(),b[s][plot].GetRMS()))
        counts = counts+1

    disp.Draw("same")    
    ROOT.DrawPrelimLabel(canvas)
    ROOT.DrawLumiLabel(canvas,'200 PU')

    leg.Draw("same")
    
    ROOT.SaveCanvas(canvas, outpath + plot)



outFile.Close()



