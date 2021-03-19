# algorithm:
# 0 = Std
# 1 = Bayes

import ROOT
from ROOT import TH1F, TCanvas, gStyle, TLegend
from ROOT import kRed, kGreen

splitting_paths_bayes = False

input_file = "../../../../dump.txt"

with open(input_file) as f:
    lines = f.readlines()

h_filtering_std = TH1F("h_filtering_std", "h_filtering_std", 50, 0, 500)
h_filtering_bay = TH1F("h_filtering_bay", "h_filtering_bay", 50, 0, 500)

h_fitting_std = TH1F("h_fitting_std", "h_fitting_std", 50, 0, 500)
h_fitting_bay = TH1F("h_fitting_bay", "h_fitting_bay", 50, 0, 500)

h_correlation_std = TH1F("h_correlation_std", "h_correlation_std", 50, 0, 500)
h_correlation_bay = TH1F("h_correlation_bay", "h_correlation_bay", 50, 0, 500)

for line in lines:
    
    # Muon paths after filtering
    if "MUON PATH" and "after filtering" in line:
        numbers = [float(s) for s in line.split() if s.isdigit()]
        if numbers[0] == 0:
            h_filtering_std.Fill(numbers[1])
        elif numbers[0] == 1:
            if splitting_paths_bayes == True:
                h_filtering_bay.Fill(3*numbers[1])
            else:
                h_filtering_bay.Fill(numbers[1])

    # Muon paths after fitting
    elif "after fitting" in line:
        numbers = [float(s) for s in line.split() if s.isdigit()]
        if numbers[0] == 0:
            h_fitting_std.Fill(numbers[1])
        elif numbers[0] == 1:
            h_fitting_bay.Fill(numbers[1])

    # Muon paths after correlation
    elif "Correlated metaPrimitive" in line:
        numbers = [float(s) for s in line.split() if s.isdigit()]
        if numbers[0] == 0:
            h_correlation_std.Fill(numbers[1])
        elif numbers[0] == 1:
            h_correlation_bay.Fill(numbers[1])


gStyle.SetOptStat(0)
gStyle.SetLegendTextSize(0.035)
ROOT.gROOT.SetBatch(True)

# Filtering
c_filtering = TCanvas("c_filtering","c_filtering", 600, 600)
c_filtering.cd()

h_filtering_std.SetTitle("Muon paths before fitting")
h_filtering_std.GetXaxis().SetTitle("Number of Muon Paths per Event")
h_filtering_std.GetYaxis().SetRangeUser(0., 1.4*max(h_filtering_std.GetMaximum(), h_filtering_bay.GetMaximum()))

h_filtering_std.SetLineColor(kRed+1)
h_filtering_bay.SetLineColor(kGreen+1)
h_filtering_std.SetLineWidth(2)
h_filtering_bay.SetLineWidth(2)

h_filtering_std.Draw()
h_filtering_bay.Draw("same")

mean_std = h_filtering_std.GetMean()
mean_bay = h_filtering_bay.GetMean()

leg_filtering = TLegend(0.20,0.73,0.50,0.88)
leg_filtering.SetLineColor(0)
leg_filtering.AddEntry(h_filtering_std, "Standard: mean = {:.0f}".format(mean_std),"l")
leg_filtering.AddEntry(h_filtering_bay, "Bayes: mean = {:.0f}".format(mean_bay),   "l")
leg_filtering.Draw()

c_filtering.Print("../../../../Groupings/StdToBayes/h_1_filtering.png")

# Fitting
c_fitting = TCanvas("c_fitting","c_fitting", 600, 600)
c_fitting.cd()

h_fitting_std.SetTitle("Muon paths after fitting")
h_fitting_std.GetXaxis().SetTitle("Number of Muon Paths per Event")
h_fitting_std.GetYaxis().SetRangeUser(0., 1.4*max(h_fitting_std.GetMaximum(), h_fitting_bay.GetMaximum()))

h_fitting_std.SetLineColor(kRed+1)
h_fitting_bay.SetLineColor(kGreen+1)
h_fitting_std.SetLineWidth(2)
h_fitting_bay.SetLineWidth(2)

h_fitting_std.Draw()
h_fitting_bay.Draw("same")

mean_std = h_fitting_std.GetMean()
mean_bay = h_fitting_bay.GetMean()

leg_fitting = TLegend(0.20,0.73,0.50,0.88)
leg_fitting.SetLineColor(0)
leg_fitting.AddEntry(h_fitting_std, "Standard: mean = {:.0f}".format(mean_std),"l")
leg_fitting.AddEntry(h_fitting_bay, "Bayes: mean = {:.0f}".format(mean_bay),   "l")
leg_fitting.Draw()

c_fitting.Print("../../../../Groupings/StdToBayes/h_2_fitting.png")

# Correlation
c_correlation = TCanvas("c_correlation","c_correlation", 600, 600)
c_correlation.cd()

h_correlation_std.SetTitle("Muon pathsa after correlation and cleaning")
h_correlation_std.GetXaxis().SetTitle("Number of Muon Paths per Event")
h_correlation_std.GetYaxis().SetRangeUser(0., 1.4*max(h_correlation_std.GetMaximum(), h_correlation_bay.GetMaximum()))

h_correlation_std.SetLineColor(kRed+1)
h_correlation_bay.SetLineColor(kGreen+1)
h_correlation_std.SetLineWidth(2)
h_correlation_bay.SetLineWidth(2)

h_correlation_std.Draw()
h_correlation_bay.Draw("same")

mean_std = h_correlation_std.GetMean()
mean_bay = h_correlation_bay.GetMean()

leg_correlation = TLegend(0.20,0.73,0.50,0.88)
leg_correlation.SetLineColor(0)
leg_correlation.AddEntry(h_correlation_std, "Standard: mean = {:.0f}".format(mean_std),"l")
leg_correlation.AddEntry(h_correlation_bay, "Bayes: mean = {:.0f}".format(mean_bay),   "l")
leg_correlation.Draw()

c_correlation.Print("../../../../Groupings/StdToBayes/h_3_correlation.png")


f.close()    
