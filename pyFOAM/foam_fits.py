from array import array
import numpy as np
import matplotlib.pyplot as plt
from ROOT import *

# ## AREA SCALING
# iteration           = np.linspace(0,1699,1700)
# # mean_area_local     = np.load('saved/experiments/local_1699_average_area.npy')[0:200]
# mean_area_quasi     = np.load('saved/experiments/quasi_global_1699_average_area_large.npy')
# # #mean_area_quasi     = np.load('saved/experiments/quasi_global_999_average_area_large.npy')
# # # mean_area_local_err = np.load('saved/experiments/local_1699_average_area_error.npy')[0:200]
# # mean_area_quasi_err = np.load('saved/experiments/quasi_global_1699_average_area_equal_boundary_error.npy')
# mean_area_quasi_err = np.load('saved/experiments/quasi_global_999_average_area_large_error.npy')
# #
# #
# # # Convert to arrays for Root to handle
# iteration_array            = array('d', iteration)
# # mean_area_local_array      = array('d', mean_area_local)
# mean_area_quasi_array      = array('d', mean_area_quasi)
# # mean_area_local_err_array  = array('d', mean_area_local_err)
# mean_area_quasi_err_array  = array('d', mean_area_quasi_err)
#
# # Setting of general plotting style:
# gStyle.SetCanvasColor(0)
# gStyle.SetFillColor(1)
# # Setting what to be shown in statistics box:
# gStyle.SetOptStat("emr")
# gStyle.SetOptFit(1111)
# gStyle.SetStatX(0.5)
# gStyle.SetStatY(1.0)

# #
# # Local graph
# c1 = TCanvas("c1", "", 100, 320, 600, 450)
# c1.SetFillColor(0)
# graph = TGraphErrors(len(iteration),iteration_array,mean_area_local_array,mean_area_local_err_array)
# graph.SetLineWidth(1)
# graph.SetMarkerStyle(20)
# graph.SetMarkerColor(4)
# graph.SetMarkerSize(0.4)
# graph.Draw("AP")
# graph.GetXaxis().SetTitle("Iteration")
# graph.GetYaxis().SetTitle("Mean cell area")
# graph.SetTitle("")
# ## Fit
# fit = TF1("fit", "[0]+[1]*x", 0, 1699)
# fit.SetParameters(0.1,1.0)
# fit.SetLineWidth(2)
# fit.SetLineColor(2)
# graph.Fit("fit", "R")
# c1.SaveAs("local_area_fit_200.pdf")
# c1.Update()
#
#
# # Quasi global graph
# c2 = TCanvas("c2", "", 0, 320, 600, 450)
# c2.SetFillColor(0)
# graph = TGraphErrors(len(iteration),iteration_array,mean_area_quasi_array,mean_area_quasi_err_array)
# graph.SetLineWidth(1)
# graph.SetMarkerStyle(20)
# graph.SetMarkerColor(4)
# graph.SetMarkerSize(0.4)
# graph.Draw("AP")
# graph.GetXaxis().SetTitle("Iteration")
# graph.GetYaxis().SetTitle("Mean cell area")
# graph.SetTitle("")
# ## Fit
# fit = TF1("fit", "[0]+[1]*x", 0, 1699)
# fit.SetParameters(0.1,1.0)
# fit.SetLineWidth(2)
# fit.SetLineColor(2)
# graph.Fit("fit", "R")
# c2.SaveAs("quasi_area_fit_1699_equal_boundary.pdf")
# c2.Update()




# iteration         = np.linspace(0,1699,1700)[20:100]
# merit_local       = np.load('saved/experiments/merit_local_1699_energy.npy')[20:100]
# merit_quasi       = np.load('saved/experiments/merit_quasi_global_1699_energy.npy')[20:100]
# merit_local_norm  = np.load('saved/experiments/merit_norm_local_1699_energy.npy')[20:100]
# merit_quasi_norm  = np.load('saved/experiments/merit_norm_quasi_global_1699_energy.npy')[20:100]
# merit_global      = np.load('saved/experiments/merit_global_99_energy.npy')[20:100]
# merit_global_norm = np.load('saved/experiments/merit_norm_global_99_energy.npy')[20:100]



# ## MERIT FUNCTIONS
#
# iteration         = np.linspace(0,1699,1700)
# merit_local       = np.load('saved/experiments/merit_local_1699_energy.npy')
# merit_quasi       = np.load('saved/experiments/merit_quasi_global_1699_energy.npy')
# merit_local_norm  = np.load('saved/experiments/merit_norm_local_1699_energy.npy')
# merit_quasi_norm  = np.load('saved/experiments/merit_norm_quasi_global_1699_energy.npy')
# merit_global      = np.load('saved/experiments/merit_global_99_energy.npy')
# merit_global_norm = np.load('saved/experiments/merit_norm_global_99_energy.npy')
#
#
# # Convert to arrays for Root to handle
# iteration_array         = array('d', iteration)
# merit_local_array       = array('d', merit_local)
# merit_quasi_array       = array('d', merit_quasi)
# merit_global_array      = array('d', merit_global)
# merit_local_norm_array  = array('d', merit_local_norm)
# merit_quasi_norm_array  = array('d', merit_quasi_norm)
# merit_global_norm_array = array('d', merit_global_norm)
#
#
# # Setting of general plotting style:
# gStyle.SetCanvasColor(0)
# gStyle.SetFillColor(1)
# # Setting what to be shown in statistics box:
# gStyle.SetOptStat("emr")
# gStyle.SetOptFit(0111)
# gStyle.SetStatX(0.9)
# gStyle.SetStatY(1.0)


# # LOCAL
# c1 = TCanvas("c1", "", 100, 400, 600, 450)
# c1.SetFillColor(0)
# graph = TGraphErrors(len(iteration),iteration_array,merit_local_norm_array)
# #graph = TGraphErrors(len(iteration),iteration_array,merit_local_array)
# graph.SetLineWidth(1)
# graph.SetMarkerStyle(20)
# graph.SetMarkerColor(4)
# graph.SetMarkerSize(0.4)
# graph.Draw("AP")
# graph.GetXaxis().SetTitle("Iteration")
# graph.GetYaxis().SetTitle("1/(2*num_{cells}) * F^{T}F")
# #graph.GetYaxis().SetTitle("1/2 * F^{T}F")
# graph.SetTitle("")
# ## Fit
# fit = TF1("fit", "[0]+[1]*x", 0, 1699)
# fit.SetParameters(0,0.1)
# fit.SetLineWidth(2)
# fit.SetLineColor(2)
# graph.Fit("fit", "R")
# c1.SaveAs("merit_norm_local_fit.pdf")
# #c1.SaveAs("merit_local_fit.pdf")
# c1.Update()
#
#
# # QUASI-GLOBAL
# c2 = TCanvas("c2", "", 100, 400, 600, 450)
# c2.SetFillColor(0)
# graph = TGraphErrors(len(iteration),iteration_array,merit_quasi_norm_array)
# #graph = TGraphErrors(len(iteration),iteration_array,merit_quasi_array)
# graph.SetLineWidth(1)
# graph.SetMarkerStyle(20)
# graph.SetMarkerColor(4)
# graph.SetMarkerSize(0.4)
# graph.Draw("AP")
# graph.GetXaxis().SetTitle("Iteration")
# graph.GetYaxis().SetTitle("1/(2*num_{cells}) * F^{T}F")
# #graph.GetYaxis().SetTitle("1/2 * F^{T}F")
# graph.SetTitle("")
# ## Fit
# fit = TH1("fit", "[0]+[1]*x", 0, 1699)
# fit.SetParameters(0,0.1)
# fit.SetLineWidth(2)
# fit.SetLineColor(2)
# graph.Fit("fit", "R")
# c2.SaveAs("merit_norm_quasi_fit.pdf")
# #c2.SaveAs("merit_quasi_fit.pdf")
# c2.Update()
#
#
#
# # GLOBAL
# c3 = TCanvas("c3", "", 100, 400, 600, 450)
# c3.SetFillColor(0)
# graph = TGraphErrors(100,iteration_array[0:100],merit_global_norm_array)
# #graph = TGraphErrors(100,iteration_array[0:100],merit_global_array)
# graph.SetLineWidth(1)
# graph.SetMarkerStyle(20)
# graph.SetMarkerColor(4)
# graph.SetMarkerSize(0.4)
# graph.Draw("AP")
# graph.GetXaxis().SetTitle("Iteration")
# graph.GetYaxis().SetTitle("1/(2*num_{cells}) * F^{T}F")
# #graph.GetYaxis().SetTitle("1/2 * F^{T}F")
# graph.SetTitle("")
# ## Fit
# fit = TF1("fit", "[0]+[1]*x", 0, 99)
# fit.SetParameters(0,0.1)
# fit.SetLineWidth(2)
# fit.SetLineColor(2)
# graph.Fit("fit", "R")
# #c3.SaveAs("merit_norm_global_fit.pdf")
# #c3.SaveAs("merit_global_fit.pdf")
# c3.Update()
#
#
#
#
# mean_val = np.mean(merit_global_norm)
#
# SSres = 0
# SStot = 0
#
# a = 0.06173
# b = -0.0003165
#
# for i in range(len(merit_global_norm)):
#     SSres += (merit_global_norm[i] - (a+b*i))**2
#     SStot += (merit_global_norm[i] - mean_val)**2
#
# print 1.0 - SSres/SStot


#
# ## PRESSIRE DISTRIBUTIONS
# iteration         = np.linspace(0,1699,1700)
# pressure_local    = np.load('saved/experiments/pressure_list_local.npy')
# pressure_quasi    = np.load('saved/experiments/pressure_list_large_quasi.npy')
# pressure_global   = np.load('saved/experiments/pressure_list_global.npy')
#
#
# # Convert to arrays for Root to handle
# iteration_array         = array('d', iteration)
# pressure_local_array    = array('d', pressure_local)
# pressure_quasi_array    = array('d', pressure_quasi)
# pressure_global_array   = array('d', pressure_global)
#
#
#
# Hist_Trans   = TH1F("Hist_Trans", "", 100, 0.0, 2.0)
#
# for x in pressure_quasi: Hist_Trans.Fill(x)
#
# c0 = TCanvas("c0", "", 20, 20, 1000, 600)
#
#
# gStyle.SetOptStat("emr")
# gStyle.SetOptFit(1111)
# gStyle.SetStatX(1.0)
# gStyle.SetStatY(1.0)
#
#
# # Setting line width and color and axis titles of histograms:
# Hist_Trans.GetXaxis().SetTitle("Pressure")
# Hist_Trans.GetYaxis().SetTitle("Counts")
# Hist_Trans.SetLineWidth(2)
# Hist_Trans.SetLineColor(kBlue)
#
# # Fitting histogram with Transformation numbers:
# fit_Trans = TF1("fit_Trans","[0]*exp(-(x-[1])**2/(2*[2]**2))", 0.0, 2.0)             # Draw function over histogram
# #fit_Trans = TF1("fit_Trans","[0]*x^[1]", 0.01, 1.0)                 # Fit function
# fit_Trans.SetParameters(1.0, 1.0, 1.0)
# fit_Trans.SetLineColor(kRed)
# Hist_Trans.Fit("fit_Trans", "RL")                                       # Fit function
# Hist_Trans.Draw()
# fit_Trans.Draw("same")
# gPad.Update()
# st1 = Hist_Trans.GetListOfFunctions().FindObject("stats")
# st1.SetX1NDC(0.12)
# st1.SetX2NDC(0.40)
# st1.SetY1NDC(0.89)
# st1.SetY2NDC(0.70)
# st1.SetLineColor(4)
#
# c0.Update()
#
# c0.SaveAs("pressure_large_quasi.pdf")



# ## MEAN PRESSURE DEVELOPMENT
# iteration           = np.linspace(0,1699,1700)
# mean_pressure_local = np.load('saved/experiments/average_pressure_list_local.npy')
# mean_pressure_quasi = np.load('saved/experiments/average_pressure_list_quasi.npy')
#
#
# # Convert to arrays for Root to handle
# iteration_array            = array('d', iteration)
# mean_pressure_local_array      = array('d', mean_pressure_local)
# mean_pressure_quasi_array      = array('d', mean_pressure_quasi)
#
#
# # Setting of general plotting style:
# gStyle.SetCanvasColor(0)
# gStyle.SetFillColor(1)
# # Setting what to be shown in statistics box:
# gStyle.SetOptStat("emr")
# gStyle.SetOptFit(1111)
# gStyle.SetStatX(0.5)
# gStyle.SetStatY(1.0)
#
# #
# # Local graph
# c1 = TCanvas("c1", "", 100, 320, 600, 450)
# c1.SetFillColor(0)
# graph = TGraphErrors(len(iteration),iteration_array,mean_pressure_quasi_array)
# graph.SetLineWidth(1)
# graph.SetMarkerStyle(20)
# graph.SetMarkerColor(4)
# graph.SetMarkerSize(0.4)
# graph.Draw("AP")
# graph.GetXaxis().SetTitle("Iteration")
# graph.GetYaxis().SetTitle("Mean pressure")
# graph.SetTitle("")
#
# c1.SaveAs("quasi_pressure.pdf")
# c1.Update()


# ## ENERGY DECAY
iteration        = np.linspace(0,1699,1700)
energy_quasi     = np.load('saved/experiments/quasi_global_1700_energy.npy')

# # Convert to arrays for Root to handle
iteration_array     = array('d', iteration)
energy_quasi_array  = array('d', energy_quasi)

# Setting of general plotting style:
gStyle.SetCanvasColor(0)
gStyle.SetFillColor(1)
# Setting what to be shown in statistics box:
gStyle.SetOptStat("emr")
gStyle.SetOptFit(1111)
gStyle.SetStatX(1.0)
gStyle.SetStatY(1.0)

#
# Local graph
c1 = TCanvas("c1", "", 100, 320, 600, 450)
c1.SetFillColor(0)
graph = TGraphErrors(len(iteration),iteration_array,energy_quasi_array)
graph.SetLineWidth(1)
graph.SetMarkerStyle(20)
graph.SetMarkerColor(4)
graph.SetMarkerSize(0.4)
graph.Draw("AP")
graph.GetXaxis().SetTitle("Iteration")
graph.GetYaxis().SetTitle("Total surface energy")
graph.SetTitle("")
## Fit
# fit = TF1("fit", "[0]*exp([1]*x)+[2]", 0, 1700)
# fit.SetParameters(1900.0,-1.0,0.1)
fit = TF1("fit", "[0]*x**(-[1])+[2]", 1, 1700)
fit.SetParameters(1900.0,1.0,500.0)
fit.SetLineWidth(2)
fit.SetLineColor(2)
graph.Fit("fit", "R")
#c1.SaveAs("quasi_energy_fit.pdf")
c1.Update()
#
#
# mean_val = np.mean(energy_quasi)
#
# SSres = 0
# SStot = 0
#
# a = 1047
# b = -0.001654
# c = 743.9
#
# for i in range(len(energy_quasi)):
#     SSres += (energy_quasi[i] - (a*exp(b*i)+c))**2
#     SStot += (energy_quasi[i] - mean_val)**2
#
# print 1.0 - SSres/SStot




# ## Cell Count
# iteration        = np.linspace(0,1700,1699)
# cell_count_quasi = np.load('saved/experiments/cell_counter_list_quasi_large.npy')
#
# # # Convert to arrays for Root to handle
# iteration_array        = array('d', iteration)
# cell_count_quasi_array = array('d', cell_count_quasi)
#
# # Setting of general plotting style:
# gStyle.SetCanvasColor(0)
# gStyle.SetFillColor(1)
# # Setting what to be shown in statistics box:
# gStyle.SetOptStat("emr")
# gStyle.SetOptFit(1111)
# gStyle.SetStatX(1.0)
# gStyle.SetStatY(1.0)
#
# #
# # Local graph
# c1 = TCanvas("c1", "", 100, 320, 600, 450)
# c1.SetFillColor(0)
# graph = TGraphErrors(len(iteration),iteration_array,cell_count_quasi_array)
# graph.SetLineWidth(1)
# graph.SetMarkerStyle(20)
# graph.SetMarkerColor(4)
# graph.SetMarkerSize(0.4)
# graph.Draw("AP")
# graph.GetXaxis().SetTitle("Iteration")
# graph.GetYaxis().SetTitle("Number of cells")
# graph.SetTitle("")
# ## Fit
# #fit = TF1("fit", "[0]*exp([1]*x)+[2]", 0, 1700)
# fit = TF1("fit", "[0]*x**(-[1])+[2]", 1, 1700)
# fit.SetParameters(1900.0,1.0,500.0)
# fit.SetLineWidth(2)
# fit.SetLineColor(2)
# graph.Fit("fit", "R")
# #c1.SaveAs("quasi_cell_count_fit.pdf")
# c1.Update()












# # ## ABOAV-WEAIRE
# n            = np.linspace(3,10,8)
# n_err        = np.zeros((8,1))
# aw_quasi     = np.load('saved/experiments/quasi_global_1000_aboav_weaire.npy')
# aw_quasi_err = np.load('saved/experiments/quasi_global_1000_aboav_weaire_error.npy')
#
#
# # # Convert to arrays for Root to handle
# n_array            = array('d', n)
# n_err_array        = array('d', n_err)
# aw_quasi_array     = array('d', aw_quasi)
# aw_quasi_err_array = array('d', aw_quasi_err)
#
#
# # Setting of general plotting style:
# gStyle.SetCanvasColor(0)
# gStyle.SetFillColor(1)
# # Setting what to be shown in statistics box:
# gStyle.SetOptStat("emr")
# gStyle.SetOptFit(1111)
# gStyle.SetStatX(0.5)
# gStyle.SetStatY(1.0)
#
# #
# # Local graph
# c1 = TCanvas("c1", "", 100, 320, 600, 450)
# c1.SetFillColor(0)
# graph = TGraphErrors(len(n_array),n_array,aw_quasi_array)#,n_err_array,aw_quasi_err_array)
# graph.SetLineWidth(1)
# graph.SetMarkerStyle(20)
# graph.SetMarkerColor(4)
# graph.SetMarkerSize(0.4)
# graph.Draw("AP")
# graph.GetXaxis().SetTitle("n")
# graph.GetYaxis().SetTitle("n*m(n)")
# graph.SetTitle("")
# ## Fit
# fit = TF1("fit", "(6-[0])*x + 6*[0] + 1.57623488629", 2, 11)
# #fit = TF1("fit", "6-[0] + (6*[0] + [1])/x", 2, 11)
# fit.SetParameters(1.0,0.0)
# fit.SetLineWidth(2)
# fit.SetLineColor(2)
# graph.Fit("fit", "R")
# c1.SaveAs("aw_1000.pdf")
# c1.Update()










raw_input('Press Enter to exit')
