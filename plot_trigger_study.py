#!/usr/bin/env python3
"""
Plot all histograms from trigger_study.C output ROOT files.

Reads ROOT files produced by trigger_study.C and creates plots for:
- Cutflow histograms
- Di-muon invariant mass distributions
- Individual L1 trigger counts
- Trigger efficiency vs kinematic variables (leading pT, subleading pT, di-pT)

Usage:
    python plot_trigger_study.py
    (configure base_dirs and regions in __main__ section)
"""

import ROOT
import os
import glob
import re
from tqdm import tqdm

# CMS Style settings
def setCMSStyle():
    """Set CMS style for plots"""
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gErrorIgnoreLevel = ROOT.kWarning

    # Canvas settings
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
    ROOT.gStyle.SetCanvasDefH(600)
    ROOT.gStyle.SetCanvasDefW(600)
    ROOT.gStyle.SetCanvasDefX(0)
    ROOT.gStyle.SetCanvasDefY(0)

    # Pad settings
    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetPadColor(ROOT.kWhite)
    ROOT.gStyle.SetPadGridX(False)
    ROOT.gStyle.SetPadGridY(False)
    ROOT.gStyle.SetGridColor(0)
    ROOT.gStyle.SetGridStyle(3)
    ROOT.gStyle.SetGridWidth(1)

    # Frame settings
    ROOT.gStyle.SetFrameBorderMode(0)
    ROOT.gStyle.SetFrameBorderSize(1)
    ROOT.gStyle.SetFrameFillColor(0)
    ROOT.gStyle.SetFrameFillStyle(0)
    ROOT.gStyle.SetFrameLineColor(1)
    ROOT.gStyle.SetFrameLineStyle(1)
    ROOT.gStyle.SetFrameLineWidth(1)

    # Margins
    ROOT.gStyle.SetPadTopMargin(0.08)
    ROOT.gStyle.SetPadBottomMargin(0.13)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadRightMargin(0.05)

    # Axis settings
    ROOT.gStyle.SetTitleSize(0.05, "XYZ")
    ROOT.gStyle.SetTitleXOffset(0.9)
    ROOT.gStyle.SetTitleYOffset(1.4)
    ROOT.gStyle.SetLabelSize(0.04, "XYZ")
    ROOT.gStyle.SetLabelOffset(0.007, "XYZ")

    # Legend settings
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetLegendFillColor(0)
    ROOT.gStyle.SetLegendFont(42)
    ROOT.gStyle.SetLegendTextSize(0.03)

    ROOT.gStyle.SetPalette(ROOT.kBird)
    ROOT.gStyle.SetNumberContours(99)

    ROOT.gROOT.ForceStyle()


def addCMSText(canvas, lumi_text="13.6 TeV", extra_text="Work in Progress"):
    """Add CMS label and luminosity text"""
    canvas.cd()

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)

    # CMS text
    latex.SetTextSize(0.05)
    latex.SetTextFont(61)
    latex.DrawLatex(0.15, 0.93, "CMS")

    # Extra text (e.g., Preliminary)
    if extra_text:
        latex.SetTextSize(0.04)
        latex.SetTextFont(52)
        latex.DrawLatex(0.25, 0.93, extra_text)

    # Luminosity/energy text
    if lumi_text:
        latex.SetTextSize(0.04)
        latex.SetTextFont(42)
        latex.DrawLatex(0.70, 0.93, lumi_text)

    return latex


def plot_cutflow_overlay(root_files, output_dir, sample_base, normalize=False):
    """Plot overlaid cutflow histograms for samples with same origin.

    Args:
        normalize: If True, normalize each histogram to its first bin
    """
    if not root_files:
        return

    suffix = "_normalized" if normalize else ""
    canvas = ROOT.TCanvas(f"c_cutflow_overlay_{sample_base}{suffix}", f"c_cutflow_overlay_{sample_base}{suffix}", 800, 800)
    canvas.cd()
    canvas.SetLogy()

    colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2,
              ROOT.kMagenta, ROOT.kOrange+7, ROOT.kCyan+1, ROOT.kViolet]

    legend = ROOT.TLegend(0.18, 0.12, 0.65, 0.60)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.022)
    legend.SetHeader(sample_base, "L")

    histograms = []
    max_val = 0

    for i, (filepath, tag) in enumerate(root_files):
        tfile = ROOT.TFile.Open(filepath)
        if not tfile or tfile.IsZombie():
            continue

        hist = tfile.Get("h_cutflow")
        if not hist:
            tfile.Close()
            continue

        hist_clone = hist.Clone(f"h_cutflow_{tag}")
        hist_clone.SetDirectory(0)
        hist_clone.SetLineColor(colors[i % len(colors)])
        hist_clone.SetLineWidth(2)
        hist_clone.GetYaxis().SetTitle("Events")

        max_val = max(max_val, hist_clone.GetMaximum())
        histograms.append(hist_clone)

        # Extract label from tag
        # For signal samples: extract mass parameters (e.g., Par-MTp1000-MS100)
        # For BkgMC samples: extract process and bin info (e.g., DYto2Mu_Bin-MLL-10to50)
        # For data samples: extract run period (e.g., Run2024F)
        label_parts = tag.split('_')
        # Remove region prefix if present
        if label_parts[0] in ['sr', 'vr']:
            label_parts = label_parts[1:]

        if any('Par-' in part for part in label_parts):
            # Signal sample: extract the parameter part
            label = [part for part in label_parts if 'Par-' in part][0]
        elif any('Bin-' in part for part in label_parts):
            # BkgMC sample: extract process name and bin info
            process_name = label_parts[0] if label_parts else ""
            bin_info = [part for part in label_parts if 'Bin-' in part]
            label = f"{process_name}_{bin_info[0]}" if bin_info else process_name
        else:
            # Data sample: extract last part (e.g., Run2024F)
            label = label_parts[-1] if len(label_parts) > 1 else tag
        legend.AddEntry(hist_clone, label, "l")

        tfile.Close()

    if not histograms:
        canvas.Close()
        return

    # Normalize to first bin if requested
    if normalize:
        for hist in histograms:
            first_bin_val = hist.GetBinContent(1)
            if first_bin_val > 0:
                hist.Scale(1.0 / first_bin_val)
        # Recalculate max_val after normalization
        max_val = max(hist.GetMaximum() for hist in histograms)

    for i, hist in enumerate(histograms):
        if normalize:
            hist.GetYaxis().SetTitle("Fraction of Events")
            hist.SetMinimum(0.0001)
            hist.SetMaximum(max_val * 10)
        else:
            hist.SetMinimum(0.1)
            hist.SetMaximum(max_val * 1000)
        hist.Draw("HIST" if i == 0 else "HIST SAME")

    legend.Draw()
    addCMSText(canvas, lumi_text="13.6 TeV", extra_text="Work in Progress")

    suffix = "_normalized" if normalize else ""
    canvas.SaveAs(os.path.join(output_dir, f"cutflow_overlay_{sample_base}{suffix}.png"))
    canvas.SaveAs(os.path.join(output_dir, f"cutflow_overlay_{sample_base}{suffix}.pdf"))
    canvas.Update()
    canvas.Close()




def plot_mass_histograms(root_file, output_dir, tag):
    """Plot di-muon invariant mass histograms overlaid."""
    tfile = ROOT.TFile.Open(root_file)
    if not tfile or tfile.IsZombie():
        return

    canvas = ROOT.TCanvas(f"c_mass_{tag}", f"c_mass_{tag}", 800, 800)
    canvas.cd()
    canvas.SetLogy()

    colors = [ROOT.kBlue, ROOT.kRed, ROOT.kGreen+2]
    labels = ["2 quality #mu", "di-p_{T} > 20 GeV"]
    hist_names = ["h_mass_dimu", "h_mass_dipt"]

    legend = ROOT.TLegend(0.55, 0.70, 0.88, 0.88)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)

    histograms = []
    fit_results = []
    fit_functions = []
    for i, hname in enumerate(hist_names):
        h = tfile.Get(hname)
        if not h:
            continue
        h_clone = h.Clone(f"{hname}_clone")
        h_clone.SetDirectory(0)
        h_clone.SetLineColor(colors[i])
        h_clone.SetLineWidth(2)
        histograms.append(h_clone)
        legend.AddEntry(h_clone, labels[i], "l")

        # Fit Gaussian to Z peak
        fit_func = ROOT.TF1(f"gaus_mass_{i}", "gaus", 70, 110)
        fit_func.SetLineColor(colors[i])
        fit_func.SetLineStyle(2)
        fit_func.SetLineWidth(2)
        fit_func.SetParameters(h_clone.GetMaximum(), 91, 3)  # Initial guess
        fit_result = h_clone.Fit(fit_func, "RQNS")  # R=range, Q=quiet, N=no draw, S=save result
        if fit_result.IsValid():
            mean = fit_func.GetParameter(1)
            sigma = fit_func.GetParameter(2)
            fit_results.append((mean, sigma, colors[i]))
            fit_functions.append(fit_func)

    if not histograms:
        tfile.Close()
        canvas.Close()
        return

    max_val = max(h.GetMaximum() for h in histograms)

    for i, h in enumerate(histograms):
        h.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
        h.GetYaxis().SetTitle("Events / bin")
        h.SetMaximum(max_val * 10)
        h.SetMinimum(0.5)
        h.Draw("HIST" if i == 0 else "HIST SAME")

    # Draw fit functions (limited range to avoid overlap with text)
    for fit_func in fit_functions:
        # Draw only the central part of the fit (above y ~ 1000)
        fit_func.SetRange(82, 100)
        fit_func.Draw("SAME")

    legend.Draw()

    # Display fit results in bottom left
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.032)
    y_pos = 0.30
    for mean, sigma, color in fit_results:
        latex.SetTextColor(color)
        latex.DrawLatex(0.18, y_pos, f"#mu = {mean:.2f} GeV, #sigma = {sigma:.2f} GeV")
        y_pos -= 0.05

    addCMSText(canvas, lumi_text="13.6 TeV", extra_text="Work in Progress")

    canvas.SaveAs(os.path.join(output_dir, f"mass_{tag}.png"))
    canvas.SaveAs(os.path.join(output_dir, f"mass_{tag}.pdf"))
    canvas.Update()
    canvas.Close()

    tfile.Close()


def plot_trig_fired_overlay(root_files, output_dir, sample_base, normalize=False):
    """Plot overlaid trigger fired histograms for different eras.

    Args:
        normalize: If True, normalize each histogram to its first bin
    """
    if not root_files:
        return

    suffix = "_normalized" if normalize else ""
    canvas = ROOT.TCanvas(f"c_trig_fired_overlay_{sample_base}{suffix}", f"c_trig_fired_overlay_{sample_base}{suffix}", 800, 800)
    canvas.cd()
    canvas.SetLogy()

    colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2,
              ROOT.kMagenta, ROOT.kOrange+7, ROOT.kCyan+1, ROOT.kViolet]

    legend = ROOT.TLegend(0.18, 0.12, 0.65, 0.60)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.022)
    legend.SetHeader(sample_base, "L")

    histograms = []
    max_val = 0

    for i, (filepath, tag) in enumerate(root_files):
        tfile = ROOT.TFile.Open(filepath)
        if not tfile or tfile.IsZombie():
            continue

        hist = tfile.Get("h_trig_fired")
        if not hist:
            tfile.Close()
            continue

        hist_clone = hist.Clone(f"h_trig_fired_{tag}")
        hist_clone.SetDirectory(0)
        hist_clone.SetLineColor(colors[i % len(colors)])
        hist_clone.SetLineWidth(2)
        hist_clone.GetXaxis().SetTitle("")
        hist_clone.GetXaxis().SetLabelSize(0.025)
        hist_clone.GetYaxis().SetTitle("Events")
        hist_clone.LabelsOption("a", "X")  # Rotate labels

        max_val = max(max_val, hist_clone.GetMaximum())
        histograms.append(hist_clone)

        # Extract label from tag
        # For signal samples: extract mass parameters (e.g., Par-MTp1000-MS100)
        # For BkgMC samples: extract process and bin info (e.g., DYto2Mu_Bin-MLL-10to50)
        # For data samples: extract run period (e.g., Run2024F)
        label_parts = tag.split('_')
        # Remove region prefix if present
        if label_parts[0] in ['sr', 'vr']:
            label_parts = label_parts[1:]

        if any('Par-' in part for part in label_parts):
            # Signal sample: extract the parameter part
            label = [part for part in label_parts if 'Par-' in part][0]
        elif any('Bin-' in part for part in label_parts):
            # BkgMC sample: extract process name and bin info
            process_name = label_parts[0] if label_parts else ""
            bin_info = [part for part in label_parts if 'Bin-' in part]
            label = f"{process_name}_{bin_info[0]}" if bin_info else process_name
        else:
            # Data sample: extract last part (e.g., Run2024F)
            label = label_parts[-1] if len(label_parts) > 1 else tag
        legend.AddEntry(hist_clone, label, "l")

        tfile.Close()

    if not histograms:
        canvas.Close()
        return

    # Normalize to first bin if requested
    if normalize:
        for hist in histograms:
            first_bin_val = hist.GetBinContent(1)
            if first_bin_val > 0:
                hist.Scale(1.0 / first_bin_val)
        # Recalculate max_val after normalization
        max_val = max(hist.GetMaximum() for hist in histograms)

    for i, hist in enumerate(histograms):
        if normalize:
            hist.GetYaxis().SetTitle("Fraction of Events")
            hist.SetMinimum(0.0001)
            hist.SetMaximum(max_val * 10)
        else:
            hist.SetMinimum(0.1)
            hist.SetMaximum(max_val * 10)
        hist.Draw("HIST" if i == 0 else "HIST SAME")

    legend.Draw()
    addCMSText(canvas, lumi_text="13.6 TeV", extra_text="Work in Progress")

    suffix = "_normalized" if normalize else ""
    canvas.SaveAs(os.path.join(output_dir, f"trig_fired_overlay_{sample_base}{suffix}.png"))
    canvas.SaveAs(os.path.join(output_dir, f"trig_fired_overlay_{sample_base}{suffix}.pdf"))
    canvas.Update()
    canvas.Close()


def plot_trigger_efficiency(root_file, output_dir, tag):
    """
    Plot trigger efficiency vs kinematic variables with all triggers overlaid.

    Creates ProfileX from 2D histograms to show efficiency vs pT.
    Different triggers shown with different colors.
    """
    tfile = ROOT.TFile.Open(root_file)
    if not tfile or tfile.IsZombie():
        return

    # All triggers to overlay
    triggers = [
        ("DoubleMu Triggers OR", "DoubleMu Triggers OR"),
        ("DoubleMu_12_5", "DoubleMu_12_5"),
        ("DoubleMu_15_7", "DoubleMu_15_7"),
        ("DoubleMu4p5er2p0_SQ_OS_Mass_Min7", "DoubleMu4p5er2p0_SQ_OS_Mass_Min7"),
        ("DoubleMu4p5er2p0_SQ_OS_Mass_7to18", "DoubleMu4p5er2p0_SQ_OS_Mass_7to18"),
        ("DoubleMu4_SQ_OS_dR_Max1p2", "DoubleMu4_SQ_OS_dR_Max1p2"),
        ("DoubleMu4p5_SQ_OS_dR_Max1p2", "DoubleMu4p5_SQ_OS_dR_Max1p2"),
    ]

    colors = [
        ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2,
        ROOT.kMagenta, ROOT.kOrange+7, ROOT.kCyan+1
    ]

    markers = [20, 21, 22, 23, 33, 34, 29]

    # Kinematic variables to plot
    kin_vars = {
        "leading": ("p_{T}^{lead} [GeV]", 0, 100),
        "subleading": ("p_{T}^{sublead} [GeV]", 0, 100),
        "dipt": ("p_{T}^{#mu#mu} [GeV]", 0, 200),
    }

    # Plot for both inclusive and quality pair versions
    for suffix, suffix_label in [("", "inclusive"), ("_2mu", "quality pair")]:
        for var, (xlabel, xmin, xmax) in kin_vars.items():
            canvas = ROOT.TCanvas(f"c_trigeff_{var}{suffix}_{tag}",
                                f"c_trigeff_{var}{suffix}_{tag}", 800, 800)
            canvas.cd()
            canvas.SetRightMargin(0.05)
            canvas.SetTopMargin(0.10)

            # Create legend at the top
            legend = ROOT.TLegend(0.18, 0.65, 0.90, 0.88)
            legend.SetFillStyle(0)
            legend.SetBorderSize(0)
            legend.SetTextSize(0.018)
            # legend.SetNColumns(2)

            # Parse tag to extract sample and era
            tag_parts = tag.split('_')
            # Remove region prefix if present
            if tag_parts[0] in ['sr', 'vr']:
                tag_parts = tag_parts[1:]
            sample_name = tag_parts[0] if len(tag_parts) > 0 else ""
            era_name = tag_parts[1] if len(tag_parts) > 1 else ""

            # Format header text
            suffix_text = "requiring a quality pair" if suffix_label == "quality pair" else "inclusive"
            legend.SetHeader(f"{sample_name}  {era_name} ({suffix_text})", "L")

            profiles = []
            first_drawn = False

            for i, (trig_name, trig_label) in enumerate(triggers):
                hname = f"h2_trigeff_{var}{suffix}_{trig_name}"
                h2 = tfile.Get(hname)
                if not h2:
                    continue

                # Create profile without rebinning
                prof = h2.ProfileX(f"prof_{hname}")
                prof.SetDirectory(0)
                prof.SetLineColor(colors[i % len(colors)])
                prof.SetLineWidth(2)
                prof.SetMarkerColor(colors[i % len(colors)])
                prof.SetMarkerStyle(markers[i % len(markers)])
                prof.SetMarkerSize(0.9)

                if not first_drawn:
                    prof.GetXaxis().SetTitle(xlabel)
                    prof.GetYaxis().SetTitle("Trigger Efficiency wrt L1 H_{T} triggers")
                    prof.GetXaxis().SetRangeUser(xmin, xmax)
                    prof.SetMinimum(0.0)
                    prof.SetMaximum(2.0)
                    prof.SetTitle("")
                    prof.Draw("E")
                    first_drawn = True
                else:
                    prof.Draw("E SAME")

                legend.AddEntry(prof, trig_label, "lpe")
                profiles.append(prof)

            if not profiles:
                canvas.Close()
                continue

            # Add horizontal line at 1.0
            line = ROOT.TLine(xmin, 1.0, xmax, 1.0)
            line.SetLineStyle(2)
            line.SetLineColor(ROOT.kGray+1)
            line.SetLineWidth(1)
            line.Draw()

            legend.Draw()

            addCMSText(canvas, lumi_text="13.6 TeV", extra_text="Work in Progress")

            varname_safe = var + suffix.replace("_", "")
            canvas.SaveAs(os.path.join(output_dir, f"trigeff_overlay_{varname_safe}_{tag}.png"))
            canvas.SaveAs(os.path.join(output_dir, f"trigeff_overlay_{varname_safe}_{tag}.pdf"))
            canvas.Update()
            canvas.Close()

    tfile.Close()


def plot_trigger_efficiency_with_fit(root_file, output_dir, tag):
    """
    Plot DoubleMu Triggers OR efficiency vs kinematic variables with turn-on function fit.

    Creates ProfileX from 2D histograms and fits with plain sigmoid function (equation 19):
    f(x) = a / (1 + exp((b - x) / c))
    where a = plateau, b = turn-on position, c = width
    """
    tfile = ROOT.TFile.Open(root_file)
    if not tfile or tfile.IsZombie():
        return

    # Only plot DoubleMu Triggers OR
    trigger_name = "DoubleMu Triggers OR"

    # Kinematic variables to plot
    kin_vars = {
        "leading": ("p_{T}^{lead} [GeV]", 0, 100),
        "subleading": ("p_{T}^{sublead} [GeV]", 0, 100),
        "dipt": ("p_{T}^{#mu#mu} [GeV]", 0, 200),
    }

    # Define turn-on function using plain sigmoid (equation 19)
    # f(x) = a / (1 + exp((b - x) / c))
    # [0] = a (plateau, efficiency at saturation)
    # [1] = b (turn-on position)
    # [2] = c (width parameter; smaller = sharper turn-on)
    turnon_func_str = "[0] / (1.0 + TMath::Exp(([1] - x) / [2]))"

    # Plot for both inclusive and quality pair versions
    for suffix, suffix_label in [("", "inclusive"), ("_2mu", "quality pair")]:
        for var, (xlabel, xmin, xmax) in kin_vars.items():
            hname = f"h2_trigeff_{var}{suffix}_{trigger_name}"
            h2 = tfile.Get(hname)
            if not h2:
                continue

            # Create profile
            prof = h2.ProfileX(f"prof_{hname}")
            prof.SetDirectory(0)

            # Rebin by factor of 4 to smooth statistical fluctuations
            prof.Rebin(4)

            # Skip if profile is empty
            if prof.GetEntries() == 0:
                continue

            canvas = ROOT.TCanvas(f"c_trigeff_fit_{var}{suffix}_{tag}",
                                f"c_trigeff_fit_{var}{suffix}_{tag}", 800, 800)
            canvas.cd()
            canvas.SetRightMargin(0.05)
            canvas.SetTopMargin(0.10)

            # Style the profile
            prof.SetLineColor(ROOT.kBlack)
            prof.SetLineWidth(2)
            prof.SetMarkerColor(ROOT.kBlack)
            prof.SetMarkerStyle(20)
            prof.SetMarkerSize(0.9)
            prof.GetXaxis().SetTitle(xlabel)
            prof.GetYaxis().SetTitle("Trigger Efficiency wrt L1 H_{T} triggers")
            prof.GetXaxis().SetRangeUser(xmin, xmax)
            prof.SetMinimum(0.0)
            prof.SetMaximum(1.4)
            prof.SetTitle("")
            prof.Draw("E")

            # Two-stage fitting approach:
            # Stage 1: Fit plateau region to get accurate plateau value
            plateau_xmin = 15.0 if var in ["leading", "subleading"] else 30.0
            plateau_func = ROOT.TF1(f"plateau_{var}{suffix}", "[0]", plateau_xmin, xmax)
            plateau_func.SetParameter(0, 0.97)
            prof.Fit(plateau_func, "RQN0")
            plateau_value = plateau_func.GetParameter(0)

            # Stage 2: Fit turn-on region with plateau fixed
            turnon_xmin = 1.0 if var in ["leading", "subleading"] else 5.0
            turnon_xmax = 100.0 if var in ["leading", "subleading"] else 40.0

            fit_func = ROOT.TF1(f"turnon_{var}{suffix}", turnon_func_str, turnon_xmin, turnon_xmax)
            fit_func.SetParName(0, "a")
            fit_func.SetParName(1, "b")
            fit_func.SetParName(2, "c")

            # Fix plateau, set initial guesses for turn-on parameters
            fit_func.FixParameter(0, plateau_value)
            fit_func.SetParameter(1, 6.0)   # b (turn-on position)
            fit_func.SetParameter(2, 2.5)   # c (width)

            # Parameter limits for turn-on parameters
            fit_func.SetParLimits(1, 4.0, 9.0)    # b (turn-on position)
            fit_func.SetParLimits(2, 2.1, 6.5)    # c (width)

            fit_func.SetLineColor(ROOT.kRed)
            fit_func.SetLineWidth(2)

            # Perform fit
            fit_result = prof.Fit(fit_func, "RQS")

            # Extend function drawing range to show extrapolation
            fit_func.SetRange(turnon_xmin, turnon_xmax)

            # Draw the fit function (over fit range only)
            fit_func.Draw("same")

            # Draw horizontal line at 1.0
            line = ROOT.TLine(xmin, 1.0, xmax, 1.0)
            line.SetLineStyle(2)
            line.SetLineColor(ROOT.kGray+1)
            line.SetLineWidth(1)
            line.Draw()

            # Create legend
            legend = ROOT.TLegend(0.18, 0.70, 0.60, 0.88)
            legend.SetFillStyle(0)
            legend.SetBorderSize(0)
            legend.SetTextSize(0.025)

            # Parse tag to extract sample and era
            tag_parts = tag.split('_')
            if tag_parts[0] in ['sr', 'vr']:
                tag_parts = tag_parts[1:]
            sample_name = tag_parts[0] if len(tag_parts) > 0 else ""
            era_name = tag_parts[1] if len(tag_parts) > 1 else ""

            # Format header text
            suffix_text = "requiring a quality pair" if suffix_label == "quality pair" else "inclusive"
            legend.SetHeader(f"{trigger_name}", "L")
            legend.AddEntry(prof, f"{sample_name} {era_name} ({suffix_text})", "lpe")

            # Add fit results to legend if fit is valid
            if fit_result and fit_result.IsValid():
                import math

                a = fit_func.GetParameter(0)
                b = fit_func.GetParameter(1)
                c = fit_func.GetParameter(2)

                # Calculate x90: point where efficiency reaches 90% of plateau
                # For sigmoid f(x) = a / (1 + exp((b - x) / c))
                # At x90: 0.9*a = a / (1 + exp((b - x90) / c))
                # Solving: x90 = b + c * ln(9)
                ln9 = math.log(9.0)
                x90 = b + c * ln9

                legend.AddEntry(fit_func, "Turn-on fit", "l")
                legend.AddEntry(ROOT.nullptr, f"plateau = {a:.3f},  x_{{90}} = {x90:.2f} GeV,  width = {c:.2f} GeV", "")
            else:
                # If fit failed, still add function to legend but without parameters
                legend.AddEntry(fit_func, "Turn-on fit (failed)", "l")

            legend.Draw()

            addCMSText(canvas, lumi_text="13.6 TeV", extra_text="Work in Progress")

            varname_safe = var + suffix.replace("_", "")
            canvas.SaveAs(os.path.join(output_dir, f"trigeff_fit_{varname_safe}_{tag}.png"))
            canvas.SaveAs(os.path.join(output_dir, f"trigeff_fit_{varname_safe}_{tag}.pdf"))
            canvas.Update()
            canvas.Close()

    tfile.Close()


def process_file(filepath, output_dir):
    """Process a single ROOT file and produce all plots."""
    # Extract tag from filename
    tag = os.path.basename(filepath).replace("trigger_study_", "").replace(".root", "")

    # Note: Individual cutflow and trig_fired plots skipped - using overlay plots instead
    plot_mass_histograms(filepath, output_dir, tag)
    plot_trigger_efficiency(filepath, output_dir, tag)
    plot_trigger_efficiency_with_fit(filepath, output_dir, tag)


def plot_all_trigger_studies(input_dir, output_dir, region=""):
    """
    Main function to create all trigger study plots for files in a directory.

    Args:
        input_dir: Directory containing trigger_study ROOT files
        output_dir: Directory to save output figures
        region: Region label (sr, vr) for output naming
    """
    # Enable batch mode
    ROOT.gROOT.SetBatch(True)

    # Set CMS style
    setCMSStyle()

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Get list of ROOT files
    root_files = glob.glob(os.path.join(input_dir, "trigger_study_*.root"))

    # Filter out job-split files (keep only merged/hadded files)
    root_files = [f for f in root_files if "_job" not in os.path.basename(f)]

    if not root_files:
        print(f"No trigger_study ROOT files found in {input_dir}")
        return

    # Check if processing background MC (group all BkgMC together)
    is_bkgmc = 'BkgMC' in input_dir

    # Group files by sample base name for cutflow overlays
    cutflow_groups = {}
    for root_file in root_files:
        tag = os.path.basename(root_file).replace("trigger_study_", "").replace(".root", "")

        # Extract base sample name (remove region prefix and mass/variant suffix)
        tag_parts = tag.split('_')
        if tag_parts[0] in ['sr', 'vr']:
            tag_parts = tag_parts[1:]

        # For BkgMC: group all samples together
        if is_bkgmc:
            base_sample = "BkgMC"
        # For signal samples: group by process name (e.g., all TpTpTo2T2STo2Mu2B regardless of mass)
        elif tag_parts and ('TpTp' in tag_parts[0] or 'Par-' in tag):
            base_sample = tag_parts[0]
        # For data samples: group by base name (e.g., ScoutingPFRun3)
        else:
            if len(tag_parts) > 1:
                base_sample = '_'.join(tag_parts[:-1])
            else:
                base_sample = tag_parts[0] if tag_parts else tag

        if base_sample not in cutflow_groups:
            cutflow_groups[base_sample] = []
        cutflow_groups[base_sample].append((root_file, tag))

    # Create overlay cutflow plots for each group (including single-file groups)
    overlay_groups = [(bs, fl) for bs, fl in cutflow_groups.items()]
    for base_sample, file_list in tqdm(overlay_groups, desc="Creating cutflow overlays", disable=len(overlay_groups)==0):
        plot_cutflow_overlay(file_list, output_dir, base_sample, normalize=False)
        plot_cutflow_overlay(file_list, output_dir, base_sample, normalize=True)

    # Create overlay trig_fired plots for each group (including single-file groups)
    for base_sample, file_list in tqdm(overlay_groups, desc="Creating trig_fired overlays", disable=len(overlay_groups)==0):
        plot_trig_fired_overlay(file_list, output_dir, base_sample, normalize=False)
        plot_trig_fired_overlay(file_list, output_dir, base_sample, normalize=True)

    # Process each file
    for root_file in tqdm(root_files, desc="Processing individual files"):
        process_file(root_file, output_dir)


if __name__ == "__main__":
    # Base directories and their labels
    base_dirs = {
        "/ceph/cms/store/user/tvami/DiMuonPlusX/Ntuples_v1.0.1_trigger_study/Data/": "Data",
        "/ceph/cms/store/user/tvami/DiMuonPlusX/Ntuples_v1.0.1_trigger_study/Signal/": "Signal",
        "/ceph/cms/store/user/tvami/DiMuonPlusX/Ntuples_v1.0.1_trigger_study/BkgMC/": "BkgMC",
    }

    # Define regions to loop over
    regions = ['sr', 'vr']

    # Loop over all combinations
    for base_dir, sample_label in base_dirs.items():
        # Extract Ntuple version from path (e.g., "v1.0.1_trigger_study" from "Ntuples_v1.0.1_trigger_study")
        version_match = re.search(r'Ntuples_(v[\d.]+(?:_trigger_study)?)', base_dir)
        version = version_match.group(1) if version_match else "unknown"
        output_base = f"figures/trigger_study_{version}"
        os.makedirs(output_base, exist_ok=True)

        for region in regions:
            input_dir = os.path.join(base_dir, region)

            # Check if directory exists
            if not os.path.exists(input_dir):
                print(f"Warning: Directory {input_dir} does not exist, skipping...")
                continue

            print(f"\n{'='*100}")
            print(f"Processing: {sample_label}/{region}")
            print(f"{'='*100}\n")

            # Create output subdirectory for this sample and region
            output_dir = os.path.join(output_base, sample_label, region)
            plot_all_trigger_studies(input_dir, output_dir, region=region)

    print("\nDone!")
