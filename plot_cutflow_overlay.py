#!/usr/bin/env python3
import ROOT
import os
import glob
import re

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

def addCMSText(canvas, lumi_text="13.6 TeV", extra_text="Preliminary"):
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

def extract_parameters(filename):
    """Extract MS (scalar mass), MTp, and decay channel from filename"""
    params = {}

    # Extract MS (scalar mass) -> stored as minp for compatibility
    match = re.search(r'-MS(\d+)', filename)
    if match:
        params['minp'] = int(match.group(1))

    # Extract MTp (top partner mass)
    match = re.search(r'-MTp(\d+)', filename)
    if match:
        params['mtp'] = int(match.group(1))

    # Extract decay channel (e.g., 2Mu2B, 2Mu2G, 2MuInv)
    match = re.search(r'TpTpTo2T2STo(\w+?)_', filename)
    if match:
        params['channel'] = match.group(1)

    return params

def get_grouping_key(params):
    """
    Create a grouping key from parameters (excluding MS which is the x-axis)
    Returns a tuple of (channel, mtp)
    """
    return (params.get('channel'), params.get('mtp'))

def format_group_label(group_key):
    """Create a readable label from the grouping key"""
    channel, mtp = group_key

    parts = []
    if channel is not None:
        parts.append(f"{channel}")
    if mtp is not None:
        parts.append(f"M_{{Tp}}={mtp} GeV")

    if not parts:
        return "Default parameters"

    return ", ".join(parts)

def plot_all_cutflow_analysis(input_dir="skimmed_volt2", hist_name="h_cutflow", output_name="cutflow_overlay.png"):
    """
    Create all cutflow plots: 1D overlay, 2D heatmap, and efficiency vs MinP
    
    Args:
        input_dir: Directory containing ROOT files
        hist_name: Name of the histogram to overlay
        output_name: Base output filename for the plots
    """
    
    # Enable batch mode
    ROOT.gROOT.SetBatch(True)
    
    # Set CMS style
    setCMSStyle()
    
    # Get list of ROOT files
    root_files = glob.glob(os.path.join(input_dir, "*.root"))
    
    if not root_files:
        print(f"No ROOT files found in {input_dir}")
        return
    
    # Sort files by MinP value
    def extract_minp(filename):
        """Extract MS value from filename for sorting"""
        match = re.search(r'-MS(\d+)', filename)
        if match:
            return int(match.group(1))
        return 0  # Default value if MS not found
    
    def apply_minp_linestyle(obj, minp):
        if minp is None:
            return
        obj.SetLineStyle(2 if minp > 9999 else 1)

    
    root_files = [f for f in root_files if "job" not in os.path.basename(f)]
    root_files = sorted(root_files, key=extract_minp)

    print(f"Found {len(root_files)} ROOT files")

    # Color palette - expanded for 30+ curves
    colors = [
        ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta,
        ROOT.kCyan+1, ROOT.kOrange, ROOT.kViolet, ROOT.kTeal-1, ROOT.kPink+1,
        ROOT.kAzure+7, ROOT.kSpring-7, ROOT.kYellow+2, ROOT.kRed+2, ROOT.kBlue+2,
        ROOT.kGreen+3, ROOT.kMagenta+2, ROOT.kCyan+3, ROOT.kOrange+7, ROOT.kViolet+2,
        ROOT.kTeal+2, ROOT.kPink-2, ROOT.kAzure-3, ROOT.kSpring+3, ROOT.kYellow-2,
        ROOT.kRed-3, ROOT.kBlue-3, ROOT.kGreen-2, ROOT.kMagenta-3, ROOT.kOrange-3,
        ROOT.kViolet-3, ROOT.kTeal-3, ROOT.kPink+3, ROOT.kAzure+2, ROOT.kSpring-2,
        ROOT.kYellow+3, ROOT.kRed+3, ROOT.kBlue+3, ROOT.kGreen-3, ROOT.kMagenta+3, ROOT.kCyan-7,
        ROOT.kOrange+3, ROOT.kViolet+3, ROOT.kTeal+3, ROOT.kPink-3, ROOT.kAzure-2, ROOT.kSpring+2,
        ROOT.kYellow-3, ROOT.kRed-2, ROOT.kBlue-2, ROOT.kGreen+1, ROOT.kMagenta-2, ROOT.kCyan+7, ROOT.kOrange-7
    ]

    # Build mappings: MS -> color
    all_minp_values = set()
    file_params = []  # Store params for each file in order
    for root_file in root_files:
        params = extract_parameters(os.path.basename(root_file))
        file_params.append(params)
        if 'minp' in params:
            all_minp_values.add(params['minp'])

    sorted_minp = sorted(all_minp_values)
    minp_to_color = {minp: colors[j % len(colors)] for j, minp in enumerate(sorted_minp)}
    
    # Extract region and sample label for unique canvas names
    path_parts = input_dir.rstrip('/').split('/')
    canvas_suffix = ""
    for part in reversed(path_parts):
        if part in ['vr', 'sr', 'Signal', 'Data', 'BkgMC']:
            canvas_suffix = part + "_" + canvas_suffix
    canvas_suffix = (canvas_suffix + hist_name).rstrip('_')
    
    histograms = []
    histogram_labels = []  # Store labels for reuse
    
    # Loop over files and get histograms
    for i, root_file in enumerate(root_files):
        try:
            tfile = ROOT.TFile.Open(root_file)
        except (OSError, RuntimeError) as e:
            print(f"Warning: Could not open {root_file}: {e}")
            continue
        
        if not tfile or tfile.IsZombie():
            print(f"Warning: Zombie file {root_file}")
            continue
        
        hist = tfile.Get(hist_name)
        if not hist:
            print(f"Warning: {hist_name} not found in {root_file}")
            tfile.Close()
            continue
        
        # Clone histogram to keep it after file closes
        hist_clone = hist.Clone(f"h_{i}")
        hist_clone.SetDirectory(0)
        
        # Extract meaningful label from filename
        filename = os.path.basename(root_file)
        label = filename.replace("skimmed_", "").replace(".root", "")
        
        # Clean up label for readability
        label = re.sub(r'^(sr|vr)_', '', label)
        label = label.replace("_TuneCP5_13p6TeV_madgraph-pythia8", "")
        label = label.replace("_TuneCP5_13p6TeV_powheg-pythia8", "")
        label = label.replace("_TuneCP5_13p6TeV_pythia8", "")
        label = label.replace("TpTpTo2T2STo", "")
        label = label.replace("_Par-", " ")
        label = label.replace("-", " ")
        label = label.replace("_", " ")
        
        # Set style: color by MS
        params = file_params[i]
        color = minp_to_color.get(params.get('minp'), colors[i % len(colors)])
        hist_clone.SetLineColor(color)
        hist_clone.SetLineWidth(2)
        hist_clone.SetMarkerColor(color)
        apply_minp_linestyle(hist_clone, params.get('minp'))
        hist_clone.SetMarkerSize(0.8)

        histograms.append((hist_clone, tfile))
        histogram_labels.append(label)
    
    if not histograms:
        print("No valid histograms found!")
        return
    
    # Extract region from input_dir for output filenames
    path_parts = input_dir.rstrip('/').split('/')
    region = ""
    for part in reversed(path_parts):
        if part in ['vr', 'sr']:
            region = part
            break

    # Create base name with region
    base_name = output_name.replace(".png", "").replace(".pdf", "")
    if region:
        base_name = f"{base_name}_{region}"
    
    # ===== Create normalized version (normalized to first bin) =====
    canvas2 = ROOT.TCanvas(f"c2_{canvas_suffix}_{region}", f"c2_{canvas_suffix}_{region}", 800, 800)
    canvas2.cd()
    
    # Create new legend for normalized plot - spanning full width
    legend2 = ROOT.TLegend(0.17, 0.5, 0.9, 0.9)
    legend2.SetHeader(f"Region: {'pre-SR' if region.upper() == 'SR' else region.upper()}")
    legend2.SetTextSize(0.03)
    legend2.SetFillStyle(0)
    legend2.SetBorderSize(0)
    legend2.SetTextSize(0.015)
    legend2.SetNColumns(3)
    
    histograms_norm = []
    max_y_norm = 0
    min_y_norm = 1.0
    
    # Normalize each histogram to its first bin
    for i, (hist, tfile) in enumerate(histograms):
        hist_norm = hist.Clone(f"h_norm_{i}")
        
        # Get first bin content
        first_bin_content = hist_norm.GetBinContent(1)
        
        if first_bin_content > 0:
            # Normalize to first bin
            hist_norm.Scale(1.0 / first_bin_content)
        
        # Track maximum and minimum (excluding zeros)
        for bin_i in range(1, hist_norm.GetNbinsX() + 1):
            bin_content = hist_norm.GetBinContent(bin_i)
            if bin_content > 0:
                if bin_content > max_y_norm:
                    max_y_norm = bin_content
                if bin_content < min_y_norm:
                    min_y_norm = bin_content
        
        histograms_norm.append(hist_norm)
    
    if not histograms_norm:
        print("No valid normalized histograms found!")
        return
    
    # Set axis range on first histogram BEFORE drawing anything
    first_hist_norm = histograms_norm[0]
    first_hist_norm.GetYaxis().SetTitle("Efficiency")
    first_hist_norm.GetXaxis().SetTitle("")
    
    # Use fixed Y-axis range for normalized plot
    first_hist_norm.SetMaximum(1e5)
    first_hist_norm.SetMinimum(1e-4)
    
    # Set log scale AFTER setting range
    canvas2.SetLogy()
    
    # Draw normalized histograms
    for i, hist_norm in enumerate(histograms_norm):
        draw_option = "HISTTEXT0" if i == 0 else "HISTTEXT0 SAME"
        hist_norm.Draw(draw_option)
        
        # Use the stored label
        label = histogram_labels[i]
        
        # Add to legend
        legend2.AddEntry(hist_norm, label, "l")
    
    # Draw legend
    legend2.Draw()
    
    # Add CMS text
    addCMSText(canvas2, lumi_text="13.6 TeV", extra_text="Work in Progress")
    
    # Save normalized canvas
    png_name_norm = base_name + "_normalized.png"
    pdf_name_norm = base_name + "_normalized.pdf"
    
    canvas2.SaveAs(png_name_norm)
    canvas2.SaveAs(pdf_name_norm)
    
    canvas2.Update()
    
    # ===== Create 2D histogram (cutflow vs samples) =====
    canvas3 = ROOT.TCanvas(f"c3_{canvas_suffix}_{region}", f"c3_{canvas_suffix}_{region}", 1200, 800)
    canvas3.cd()
    canvas3.SetRightMargin(0.15)
    canvas3.SetLeftMargin(0.2)
    canvas3.SetBottomMargin(0.07)
    
    # Get number of bins from first histogram
    n_bins = histograms_norm[0].GetNbinsX()
    n_samples = len(histograms_norm)
    
    # Create 2D histogram
    h2d = ROOT.TH2F("h2d_cutflow", "Cutflow vs Samples", 
                    n_bins, 0.5, n_bins + 0.5,
                    n_samples, 0.5, n_samples + 0.5)
    
    # Set X-axis bin labels (cutflow steps)
    for bin_x in range(1, n_bins + 1):
        h2d.GetXaxis().SetBinLabel(bin_x, histograms_norm[0].GetXaxis().GetBinLabel(bin_x))

    # Fill the 2D histogram and set Y-axis bin labels
    for i, hist_norm in enumerate(histograms_norm):
        # Use the stored label (already sorted by MinP)
        label = histogram_labels[i]
        
        # Set Y-axis label for this sample
        h2d.GetYaxis().SetBinLabel(i + 1, label)

        # Fill each bin
        for bin_x in range(1, n_bins + 1):
            yield_value = hist_norm.GetBinContent(bin_x)
            h2d.SetBinContent(bin_x, i + 1, yield_value)
    
    # Set axis titles
    h2d.GetZaxis().SetTitle("Efficiency")
    
    # Set label size for readability
    h2d.GetYaxis().SetLabelSize(0.02)
    h2d.GetXaxis().SetLabelSize(0.04)
    
    # Set Z-axis range
    h2d.SetMinimum(1e-4)
    h2d.SetMaximum(1.0)
    
    # Draw with COLZ option for color scale, and TEXT for bin values
    ROOT.gStyle.SetPalette(ROOT.kBird)
    ROOT.gStyle.SetPaintTextFormat("4.3f")
    h2d.Draw("COLZ TEXT")
    
    # Set log scale for Z-axis
    canvas3.SetLogz()
    
    # Add CMS text
    addCMSText(canvas3, lumi_text="13.6 TeV", extra_text="Work in Progress")
    
    # Save 2D plot
    png_name_2d = base_name + "_2d_cutflow.png"
    pdf_name_2d = base_name + "_2d_cutflow.pdf"
    
    canvas3.SaveAs(png_name_2d)
    canvas3.SaveAs(pdf_name_2d)
    
    canvas3.Update()
    
    # ===== Create efficiency vs MinP plots =====
    # Skip for BkgMC (MinP is not meaningful for background)

    # Dictionary to organize data by parameter combinations
    data_by_group = {}

    # Loop over files again and extract efficiency data
    for root_file in root_files:
        # Extract parameters from filename
        filename = os.path.basename(root_file)
        params = extract_parameters(filename)
        
        if 'minp' not in params:
            continue
        
        minp = params['minp']
        group_key = get_grouping_key(params)
        
        # Open ROOT file and get histogram
        try:
            tfile = ROOT.TFile.Open(root_file)
        except (OSError, RuntimeError) as e:
            continue
        
        if not tfile or tfile.IsZombie():
            continue
        
        hist = tfile.Get(hist_name)
        if not hist:
            tfile.Close()
            continue
        
        # Get efficiency from the last cut step
        n_bins = hist.GetNbinsX()
        first_bin = hist.GetBinContent(1)
        last_bin = hist.GetBinContent(n_bins)
        
        if first_bin > 0:
            efficiency = last_bin / first_bin
        else:
            efficiency = 0

        if (efficiency < 0.01 and region.upper() == "SR" and "BkgMC" not in input_dir):
            print(f"Warning: Efficiency for {root_file} is very low ({efficiency:.4f})")
        
        # Store data
        if group_key not in data_by_group:
            data_by_group[group_key] = []
        
        data_by_group[group_key].append((minp, efficiency))
        
        tfile.Close()
    
    if data_by_group and "BkgMC" not in input_dir:
        # Print grouping information
        # print(f"Found {len(data_by_group)} different parameter combinations:")
        # for group_key in data_by_group:
        #     label = format_group_label(group_key)
        #     print(f"  - {label}: {len(data_by_group[group_key])} points")
        
        # Sort data for each group
        for group_key in data_by_group:
            data_by_group[group_key].sort(key=lambda x: x[0])
        
        # Create canvas for efficiency vs MinP
        canvas_eff = ROOT.TCanvas("c_eff_vs_ms"+canvas_suffix, "c_eff_vs_ms"+canvas_suffix, 800, 800)
        canvas_eff.cd()
        canvas_eff.SetRightMargin(0.1)
        
        # Marker styles
        markers = [20, 21, 22, 23, 29, 33, 34, 47, 43, 45]
        
        # Create TGraphs for each parameter combination
        graphs = []
        legend_eff = ROOT.TLegend(0.32, 0.5, 0.83, 0.9)
        legend_eff.SetHeader(f"Region: {'pre-SR' if region.upper() == 'SR' else region.upper()}")
        legend_eff.SetTextSize(0.03)
        legend_eff.SetFillStyle(0)
        legend_eff.SetBorderSize(0)
        legend_eff.SetTextSize(0.03)
        
        # Sort groups for consistent ordering
        def sort_key(group_tuple):
            channel, mtp = group_tuple
            return (channel if channel is not None else "",
                    mtp if mtp is not None else -1)

        sorted_groups = sorted(data_by_group.keys(), key=sort_key)
        
        for i, group_key in enumerate(sorted_groups):
            data = data_by_group[group_key]
            n_points = len(data)
            
            # Create arrays for TGraph
            import array
            x_array = array.array('d', [point[0] for point in data])
            y_array = array.array('d', [point[1] for point in data])
            
            graph = ROOT.TGraph(n_points, x_array, y_array)
            
            # Set style
            color = colors[i % len(colors)]
            marker = markers[i % len(markers)]
            
            graph.SetLineColor(color)
            graph.SetLineWidth(2)
            graph.SetMarkerColor(color)
            graph.SetMarkerStyle(marker)
            graph.SetMarkerSize(1.2)
            
            # Create label
            label = format_group_label(group_key)

            legend_eff.AddEntry(graph, label, "lp")
            graphs.append(graph)
        
        if graphs:
            # Draw graphs
            first_graph = graphs[0]
            first_graph.SetTitle("")
            first_graph.GetXaxis().SetTitle("M_{S} [GeV]")
            first_graph.GetYaxis().SetTitle("Pre-selection efficiency")
            first_graph.GetYaxis().SetTitleOffset(1.05)
            
            # Set axis ranges
            first_graph.GetYaxis().SetRangeUser(0, 0.6)
            first_graph.GetXaxis().SetMaxDigits(3)  # Use scientific notation for large numbers
            first_graph.GetXaxis().SetNdivisions(510) # More divisions for better tick marks
            
            first_graph.Draw("ALP")
            
            for graph in graphs[1:]:
                graph.Draw("LP SAME")
            
            # Draw legend
            legend_eff.Draw()
            
            # Add CMS text
            addCMSText(canvas_eff, lumi_text="13.6 TeV", extra_text="Work in Progress")
            
            # Save canvas
            png_name_eff = base_name.replace("_overlay", "_efficiency_vs_ms") + ".png"
            pdf_name_eff = base_name.replace("_overlay", "_efficiency_vs_ms") + ".pdf"
            
            canvas_eff.SaveAs(png_name_eff)
            canvas_eff.SaveAs(pdf_name_eff)
            
            canvas_eff.Update()

    # ===== Create trigger efficiency vs MinP plots =====
    # print("\nCreating trigger efficiency vs MinP plots...")

    # Dictionary to organize trigger efficiency data by parameter combinations
    trigger_data_by_group = {}

    # Loop over files again and extract trigger efficiency data
    for root_file in root_files:
        # Extract parameters from filename
        filename = os.path.basename(root_file)
        params = extract_parameters(filename)

        if 'minp' not in params:
            continue

        minp = params['minp']
        group_key = get_grouping_key(params)

        # Open ROOT file and get histogram
        try:
            tfile = ROOT.TFile.Open(root_file)
        except (OSError, RuntimeError) as e:
            continue

        if not tfile or tfile.IsZombie():
            continue

        hist = tfile.Get(hist_name)
        if not hist:
            tfile.Close()
            continue

        # Get trigger efficiency (2nd bin / 1st bin)
        first_bin = hist.GetBinContent(1)  # All events
        second_bin = hist.GetBinContent(2)  # After trigger

        if first_bin > 0:
            trigger_efficiency = second_bin / first_bin
        else:
            trigger_efficiency = 0

        # Store data
        if group_key not in trigger_data_by_group:
            trigger_data_by_group[group_key] = []

        trigger_data_by_group[group_key].append((minp, trigger_efficiency))

        tfile.Close()

    if trigger_data_by_group and "BkgMC" not in input_dir:
        # Sort data for each group
        for group_key in trigger_data_by_group:
            trigger_data_by_group[group_key].sort(key=lambda x: x[0])

        # Create canvas for trigger efficiency vs MinP
        canvas_trig = ROOT.TCanvas("c_trig_eff_vs_ms"+canvas_suffix, "c_trig_eff_vs_ms"+canvas_suffix, 800, 800)
        canvas_trig.cd()
        canvas_trig.SetRightMargin(0.1)


        # Marker styles
        markers = [20, 21, 22, 23, 29, 33, 34, 47, 43, 45]

        # Create TGraphs for each parameter combination
        trig_graphs = []
        legend_trig = ROOT.TLegend(0.32, 0.15, 0.83, 0.55)
        legend_trig.SetHeader(f"Region: {'pre-SR' if region.upper() == 'SR' else region.upper()}")
        legend_trig.SetTextSize(0.03)
        legend_trig.SetFillStyle(0)
        legend_trig.SetBorderSize(0)
        legend_trig.SetTextSize(0.03)

        # Sort groups for consistent ordering
        def sort_key_trig(group_tuple):
            channel, mtp = group_tuple
            return (channel if channel is not None else "",
                    mtp if mtp is not None else -1)

        sorted_groups_trig = sorted(trigger_data_by_group.keys(), key=sort_key_trig)

        for i, group_key in enumerate(sorted_groups_trig):
            data = trigger_data_by_group[group_key]
            n_points = len(data)

            # Create arrays for TGraph
            import array
            x_array = array.array('d', [point[0] for point in data])
            y_array = array.array('d', [point[1] for point in data])

            graph = ROOT.TGraph(n_points, x_array, y_array)

            # Set style
            color = colors[i % len(colors)]
            marker = markers[i % len(markers)]

            graph.SetLineColor(color)
            graph.SetLineWidth(2)
            graph.SetMarkerColor(color)
            graph.SetMarkerStyle(marker)
            graph.SetMarkerSize(1.2)

            # Create label
            label = format_group_label(group_key)

            legend_trig.AddEntry(graph, label, "lp")
            trig_graphs.append(graph)

        if trig_graphs:
            # Draw graphs
            first_graph_trig = trig_graphs[0]
            first_graph_trig.SetTitle("")
            first_graph_trig.GetXaxis().SetTitle("M_{S} [GeV]")
            first_graph_trig.GetYaxis().SetTitle("Trigger efficiency")
            first_graph_trig.GetYaxis().SetTitleOffset(1.05)

            # Set axis ranges
            first_graph_trig.GetYaxis().SetRangeUser(0, 1.0)
            first_graph_trig.GetXaxis().SetMaxDigits(3)  # Use scientific notation for large numbers
            first_graph_trig.GetXaxis().SetNdivisions(510) # More divisions for better tick marks

            first_graph_trig.Draw("ALP")

            for graph in trig_graphs[1:]:
                graph.Draw("LP SAME")

            # Draw legend
            legend_trig.Draw()

            # Add CMS text
            addCMSText(canvas_trig, lumi_text="13.6 TeV", extra_text="Work in Progress")

            # Save canvas
            png_name_trig = base_name.replace("_overlay", "_trigger_efficiency_vs_ms") + ".png"
            pdf_name_trig = base_name.replace("_overlay", "_trigger_efficiency_vs_ms") + ".pdf"

            canvas_trig.SaveAs(png_name_trig)
            canvas_trig.SaveAs(pdf_name_trig)

            canvas_trig.Update()

    # ===== Create mean object count at last cutflow bin vs MinP plots =====
    objcount_data_by_group = {}

    for root_file in root_files:
        filename = os.path.basename(root_file)
        params = extract_parameters(filename)

        if 'minp' not in params:
            continue

        minp = params['minp']
        group_key = get_grouping_key(params)

        try:
            tfile = ROOT.TFile.Open(root_file)
        except (OSError, RuntimeError) as e:
            continue

        if not tfile or tfile.IsZombie():
            continue

        h2 = tfile.Get("h2_cutflow_nobj")
        if not h2:
            tfile.Close()
            continue

        # ProfileX gives mean N_obj at each cutflow step
        prof = h2.ProfileX(f"prof_objcount_{minp}")

        # Get mean N_obj at last bin
        last_bin = prof.GetNbinsX()
        mean_nobj = prof.GetBinContent(last_bin)

        if group_key not in objcount_data_by_group:
            objcount_data_by_group[group_key] = []

        objcount_data_by_group[group_key].append((minp, mean_nobj))

        tfile.Close()

    if objcount_data_by_group and "BkgMC" not in input_dir:
        for group_key in objcount_data_by_group:
            objcount_data_by_group[group_key].sort(key=lambda x: x[0])

        canvas_objcount = ROOT.TCanvas("c_objcount_vs_ms"+canvas_suffix, "c_objcount_vs_ms"+canvas_suffix, 800, 800)
        canvas_objcount.cd()
        canvas_objcount.SetRightMargin(0.1)

        markers = [20, 21, 22, 23, 29, 33, 34, 47, 43, 45]

        objcount_graphs = []
        legend_objcount = ROOT.TLegend(0.32, 0.5, 0.83, 0.9)
        legend_objcount.SetHeader(f"Region: {'pre-SR' if region.upper() == 'SR' else region.upper()}")
        legend_objcount.SetTextSize(0.03)
        legend_objcount.SetFillStyle(0)
        legend_objcount.SetBorderSize(0)
        legend_objcount.SetTextSize(0.03)

        def sort_key_obj(group_tuple):
            channel, mtp = group_tuple
            return (channel if channel is not None else "",
                    mtp if mtp is not None else -1)

        sorted_groups_obj = sorted(objcount_data_by_group.keys(), key=sort_key_obj)

        for i, group_key in enumerate(sorted_groups_obj):
            data = objcount_data_by_group[group_key]
            n_points = len(data)

            import array
            x_array = array.array('d', [point[0] for point in data])
            y_array = array.array('d', [point[1] for point in data])

            graph = ROOT.TGraph(n_points, x_array, y_array)

            color = colors[i % len(colors)]
            marker = markers[i % len(markers)]

            graph.SetLineColor(color)
            graph.SetLineWidth(2)
            graph.SetMarkerColor(color)
            graph.SetMarkerStyle(marker)
            graph.SetMarkerSize(1.2)

            label = format_group_label(group_key)
            legend_objcount.AddEntry(graph, label, "lp")
            objcount_graphs.append(graph)

        if objcount_graphs:
            first_graph_obj = objcount_graphs[0]
            first_graph_obj.SetTitle("")
            first_graph_obj.GetXaxis().SetTitle("M_{S} [GeV]")
            first_graph_obj.GetYaxis().SetTitle("<N_{obj}> at last cutflow step")
            first_graph_obj.GetYaxis().SetTitleOffset(1.05)

            all_y = [pt[1] for data in objcount_data_by_group.values() for pt in data]
            max_y = max(all_y) if all_y else 5
            first_graph_obj.GetYaxis().SetRangeUser(0, max_y * 3.0)
            first_graph_obj.GetXaxis().SetMaxDigits(3)
            first_graph_obj.GetXaxis().SetNdivisions(510)

            first_graph_obj.Draw("ALP")

            for graph in objcount_graphs[1:]:
                graph.Draw("LP SAME")

            legend_objcount.Draw()
            addCMSText(canvas_objcount, lumi_text="13.6 TeV", extra_text="Work in Progress")

            png_name_objcount = base_name.replace("_overlay", "_objcount_vs_ms") + ".png"
            pdf_name_objcount = base_name.replace("_overlay", "_objcount_vs_ms") + ".pdf"

            canvas_objcount.SaveAs(png_name_objcount)
            canvas_objcount.SaveAs(pdf_name_objcount)
            canvas_objcount.Update()

    # ===== Create overlay plots for final kinematic distributions =====
    # print("\nCreating final kinematic distribution overlays...")
    
    # Define histograms to overlay
    kinematic_hists = {
        # 'h_pt_final': {'title': 'p_{T} [GeV]', 'log_y': True},
        # 'h_pt_trigger': {'title': 'p_{T} [GeV]', 'log_y': True},
        # 'h_pt_pretrigger': {'title': 'p_{T} [GeV]', 'log_y': True},
        # 'h_eta_final': {'title': '#eta', 'log_y': True},
        # 'h_eta_trigger': {'title': '#eta', 'log_y': True},
        # 'h_eta_pretrigger': {'title': '#eta', 'log_y': True},
        # 'h_phi_final': {'title': '#phi', 'log_y': True},
        # 'h_phi_trigger': {'title': '#phi', 'log_y': True},
        # 'h_phi_pretrigger': {'title': '#phi', 'log_y': True},
        # 'h_ntrack_final': {'title': 'N_{tracks}', 'log_y': True},
        # 'h_ntrack_trigger': {'title': 'N_{tracks}', 'log_y': True},
        # 'h_ntrack_nminus1': {'title': 'N_{tracks}', 'log_y': True},
        # 'h_nhits_final': {'title': 'N_{valid hits}', 'log_y': True},
        # 'h_nhits_trigger': {'title': 'N_{valid hits}', 'log_y': True},
        # 'h_nhits_nminus1': {'title': 'N_{valid hits}', 'log_y': True},
        # 'h_chi2ndof_final': {'title': '#chi^{2}/ndof', 'log_y': True},
        # 'h_chi2ndof_trigger': {'title': '#chi^{2}/ndof', 'log_y': True},
        # 'h_chi2ndof_nminus1': {'title': '#chi^{2}/ndof', 'log_y': True}
    }
    
    for hist_name_kin, hist_config in kinematic_hists.items():
        # print(f"  Processing {hist_name_kin}...")
        
        # Collect histograms from all files
        kin_histograms = []
        kin_labels = []
        
        for i, root_file in enumerate(root_files):
            try:
                tfile = ROOT.TFile.Open(root_file)
            except (OSError, RuntimeError) as e:
                continue
            
            if not tfile or tfile.IsZombie():
                continue
            
            hist_kin = tfile.Get(hist_name_kin)
            if not hist_kin:
                tfile.Close()
                continue
            
            # Clone histogram
            hist_kin_clone = hist_kin.Clone(f"{hist_name_kin}_{i}")
            hist_kin_clone.SetDirectory(0)
            
            # Normalize to unit area if histogram is not empty
            integral = hist_kin_clone.Integral()
            if integral > 0:
                hist_kin_clone.Scale(1.0 / integral)

            # Set style: color by MS
            params = file_params[i]
            color = minp_to_color.get(params.get('minp'), colors[i % len(colors)])
            hist_kin_clone.SetLineColor(color)
            hist_kin_clone.SetLineWidth(2)
            hist_kin_clone.SetMarkerColor(color)
            apply_minp_linestyle(hist_kin_clone, params.get('minp'))
            hist_kin_clone.SetMarkerSize(0.8)
            
            # Use the same labels as before
            kin_histograms.append((hist_kin_clone, tfile))
            kin_labels.append(histogram_labels[i] if i < len(histogram_labels) else "")
        
        if not kin_histograms:
            print(f"  Warning: No valid histograms found for {hist_name_kin}")
            continue
        
        # Create canvas
        canvas_kin = ROOT.TCanvas(f"c_{hist_name_kin}_{canvas_suffix}_{i}", f"c_{hist_name_kin}_{canvas_suffix}_{i}", 800, 800)
        canvas_kin.cd()
        
        # Set log scale if requested
        if hist_config['log_y']:
            canvas_kin.SetLogy()
        
        # Create legend
        legend_kin = ROOT.TLegend(0.17, 0.5, 0.9, 0.9)
        legend_kin.SetHeader(f"Region: {'pre-SR' if region.upper() == 'SR' else region.upper()}")
        legend_kin.SetTextSize(0.03)
        legend_kin.SetFillStyle(0)
        legend_kin.SetBorderSize(0)
        legend_kin.SetTextSize(0.03)
        legend_kin.SetTextSize(0.015)
        legend_kin.SetNColumns(3)
        
        # Find maximum for y-axis range
        max_val = 0
        for hist_kin_clone, _ in kin_histograms:
            hist_max = hist_kin_clone.GetMaximum()
            if hist_max > max_val:
                max_val = hist_max
        
        # Draw histograms
        for i, (hist_kin_clone, tfile) in enumerate(kin_histograms):
            # Set axis titles and range on first histogram
            if i == 0:
                hist_kin_clone.GetXaxis().SetTitle(hist_config['title'])
                hist_kin_clone.GetYaxis().SetTitle("Normalized Events")
                if hist_config['log_y']:
                    hist_kin_clone.SetMaximum(max_val * 100000)
                    hist_kin_clone.SetMinimum(max_val * 1e-4)
                else:
                    hist_kin_clone.SetMaximum(max_val * 1.2)
                
                draw_option = "HIST"
            else:
                draw_option = "HIST SAME"
            
            hist_kin_clone.Draw(draw_option)
            legend_kin.AddEntry(hist_kin_clone, kin_labels[i], "l")
        
        # Draw legend
        legend_kin.Draw()
        
        # Add CMS text
        addCMSText(canvas_kin, lumi_text="13.6 TeV", extra_text="Work in Progress")
        
        # Save canvas
        png_name_kin = base_name.replace("_overlay", f"_{hist_name_kin}") + ".png"
        pdf_name_kin = base_name.replace("_overlay", f"_{hist_name_kin}") + ".pdf"
        
        canvas_kin.SaveAs(png_name_kin)
        canvas_kin.SaveAs(pdf_name_kin)
        canvas_kin.Update()
        
        # Cleanup
        for hist_kin_clone, tfile in kin_histograms:
            tfile.Close()
    
    # ===== Create ProfileX overlays for 2D trigger histograms =====
    # print("\nCreating ProfileX overlays for 2D trigger histograms...")

    profile_hists = {
        # 'h_eta_vs_trigger': {'x_title': '#eta', 'y_title': '<Trigger>', 'log_y': False, 'rebin': 20},
        # 'h_pt_vs_trigger':  {'x_title': 'p_{T} [GeV]', 'y_title': '<Trigger>', 'log_y': False, 'rebin': 20},
        # 'h_phi_vs_trigger': {'x_title': '#phi', 'y_title': '<Trigger>', 'log_y': False, 'rebin': 20},
        # 'h_eta_vs_trigger_obj_eta': {'x_title': '#eta', 'y_title': '<Trigger> (obj found, |#eta|<0.9)', 'log_y': False, 'rebin': 20},
        # 'h_pt_vs_trigger_obj_eta':  {'x_title': 'p_{T} [GeV]', 'y_title': '<Trigger> (obj found, |#eta|<0.9)', 'log_y': False, 'rebin': 20},
        # 'h_phi_vs_trigger_obj_eta': {'x_title': '#phi', 'y_title': '<Trigger> (obj found, |#eta|<0.9)', 'log_y': False, 'rebin': 20},
        'h2_trigeff_leading':    {'x_title': 'p_{T}^{lead} [GeV]',    'y_title': 'L1 Trigger OR efficiency', 'log_y': False, 'rebin': 5},
        'h2_trigeff_subleading': {'x_title': 'p_{T}^{sublead} [GeV]', 'y_title': 'L1 Trigger OR efficiency', 'log_y': False, 'rebin': 5},
        'h2_trigeff_dipt':       {'x_title': 'p_{T}^{#mu#mu} [GeV]',  'y_title': 'L1 Trigger OR efficiency', 'log_y': False, 'rebin': 5},
    }

    for hist_name_2d, prof_config in profile_hists.items():
        # print(f"  Processing ProfileX of {hist_name_2d}...")

        prof_histograms = []
        prof_labels = []

        for i, root_file in enumerate(root_files):
            try:
                tfile = ROOT.TFile.Open(root_file)
            except (OSError, RuntimeError):
                continue

            if not tfile or tfile.IsZombie():
                continue

            h2 = tfile.Get(hist_name_2d)
            if not h2:
                tfile.Close()
                continue

            # Correct Y-axis: shift entries from [0,2] to [-0.5,1.5] binning
            # Move bin contents from Y-bin 1 (center 0.5) to Y-bin that maps to 0,
            # and Y-bin 2 (center 1.5) to Y-bin that maps to 1.
            # Recreate with corrected Y-axis
            nx = h2.GetNbinsX()
            xmin = h2.GetXaxis().GetXmin()
            xmax = h2.GetXaxis().GetXmax()
            h2_corr = ROOT.TH2F(f"h2_corr_{hist_name_2d}_{i}", "",
                                 nx, xmin, xmax, 2, -0.5, 1.5)
            for bx in range(1, nx + 1):
                for by in range(1, 3):
                    h2_corr.SetBinContent(bx, by, h2.GetBinContent(bx, by))
                    h2_corr.SetBinError(bx, by, h2.GetBinError(bx, by))

            # Rebin X-axis before profiling
            h2_corr.RebinX(prof_config.get('rebin', 20))
            # Create ProfileX
            prof = h2_corr.ProfileX(f"prof_{hist_name_2d}_{i}")
            prof.SetDirectory(0)

            # Set style: color by MS
            params = file_params[i]
            color = minp_to_color.get(params.get('minp'), colors[i % len(colors)])
            prof.SetLineColor(color)
            prof.SetLineWidth(2)
            prof.SetMarkerColor(color)
            apply_minp_linestyle(prof, params.get('minp'))
            prof.SetMarkerSize(0.8)

            prof_histograms.append(prof)
            prof_labels.append(histogram_labels[i] if i < len(histogram_labels) else "")

            tfile.Close()

        if not prof_histograms:
            print(f"  Warning: No valid 2D histograms found for {hist_name_2d}")
            continue

        # Create canvas
        canvas_prof = ROOT.TCanvas(f"c_prof_{hist_name_2d}_{canvas_suffix}", f"c_prof_{hist_name_2d}_{canvas_suffix}", 800, 800)
        canvas_prof.cd()

        if prof_config['log_y']:
            canvas_prof.SetLogy()

        # Create legend
        legend_prof = ROOT.TLegend(0.17, 0.65, 0.9, 0.9)
        legend_prof.SetHeader(f"Region: {'pre-SR' if region.upper() == 'SR' else region.upper()}")
        legend_prof.SetFillStyle(0)
        legend_prof.SetBorderSize(0)
        legend_prof.SetTextSize(0.015)
        legend_prof.SetNColumns(3)

        # Find y-axis range
        max_val = 0
        min_val = 1e99
        for prof in prof_histograms:
            for b in range(1, prof.GetNbinsX() + 1):
                val = prof.GetBinContent(b)
                if val > max_val:
                    max_val = val
                if val > 0 and val < min_val:
                    min_val = val

        # Draw profiles
        for i, prof in enumerate(prof_histograms):
            if i == 0:
                prof.GetXaxis().SetTitle(prof_config['x_title'])
                prof.GetYaxis().SetTitle(prof_config['y_title'])
                prof.SetMaximum(2.0)
                prof.SetMinimum(0)
                draw_option = "E"
            else:
                draw_option = "E SAME"

            prof.Draw(draw_option)
            legend_prof.AddEntry(prof, prof_labels[i], "lp")

        legend_prof.Draw()
        addCMSText(canvas_prof, lumi_text="13.6 TeV", extra_text="Work in Progress")

        # Save
        png_name_prof = base_name.replace("_overlay", f"_profileX_{hist_name_2d}") + ".png"
        pdf_name_prof = base_name.replace("_overlay", f"_profileX_{hist_name_2d}") + ".pdf"

        canvas_prof.SaveAs(png_name_prof)
        canvas_prof.SaveAs(pdf_name_prof)
        canvas_prof.Update()

    # Cleanup
    for hist, tfile in histograms:
        tfile.Close()

if __name__ == "__main__":
    # Base directories and their labels
    base_dirs = {
        # "/ceph/cms/store/user/tvami/DiMuonPlusX/Ntuples_v1.0.1/Data/": "Data",
        # "/ceph/cms/store/user/tvami/DiMuonPlusX/Ntuples_v1.0.1/Signal/": "Signal",
        "/ceph/cms/store/user/tvami/DiMuonPlusX/Ntuples_v1.0.1/BkgMC/": "BkgMC",
    }

    # Define regions to loop over
    regions = ['sr', 'vr']

    # Loop over all combinations
    for base_dir, sample_label in base_dirs.items():
        # Extract Ntuple version from path (e.g., "v1.0.1" from "Ntuples_v1.0.1")
        version_match = re.search(r'Ntuples_(v[\d.]+)', base_dir)
        version = version_match.group(1) if version_match else "unknown"
        output_dir = f"figures/cutflow_plots_{version}"
        os.makedirs(output_dir, exist_ok=True)

        for region in regions:
            input_dir = os.path.join(base_dir, region)

            # Check if directory exists
            if not os.path.exists(input_dir):
                print(f"Warning: Directory {input_dir} does not exist, skipping...")
                continue

            print(f"\n{'='*100}")
            print(f"Processing: {sample_label}/{region}")
            print(f"{'='*100}\n")

            # Run the plotting function (creates all three types of plots)
            plot_all_cutflow_analysis(
                input_dir=input_dir,
                hist_name="h_cutflow",
                output_name=os.path.join(output_dir, f"{sample_label}_cutflow_overlay.png")
            )
