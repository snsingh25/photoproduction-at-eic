#!/usr/bin/env python3
"""
Differential Jet Shapes Analysis - Python Version
Creates differential jet shape plots (rho(r)) for QQ and GG jets
Replicates functionality of plot_diffshapes.C
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, MultipleLocator
import numpy as np
from scipy.interpolate import make_interp_spline

def plot_diffshapes():
    """
    Create differential jet shape plots for QQ and GG jets
    across different eta bins and center-of-mass energies.
    """
    
    # --- Plotting Style Configuration ---
    plt.style.use('default')
    plt.rcParams.update({
        'font.size': 24,
        'font.family': 'serif',
        'font.serif': ['Computer Modern'], # Matches ROOT's Times/Helvetica look when Tex is on
        'text.usetex': True,
        'axes.linewidth': 1.2,
        'axes.spines.top': False,   # Clean look
        'axes.spines.right': False, # Clean look
        'xtick.major.size': 6,
        'xtick.minor.size': 3,
        'ytick.major.size': 6,
        'ytick.minor.size': 3,
        'lines.linewidth': 2,
        'legend.frameon': False,
        'figure.facecolor': 'white',
        'axes.facecolor': 'white'
    })
    
    # Data for r values (x-axis)
    r_values = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    
    def create_plot(qq_64, qq_105, qq_141, qq_300, 
                   gg_64, gg_105, gg_141, gg_300, 
                   eta_label, file_name):
        """
        Create and save individual plot for each eta bin
        """
        
        # Create figure and axis
        fig, ax = plt.subplots(figsize=(8, 8))
        
        # Define colors
        qq_color = '#d62728' # Standard MPL Red (close to kRed+1)
        gg_color = '#1f77b4' # Standard MPL Blue (close to kBlue+1)
        
        # --- Data Smoothing (to mimic ROOT "C" draw option) ---
        # We create a dense x-axis (300 points) and interpolate the y-values
        r_smooth = np.linspace(r_values.min(), r_values.max(), 300)
        
        # Helper to smooth data
        def get_smooth(y_data):
            return make_interp_spline(r_values, y_data, k=3)(r_smooth)

        # Plot QQ jets (Red)
        ax.plot(r_smooth, get_smooth(qq_64), color=qq_color, linestyle=':', linewidth=2, label='ep 64 GeV')
        ax.plot(r_smooth, get_smooth(qq_105), color=qq_color, linestyle='--', linewidth=2, label='ep 105 GeV')
        ax.plot(r_smooth, get_smooth(qq_141), color=qq_color, linestyle='-.', linewidth=2, label='ep 141 GeV')
        ax.plot(r_smooth, get_smooth(qq_300), color=qq_color, linestyle='-', linewidth=2, label='ep 300 GeV')
        
        # Plot GG jets (Blue)
        ax.plot(r_smooth, get_smooth(gg_64), color=gg_color, linestyle=':', linewidth=2)
        ax.plot(r_smooth, get_smooth(gg_105), color=gg_color, linestyle='--', linewidth=2)
        ax.plot(r_smooth, get_smooth(gg_141), color=gg_color, linestyle='-.', linewidth=2)
        ax.plot(r_smooth, get_smooth(gg_300), color=gg_color, linestyle='-', linewidth=2)
        
        # --- Axis Properties ---
        ax.set_xlim(0.05, 1.0)
        ax.set_ylim(0.0, 4.0) # Linear scale 0 to 4
        
        ax.set_xlabel(r'$r$', fontsize=28)
        ax.set_ylabel(r'$\rho(r)$', fontsize=28, labelpad=10)
        
        # Center axis labels (approximate centering logic)
        ax.xaxis.set_label_coords(0.5, -0.08)
        ax.yaxis.set_label_coords(-0.08, 0.5)
        
        # Set tick parameters
        ax.tick_params(axis='both', which='major', labelsize=20)
        ax.xaxis.set_minor_locator(MultipleLocator(0.05)) # Minor ticks every 0.05
        ax.yaxis.set_minor_locator(MultipleLocator(0.2))  # Minor ticks every 0.2
        ax.minorticks_on()
        
        # --- Legend Construction ---
        # We create black dummy lines for the legend so it only shows line styles
        legend_elements = [
            plt.Line2D([0], [0], color='black', linestyle=':', linewidth=2, label='64 GeV'),
            plt.Line2D([0], [0], color='black', linestyle='--', linewidth=2, label='105 GeV'),
            plt.Line2D([0], [0], color='black', linestyle='-.', linewidth=2, label='141 GeV'),
            plt.Line2D([0], [0], color='black', linestyle='-', linewidth=2, label='300 GeV')
        ]
        
        # Position legend (Top Right in C++ code, but usually differential plots are empty in top right)
        ax.legend(handles=legend_elements, loc='upper right', frameon=False, fontsize=22)
        
        # --- Text Annotations ---
        # Add "QQ / GG" colored text
        # Positioned relative to axes (0,0 is bottom-left, 1,1 is top-right)
        ax.text(0.60, 0.60, 'QQ', transform=ax.transAxes, fontsize=24, color=qq_color, weight='bold')
        ax.text(0.70, 0.60, ' / ', transform=ax.transAxes, fontsize=24, color='black', weight='bold')  
        ax.text(0.75, 0.60, 'GG', transform=ax.transAxes, fontsize=24, color=gg_color, weight='bold')
        
        # Add Eta Label
        # Processing string to ensure proper LaTeX rendering if raw string is passed
        if "eta" in eta_label and "\\" not in eta_label:
             # Basic fix if user passes simple text, but our inputs below are raw strings
             pass
        ax.text(0.15, 0.85, eta_label, transform=ax.transAxes, fontsize=24)
        
        # Adjust layout and save
        plt.tight_layout()
        plt.savefig(file_name, dpi=600, bbox_inches='tight', facecolor='white')
        print(f"Plot created for {eta_label} -> {file_name}")
        plt.close()

    # --- Data Definitions (Transcribed from C++) ---
    
    # ETA BIN 1: (-1 < eta < 0)
    qq_64_eta1 = np.array([2.7067, 2.1164, 1.3982, 0.8514, 0.5688, 0.396, 0.3346, 0.2481, 0.1452, 0.0181])
    qq_105_eta1 = np.array([3.2144, 2.4553, 1.5482, 0.989, 0.6732, 0.478, 0.3704, 0.2872, 0.1802, 0.0227])
    qq_141_eta1 = np.array([3.1969, 2.3892, 1.4322, 0.8947, 0.6225, 0.4458, 0.3458, 0.2681, 0.1735, 0.0244])
    qq_300_eta1 = np.array([3.1127, 2.2566, 1.3594, 0.871, 0.6213, 0.4571, 0.3672, 0.3044, 0.2078, 0.0335])

    gg_64_eta1 = np.array([1.2393, 1.7212, 1.4479, 1.1181, 0.9598, 0.9094, 0.5361, 0.6709, 0.2544, 0.0162])
    gg_105_eta1 = np.array([1.1629, 1.7121, 1.5188, 1.2898, 0.9642, 0.7612, 0.7216, 0.5244, 0.2693, 0.0227])
    gg_141_eta1 = np.array([1.2213, 1.5389, 1.477, 1.2965, 1.0624, 0.8243, 0.6386, 0.5299, 0.3096, 0.0457])
    gg_300_eta1 = np.array([1.4184, 1.7653, 1.5791, 1.299, 1.0366, 0.8584, 0.7193, 0.5869, 0.3915, 0.064])

    # ETA BIN 2: (0 < eta < 1)
    qq_64_eta2 = np.array([3.2144, 2.4553, 1.5482, 0.989, 0.6732, 0.478, 0.3704, 0.2872, 0.1802, 0.0227])
    qq_105_eta2 = np.array([3.0471, 2.4194, 1.4517, 0.9281, 0.6276, 0.4586, 0.3492, 0.2701, 0.1728, 0.0229])
    qq_141_eta2 = np.array([3.0451, 2.3809, 1.4758, 0.9425, 0.673, 0.4836, 0.3814, 0.3, 0.1988, 0.0283])
    qq_300_eta2 = np.array([3.0955, 2.3324, 1.4234, 0.8844, 0.6038, 0.447, 0.3407, 0.2733, 0.1797, 0.0235])

    gg_64_eta2 = np.array([1.4013, 1.8083, 1.6386, 1.2993, 1.0333, 0.822, 0.6873, 0.5343, 0.3137, 0.0494])
    gg_105_eta2 = np.array([1.3564, 1.8137, 1.5139, 1.3047, 1.0164, 0.8221, 0.6815, 0.5232, 0.3261, 0.0516])
    gg_141_eta2 = np.array([1.2158, 1.5969, 1.5567, 1.3275, 1.0829, 0.9092, 0.7596, 0.6151, 0.3974, 0.0674])
    gg_300_eta2 = np.array([1.4663, 1.7894, 1.6186, 1.3163, 1.057, 0.8749, 0.7399, 0.6255, 0.4401, 0.0749])

    # ETA BIN 3: (1 < eta < 1.5)
    qq_64_eta3 = np.array([3.1969, 2.3892, 1.4322, 0.8947, 0.6225, 0.4458, 0.3458, 0.2681, 0.1735, 0.0244])
    qq_105_eta3 = np.array([3.1965, 2.3338, 1.4193, 0.8811, 0.6035, 0.4455, 0.348, 0.2743, 0.179, 0.026])
    qq_141_eta3 = np.array([3.074, 2.2656, 1.4013, 0.9164, 0.6398, 0.4709, 0.3716, 0.3072, 0.205, 0.0308])
    qq_300_eta3 = np.array([3.0471, 2.4194, 1.3875, 0.9281, 0.6276, 0.4586, 0.3492, 0.2701, 0.1728, 0.0229])

    gg_64_eta3 = np.array([1.5157, 1.9342, 1.7226, 1.3443, 1.1092, 0.8719, 0.7479, 0.5836, 0.3877, 0.0596])
    gg_105_eta3 = np.array([1.4394, 1.8063, 1.6565, 1.3035, 1.0745, 0.8891, 0.7483, 0.5888, 0.4094, 0.0607])
    gg_141_eta3 = np.array([1.3305, 1.67, 1.5658, 1.3375, 1.1313, 0.9602, 0.7895, 0.6468, 0.4453, 0.0725])
    gg_300_eta3 = np.array([1.499, 1.8235, 1.6476, 1.3386, 1.0851, 0.9005, 0.7639, 0.6471, 0.4643, 0.0799])

    # ETA BIN 4: (1.5 < eta < 2)
    qq_64_eta4 = np.array([3.0729, 2.2427, 1.3847, 0.8767, 0.6124, 0.4351, 0.328, 0.2706, 0.1829, 0.0241])
    qq_105_eta4 = np.array([3.2012, 2.3082, 1.4196, 0.9112, 0.6204, 0.46, 0.3616, 0.2799, 0.1963, 0.0282])
    qq_141_eta4 = np.array([2.149, 2.0844, 1.5643, 1.1801, 0.8991, 0.7025, 0.5879, 0.4761, 0.3349, 0.0515])
    qq_300_eta4 = np.array([3.1127, 2.2566, 1.3594, 0.871, 0.6213, 0.4571, 0.357, 0.3044, 0.2078, 0.0335])

    gg_64_eta4 = np.array([1.4914, 1.7605, 1.5187, 1.4257, 1.1041, 0.8632, 0.7323, 0.5641, 0.3963, 0.0613])
    gg_105_eta4 = np.array([1.4394, 1.8063, 1.6565, 1.3035, 1.0745, 0.8891, 0.7483, 0.5888, 0.4094, 0.0607])
    gg_141_eta4 = np.array([1.3708, 1.7405, 1.5921, 1.3568, 1.168, 0.9773, 0.8288, 0.676, 0.4751, 0.0811])
    gg_300_eta4 = np.array([1.4022, 1.7214, 1.5166, 1.225, 0.9666, 0.775, 0.6535, 0.5014, 0.328, 0.0502])
    
    # --- Generate Plots ---
    create_plot(qq_64_eta1, qq_105_eta1, qq_141_eta1, qq_300_eta1,
                gg_64_eta1, gg_105_eta1, gg_141_eta1, gg_300_eta1,
                r"$-1 < \eta < 0$", "diffjetshapesminus1to0.pdf")
    
    create_plot(qq_64_eta2, qq_105_eta2, qq_141_eta2, qq_300_eta2,
                gg_64_eta2, gg_105_eta2, gg_141_eta2, gg_300_eta2,
                r"$0 < \eta < 1$", "diffjetshapes0to1.pdf")
    
    create_plot(qq_64_eta3, qq_105_eta3, qq_141_eta3, qq_300_eta3,
                gg_64_eta3, gg_105_eta3, gg_141_eta3, gg_300_eta3,
                r"$1 < \eta < 1.5$", "diffjetshapes1to1p5.pdf")
    
    create_plot(qq_64_eta4, qq_105_eta4, qq_141_eta4, qq_300_eta4,
                gg_64_eta4, gg_105_eta4, gg_141_eta4, gg_300_eta4,
                r"$1.5 < \eta < 2$", "diffjetshapes1p5tominus2.pdf")
    
    print("\nAll differential plots created successfully!")

if __name__ == "__main__":
    plot_diffshapes()