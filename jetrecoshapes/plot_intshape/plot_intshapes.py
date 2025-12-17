#!/usr/bin/env python3
"""
Integrated Jet Shapes Analysis - Python Version
Creates integrated jet shape plots (Psi(r)) for QQ and GG jets
Replicates functionality with updated Tufte/ROOT-like styling
"""

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, MultipleLocator, LogLocator
import numpy as np
from scipy.interpolate import make_interp_spline

def plot_intshapes():
    """
    Create integrated jet shape plots for QQ and GG jets
    across different eta bins and center-of-mass energies
    """
    
    # --- Plotting Style Configuration ---
    plt.style.use('default')
    plt.rcParams.update({
        'font.size': 24,
        'font.family': 'serif',
        'font.serif': ['Computer Modern'],
        'text.usetex': True,
        'axes.linewidth': 1.2,
        'axes.spines.top': False,
        'axes.spines.right': False,
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
        
        # Define colors (Matched to Differential Script)
        qq_color = '#d62728' # Standard MPL Red
        gg_color = '#1f77b4' # Standard MPL Blue
        
        # --- Data Smoothing ---
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
        
        # Log Scale logic for Integrated Shapes
        ax.set_yscale('log')
        ax.set_yticks([0.05, 0.1, 1.0])
        ax.set_yticklabels(['0.05', '0.1', '1.0'])
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.yaxis.get_major_formatter().set_scientific(False)
        
        ax.set_xlabel(r'$r$', fontsize=28)
        ax.set_ylabel(r'$\Psi(r)$', fontsize=28, labelpad=10)
        
        # Center axis labels
        ax.xaxis.set_label_coords(0.5, -0.08)
        ax.yaxis.set_label_coords(-0.08, 0.5)
        
        # Set tick parameters
        ax.tick_params(axis='both', which='major', labelsize=20)
        ax.tick_params(axis='both', which='minor', labelsize=10)
        ax.minorticks_on()
        
        # --- Legend Construction ---
        legend_elements = [
            plt.Line2D([0], [0], color='black', linestyle=':', linewidth=2, label='64 GeV'),
            plt.Line2D([0], [0], color='black', linestyle='--', linewidth=2, label='105 GeV'),
            plt.Line2D([0], [0], color='black', linestyle='-.', linewidth=2, label='141 GeV'),
            plt.Line2D([0], [0], color='black', linestyle='-', linewidth=2, label='300 GeV')
        ]
        
        # Position legend (Lower Right is best for Integrated shapes as curve goes Top-Left)
        legend = ax.legend(handles=legend_elements, loc='lower right', bbox_to_anchor=(0.95, 0.15), fontsize=22)
        
        # --- Text Annotations ---
        # Add "QQ / GG" colored text
        ax.text(0.60, 0.50, 'QQ', transform=ax.transAxes, fontsize=24, color=qq_color, weight='bold')
        ax.text(0.70, 0.50, ' / ', transform=ax.transAxes, fontsize=24, color='black', weight='bold')  
        ax.text(0.75, 0.50, 'GG', transform=ax.transAxes, fontsize=24, color=gg_color, weight='bold')
        
        # Add Eta Label
        ax.text(0.65, 0.1, eta_label, transform=ax.transAxes, fontsize=24)        
        
        # Adjust layout
        plt.tight_layout()
        
        # Save plot
        plt.savefig(file_name, dpi=600, bbox_inches='tight', facecolor='white')
        print(f"Plot created for {eta_label} -> {file_name}")
        plt.close()
    
    # --- Data Definitions ---
    
    # ETA BIN 1: (-1 < eta < 0)
    qq_64_eta1 = np.array([0.1904, 0.4473, 0.6174, 0.7296, 0.7957, 0.8445, 0.8794, 0.9089, 0.9301, 0.9373])
    qq_105_eta1 = np.array([0.1915, 0.4444, 0.6182, 0.7205, 0.7888, 0.8354, 0.8696, 0.8989, 0.9201, 0.9311])
    qq_141_eta1 = np.array([0.2025, 0.4561, 0.6217, 0.7234, 0.7902, 0.8338, 0.8674, 0.9025, 0.9245, 0.9377])
    qq_300_eta1 = np.array([0.2182, 0.4756, 0.6244, 0.7405, 0.8175, 0.8755, 0.8957, 0.9478, 0.9721, 0.9854])

    gg_64_eta1 = np.array([0.0754, 0.2284, 0.3793, 0.5217, 0.6335, 0.7193, 0.7915, 0.8542, 0.9024, 0.9131])
    gg_105_eta1 = np.array([0.0587, 0.2226, 0.3755, 0.5146, 0.6219, 0.7008, 0.7678, 0.8325, 0.8787, 0.8948])
    gg_141_eta1 = np.array([0.0729, 0.2162, 0.3626, 0.4982, 0.6073, 0.6984, 0.7669, 0.8243, 0.8676, 0.8922])
    gg_300_eta1 = np.array([0.0837, 0.2442, 0.4, 0.531, 0.6323, 0.7137, 0.781, 0.835, 0.8778, 0.9023])

    # ETA BIN 2: (0 < eta < 1)
    qq_64_eta2 = np.array([0.18, 0.4307, 0.6062, 0.7189, 0.7929, 0.8448, 0.8818, 0.9105, 0.932, 0.9426])
    qq_105_eta2 = np.array([0.18, 0.435, 0.613, 0.7246, 0.7968, 0.8482, 0.8852, 0.9136, 0.9344, 0.9449])
    qq_141_eta2 = np.array([0.1751, 0.4236, 0.5973, 0.7074, 0.7811, 0.8334, 0.8736, 0.9046, 0.9273, 0.9394])
    qq_300_eta2 = np.array([0.1913, 0.4452, 0.6171, 0.7244, 0.7942, 0.8432, 0.8802, 0.9094, 0.9319, 0.9436])

    gg_64_eta2 = np.array([0.0632, 0.2148, 0.379, 0.5192, 0.6308, 0.7222, 0.7941, 0.8523, 0.8936, 0.9174])
    gg_105_eta2 = np.array([0.068, 0.2229, 0.3782, 0.5148, 0.6282, 0.7205, 0.7947, 0.8545, 0.8976, 0.9197])
    gg_141_eta2 = np.array([0.0567, 0.1891, 0.3363, 0.4707, 0.5833, 0.6784, 0.759, 0.8242, 0.8725, 0.9004])
    gg_300_eta2 = np.array([0.0711, 0.2179, 0.371, 0.5056, 0.6141, 0.7026, 0.7756, 0.8362, 0.8825, 0.9085])

    # ETA BIN 3: (1 < eta < 1.5)
    qq_64_eta3 = np.array([0.1926, 0.4548, 0.6367, 0.7522, 0.8264, 0.8781, 0.9157, 0.9451, 0.9653, 0.9762])
    qq_105_eta3 = np.array([0.1848, 0.4474, 0.6304, 0.7453, 0.821, 0.8736, 0.912, 0.9421, 0.9628, 0.9737])
    qq_141_eta3 = np.array([0.1786, 0.4313, 0.6121, 0.7382, 0.8063, 0.8623, 0.9044, 0.9368, 0.9602, 0.9729])
    qq_300_eta3 = np.array([0.1868, 0.4456, 0.6251, 0.7433, 0.8214, 0.8759, 0.9177, 0.95, 0.9737, 0.9869])

    gg_64_eta3 = np.array([0.0791, 0.2335, 0.3923, 0.5387, 0.6585, 0.7504, 0.8277, 0.8875, 0.9367, 0.9628])
    gg_105_eta3 = np.array([0.0682, 0.215, 0.3741, 0.5204, 0.6393, 0.7363, 0.8177, 0.8823, 0.9303, 0.9568])
    gg_141_eta3 = np.array([0.0613, 0.1989, 0.3496, 0.49, 0.613, 0.7145, 0.8002, 0.8702, 0.9209, 0.9512])
    gg_300_eta3 = np.array([0.0683, 0.2177, 0.3781, 0.5204, 0.6375, 0.7334, 0.8129, 0.8772, 0.9262, 0.9557])

    # ETA BIN 4: (1.5 < eta < 2)
    qq_64_eta4 = np.array([0.1937, 0.4533, 0.6353, 0.753, 0.8293, 0.883, 0.9211, 0.9501, 0.9717, 0.983])
    qq_105_eta4 = np.array([0.188, 0.4526, 0.6334, 0.7503, 0.8287, 0.8823, 0.9219, 0.9531, 0.9749, 0.9869])
    qq_141_eta4 = np.array([0.188, 0.4445, 0.6244, 0.7405, 0.8018, 0.8473, 0.8882, 0.9092, 0.9309, 0.9438])
    qq_300_eta4 = np.array([0.186, 0.4447, 0.6258, 0.7419, 0.8163, 0.8689, 0.909, 0.9401, 0.9631, 0.9751])

    gg_64_eta4 = np.array([0.0694, 0.2216, 0.3861, 0.5318, 0.6527, 0.7591, 0.8361, 0.8987, 0.9392, 0.9704])
    gg_105_eta4 = np.array([0.068, 0.215, 0.3741, 0.5204, 0.6393, 0.7363, 0.8177, 0.8823, 0.9303, 0.9568])
    gg_141_eta4 = np.array([0.0653, 0.2073, 0.3604, 0.5052, 0.6239, 0.7285, 0.8161, 0.8871, 0.9394, 0.9726])
    gg_300_eta4 = np.array([0.0694, 0.2192, 0.3823, 0.5268, 0.6474, 0.7405, 0.8292, 0.8952, 0.9444, 0.9753])
    
    # Create plots for all eta bins
    create_plot(qq_64_eta1, qq_105_eta1, qq_141_eta1, qq_300_eta1,
                gg_64_eta1, gg_105_eta1, gg_141_eta1, gg_300_eta1,
                r"$-1 < \eta < 0$", "injetshapesminus1to0.pdf")
    
    create_plot(qq_64_eta2, qq_105_eta2, qq_141_eta2, qq_300_eta2,
                gg_64_eta2, gg_105_eta2, gg_141_eta2, gg_300_eta2,
                r"$0 < \eta < 1$", "injetshapes0to1.pdf")
    
    create_plot(qq_64_eta3, qq_105_eta3, qq_141_eta3, qq_300_eta3,
                gg_64_eta3, gg_105_eta3, gg_141_eta3, gg_300_eta3,
                r"$1 < \eta < 1.5$", "injetshapes1to1p5.pdf")
    
    create_plot(qq_64_eta4, qq_105_eta4, qq_141_eta4, qq_300_eta4,
                gg_64_eta4, gg_105_eta4, gg_141_eta4, gg_300_eta4,
                r"$1.5 < \eta < 2$", "injetshapes1p5to2.pdf")
    
    print("\nAll integrated plots created successfully!")

if __name__ == "__main__":
    plot_intshapes()