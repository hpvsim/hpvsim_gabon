"""
Plot residual cervical cancer burden under combined screening and vaccination scenarios
"""


import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import sciris as sc
import numpy as np
import utils as ut


def plot_fig1(filestem=''):
    """ Plot ASR cancer incidence over time (panel A) and cumulative cancers 2025-2100 (panel B) """
    ut.set_font(20)
    fig = plt.figure(layout="tight", figsize=(20, 6))
    gs = fig.add_gridspec(1, 2)  # 1 row, 2 columns

    # Load Gabon scenario data
    msim_dict = sc.loadobj(f'results/scens_gabon{filestem}.obj')

    # What to plot
    start_year = 2016
    end_year = 2100
    vax_start_year = 2025  # Vaccination scenarios start in 2025; used as the cutoff for cumulative cancers averted
    ymax = 35  # y-axis max for ASR incidence panel, chosen to give headroom above the baseline curve
    start_idx = sc.findinds(msim_dict['Baseline'].year, start_year)[0]
    end_idx = sc.findinds(msim_dict['Baseline'].year, end_year)[0]
    vax_start_idx = sc.findinds(msim_dict['Baseline'].year, vax_start_year)[0]

    # Define screening levels and vaccination status
    screening_levels = ['10%', '40%', '90%']
    vax_status = ['No vaccination', '90% vax coverage']

    # Create color palette for screening levels
    screening_colors = sc.vectocolor(len(screening_levels)).tolist()

    # Line styles for vaccination status
    line_styles = {
        'No vaccination': '-',
        '90% vax coverage': '--'
    }

    ######################################################
    # Left Panel: Time series of ASR cancer incidence
    ######################################################
    ax = fig.add_subplot(gs[0])

    # Plot baseline
    ax = ut.plot_single(ax, msim_dict['Baseline'], 'asr_cancer_incidence', start_idx, end_idx,
                        color='k', label='Baseline')

    # Plot each combination
    for screen_idx, screen_level in enumerate(screening_levels):
        for vax in vax_status:
            scen_key = f'Screen {screen_level} + {vax}'
            ls = line_styles[vax]
            label = f'{screen_level} screening' if vax == 'No vaccination' else ''
            ax = ut.plot_single(ax, msim_dict[scen_key], 'asr_cancer_incidence', start_idx, end_idx,
                               color=screening_colors[screen_idx], ls=ls, label=label)

    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('ASR cervical cancer incidence, 2025-2100\nScreening and prophylactic vaccination in Gabon')

    # Create legends
    # Screening level legend
    screen_handles = [Patch(facecolor=screening_colors[i], label=f'{screening_levels[i]} screening')
                     for i in range(len(screening_levels))]
    legend1 = ax.legend(handles=screen_handles, title='',
                       loc='lower left', frameon=False)
    ax.add_artist(legend1)

    # Vaccination legend
    vax_handles = [Line2D([0], [0], color='k', linestyle='-', lw=2, label='No vaccination'),
                   Line2D([0], [0], color='k', linestyle='--', lw=2, label='90% vax coverage')]
    ax.legend(handles=vax_handles, title='', loc='lower left', bbox_to_anchor=(0.3, 0), frameon=False)

    # Add panel label
    ax.text(-0.1, 1.05, 'A', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    ######################################################
    # Right Panel: Cumulative cancers
    ######################################################
    ax = fig.add_subplot(gs[1])

    # Set up grouped bars; offsets generalize to however many vax_status entries there are
    # (for n_vax=2 this reproduces the original +/-bar_width/2 spacing exactly)
    bar_width = 0.35
    x_base = np.arange(len(screening_levels))
    n_vax = len(vax_status)
    offsets = [(i - (n_vax-1)/2) * bar_width for i in range(n_vax)]

    # Colors for vaccination status (extend this list if more vax_status entries are added)
    vax_colors = ['gray', 'lightblue']
    assert n_vax <= len(vax_colors), f'Need at least {n_vax} vax_colors, only have {len(vax_colors)}'

    for vax_idx, vax in enumerate(vax_status):
        cum_cancers, cum_low, cum_high = [], [], []

        for screen_level in screening_levels:
            scen_key = f'Screen {screen_level} + {vax}'
            res = msim_dict[scen_key]['cancers']
            val = res.values[vax_start_idx:].sum()
            cum_cancers.append(val)
            cum_low.append(res.low[vax_start_idx:].sum())
            cum_high.append(res.high[vax_start_idx:].sum())
            print(f'{scen_key}: {val:.0f} cancers ({cum_low[-1]:.0f}, {cum_high[-1]:.0f})')

        yerr = np.array([np.array(cum_cancers) - np.array(cum_low),
                         np.array(cum_high) - np.array(cum_cancers)])
        bars = ax.bar(x_base + offsets[vax_idx], cum_cancers, width=bar_width,
                      color=vax_colors[vax_idx], label=vax,
                      yerr=yerr, capsize=5, error_kw={'ecolor': 'black', 'lw': 1.2})

        # Add value labels above the error bar top
        for bar, hi in zip(bars, cum_high):
            ax.text(bar.get_x() + bar.get_width()/2., hi,
                    f'{int(bar.get_height()):,}',
                    ha='center', va='bottom', fontsize=16)

    ax.set_xticks(x_base)
    ax.set_xticklabels([f'{level}' for level in screening_levels])
    ax.set_xlabel('Screening coverage')
    ax.set_title('Cumulative cancers\n2025-2100')
    sc.SIticks()
    ax.set_ylim([0, 70e3])
    ax.legend(title='', loc='upper right', frameon=False)

    # Add panel label
    ax.text(-0.1, 1.05, 'B', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    fig.tight_layout()
    fig_name = f'figures/gabon_vax_screening{filestem}.png'
    sc.savefig(fig_name, dpi=100)

    return msim_dict


# %% Run as a script
if __name__ == '__main__':

    msim_dict = plot_fig1()

    mbase = msim_dict['Screen 10% + 90% vax coverage']
    mno = msim_dict['Screen 10% + No vaccination']
    start_year = 2016
    end_year = 2100
    start_idx = sc.findinds(mbase.year, start_year)[0]
    end_idx = sc.findinds(mbase.year, end_year)[0]
    vax_start_idx = sc.findinds(mbase.year, 2025)[0]

    print(f'Cancers in 2025: {mbase.cancers[start_idx]} ({mbase.cancers.low[start_idx]}, {mbase.cancers.high[start_idx]})')
    print(f'Cancers in 2100: {mbase.cancers[end_idx]} ({mbase.cancers.low[end_idx]}, {mbase.cancers.high[end_idx]})')
    print(f'Cancers averted: {mno.cancers[vax_start_idx:].sum()-mbase.cancers[vax_start_idx:].sum()}')

