"""
Plot residual burden
"""


import pylab as pl
import sciris as sc
import numpy as np
import utils as ut

 
def plot_fig1():
    ut.set_font(20)
    fig = pl.figure(layout="tight", figsize=(14, 6))
    gs = fig.add_gridspec(1, 2)  # 1 row, 2 columns

    # Load Gabon scenario data
    msim_dict = sc.loadobj('results/scens_gabon.obj')

    # What to plot
    start_year = 2016
    end_year = 2100
    ymax = 25
    si = sc.findinds(msim_dict['Baseline'].year, start_year)[0]
    ei = sc.findinds(msim_dict['Baseline'].year, end_year)[0]
    fi = sc.findinds(msim_dict['Baseline'].year, 2025)[0]

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
    ax = ut.plot_single(ax, msim_dict['Baseline'], 'asr_cancer_incidence', si, ei,
                        color='k', label='Baseline')

    # Plot each combination
    for screen_idx, screen_level in enumerate(screening_levels):
        for vax in vax_status:
            scen_key = f'Screen {screen_level} + {vax}'
            ls = line_styles[vax]
            label = f'{screen_level} screening' if vax == 'No vaccination' else ''
            ax = ut.plot_single(ax, msim_dict[scen_key], 'asr_cancer_incidence', si, ei,
                               color=screening_colors[screen_idx], ls=ls, label=label)

    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('ASR cervical cancer incidence, 2025-2100\nScreening and prophylactic vaccination in Gabon')

    # Create legends
    # Screening level legend
    from matplotlib.patches import Patch
    screen_handles = [Patch(facecolor=screening_colors[i], label=f'{screening_levels[i]} screening')
                     for i in range(len(screening_levels))]
    legend1 = ax.legend(handles=screen_handles, title='Screening coverage',
                       loc='upper right', bbox_to_anchor=(1, 0.8), frameon=False)
    ax.add_artist(legend1)

    # Vaccination legend
    from matplotlib.lines import Line2D
    vax_handles = [Line2D([0], [0], color='k', linestyle='-', lw=2, label='No vaccination'),
                   Line2D([0], [0], color='k', linestyle='--', lw=2, label='90% vax coverage')]
    ax.legend(handles=vax_handles, title='Vaccination', loc='upper right', frameon=False)

    # Add panel label
    ax.text(-0.1, 1.05, 'A', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    ######################################################
    # Right Panel: Cumulative cancers
    ######################################################
    ax = fig.add_subplot(gs[1])

    # Set up grouped bars
    bar_width = 0.35
    x_base = np.arange(len(screening_levels))
    offsets = [-bar_width/2, bar_width/2]

    # Colors for vaccination status
    vax_colors = ['gray', 'lightblue']

    for vax_idx, vax in enumerate(vax_status):
        cum_cancers = []

        for screen_level in screening_levels:
            scen_key = f'Screen {screen_level} + {vax}'
            val = msim_dict[scen_key]['cancers'].values[fi:].sum()
            cum_cancers.append(val)
            print(f'{scen_key}: {val} cancers')

        bars = ax.bar(x_base + offsets[vax_idx], cum_cancers, width=bar_width,
                      color=vax_colors[vax_idx], label=vax)

        # Add value labels on bars
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height):,}',
                    ha='center', va='bottom', fontsize=10)

    ax.set_xticks(x_base)
    ax.set_xticklabels([f'{level}' for level in screening_levels])
    ax.set_xlabel('Screening coverage')
    ax.set_title('Cumulative cancers\n2025-2100')
    sc.SIticks()
    ax.legend(title='Vaccination', loc='upper right', frameon=False)

    # Add panel label
    ax.text(-0.1, 1.05, 'B', transform=ax.transAxes, fontsize=24, fontweight='bold', va='top')

    fig.tight_layout()
    fig_name = 'figures/gabon_vax_screening.png'
    sc.savefig(fig_name, dpi=100)

    return


# %% Run as a script
if __name__ == '__main__':

    location = 'gabon'
    plot_fig1()

    msim_dict = sc.loadobj('results/scens_gabon.obj')
    mbase = msim_dict['Screen 10% + 90% vax coverage']
    mno = msim_dict['Screen 10% + No vaccination']
    start_year = 2016
    end_year = 2100
    si = sc.findinds(mbase.year, start_year)[0]
    ei = sc.findinds(mbase.year, end_year)[0]
    fi = sc.findinds(mbase.year, 2025)[0]
 
    # elim_year = sc.findfirst(msim_dict['60% routine coverage']['asr_cancer_incidence'][si:]<4, die=False)+si

    # print(f'Elim year: {mbase.year[elim_year]}')
    print(f'Cancers in 2025: {mbase.cancers[si]} ({mbase.cancers.low[si]}, {mbase.cancers.high[si]})')
    print(f'Cancers in 2100: {mbase.cancers[ei]} ({mbase.cancers.low[ei]}, {mbase.cancers.high[ei]})')
    print(f'Cancers averted: {mno.cancers[fi:].sum()-mbase.cancers[fi:].sum()}')

