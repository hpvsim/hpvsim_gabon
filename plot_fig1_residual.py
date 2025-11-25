"""
Plot residual burden
"""


import pylab as pl
import sciris as sc
import numpy as np
import utils as ut

 
def plot_fig1():
    ut.set_font(20)
    fig = pl.figure(layout="tight", figsize=(12, 5))

    gs = fig.add_gridspec(1, 3)
    ax = fig.add_subplot(gs[:2])

    resname = 'asr_cancer_incidence'

    # What to plot
    start_year = 2016
    end_year = 2100
    ymax = 30

    msim_dict = sc.loadobj(f'results/vx_scens_{location}.obj')
    mbase = msim_dict['90% routine coverage']
    si = sc.findinds(mbase.year, start_year)[0]
    ei = sc.findinds(mbase.year, end_year)[0]

    this_dict = {
        'No vaccination': msim_dict['No vaccination'],
        'Scale-up to 90% coverage': msim_dict['90% routine coverage'],
    }
    colors = sc.gridcolors(len(this_dict))
    for idx, (slabel, mres) in enumerate(this_dict.items()):
        # ls = ':' if slabel == 'No interventions' else '-'
        ax = ut.plot_single(ax, mres, resname, si, ei, color=colors[idx], label=slabel)
    ax.set_ylim(bottom=0, top=ymax)
    ax.set_title('ASR cervical cancer incidence, 2025-2100')
    ax.legend(loc="upper left", frameon=False, bbox_to_anchor=(0, 0.6))

    # Assemble cumulative results
    cum_res = dict()
    resname = 'cancers'
    fi = sc.findinds(mbase.year, 2025)[0]
    for sname, scen in this_dict.items():
        cum_res[sname] = scen[resname].values[fi:].sum()
        lb, ub = scen[resname].low[fi:].sum(), scen[resname].high[fi:].sum()
        print(f'{sname}: {cum_res[sname]} ({lb}, {ub}) cancers')
    # Make a bar plot
    ax = fig.add_subplot(gs[2])
    bars = list(cum_res.values())
    labels = list(cum_res.keys())
    x = np.arange(len(labels))
    ax.bar(x, bars, color='k')
    ax.set_xticks(x)
    ax.set_xticklabels(['No\nvax', 'Scale-up'])
    ax.set_title('Cumulative cancers')
    sc.SIticks()

    fig.tight_layout()
    fig_name = 'figures/fig1_residual.png'
    sc.savefig(fig_name, dpi=100)

    return


# %% Run as a script
if __name__ == '__main__':

    location = 'gabon'
    plot_fig1()

    msim_dict = sc.loadobj(f'results/vx_scens_{location}.obj')
    mbase = msim_dict['90% routine coverage']
    mno = msim_dict['No vaccination']
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

