#===============================================================================
# fgwasplot.py
#===============================================================================

"""Plotting function for fgwas .params files"""



# Imports ======================================================================

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns




# Functions ====================================================================

def plot_params(params_path, pdf_path, title='', lower_bound=-10,
                upper_bound=10):
    df = pd.read_csv(params_path, sep=' ', header=None, skiprows=1)
    df.columns = ['Annotation', 'CI_95L', 'Estimate', 'CI_95H']
    df['CI_95L'] = df['CI_95L'].astype(str).str.replace('<','').astype(float)
    df['CI_95H'] = (df['CI_95H'].astype(str).str.replace('>', '')
                    .str.replace('fail', str(upper_bound)).astype(float))
    df['Annotation'] = df['Annotation'].replace(regex='_ln',value='')

    plt.style.use('seaborn-notebook')

    y_count = range(1,len(df.index) + 1)
    studies = list(reversed(df['Annotation'].tolist()))
    estimates = list(reversed(df['Estimate'].tolist()))
    estimates_fixed = [lower_bound if est < lower_bound else upper_bound if
                       est > upper_bound else est for est in estimates]
    lower_ci = list(reversed(df['CI_95L'].tolist()))
    upper_ci = list(reversed(df['CI_95H'].tolist()))
    colors = ['red' if l_ci > 0 and u_ci < upper_bound else 'lightgrey' for
              l_ci, u_ci in zip(lower_ci, upper_ci)]
            
    fig, ax1 = plt.subplots(1, 1, figsize=(3, len(df) / 2))
    ax1.scatter(estimates_fixed, y_count, c=colors)
    ax1.hlines(y_count, lower_ci, upper_ci, colors=colors, linewidth=2.2)

    buffer_zone = 0.5
    for i in range(len(y_count)):
        if lower_ci[i] < lower_bound and upper_ci[i] < upper_bound:
            ax1.annotate(s='', xy=(lower_bound - buffer_zone, y_count[i]), 
                         xytext=(upper_ci[i], y_count[i]),
                         arrowprops=dict(width=2, headwidth=7, headlength=7,
                                         facecolor='lightgrey', 
                                         edgecolor='lightgrey'))
        elif lower_ci[i] > lower_bound and upper_ci[i] > upper_bound:
            ax1.annotate(s='', xy=(upper_bound + buffer_zone, y_count[i]), 
                         xytext=(lower_ci[i], y_count[i]),
                         arrowprops=dict(width=2, headwidth=7, headlength=7,
                                         facecolor='lightgrey', 
                                         edgecolor='lightgrey'))
        elif lower_ci[i] < lower_bound and upper_ci[i] > upper_bound:
            ax1.annotate(s='', xy=(upper_bound + buffer_zone, y_count[i]), 
                         xytext=(0, y_count[i]),
                         arrowprops=dict(width=2, headwidth=7, headlength=7,
                                         facecolor='lightgrey', 
                                         edgecolor='lightgrey'))        
            ax1.annotate(s='', xy=(lower_bound - buffer_zone, y_count[i]), 
                         xytext=(0, y_count[i]),
                         arrowprops=dict(width=2, headwidth=7, headlength=7,
                                         facecolor='lightgrey', 
                                         edgecolor='lightgrey'))
        
    ax1.set_xlim(lower_bound - 0.5, upper_bound + buffer_zone)
    ax1.xaxis.set_tick_params(labelsize=14, length=5, width=1.5)
    ax1.yaxis.set_ticks(y_count)
    ax1.yaxis.set_tick_params(labelsize=14, length=5, width=1.5)
    ax1.set_yticklabels(studies)
    for axis in ['bottom','left']:
        ax1.spines[axis].set_linewidth(1.5)

    ax1.vlines(0, min(y_count) - 0.5, max(y_count) + 0.5, color='black',
               linestyle='dotted', lw=1.5)
    ax1.set_title(title, fontsize=20)
    ax1.set_xlabel('ln(enrichment)', fontsize=18)
    sns.despine(ax=ax1, offset=10, trim=True)
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight')




# Exceptions ===================================================================

class Error(Exception):
   '''Base class for other exceptions'''
   pass


class BadArgumentError(Error):
    '''Bad argument error'''
    pass
