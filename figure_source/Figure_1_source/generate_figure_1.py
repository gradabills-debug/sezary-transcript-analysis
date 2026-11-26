"""Render Figure 1 from audited counts; no clustering or statistical inference."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter

HERE=Path(__file__).resolve().parent

def load_data():
    return pd.read_csv(HERE.parents[1]/"results"/"Figure_1_source_data.csv")

def configure_font(font_dir):
    if font_dir:
        for p in Path(font_dir).glob('*.ttf'):
            font_manager.fontManager.addfont(str(p))
    path=font_manager.findfont(font_manager.FontProperties(family='Arial'),fallback_to_default=False)
    assert font_manager.FontProperties(fname=path).get_name()=='Arial'
    plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.labelsize':9.5,'axes.titlesize':10.5,'xtick.labelsize':9,'ytick.labelsize':9.5,'axes.linewidth':.7,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white','figure.facecolor':'white','mathtext.fontset':'custom','mathtext.rm':'Arial','mathtext.it':'Arial:italic','mathtext.bf':'Arial:bold'})
    return path

def finish_axis(ax):
    ax.spines[['top','right','left']].set_visible(False)
    ax.tick_params(axis='y',length=0,pad=10)
    ax.tick_params(axis='x',length=3,width=.65,pad=4)
    ax.set_axisbelow(True)
    ax.grid(axis='x',color='#E4E4E4',linewidth=.5)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--font-dir',default=str(HERE.parent/'figure_fonts'))
    args=parser.parse_args()
    font_path=configure_font(args.font_dir)
    data=load_data()
    a=data[data.panel=='a'].copy();b=data[data.panel=='b'].copy()
    assert a.selected_cells.sum()==7929 and a.paired_cells.sum()==4657
    assert b[b.gex_accession!='GSM8200498'].clonotype_cells.sum()==2084
    assert b.clonotype_cells.tolist()==[773,1311,238,236]
    assert np.allclose(b.percentage,[100*773/938,100*1311/1741,100*238/1978,100*236/1978])

    fig=plt.figure(figsize=(7.2,6.2),dpi=180)
    fig.text(.045,.97,'Clonotype recovery and abundance in Sézary syndrome',fontsize=12,fontweight='bold',ha='left',va='top')
    fig.text(.045,.87,'a',fontsize=14,fontweight='bold',ha='left',va='center')
    fig.text(.235,.87,'Paired TCR recovery',fontsize=10.5,fontweight='bold',ha='left',va='center')
    ax=fig.add_axes([.235,.575,.61,.225])
    y=np.arange(len(a))
    ax.barh(y,a.paired_cells,height=.60,color='#1F4E79',edgecolor='#163A5C',linewidth=.6,zorder=3)
    ax.barh(y,a.unpaired_cells,left=a.paired_cells,height=.60,color='#E6E6E6',edgecolor='#999999',linewidth=.6,zorder=3)
    for i,(_,r) in enumerate(a.iterrows()):
        ax.text(r.paired_cells/2,i,f'{int(r.paired_cells):,} ({r.percentage:.1f}%)',color='white',ha='center',va='center',fontsize=8.5,zorder=4)
        ax.text(r.paired_cells+r.unpaired_cells/2,i,f'{int(r.unpaired_cells):,}',color='black',ha='center',va='center',fontsize=9,zorder=4)
        ax.text(1.105,i,f'{int(r.selected_cells):,}',transform=ax.get_yaxis_transform(),ha='center',va='center',fontsize=9.5,clip_on=False)
    ax.text(1.105,1.10,'Total cells',transform=ax.transAxes,ha='center',va='bottom',fontsize=8.5,clip_on=False)
    ax.set_yticks(y,a.label)
    ax.set_ylim(2.65,-.65)
    ax.set_xlim(0,3500)
    ax.set_xticks([0,1000,2000,3000])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{x:,.0f}'))
    ax.set_xlabel(r'$\mathit{CD3D}$-selected cells',labelpad=7)
    ax.axhline(1.5,color='#BDBDBD',linestyle=(0,(2,3)),linewidth=.6,zorder=1)
    finish_axis(ax)
    handles=[Patch(facecolor='#1F4E79',edgecolor='#163A5C',label='Paired TCR recovered'),Patch(facecolor='#E6E6E6',edgecolor='#999999',label='No paired TCR recovered')]
    fig.legend(handles=handles,loc='center left',bbox_to_anchor=(.235,.833),ncol=2,frameon=False,fontsize=8.3,handlelength=1.2,handleheight=.9,columnspacing=1.2,borderaxespad=0)

    fig.text(.045,.465,'b',fontsize=14,fontweight='bold',ha='left',va='center')
    fig.text(.235,.465,'Leading paired TCR clonotypes',fontsize=10.5,fontweight='bold',ha='left',va='center')
    bx=fig.add_axes([.235,.16,.61,.255])
    z=np.arange(len(b))
    bx.barh(z,b.percentage,height=.58,color='#1F4E79',edgecolor='#163A5C',linewidth=.6,zorder=3)
    for i,(_,r) in enumerate(b.iterrows()):
        bx.text(r.percentage+1.5,i,f'{r.percentage:.1f}%',ha='left',va='center',fontsize=9,zorder=4)
        bx.text(1.105,i,f'{int(r.numerator):,}/{int(r.denominator):,}',transform=bx.get_yaxis_transform(),ha='center',va='center',fontsize=9,clip_on=False)
    bx.text(1.105,1.08,'Clone / paired',transform=bx.transAxes,ha='center',va='bottom',fontsize=8.5,clip_on=False)
    bx.set_yticks(z,b.label)
    bx.set_ylim(3.65,-.65)
    bx.set_xlim(0,100)
    bx.set_xticks([0,20,40,60,80,100])
    bx.set_xlabel('Cells with paired TCR recovery (%)',labelpad=7)
    bx.axhline(1.5,color='#BDBDBD',linestyle=(0,(2,3)),linewidth=.6,zorder=1)
    finish_axis(bx)
    fig.text(.045,.055,'HLA08 reactions 1 and 2 are from one specimen; HLA19 is a second patient.',ha='left',va='center',fontsize=8.7)

    for ext in ['png','pdf','svg','tiff']:
        options={'dpi':600,'facecolor':'white'}
        if ext=='tiff':options['pil_kwargs']={'compression':'tiff_lzw'}
        fig.savefig(HERE/f'Figure_1.{ext}',**options)
    fig.savefig(HERE/'Figure_1_preview.png',dpi=180,facecolor='white')
    plt.close(fig)
    legend=('Figure 1. Clonotype recovery and abundance in Sézary syndrome. '
        '(a) Paired-TCR recovery among CD3D-selected cells in three libraries from two patients. Navy segments show cells with recovered productive, high-confidence TRA and TRB sequences; light gray segments show cells without a recovered pair. Segment labels give cell counts, with paired-recovery percentages in parentheses; right-hand labels show total selected cells. '
        '(b) Abundance of the leading HLA08 clonotype in each reaction and the two leading HLA19 clonotypes (pairs 1 and 2). Percentages use cells with paired-TCR recovery within the corresponding library as the denominator; right-hand labels give numerator/denominator. HLA08 reactions originated from the same specimen at one sampling time. Comparisons are descriptive.')
    (HERE/'Figure_1_legend.txt').write_text(legend+'\n')
    readme=('Figure 1 uses reconstructed processed data from GEO GSE263881.\n'
        'The source CSV supplies all plotted counts, denominators, and paired CDR3 sequence identities.\n'
        'The HLA08 paired sequence is shared between reactions. HLA19 pair 1 and pair 2 are different paired-TCR clonotypes; they are not assumed to be one malignant clone.\n'
        'All comparisons are descriptive. No statistical error bars or significance tests are presented.\n'
        'The PDF is vector artwork with embedded Arial subsets. The PNG and LZW-compressed TIFF are 600 dpi; the SVG preserves editable text and requires Arial for matching typography.\n'
        'The bars use navy (#1F4E79), with light gray (#E6E6E6) for cells without paired-TCR recovery.\n'
        'To reproduce, install numpy, pandas, and matplotlib and run generate_figure_1.py --font-dir /path/to/Arial/fonts. Arial font files are not included. The source CSV is used directly.\n')
    (HERE/'README.txt').write_text(readme)
    print(json.dumps({'font':font_path,'legend_words':len(legend.split()),'panel_a':a[['label','selected_cells','paired_cells','percentage']].to_dict('records'),'panel_b':b[['label','numerator','denominator','percentage']].to_dict('records'),'outputs':[{ 'file':p.name,'bytes':p.stat().st_size} for p in HERE.iterdir() if p.suffix in ['.png','.pdf','.svg','.tiff','.zip']]},indent=2))

if __name__=='__main__':main()
