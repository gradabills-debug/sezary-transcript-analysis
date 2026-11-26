"""Figure 2 from audited cell-level transcript counts; descriptive only."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

HERE=Path(__file__).resolve().parent
GENES=['PTGDR2','HPGDS','CRLF2','KLRB1']
JOINT_GENES=['PTGDR2','HPGDS','CRLF2']

def load_data():
    return pd.read_csv(HERE.parents[1]/"results"/"Figure_2_source_data.csv")

def configure_font(font_dir):
    if font_dir:
        for file in Path(font_dir).glob('*.ttf'):font_manager.fontManager.addfont(str(file))
    path=font_manager.findfont(font_manager.FontProperties(family='Arial'),fallback_to_default=False)
    assert font_manager.FontProperties(fname=path).get_name()=='Arial'
    plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.labelsize':9.5,'axes.titlesize':10.5,'xtick.labelsize':9,'ytick.labelsize':9.5,'axes.linewidth':.7,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white','savefig.facecolor':'white'})
    return path

def text_color(rgba):
    rgb=np.array(rgba[:3]);linear=np.where(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)**2.4)
    luminance=float(np.dot(linear,[.2126,.7152,.0722]))
    white=1.05/(luminance+.05);black=(luminance+.05)/.05
    return 'white' if white>black else 'black'

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--font-dir',default=str(HERE.parent/'figure_fonts'));args=parser.parse_args()
    font=configure_font(args.font_dir)
    data=load_data()
    assert data.n_cells.tolist()==[773,1311,595,238,236]
    assert data.joint_positive.tolist()==[118,179,0,1,1]
    for gene in GENES:
        assert data[gene+'_positive'].between(0,data.n_cells).all()
        assert np.allclose(data[gene+'_pct'],100*data[gene+'_positive']/data.n_cells)
    assert (data.joint_positive<=data[[g+'_positive' for g in JOINT_GENES]].min(axis=1)).all()
    matrix=data[[g+'_pct' for g in GENES]].to_numpy()

    fig=plt.figure(figsize=(7.2,7.0),dpi=180)
    fig.text(.045,.973,'Transcript detection within expanded Sézary syndrome clonotypes',fontsize=12,fontweight='bold',ha='left',va='top')
    fig.text(.045,.895,'a',fontsize=14,fontweight='bold',ha='left',va='center')
    fig.text(.30,.895,'Individual transcript detection',fontsize=10.5,fontweight='bold',ha='left',va='center')
    ax=fig.add_axes([.30,.535,.53,.30])
    mesh=ax.pcolormesh(np.arange(5),np.arange(6),matrix,cmap='Reds',vmin=0,vmax=100,edgecolors='white',linewidth=.8,rasterized=False)
    ax.set_xlim(0,4);ax.set_ylim(5,0)
    ax.set_xticks(np.arange(4)+.5,GENES,fontstyle='italic',fontsize=10)
    ax.xaxis.tick_top();ax.tick_params(axis='x',length=0,pad=8)
    ax.set_yticks(np.arange(5)+.5,data.plot_label);ax.tick_params(axis='y',length=0,pad=10)
    for spine in ax.spines.values():spine.set_visible(False)
    cmap=plt.get_cmap('Reds')
    for i in range(5):
        for j in range(4):
            value=matrix[i,j]
            ax.text(j+.5,i+.5,f'{value:.1f}%',ha='center',va='center',fontsize=10,color=text_color(cmap(value/100)))
        ax.text(1.14,i+.5,f'{int(data.iloc[i].n_cells):,}',transform=ax.get_yaxis_transform(),ha='center',va='center',fontsize=9.5,clip_on=False)
    ax.text(1.14,1.045,'Cells',transform=ax.transAxes,ha='center',va='bottom',fontsize=9.5,fontweight='bold',clip_on=False)
    ax.axhline(3,color='#8A8A8A',linewidth=.6,linestyle=(0,(2,3)))
    cax=fig.add_axes([.55,.487,.28,.013])
    cb=fig.colorbar(mesh,cax=cax,orientation='horizontal',ticks=[0,50,100])
    cb.solids.set_rasterized(False)
    cb.outline.set_visible(False);cb.ax.tick_params(length=2,width=.5,pad=2,labelsize=8)
    cb.set_label('Transcript-positive cells (%)',fontsize=8.5,labelpad=3)

    fig.text(.045,.408,'b',fontsize=14,fontweight='bold',ha='left',va='center')
    fig.text(.30,.408,'Joint detection of three transcripts',fontsize=10.5,fontweight='bold',ha='left',va='center')
    fig.text(.30,.386,'PTGDR2, HPGDS, and CRLF2',fontsize=8.5,fontstyle='italic',ha='left',va='center')
    bx=fig.add_axes([.30,.145,.53,.225])
    y=np.arange(5)
    bx.barh(y,data.joint_pct,height=.58,color='#1F4E79',edgecolor='#163A5C',linewidth=.6,zorder=3)
    for i,(_,row) in enumerate(data.iterrows()):
        bx.text(row.joint_pct+.35,i,f'{row.joint_pct:.1f}%',ha='left',va='center',fontsize=9)
        bx.text(1.14,i,f'{int(row.joint_positive):,}/{int(row.n_cells):,}',transform=bx.get_yaxis_transform(),ha='center',va='center',fontsize=9,clip_on=False)
    bx.text(1.14,1.08,'Joint / total',transform=bx.transAxes,ha='center',va='bottom',fontsize=8.5,clip_on=False)
    bx.set_yticks(y,data.plot_label);bx.set_ylim(4.6,-.6)
    bx.set_xlim(0,20);bx.set_xticks([0,5,10,15,20])
    bx.set_xlabel('Cells with joint transcript detection (%)',labelpad=7)
    bx.tick_params(axis='y',length=0,pad=10);bx.tick_params(axis='x',length=3,width=.65,pad=4)
    bx.spines[['top','right','left']].set_visible(False)
    bx.set_axisbelow(True);bx.grid(axis='x',color='#E4E4E4',linewidth=.5)
    bx.axhline(2.5,color='#BDBDBD',linestyle=(0,(2,3)),linewidth=.6,zorder=1)
    fig.text(.045,.045,'HLA08 reactions 1 and 2 are from one specimen; HLA19 is a second patient.',ha='left',va='center',fontsize=8.7)
    for ext in ['png','pdf','svg','tiff']:
        options={'dpi':600,'facecolor':'white'}
        if ext=='tiff':options['pil_kwargs']={'compression':'tiff_lzw'}
        fig.savefig(HERE/f'Figure_2.{ext}',**options)
    fig.savefig(HERE/'Figure_2_preview.png',dpi=180,facecolor='white')
    plt.close(fig)
    legend=('Figure 2. Transcript detection within expanded Sézary syndrome clonotypes. '
        '(a) Percentage of cells with at least one UMI for each transcript. Darker red indicates higher detection. Numbers are detection percentages; cell counts are shown at right. '
        '(b) Joint detection of PTGDR2, HPGDS, and CRLF2 in the same cell, requiring at least one UMI for each transcript. Labels give percentages and positive/total cell counts. '
        'The HLA08 clonotype is shown separately for two reactions from one specimen. Other HLA08 cells comprise recovered paired-TCR clonotypes other than the target; cells lacking paired-TCR recovery were excluded. The two HLA19 clonotypes (pairs 1 and 2) were analyzed separately. Comparisons are descriptive; the comparator was not restricted to CD4 T cells.')
    (HERE/'Figure_2_legend.txt').write_text(legend+'\n')
    readme=('Source: public processed gene-expression and paired TCR data from GEO GSE263881, as reconstructed in CTCL_public_data_audit.zip.\n'
        'The source CSV gives raw transcript-positive counts, denominators, percentages, and clonotype sequence identities for every plotted population. Joint detection is computed per cell, not inferred from marginal percentages.\n'
        'All plotted comparisons are descriptive. Transcript detection does not establish protein expression, functional signaling, malignant identity, or a definitive Th2A phenotype. Technical reactions do not establish patient-level replication.\n'
        'The PDF is vector artwork with embedded Arial subsets. PNG and LZW-compressed TIFF are 600 dpi. SVG preserves editable text and requires Arial for matching typography.\n'
        'The heatmap uses the sequential Reds palette with a fixed 0-100% range. Bars use navy (#1F4E79), matching Figure 1.\n'
        'Reproduce with Python, numpy, pandas, and matplotlib: run generate_figure_2.py --font-dir /path/to/Arial/fonts. Arial font files are not included. The source CSV is used directly.\n')
    (HERE/'README.txt').write_text(readme)
    print(json.dumps({'font':font,'legend_words':len(legend.split()),'data':data[['population','n_cells','joint_positive','joint_pct']].to_dict('records'),'files':[{'name':p.name,'bytes':p.stat().st_size} for p in HERE.iterdir() if p.suffix in ['.png','.pdf','.tiff','.svg','.zip']]},indent=2))

if __name__=='__main__':main()
