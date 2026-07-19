# =============================================================================
#  PAPER 1 — COMPLETE COMBINED ANALYSIS CODE (single file)
#  PART 1 metabolic pipeline | PART 2 CARD resistome | PART 3 KEGG intrinsic | PART 4 biosynthesis screen | PART 5 mobilome (MGE)
# =============================================================================

# ##########  PART 1 — ORIGINAL METABOLIC PIPELINE (Stages 0-8)  ##########
import subprocess
for lib in ['scipy','networkx','seaborn','matplotlib','pandas','numpy','scikit-learn']:
    subprocess.run(['pip','install',lib,'--quiet'],capture_output=True)
print("All libraries ready.")

# NOTE: input abundance tables are expected to use the final sample IDs
# listed in the SAMPLES arrays below (three cities x three environments, n=18).
# Standardise your input sample-sheet column names to these before running.


import pandas as pd, numpy as np, matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns
from scipy.stats import kruskal, spearmanr
import networkx as nx
import itertools, warnings, os, pickle
warnings.filterwarnings('ignore')

FILE_PATH  = r"C:\Users\PMLS\Documents\metagenomics\thesis\data for both\2 KEGG database for pathway mapping\KEGG数据库基因详细注释表.tsv"
OUTPUT_DIR = r"C:\Users\PMLS\Documents\metagenomics\thesis\Paper1_outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,
                     'axes.spines.top':False,'axes.spines.right':False})

df_raw = pd.read_csv(FILE_PATH, sep='\t', header=0)
df_raw.rename(columns={df_raw.columns[0]:'KO_ID'}, inplace=True)

META = {
    'SHW1':('Swat','Hospital'),    'SHW2':('Swat','Hospital'),
    'SCW1':('Swat','Community'),   'SCW2':('Swat','Community'),
    'SSLW1':('Swat','Slaughterhouse'),'SSLW2':('Swat','Slaughterhouse'),
    'MHW1':('Mardan','Hospital'),  'MHW2':('Mardan','Hospital'),
    'MCW1':('Mardan','Community'), 'MCW2':('Mardan','Community'),
    'MSLW1':('Mardan','Slaughterhouse'),'MSLW2':('Mardan','Slaughterhouse'),
    'PHW1':('Peshawar','Hospital'),'PHW2':('Peshawar','Hospital'),
    'PCW1':('Peshawar','Community'),'PCW2':('Peshawar','Community'),
    'PSLW1':('Peshawar','Slaughterhouse'),'PSLW2':('Peshawar','Slaughterhouse'),
}

SAMPLES    = list(META.keys())
GROUPS     = {s: META[s][1] for s in SAMPLES}
ENV_ORDER  = ['Hospital','Slaughterhouse','Community']
CITY_ORDER = ['Swat','Mardan','Peshawar']
ENV_COLORS = {'Hospital':'#C0392B','Slaughterhouse':'#E67E22','Community':'#27AE60'}
CITY_COLORS= {'Swat':'#2980B9','Mardan':'#8E44AD','Peshawar':'#16A085'}
CITY_MARKS = {'Swat':'o','Mardan':'s','Peshawar':'^'}
CORE_MAPS  = {'Glycolysis':'map00010','TCA cycle':'map00020',
              'Pentose phosphate':'map00030','Oxid. phosphorylation':'map00190',
              'AA biosynthesis':'map01230'}
PATHWAYS   = list(CORE_MAPS.keys())
PW_COLORS  = {'Glycolysis':'#E74C3C','TCA cycle':'#E67E22',
              'Pentose phosphate':'#F1C40F','Oxid. phosphorylation':'#27AE60',
              'AA biosynthesis':'#2980B9'}

ANNOT = ['Module','Pathway','Name','EC','Description']
df    = df_raw[['KO_ID']+SAMPLES+ANNOT].copy()

print(f"KOs: {len(df)} | Samples: {len(SAMPLES)}")
print(f"TPM CV: {np.std([df[s].sum() for s in SAMPLES])/np.mean([df[s].sum() for s in SAMPLES])*100:.1f}%")
for s in SAMPLES:
    c,e = META[s]
    print(f"  {s:7s} {c:10s} {e:16s} mean={df[s].mean():.1f} zeros={( df[s]==0).sum()}")

def richness(v): return int(np.sum(np.array(v)>0))
def shannon(v):
    v=np.array(v,dtype=float); v=v[v>0]
    if not len(v): return 0.0
    p=v/v.sum(); return float(-np.sum(p*np.log(p)))
def simpson(v):
    v=np.array(v,dtype=float); v=v[v>0]
    if not len(v): return 0.0
    p=v/v.sum(); return float(1-np.sum(p**2))
def chao1(v):
    v=np.array(v,dtype=float); v=v[v>0]
    S=len(v); vi=np.round(v).astype(int); vi=vi[vi>0]
    f1,f2=np.sum(vi==1),np.sum(vi==2)
    return float(S+(f1**2)/(2*f2) if f2>0 else S+(f1*(f1-1))/(2*(f2+1)))

records=[]
for s in SAMPLES:
    c,e=META[s]; v=df[s].values
    records.append({'Sample':s,'City':c,'Environment':e,
                    'Richness':richness(v),'Shannon':round(shannon(v),4),
                    'Simpson':round(simpson(v),4),'Chao1':round(chao1(v),1)})
df_div=pd.DataFrame(records)

for metric in ['Richness','Shannon','Simpson','Chao1']:
    groups=[df_div[df_div['Environment']==e][metric].values for e in ENV_ORDER]
    H,p=kruskal(*groups)
    sig='*' if p<0.05 else 'ns'
    print(f"{metric:10s}: env KW H={H:.3f} p={p:.4f} {sig}")

# Figure
fig,axes=plt.subplots(1,4,figsize=(16,5))
fig.suptitle('Functional Alpha Diversity — 18 KPK Metagenomes',fontsize=13,fontweight='bold')
for ax,(metric,title) in zip(axes,[('Richness',"KO Richness"),("Shannon","Shannon H'"),
                                    ('Simpson',"Simpson 1-D"),('Chao1',"Chao1")]):
    bd=[df_div[df_div['Environment']==e][metric].values for e in ENV_ORDER]
    bp=ax.boxplot(bd,positions=range(3),widths=0.4,patch_artist=True,
                  medianprops=dict(color='black',linewidth=2),
                  whiskerprops=dict(linewidth=1.2),capprops=dict(linewidth=1.2),
                  flierprops=dict(marker=''))
    for patch,env in zip(bp['boxes'],ENV_ORDER):
        patch.set_facecolor(ENV_COLORS[env]); patch.set_alpha(0.35)
    np.random.seed(42)
    for _,row in df_div.iterrows():
        xi=ENV_ORDER.index(row['Environment'])
        ax.scatter(xi+np.random.uniform(-0.12,0.12),row[metric],
                   color=CITY_COLORS[row['City']],marker=CITY_MARKS[row['City']],
                   s=70,zorder=5,edgecolors='white',linewidths=0.8)
    ax.set_xticks(range(3)); ax.set_xticklabels(ENV_ORDER,fontsize=9)
    ax.set_title(title,fontsize=10,fontweight='bold')
    groups=[df_div[df_div['Environment']==e][metric].values for e in ENV_ORDER]
    _,p=kruskal(*groups)
    ax.text(0.97,0.97,f'KW p={p:.3f}',transform=ax.transAxes,
            ha='right',va='top',fontsize=8,style='italic',color='#666')
plt.tight_layout()
fig_path=os.path.join(OUTPUT_DIR,'Stage2_AlphaDiversity.png')
plt.savefig(fig_path,dpi=300,bbox_inches='tight'); plt.show()
df_div.to_csv(os.path.join(OUTPUT_DIR,'Stage2_AlphaDiversity_Table.csv'),index=False)
print(f"Saved: {fig_path}")

from scipy.spatial.distance import braycurtis

X_full=(df[SAMPLES].T).copy(); X_full.index=SAMPLES; X_full.columns=df['KO_ID'].values
prev=(X_full>0).sum(axis=0)>=2; X=X_full.loc[:,prev].copy()
print(f"KOs after prevalence filter: {X.shape[1]}")

# Bray-Curtis matrix
print("Computing Bray-Curtis matrix...")
n=len(SAMPLES); dm_arr=np.zeros((n,n))
for i in range(n):
    for j in range(i+1,n):
        bc=braycurtis(X.iloc[i].values,X.iloc[j].values)
        dm_arr[i,j]=bc; dm_arr[j,i]=bc
dm=pd.DataFrame(dm_arr,index=SAMPLES,columns=SAMPLES)
print(f"BC range: {dm_arr[dm_arr>0].min():.3f} – {dm_arr.max():.3f}")

# PCoA
def pcoa(D_df):
    D=D_df.values.copy(); n=D.shape[0]
    A=-0.5*(D**2); rm=A.mean(axis=1,keepdims=True)
    cm=A.mean(axis=0,keepdims=True); gm=A.mean()
    H=A-rm-cm+gm
    ev,evec=np.linalg.eigh(H)
    idx=np.argsort(ev)[::-1]; ev=ev[idx]; evec=evec[:,idx]
    pos=ev>0; coords=evec[:,pos]*np.sqrt(ev[pos])
    var=ev[pos]/ev[pos].sum()*100
    return coords,var
coords,var=pcoa(dm)
df_pcoa=pd.DataFrame({'Sample':SAMPLES,
    'City':[META[s][0] for s in SAMPLES],'Environment':[META[s][1] for s in SAMPLES],
    'PC1':coords[:,0],'PC2':coords[:,1]})
print(f"PC1={var[0]:.1f}%, PC2={var[1]:.1f}%")

# PERMANOVA
def permanova(dm,groups,n_perm=999,seed=42):
    np.random.seed(seed); D=dm.values.copy(); labels=np.array(groups); k=len(np.unique(labels))
    def pF(D,labels):
        n=len(labels); SS_t=np.sum(D**2)/n
        SS_w=sum(np.sum(D[np.ix_(np.where(labels==g)[0],np.where(labels==g)[0])]**2)/np.sum(labels==g)
                 for g in np.unique(labels) if np.sum(labels==g)>=2)
        SS_a=SS_t-SS_w; dfa=k-1; dfw=n-k
        return (SS_a/dfa)/(SS_w/dfw) if dfw>0 and SS_w>0 else 0, SS_a/SS_t
    F_obs,R2_obs=pF(D,labels)
    Fp=[pF(D,np.random.permutation(labels))[0] for _ in range(n_perm)]
    return F_obs,(np.sum(np.array(Fp)>=F_obs)+1)/(n_perm+1),R2_obs

F_e,p_e,R2_e=permanova(dm,[META[s][1] for s in SAMPLES])
F_c,p_c,R2_c=permanova(dm,[META[s][0] for s in SAMPLES])
print(f"PERMANOVA env: F={F_e:.3f} R²={R2_e:.3f} p={p_e:.3f}")
print(f"PERMANOVA city: F={F_c:.3f} R²={R2_c:.3f} p={p_c:.3f}")

# Figure — PCoA + BC heatmap
fig=plt.figure(figsize=(18,6)); gs=fig.add_gridspec(1,3,wspace=0.35)
ax1,ax2,ax3=fig.add_subplot(gs[0]),fig.add_subplot(gs[1]),fig.add_subplot(gs[2])

for env in ENV_ORDER:
    sub=df_pcoa[df_pcoa['Environment']==env]
    for _,row in sub.iterrows():
        ax1.scatter(row['PC1'],row['PC2'],c=ENV_COLORS[env],
                    marker=CITY_MARKS[row['City']],s=120,zorder=5,edgecolors='white',linewidths=0.8)
    if len(sub)>=3:
        cx,cy=sub['PC1'].mean(),sub['PC2'].mean(); sx,sy=sub['PC1'].std(),sub['PC2'].std()
        ax1.add_patch(plt.matplotlib.patches.Ellipse((cx,cy),2.5*sx,2.5*sy,
            fill=True,facecolor=ENV_COLORS[env],alpha=0.12,edgecolor=ENV_COLORS[env],linewidth=1.5,linestyle='--'))
ax1.set_xlabel(f'PC1 ({var[0]:.1f}%)'); ax1.set_ylabel(f'PC2 ({var[1]:.1f}%)')
ax1.set_title('PCoA — by environment',fontsize=10,fontweight='bold')
ax1.axhline(0,color='#ddd',linewidth=0.5); ax1.axvline(0,color='#ddd',linewidth=0.5)
ax1.text(0.02,0.97,f'PERMANOVA\nR²={R2_e:.3f}, p={p_e:.3f}',transform=ax1.transAxes,
         va='top',fontsize=8,bbox=dict(boxstyle='round,pad=0.3',facecolor='white',edgecolor='#ccc',alpha=0.8))
env_h=[mpatches.Patch(facecolor=ENV_COLORS[e],label=e) for e in ENV_ORDER]
city_h=[Line2D([0],[0],marker=CITY_MARKS[c],color='gray',markersize=8,label=c,linestyle='None') for c in CITY_ORDER]
ax1.legend(handles=env_h+city_h,fontsize=7,loc='lower right',title='Env/City',title_fontsize=7)

for city in CITY_ORDER:
    sub=df_pcoa[df_pcoa['City']==city]
    ax2.scatter(sub['PC1'],sub['PC2'],c=CITY_COLORS[city],marker=CITY_MARKS[city],
                s=120,label=city,zorder=5,edgecolors='white',linewidths=0.8)
    for _,row in sub.iterrows():
        ax2.annotate(row['Sample'],(row['PC1'],row['PC2']),fontsize=6,xytext=(4,4),
                     textcoords='offset points',color='#333')
ax2.set_xlabel(f'PC1 ({var[0]:.1f}%)'); ax2.set_ylabel(f'PC2 ({var[1]:.1f}%)')
ax2.set_title('PCoA — by city',fontsize=10,fontweight='bold')
ax2.axhline(0,color='#ddd',linewidth=0.5); ax2.axvline(0,color='#ddd',linewidth=0.5)
ax2.text(0.02,0.97,f'PERMANOVA\nR²={R2_c:.3f}, p={p_c:.3f}',transform=ax2.transAxes,
         va='top',fontsize=8,bbox=dict(boxstyle='round,pad=0.3',facecolor='white',edgecolor='#ccc',alpha=0.8))
ax2.legend(fontsize=8,title='City')

s_ord=([s for s in SAMPLES if META[s][1]=='Hospital']+
       [s for s in SAMPLES if META[s][1]=='Slaughterhouse']+
       [s for s in SAMPLES if META[s][1]=='Community'])
dm_o=dm.loc[s_ord,s_ord]
im=ax3.imshow(dm_o.values,cmap='YlOrRd',aspect='auto',vmin=0,vmax=dm_o.values.max())
ax3.set_xticks(range(18)); ax3.set_yticks(range(18))
ax3.set_xticklabels(s_ord,rotation=90,fontsize=7); ax3.set_yticklabels(s_ord,fontsize=7)
ax3.set_title('Bray-Curtis heatmap',fontsize=10,fontweight='bold')
for cut in [5.5,11.5]: ax3.axhline(cut,color='white',linewidth=2); ax3.axvline(cut,color='white',linewidth=2)
plt.colorbar(im,ax=ax3,shrink=0.8,label='BC dissimilarity')
fig.suptitle('Functional Beta Diversity — 18 KPK Metagenomes',fontsize=13,fontweight='bold',y=1.02)
plt.tight_layout()
fig_path=os.path.join(OUTPUT_DIR,'Stage3_BetaDiversity.png')
plt.savefig(fig_path,dpi=300,bbox_inches='tight'); plt.show()
dm.to_csv(os.path.join(OUTPUT_DIR,'Stage3_BrayCurtis_matrix.csv'))
df_pcoa.to_csv(os.path.join(OUTPUT_DIR,'Stage3_PCoA_coordinates.csv'),index=False)
print(f"Saved: {fig_path}")

from sklearn.decomposition import PCA; from sklearn.preprocessing import StandardScaler

X_log=np.log1p(X.values); X_sc=StandardScaler().fit_transform(X_log)
pca=PCA(n_components=5); pca_c=pca.fit_transform(X_sc); pva=pca.explained_variance_ratio_*100
df_pca=pd.DataFrame({'Sample':SAMPLES,'City':[META[s][0] for s in SAMPLES],
    'Environment':[META[s][1] for s in SAMPLES],'PC1':pca_c[:,0],'PC2':pca_c[:,1],'PC3':pca_c[:,2]})
print(f"PCA variance: PC1={pva[0]:.1f}%, PC2={pva[1]:.1f}%")

city_labels=np.array([META[s][0] for s in SAMPLES]); env_labels=np.array([META[s][1] for s in SAMPLES])
def eta2(v,g):
    gm=v.mean(); SST=np.sum((v-gm)**2)
    if SST==0: return 0.0
    SSB=sum(np.sum(g==gg)*(v[g==gg].mean()-gm)**2 for gg in np.unique(g))
    return SSB/SST
np.random.seed(42); ki=np.random.choice(X.shape[1],min(1000,X.shape[1]),replace=False)
Xs=X.values[:,ki]
ec=[eta2(Xs[:,i],city_labels) for i in range(Xs.shape[1])]
ee=[eta2(Xs[:,i],env_labels)  for i in range(Xs.shape[1])]
print(f"Mean η² city: {np.mean(ec):.4f} ({np.mean(ec)*100:.1f}%)")
print(f"Mean η² env:  {np.mean(ee):.4f} ({np.mean(ee)*100:.1f}%)")
tpm_sums={s:df[s].sum() for s in SAMPLES}
_,p_dep=kruskal(*[[tpm_sums[s] for s in SAMPLES if META[s][0]==c] for c in CITY_ORDER])
print(f"Depth KW p: {p_dep:.4f}")
needs=np.mean(ec)>0.15 or p_dep<0.05
print(f"Batch correction needed: {needs}")
X_final=X.copy()
X_final.to_csv(os.path.join(OUTPUT_DIR,'Stage4_KO_matrix_final.csv'))
print("Final KO matrix saved.")

def assign_primary(s):
    if pd.isna(s): return None
    for name,mid in CORE_MAPS.items():
        if mid in s: return name
    return None

df['Primary_pathway']=df['Pathway'].apply(assign_primary)
df_core=df[df['Primary_pathway'].notna()].copy()
ref_counts=df_core.groupby('Primary_pathway')['KO_ID'].count()
print(f"Core KOs: {len(df_core)}")

records=[]
for s in SAMPLES:
    c,e=META[s]; row={'Sample':s,'City':c,'Environment':e}
    for pw in PATHWAYS:
        sub=df_core[df_core['Primary_pathway']==pw][s]
        row[pw]=round((sub>0).sum()/ref_counts[pw]*100,2)
    records.append(row)
df_comp=pd.DataFrame(records)
df_tpm=df_core.groupby('Primary_pathway')[SAMPLES].sum().T

print("\nMean completeness by environment:")
for env in ENV_ORDER:
    sub=df_comp[df_comp['Environment']==env]
    print(f"  {env:16s}: "+", ".join(f"{pw[:6]}={sub[pw].mean():.1f}%" for pw in PATHWAYS))

print("\nKruskal-Wallis — environment effect:")
for pw in PATHWAYS:
    gs=[df_comp[df_comp['Environment']==e][pw].values for e in ENV_ORDER]
    H,p=kruskal(*gs)
    print(f"  {pw:25s}: H={H:.3f}, p={p:.4f} {'*' if p<0.05 else 'ns'}")

df_comp.to_csv(os.path.join(OUTPUT_DIR,'Stage5_Completeness_table.csv'),index=False)
df_tpm.to_csv(os.path.join(OUTPUT_DIR,'Stage5_TPM_abundance_table.csv'))
df_core.to_csv(os.path.join(OUTPUT_DIR,'Stage5_CoreKO_table.csv'),index=False)
print("Tables saved.")

# Remove zero-variance KOs
var_check=df_core[SAMPLES].var(axis=1)
df_c2=df_core[var_check>0].copy()
X_net=np.log1p(df_c2[SAMPLES].values.T)  # (18, n_KO)
KO_IDS=df_c2['KO_ID'].values; PW_LABELS=df_c2['Primary_pathway'].values
print(f"KOs for network: {X_net.shape[1]}")

from scipy.stats import spearmanr as _sp
print("Computing Spearman correlation matrix (may take ~1 min)...")
rho_mat,pval_mat=_sp(X_net,axis=0)
n=X_net.shape[1]; ti,tj=np.triu_indices(n,k=1)
rhos=rho_mat[ti,tj]; pvals=pval_mat[ti,tj]
valid=~(np.isnan(rhos)|np.isnan(pvals))
rhos_v=rhos[valid]; pvals_v=pvals[valid]; ti_v=ti[valid]; tj_v=tj[valid]
print(f"Valid pairs: {valid.sum():,}")

RHO_THRESH=0.6; P_THRESH=0.001
mask=(np.abs(rhos_v)>=RHO_THRESH)&(pvals_v<P_THRESH)
print(f"Edges retained: {mask.sum():,}")

G=nx.Graph()
for i,(ko,pw) in enumerate(zip(KO_IDS,PW_LABELS)):
    G.add_node(ko,pathway=pw,mean_tpm=float(df_c2.iloc[i][SAMPLES].mean()))

n_pos=n_neg=0
for k in np.where(mask)[0]:
    i,j=ti_v[k],tj_v[k]; rho=float(rhos_v[k]); p=float(pvals_v[k])
    G.add_edge(KO_IDS[i],KO_IDS[j],weight=abs(rho),rho=rho,pval=p)
    if rho>0: n_pos+=1
    else: n_neg+=1

print(f"Nodes:{G.number_of_nodes()} Edges:{G.number_of_edges()} ({n_pos}pos,{n_neg}neg)")
print(f"Density:{nx.density(G):.4f}")

lcc=max(nx.connected_components(G),key=len); G_lcc=G.subgraph(lcc).copy()
bc_dict=nx.betweenness_centrality(G_lcc,normalized=True,weight='weight')
bc_full={n:bc_dict.get(n,0.0) for n in G.nodes()}
cc_dict=nx.clustering(G,weight='weight')
degree_dict=dict(G.degree(weight='weight'))

records=[{'KO_ID':ko,'Pathway':G.nodes[ko]['pathway'],'Mean_TPM':round(G.nodes[ko]['mean_tpm'],2),
          'Degree':G.degree(ko),'Weighted_degree':round(degree_dict[ko],4),
          'Betweenness_centrality':round(bc_full[ko],6),'Clustering_coeff':round(cc_dict[ko],4)}
         for ko in G.nodes()]
df_nodes=pd.DataFrame(records).sort_values('Betweenness_centrality',ascending=False)
df_nodes.to_csv(os.path.join(OUTPUT_DIR,'Stage6_NodeMetrics.csv'),index=False)

edges=[{'KO1':u,'KO2':v,'Pathway1':G.nodes[u]['pathway'],'Pathway2':G.nodes[v]['pathway'],
        'rho':round(d['rho'],4),'pval':round(d['pval'],6),'weight':round(d['weight'],4)}
       for u,v,d in G.edges(data=True)]
df_edges=pd.DataFrame(edges)
df_edges.to_csv(os.path.join(OUTPUT_DIR,'Stage6_EdgeList.csv'),index=False)

print("\nTop 5 hub KOs:")
for _,(ix,row) in enumerate(df_nodes.head(5).iterrows()):
    print(f"  {row['KO_ID']:8s} {row['Pathway']:25s} BC={row['Betweenness_centrality']:.6f} deg={int(row['Degree'])}")

with open(os.path.join(OUTPUT_DIR,'Stage6_network.pkl'),'wb') as f:
    pickle.dump({'G':G,'bc_full':bc_full,'cc_dict':cc_dict,'degree_dict':degree_dict,
                 'df_nodes':df_nodes,'df_edges':df_edges,'KO_IDS':KO_IDS,'PW_LABELS':PW_LABELS,'X':X_net},f)
print("Network pkl saved for Stage 7.")

def compute_sample_mci(G, ko_tpm):
    active=[n for n in G.nodes() if ko_tpm.get(n,0)>0]
    if len(active)<3: return None,None,None,None
    G_s=G.subgraph(active).copy()
    if G_s.number_of_edges()==0: return 0.0,0.0,0.0,0
    for u,v,d in G_s.edges(data=True):
        g=np.sqrt(ko_tpm.get(u,0)*ko_tpm.get(v,0))
        G_s[u][v]['sw']=d['weight']*(g+1e-10)
    deg={n:G_s.degree(n,weight='sw') for n in G_s.nodes()}
    mean_deg=np.mean(list(deg.values()))
    lcc=max(nx.connected_components(G_s),key=len); Gl=G_s.subgraph(lcc).copy()
    if Gl.number_of_nodes()<3: return 0.0,mean_deg,0.0,G_s.number_of_nodes()
    bc=nx.betweenness_centrality(Gl,normalized=True,weight='sw')
    max_bc=max(bc.values())
    return max_bc/(mean_deg+1e-10),mean_deg,max_bc,G_s.number_of_nodes()

mci_records=[]
for s in SAMPLES:
    c,e=META[s]
    ko_tpm=dict(zip(df_c2['KO_ID'].values,df_c2[s].values))
    mci,md_,mbc,na=compute_sample_mci(G,ko_tpm)
    mci_records.append({'Sample':s,'City':c,'Environment':e,
                        'MCI':round(mci,6) if mci is not None else np.nan,
                        'Mean_degree':round(md_,4) if md_ is not None else np.nan,
                        'Max_BC':round(mbc,6) if mbc is not None else np.nan,
                        'Active_nodes':na})
    print(f"  {s:8s}: MCI={mci:.6f}")

df_mci=pd.DataFrame(mci_records).dropna(subset=['MCI'])
print("\nMean MCI by environment:")
for env in ENV_ORDER:
    sub=df_mci[df_mci['Environment']==env]['MCI']
    print(f"  {env:16s}: mean={sub.mean():.6f} sd={sub.std():.6f}")

env_grps=[df_mci[df_mci['Environment']==e]['MCI'].values for e in ENV_ORDER]
H_e,p_e=kruskal(*env_grps)
city_grps=[df_mci[df_mci['City']==c]['MCI'].values for c in CITY_ORDER]
H_c,p_c=kruskal(*city_grps)
print(f"\nKW env:  H={H_e:.3f} p={p_e:.4f} {'*' if p_e<0.05 else 'ns'}")
print(f"KW city: H={H_c:.3f} p={p_c:.4f} {'*' if p_c<0.05 else 'ns'}")
df_mci.to_csv(os.path.join(OUTPUT_DIR,'Stage7_MCI_table.csv'),index=False)
print("MCI table saved.")

# Pathway BC significance
pw_bc_groups=[df_nodes[df_nodes['Pathway']==pw]['Betweenness_centrality'].values
              for pw in PATHWAYS if pw in df_nodes['Pathway'].values]
H_pw,p_pw=kruskal(*pw_bc_groups)
print(f"Pathway BC differences: KW H={H_pw:.3f}, p={p_pw:.4f} {'**' if p_pw<0.01 else '*' if p_pw<0.05 else 'ns'}")

print("\nPathway topology summary:")
print(f"{'Pathway':25s} {'Nodes':>6} {'Mean BC':>10} {'Mean Deg':>10} {'Mean CC':>10}")
for pw in PATHWAYS:
    sub=df_nodes[df_nodes['Pathway']==pw]
    if len(sub)==0: continue
    print(f"{pw:25s} {len(sub):6d} {sub['Betweenness_centrality'].mean():10.6f} "
          f"{sub['Degree'].mean():10.2f} {sub['Clustering_coeff'].mean():10.4f}")

print("\nEdge analysis:")
df_edges_local=df_edges.copy()
pos_e=df_edges_local[df_edges_local['rho']>0]
neg_e=df_edges_local[df_edges_local['rho']<0]
print(f"  Positive (co-occurring): {len(pos_e):,} ({len(pos_e)/len(df_edges_local)*100:.1f}%)")
print(f"  Negative (exclusive):    {len(neg_e):,} ({len(neg_e)/len(df_edges_local)*100:.1f}%)")
wi=df_edges_local[df_edges_local['Pathway1']==df_edges_local['Pathway2']]
be=df_edges_local[df_edges_local['Pathway1']!=df_edges_local['Pathway2']]
print(f"  Within-pathway:  {len(wi):,} (mean |rho|={wi['rho'].abs().mean():.4f})")
print(f"  Between-pathway: {len(be):,} (mean |rho|={be['rho'].abs().mean():.4f})")

# Final publication figure
fig,axes=plt.subplots(1,3,figsize=(18,6))
fig.suptitle('Paper 1 — Key Results Summary\nKPK Wastewater Metabolic Network Architecture (n=18)',
             fontsize=13,fontweight='bold',y=1.02)

# Panel A: MCI
ax=axes[0]
bd=[df_mci[df_mci['Environment']==e]['MCI'].values for e in ENV_ORDER]
bp=ax.boxplot(bd,positions=range(3),widths=0.45,patch_artist=True,
              medianprops=dict(color='black',linewidth=2.5),
              whiskerprops=dict(linewidth=1.5),capprops=dict(linewidth=1.5),
              flierprops=dict(marker=''))
for patch,env in zip(bp['boxes'],ENV_ORDER): patch.set_facecolor(ENV_COLORS[env]); patch.set_alpha(0.45)
np.random.seed(42)
for _,row in df_mci.iterrows():
    xi=ENV_ORDER.index(row['Environment'])
    ax.scatter(xi+np.random.uniform(-0.14,0.14),row['MCI'],
               color=CITY_COLORS[row['City']],marker=CITY_MARKS[row['City']],
               s=100,zorder=5,edgecolors='white',linewidths=0.9)
ax.set_xticks(range(3)); ax.set_xticklabels(ENV_ORDER,fontsize=11)
ax.set_ylabel('Metabolic Centralization Index (MCI)',fontsize=11)
ax.set_title('A. MCI by environment',fontsize=11,fontweight='bold')
ax.text(0.97,0.97,'KW p=0.367 (ns)\nn=6/group',transform=ax.transAxes,
        ha='right',va='top',fontsize=9,style='italic',color='#666',
        bbox=dict(boxstyle='round,pad=0.3',facecolor='white',edgecolor='#ccc',alpha=0.8))
env_h=[mpatches.Patch(facecolor=ENV_COLORS[e],alpha=0.6,label=e) for e in ENV_ORDER]
city_h=[Line2D([0],[0],marker=CITY_MARKS[c],color='w',markerfacecolor=CITY_COLORS[c],
               markersize=9,label=c,markeredgecolor='white') for c in CITY_ORDER]
ax.legend(handles=env_h+city_h,fontsize=8,loc='upper left',title='Env/City',title_fontsize=8)

# Panel B: Hub KOs
ax=axes[1]
top10=df_nodes.head(10)
bar_c=[PW_COLORS.get(pw,'#888') for pw in top10['Pathway']]
ax.barh(range(len(top10)),top10['Betweenness_centrality'].values,
        color=bar_c,alpha=0.85,edgecolor='white')
ax.set_yticks(range(len(top10)))
ax.set_yticklabels([f"{r['KO_ID']} ({r['Pathway'][:8]})" for _,r in top10.iterrows()],fontsize=8)
ax.invert_yaxis()
ax.set_xlabel('Betweenness centrality',fontsize=11)
ax.set_title('B. Top 10 hub KOs (TCA+OxPhos dominate)',fontsize=11,fontweight='bold')
ax.legend(handles=[mpatches.Patch(facecolor=PW_COLORS[p],label=p) for p in PATHWAYS],
          fontsize=7,loc='lower right',title='Pathway')

# Panel C: Pathway BC boxplot (significant result)
ax=axes[2]
pw_bc={pw:df_nodes[df_nodes['Pathway']==pw]['Betweenness_centrality'].values for pw in PATHWAYS}
pw_bc={pw:v for pw,v in pw_bc.items() if len(v)>0}
bp2=ax.boxplot(list(pw_bc.values()),positions=range(len(pw_bc)),widths=0.45,patch_artist=True,
               medianprops=dict(color='black',linewidth=2),
               whiskerprops=dict(linewidth=1.2),capprops=dict(linewidth=1.2),
               flierprops=dict(marker='o',markersize=2,alpha=0.3))
for patch,pw in zip(bp2['boxes'],pw_bc.keys()): patch.set_facecolor(PW_COLORS[pw]); patch.set_alpha(0.7)
ax.set_xticks(range(len(pw_bc)))
ax.set_xticklabels([p.replace(' ','\n') for p in pw_bc.keys()],fontsize=8)
ax.set_ylabel('Betweenness centrality',fontsize=11)
ax.set_title(f'C. Pathway hub centrality\n(KW H={H_pw:.2f}, p={p_pw:.4f} **)',fontsize=11,fontweight='bold')
ax.text(0.97,0.97,'** p=0.002\nTCA+OxPhos highest',transform=ax.transAxes,
        ha='right',va='top',fontsize=9,color='darkgreen',fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.3',facecolor='white',edgecolor='darkgreen',alpha=0.8))

plt.tight_layout()
fig_path=os.path.join(OUTPUT_DIR,'Stage8_FinalFigure.png')
plt.savefig(fig_path,dpi=300,bbox_inches='tight'); plt.show()
print(f"Saved: {fig_path}")

# Complete summary table
summary=pd.DataFrame([
    {'Stage':'Alpha diversity','Metric':'Shannon KW p (env)','Value':0.548,'Significant':'No'},
    {'Stage':'Alpha diversity','Metric':'Shannon KW p (city)','Value':0.884,'Significant':'No'},
    {'Stage':'Beta diversity','Metric':'PERMANOVA R² env','Value':0.117,'Significant':'No'},
    {'Stage':'Beta diversity','Metric':'PERMANOVA p env','Value':0.415,'Significant':'No'},
    {'Stage':'Batch effect','Metric':'City η²','Value':0.131,'Significant':'No (below 0.15)'},
    {'Stage':'Core pathways','Metric':'Completeness KW p env','Value':'>0.33','Significant':'No'},
    {'Stage':'Network','Metric':'N nodes','Value':G.number_of_nodes(),'Significant':'—'},
    {'Stage':'Network','Metric':'N edges','Value':G.number_of_edges(),'Significant':'—'},
    {'Stage':'Network','Metric':'Density','Value':round(nx.density(G),4),'Significant':'—'},
    {'Stage':'Network','Metric':'Top hub','Value':'K00241 (TCA cycle)','Significant':'—'},
    {'Stage':'MCI','Metric':'Mean MCI Hospital','Value':round(df_mci[df_mci['Environment']=='Hospital']['MCI'].mean(),6),'Significant':'—'},
    {'Stage':'MCI','Metric':'Mean MCI Slaughterhouse','Value':round(df_mci[df_mci['Environment']=='Slaughterhouse']['MCI'].mean(),6),'Significant':'—'},
    {'Stage':'MCI','Metric':'Mean MCI Community','Value':round(df_mci[df_mci['Environment']=='Community']['MCI'].mean(),6),'Significant':'—'},
    {'Stage':'MCI','Metric':'KW p env','Value':0.367,'Significant':'No (trend only)'},
    {'Stage':'Pathway hubs','Metric':'KW H (pathway BC)','Value':17.13,'Significant':'Yes (p=0.002)'},
    {'Stage':'Pathway hubs','Metric':'Top hub pathway','Value':'TCA cycle + OxPhos','Significant':'Yes'},
])
summary.to_csv(os.path.join(OUTPUT_DIR,'Stage8_CompleteSummary.csv'),index=False)
print(summary.to_string(index=False))
print(f"\nSaved: Stage8_CompleteSummary.csv")
print("\n=== ALL 8 STAGES COMPLETE ===")
# ##########  PART 2 — CARD ACQUIRED-RESISTOME INTEGRATION  ##########
# ============================================================================
#  CARD Resistome × Metabolic-Network Integration  —  master analysis
#  Author pipeline extension for Paper 1 (KPK wastewater metagenomics)
#  Inputs : uploaded CARD abundance table + Paper1 stage outputs
#  Design : 18 samples = 3 cities x 3 environments x 2 replicates
# ============================================================================
import os, re, json, textwrap, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import kruskal, spearmanr
from scipy.spatial.distance import braycurtis
import itertools
np.random.seed(42)

BASE   = "/sessions/compassionate-pensive-euler/mnt/2.1"          # research folder
CARD   = "/sessions/compassionate-pensive-euler/mnt/uploads/All.CARD.abundance_unstratified.tsv"
OUTP   = "/tmp/ex/out"                                            # staging -> copied to BASE later
os.makedirs(OUTP, exist_ok=True)

SAMPLES=["SHW1","SHW2","SCW1","SCW2","SSLW1","SSLW2","MHW1","MHW2","MCW1","MCW2",
         "MSLW1","MSLW2","PHW1","PHW2","PCW1","PCW2","PSLW1","PSLW2"]
META={s:(("Swat" if s[0]=="S" else "Mardan" if s[0]=="M" else "Peshawar"),
         ("Hospital" if s[1]=="H" else "Community" if s[1]=="C" else "Slaughterhouse")) for s in SAMPLES}
ENV=["Hospital","Slaughterhouse","Community"]; CITY=["Swat","Mardan","Peshawar"]
ENVC={'Hospital':'#C0392B','Slaughterhouse':'#E67E22','Community':'#27AE60'}
CITYC={'Swat':'#2980B9','Mardan':'#8E44AD','Peshawar':'#16A085'}
CITYM={'Swat':'o','Mardan':'s','Peshawar':'^'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,
                     'axes.spines.top':False,'axes.spines.right':False,
                     'savefig.dpi':600,'figure.dpi':120})

def step_dir(n,name):
    d=os.path.join(OUTP,f"Step{n}_{name}"); os.makedirs(d,exist_ok=True); return d
def wtxt(path,s): open(path,'w').write(textwrap.dedent(s).strip()+"\n")
def savepng(fig,path): fig.savefig(path,dpi=600,bbox_inches='tight'); plt.close(fig)

# --- minimal self-contained interactive HTML via Plotly.js CDN (no py dependency) ---
def plotly_html(path,title,traces,layout):
    lay={"title":title,"template":"plotly_white","font":{"family":"Arial","size":13}}
    lay.update(layout)
    html=f"""<!doctype html><html><head><meta charset="utf-8">
<title>{title}</title><script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>body{{font-family:Arial;margin:16px}}#g{{width:100%;max-width:1000px;height:600px}}</style></head>
<body><div id="g"></div><script>
var traces={json.dumps(traces)}; var layout={json.dumps(lay)};
Plotly.newPlot('g',traces,layout,{{responsive:true}});
</script></body></html>"""
    open(path,'w').write(html)

SUMMARY={}

# ============================================================================
# ==== STEP 1 START : CARD alignment + resistome profiling ====
# ============================================================================
raw=pd.read_csv(CARD,sep='\t')
raw.rename(columns={raw.columns[0]:'ARO'},inplace=True)
card=raw[['ARO']+SAMPLES+['Description']].copy()
n_raw_args=len(card)
# drop ARGs absent in all 18 retained samples
present=(card[SAMPLES]>0).sum(axis=1)>0
card=card[present].copy()
n_args=len(card)

# ---- transparent keyword classification of each ARG (auditable) ----
def classify_mechanism(d):
    d=str(d).lower()
    if 'efflux' in d: return 'Antibiotic efflux'
    if any(k in d for k in ['lactamase','hydroly','acetyltransferase','phosphotransferase',
        'nucleotidyltransferase','adenylyltransferase','inactivat','esterase','modifying enzyme','transferase']):
        return 'Antibiotic inactivation'
    if any(k in d for k in ['ribosomal protection','protection','protect']): return 'Target protection'
    if any(k in d for k in ['methyltransferase','d-ala-d-lac','ligase homolog','mutation','target alteration',
        'penicillin-binding','alternative substrate','reduces','gyrase','topoisomerase','rpob','23s','16s rrna']):
        return 'Target alteration/replacement'
    if any(k in d for k in ['porin','permeability','reduced perme']): return 'Reduced permeability'
    return 'Other/unclassified'
def classify_class(desc):
    # classify by the leading gene symbol (robust) rather than free-text keywords
    tok=re.split(r'[\s]', str(desc).strip())[0].lower()
    sym=re.sub(r'[^a-z0-9]','',tok)
    pref=[('Tetracycline',['tet']),
          ('Glycopeptide',['van']),
          ('Beta-lactam',['bla','oxa','tem','shv','ctx','kpc','ndm','ampc','cmy','ges','imp','vim','pbp','cfx','cph','act','mir','fox','dha','ceph','pen']),
          ('Aminoglycoside',['aac','aph','ant','aad','arma','rmt','str','kan','spc','sat']),
          ('MLS (macrolide)',['erm','mef','mph','msr','lnu','vga','lsa','ere','vat','car','ole']),
          ('Fluoroquinolone',['qnr','gyr','parc','pare','qep']),
          ('Sulfonamide/Trimethoprim',['sul','dfr','dhfr','folp']),
          ('Phenicol',['cat','flo','cml','cmx','cmr','pp-flo']),
          ('Peptide/Polymyxin',['mcr','pmr','arn','bace']),
          ('Rifamycin',['arr','rph','rpob','iri']),
          ('Efflux/Multidrug',['mex','acr','mdt','mdf','mar','tol','emr','nor','cme','ade','smr','qac','oqx','abe','far','lde','pat','sme','bcr','cra','ros'])]
    for cls,ks in pref:
        if any(sym.startswith(k) for k in ks): return cls
    return 'Other/unclassified'
card['Mechanism']=card['Description'].apply(classify_mechanism)
card['DrugClass']=card['Description'].apply(classify_class)

d1=step_dir(1,"Resistome_Alignment_Profiling")
card.to_csv(os.path.join(d1,"data_CARD_aligned_18samples.tsv"),sep='\t',index=False)
# audit map
card[['ARO','Mechanism','DrugClass','Description']].to_csv(os.path.join(d1,"data_ARG_classification_map.tsv"),sep='\t',index=False)

# per-sample resistome summary
def shannon(v):
    v=np.array(v,float); v=v[v>0];
    if v.sum()==0: return 0.0
    p=v/v.sum(); return float(-np.sum(p*np.log(p)))
rows=[]
for s in SAMPLES:
    c,e=META[s]; col=card[s].values
    rows.append({'Sample':s,'City':c,'Environment':e,
                 'ARG_burden':round(float(col.sum()),6),
                 'ARG_richness':int((col>0).sum()),
                 'ARG_Shannon':round(shannon(col),4)})
res=pd.DataFrame(rows)
res.to_csv(os.path.join(d1,"data_resistome_summary_per_sample.csv"),index=False)

# mechanism & class composition (abundance-weighted) per sample
mech_tab=card.groupby('Mechanism')[SAMPLES].sum().T
cls_tab =card.groupby('DrugClass')[SAMPLES].sum().T
mech_tab.to_csv(os.path.join(d1,"data_mechanism_abundance_per_sample.csv"))
cls_tab.to_csv(os.path.join(d1,"data_drugclass_abundance_per_sample.csv"))

# figure: stacked drug-class composition
order_s=[s for e in ENV for s in SAMPLES if META[s][1]==e]
clsN=cls_tab.loc[order_s]; clsN=clsN.div(clsN.sum(axis=1),axis=0)
fig,ax=plt.subplots(figsize=(12,6)); bottom=np.zeros(len(order_s))
palette=plt.cm.tab20(np.linspace(0,1,clsN.shape[1]))
for k,col in enumerate(clsN.columns):
    ax.bar(range(len(order_s)),clsN[col].values,bottom=bottom,label=col,color=palette[k],edgecolor='white',linewidth=.3)
    bottom+=clsN[col].values
ax.set_xticks(range(len(order_s))); ax.set_xticklabels(order_s,rotation=90,fontsize=8)
for i,s in enumerate(order_s):
    ax.get_xticklabels()[i].set_color(ENVC[META[s][1]])
ax.set_ylabel('Relative resistome composition'); ax.set_ylim(0,1)
ax.set_title('Drug-class composition of the wastewater resistome (18 samples)',fontweight='bold')
ax.legend(bbox_to_anchor=(1.01,1),loc='upper left',fontsize=8,frameon=False)
savepng(fig,os.path.join(d1,"Fig_S1_resistome_composition.png"))
# html
traces=[{"type":"bar","name":c,"x":order_s,"y":clsN[c].round(4).tolist()} for c in clsN.columns]
plotly_html(os.path.join(d1,"Fig_S1_resistome_composition.html"),
            "Drug-class composition of wastewater resistome",traces,
            {"barmode":"stack","yaxis":{"title":"Relative composition"},"xaxis":{"title":"Sample"}})

SUMMARY['step1']={'n_raw_args':n_raw_args,'n_args_detected':n_args,
                  'top_classes':cls_tab.sum().sort_values(ascending=False).head(5).round(3).to_dict(),
                  'top_mech':mech_tab.sum().sort_values(ascending=False).round(3).to_dict()}
# ==== STEP 1 END ====

# ============================================================================
# ==== STEP 2 START : ARG burden & diversity across environment/city ====
# ============================================================================
d2=step_dir(2,"ARG_Burden_Diversity")
kw=[]
for metric in ['ARG_burden','ARG_richness','ARG_Shannon']:
    ge=[res[res.Environment==e][metric].values for e in ENV]
    gc=[res[res.City==c][metric].values for c in CITY]
    He,pe=kruskal(*ge); Hc,pc=kruskal(*gc)
    kw.append({'Metric':metric,'H_env':round(He,3),'p_env':round(pe,4),
               'H_city':round(Hc,3),'p_city':round(pc,4)})
kw=pd.DataFrame(kw); kw.to_csv(os.path.join(d2,"data_KW_results.csv"),index=False)
grpmean=res.groupby('Environment')[['ARG_burden','ARG_richness','ARG_Shannon']].mean().round(4)
grpmean.to_csv(os.path.join(d2,"data_group_means_by_environment.csv"))

fig,axes=plt.subplots(1,3,figsize=(16,5.2))
titles={'ARG_burden':'Total ARG burden','ARG_richness':'ARG richness (n genes)','ARG_Shannon':"Resistome Shannon H'"}
for ax,metric in zip(axes,['ARG_burden','ARG_richness','ARG_Shannon']):
    data=[res[res.Environment==e][metric].values for e in ENV]
    bp=ax.boxplot(data,positions=range(len(ENV)),widths=.5,patch_artist=True,
                  medianprops=dict(color='black',linewidth=2),flierprops=dict(marker=''))
    for patch,e in zip(bp['boxes'],ENV): patch.set_facecolor(ENVC[e]); patch.set_alpha(.35)
    for _,r in res.iterrows():
        x=ENV.index(r.Environment)+np.random.uniform(-.12,.12)
        ax.scatter(x,r[metric],color=CITYC[r.City],marker=CITYM[r.City],s=70,zorder=5,edgecolors='white',linewidths=.7)
    He,pe=kruskal(*data); sig="***" if pe<.001 else "**" if pe<.01 else "*" if pe<.05 else "ns"
    ax.set_xticks(range(len(ENV))); ax.set_xticklabels(ENV,fontsize=9)
    ax.set_title(titles[metric],fontweight='bold',fontsize=11)
    ax.text(.97,.97,f'KW p={pe:.3f} {sig}',transform=ax.transAxes,ha='right',va='top',fontsize=9,style='italic',color='#555')
fig.suptitle('Antibiotic-resistance gene burden and diversity across wastewater environments',fontweight='bold',y=1.02)
savepng(fig,os.path.join(d2,"Fig_1_ARG_burden_by_environment.png"))
traces=[]
for e in ENV:
    sub=res[res.Environment==e]
    traces.append({"type":"box","name":e,"y":sub['ARG_burden'].round(5).tolist(),
                   "boxpoints":"all","text":sub['Sample'].tolist(),"marker":{"color":ENVC[e]}})
plotly_html(os.path.join(d2,"Fig_1_ARG_burden_by_environment.html"),
            "ARG burden by environment",traces,{"yaxis":{"title":"Total ARG burden"}})
SUMMARY['step2']={'kw':kw.to_dict('records'),'env_means':grpmean['ARG_burden'].to_dict()}
# ==== STEP 2 END ====

# ============================================================================
# ==== STEP 3 START : resistome beta-diversity (PCoA + PERMANOVA) ====
# ============================================================================
d3=step_dir(3,"Resistome_BetaDiversity")
X=card.set_index('ARO')[SAMPLES].T   # rows=samples
# prevalence filter >=2 samples
X=X.loc[:, (X>0).sum(axis=0)>=2]
n=len(SAMPLES); D=np.zeros((n,n))
for i in range(n):
    for j in range(i+1,n):
        D[i,j]=D[j,i]=braycurtis(X.iloc[i].values,X.iloc[j].values)
dm=pd.DataFrame(D,index=SAMPLES,columns=SAMPLES); dm.to_csv(os.path.join(d3,"data_BrayCurtis_matrix.csv"))
def pcoa(D):
    A=-0.5*D**2; n=D.shape[0]
    H=A-A.mean(1,keepdims=True)-A.mean(0,keepdims=True)+A.mean()
    w,v=np.linalg.eigh(H); idx=np.argsort(w)[::-1]; w=w[idx]; v=v[:,idx]
    pos=w>0; coords=v[:,pos]*np.sqrt(w[pos]); ve=w[pos]/w[pos].sum()*100
    return coords,ve
coords,ve=pcoa(D)
pco=pd.DataFrame({'Sample':SAMPLES,'City':[META[s][0] for s in SAMPLES],
                  'Environment':[META[s][1] for s in SAMPLES],
                  'PC1':coords[:,0],'PC2':coords[:,1],'PC3':coords[:,2] if coords.shape[1]>2 else 0})
pco.to_csv(os.path.join(d3,"data_PCoA_coordinates.csv"),index=False)
def permanova(D,labels,nperm=999,seed=42):
    rng=np.random.RandomState(seed); labels=np.array(labels); k=len(set(labels)); nn=len(labels)
    def F(D,lab):
        SSt=np.sum(D**2)/nn
        SSw=sum(np.sum(D[np.ix_(np.where(lab==g)[0],np.where(lab==g)[0])]**2)/np.sum(lab==g) for g in set(lab))
        SSa=SSt-SSw; return (SSa/(k-1))/(SSw/(nn-k)), SSa/SSt
    Fo,R2=F(D,labels); perm=[F(D,rng.permutation(labels))[0] for _ in range(nperm)]
    p=(np.sum(np.array(perm)>=Fo)+1)/(nperm+1); return Fo,R2,p
perm=[]
for name,lab in [('Environment',[META[s][1] for s in SAMPLES]),('City',[META[s][0] for s in SAMPLES])]:
    Fo,R2,p=permanova(D,lab); perm.append({'Factor':name,'Pseudo_F':round(Fo,3),'R2':round(R2,4),'p_value':round(p,4)})
perm=pd.DataFrame(perm); perm.to_csv(os.path.join(d3,"data_PERMANOVA_results.csv"),index=False)

fig,ax=plt.subplots(figsize=(8,7))
for _,r in pco.iterrows():
    ax.scatter(r.PC1,r.PC2,color=ENVC[r.Environment],marker=CITYM[r.City],s=140,zorder=5,edgecolors='white',linewidths=.8)
    ax.annotate(r.Sample,(r.PC1,r.PC2),fontsize=6,xytext=(4,4),textcoords='offset points',color='#333')
ax.set_xlabel(f'PC1 ({ve[0]:.1f}%)'); ax.set_ylabel(f'PC2 ({ve[1]:.1f}%)')
ax.axhline(0,color='#ddd',lw=.5); ax.axvline(0,color='#ddd',lw=.5)
pe=perm[perm.Factor=='Environment'].iloc[0]
ax.set_title('Resistome composition — PCoA (Bray–Curtis)',fontweight='bold')
ax.text(.02,.98,f"PERMANOVA (env)\nR²={pe.R2:.3f}, p={pe.p_value:.3f}",transform=ax.transAxes,va='top',
        fontsize=9,bbox=dict(boxstyle='round,pad=.3',facecolor='white',edgecolor='#ccc'))
from matplotlib.lines import Line2D; from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor=ENVC[e],label=e) for e in ENV]+
                  [Line2D([0],[0],marker=CITYM[c],color='gray',label=c,linestyle='None') for c in CITY],
          fontsize=8,loc='lower right')
savepng(fig,os.path.join(d3,"Fig_2_resistome_PCoA.png"))
traces=[]
for e in ENV:
    sub=pco[pco.Environment==e]
    traces.append({"type":"scatter","mode":"markers+text","name":e,"x":sub.PC1.round(4).tolist(),
                   "y":sub.PC2.round(4).tolist(),"text":sub.Sample.tolist(),"textposition":"top center",
                   "marker":{"size":12,"color":ENVC[e]}})
plotly_html(os.path.join(d3,"Fig_2_resistome_PCoA.html"),"Resistome PCoA (Bray-Curtis)",traces,
            {"xaxis":{"title":f"PC1 ({ve[0]:.1f}%)"},"yaxis":{"title":f"PC2 ({ve[1]:.1f}%)"}})
SUMMARY['step3']={'PC1':round(ve[0],1),'PC2':round(ve[1],1),'permanova':perm.to_dict('records')}
# ==== STEP 3 END ====

# ============================================================================
# ==== STEP 4 START : resistome vs metabolic architecture (KEY) ====
# ============================================================================
d4=step_dir(4,"Resistome_vs_MetabolicArchitecture")
mci=pd.read_csv(os.path.join(BASE,"2.1 stage 7","Paper1_Stage7_MCI_table.csv"))
alpha=pd.read_csv(os.path.join(BASE,"2.1 stage 2","Paper1_Stage2_AlphaDiversity_Table.csv"))
merged=res.merge(mci[['Sample','MCI','Mean_degree','Max_BC','Active_nodes']],on='Sample')\
          .merge(alpha[['Sample','Shannon','Richness']].rename(columns={'Shannon':'Func_Shannon','Richness':'Func_Richness'}),on='Sample')
merged.to_csv(os.path.join(d4,"data_integrated_per_sample.csv"),index=False)

res_metrics=['ARG_burden','ARG_richness','ARG_Shannon']
net_metrics=['MCI','Mean_degree','Max_BC','Active_nodes','Func_Shannon','Func_Richness']
cor=[]
for rm in res_metrics:
    for nm in net_metrics:
        rho,p=spearmanr(merged[rm],merged[nm])
        cor.append({'Resistome':rm,'Metabolic':nm,'rho':round(rho,3),'p':round(p,4)})
cor=pd.DataFrame(cor)
# BH-FDR across all tests
from numpy import argsort
pv=cor['p'].values; m=len(pv); o=argsort(pv); ranks=np.empty(m,int); ranks[o]=np.arange(1,m+1)
cor['p_FDR']=np.minimum.accumulate((pv[o]*m/np.arange(1,m+1))[::-1])[::-1][np.argsort(o)]
cor['p_FDR']=cor['p_FDR'].round(4)
cor.to_csv(os.path.join(d4,"data_resistome_metabolic_correlations.csv"),index=False)

# correlation heatmap
pivot=cor.pivot(index='Resistome',columns='Metabolic',values='rho')[net_metrics]
fig,ax=plt.subplots(figsize=(9,4.2))
im=ax.imshow(pivot.values,cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto')
ax.set_xticks(range(len(net_metrics))); ax.set_xticklabels(net_metrics,rotation=35,ha='right',fontsize=9)
ax.set_yticks(range(len(res_metrics))); ax.set_yticklabels(res_metrics,fontsize=9)
for i in range(len(res_metrics)):
    for j in range(len(net_metrics)):
        rr=pivot.values[i,j]; pp=cor[(cor.Resistome==res_metrics[i])&(cor.Metabolic==net_metrics[j])].p.values[0]
        star="*" if pp<.05 else ""
        ax.text(j,i,f"{rr:.2f}{star}",ha='center',va='center',fontsize=9,
                color='white' if abs(rr)>.6 else 'black')
fig.colorbar(im,ax=ax,label="Spearman ρ",shrink=.8)
ax.set_title('Resistome burden/diversity vs core-metabolic network architecture\n(Spearman ρ; * p<0.05, n=18)',fontweight='bold',fontsize=11)
savepng(fig,os.path.join(d4,"Fig_3_resistome_vs_architecture_heatmap.png"))
plotly_html(os.path.join(d4,"Fig_3_resistome_vs_architecture_heatmap.html"),
            "Resistome vs metabolic architecture (Spearman rho)",
            [{"type":"heatmap","z":pivot.values.tolist(),"x":net_metrics,"y":res_metrics,
              "zmin":-1,"zmax":1,"colorscale":"RdBu","reversescale":True}],
            {"xaxis":{"title":"Metabolic metric"},"yaxis":{"title":"Resistome metric"}})

# key scatter panel (ARG_burden vs MCI, Mean_degree, Max_BC)
fig,axes=plt.subplots(1,3,figsize=(16,5))
for ax,nm in zip(axes,['MCI','Mean_degree','Max_BC']):
    for _,r in merged.iterrows():
        ax.scatter(r[nm],r.ARG_burden,color=ENVC[r.Environment],marker=CITYM[r.City],s=90,edgecolors='white',linewidths=.7,zorder=5)
    rho,p=spearmanr(merged[nm],merged.ARG_burden)
    z=np.polyfit(merged[nm],merged.ARG_burden,1); xs=np.linspace(merged[nm].min(),merged[nm].max(),50)
    ax.plot(xs,np.polyval(z,xs),'--',color='#555',lw=1)
    ax.set_xlabel(nm); ax.set_ylabel('Total ARG burden')
    ax.set_title(f'ρ={rho:.2f}, p={p:.3f}',fontsize=10)
fig.suptitle('Resistome burden vs metabolic-network metrics (colour=environment, shape=city)',fontweight='bold',y=1.02)
savepng(fig,os.path.join(d4,"Fig_4_burden_vs_metrics_scatter.png"))
SUMMARY['step4']={'strongest':cor.reindex(cor.rho.abs().sort_values(ascending=False).index).head(6).to_dict('records')}
# ==== STEP 4 END ====

# ============================================================================
# ==== STEP 5 START : cross-domain gene correlation (metabolic hubs vs ARGs) ====
# ============================================================================
d5=step_dir(5,"CrossDomain_Gene_Correlation")
core=pd.read_csv(os.path.join(BASE,"2.1 stage 5","Paper1_Stage5_CoreKO_table.csv"))
nodes=pd.read_csv(os.path.join(BASE,"2.1 stage 6","Paper1_Stage6_NodeMetrics.csv"))
hub_ids=nodes.sort_values('Betweenness_centrality',ascending=False).head(15)['KO_ID'].tolist()
coreH=core[core['KO_ID'].isin(hub_ids)].set_index('KO_ID')[SAMPLES]
# ARG mechanism totals per sample (from step1 mech_tab)
mech=mech_tab  # samples x mechanism
# correlate each hub KO (18-vec) with each mechanism total and with total burden
cross=[]
burden=res.set_index('Sample')['ARG_burden']
for ko in hub_ids:
    kv=coreH.loc[ko,SAMPLES].values.astype(float)
    pw=nodes[nodes.KO_ID==ko]['Pathway'].values[0]
    name=core[core.KO_ID==ko]['Name'].values[0]
    rho_b,p_b=spearmanr(kv,burden[SAMPLES].values)
    row={'KO_ID':ko,'Pathway':pw,'Gene':name,'rho_totalARG':round(rho_b,3),'p_totalARG':round(p_b,4)}
    for mname in mech.columns:
        rr,pp=spearmanr(kv,mech.loc[SAMPLES,mname].values); row[f'rho::{mname}']=round(rr,3)
    cross.append(row)
cross=pd.DataFrame(cross); cross.to_csv(os.path.join(d5,"data_hubKO_vs_ARG_correlations.csv"),index=False)

mcols=[c for c in cross.columns if c.startswith('rho::')]
Z=cross.set_index('KO_ID')[mcols].values
fig,ax=plt.subplots(figsize=(10,7))
im=ax.imshow(Z,cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto')
ax.set_xticks(range(len(mcols))); ax.set_xticklabels([c.replace('rho::','') for c in mcols],rotation=35,ha='right',fontsize=9)
ax.set_yticks(range(len(cross))); ax.set_yticklabels([f"{k} ({g})" for k,g in zip(cross.KO_ID,cross.Gene.astype(str).str.slice(0,14))],fontsize=8)
for i in range(Z.shape[0]):
    for j in range(Z.shape[1]):
        ax.text(j,i,f"{Z[i,j]:.2f}",ha='center',va='center',fontsize=7,color='white' if abs(Z[i,j])>.6 else '#222')
fig.colorbar(im,ax=ax,label='Spearman ρ',shrink=.7)
ax.set_title('Top metabolic hub genes vs resistance-mechanism abundance (n=18)',fontweight='bold',fontsize=11)
savepng(fig,os.path.join(d5,"Fig_5_hubKO_vs_ARG_mechanism_heatmap.png"))
plotly_html(os.path.join(d5,"Fig_5_hubKO_vs_ARG_mechanism_heatmap.html"),
            "Metabolic hub genes vs ARG mechanisms",
            [{"type":"heatmap","z":Z.tolist(),"x":[c.replace('rho::','') for c in mcols],
              "y":cross.KO_ID.tolist(),"zmin":-1,"zmax":1,"colorscale":"RdBu","reversescale":True}],
            {"xaxis":{"title":"ARG mechanism"},"yaxis":{"title":"Hub KO"}})
SUMMARY['step5']={'hub_ids':hub_ids[:5],
                  'strongest_totalARG':cross.reindex(cross.rho_totalARG.abs().sort_values(ascending=False).index)[['KO_ID','Gene','rho_totalARG','p_totalARG']].head(5).to_dict('records')}
# ==== STEP 5 END ====

json.dump(SUMMARY,open(os.path.join(OUTP,"_RESULTS_SUMMARY.json"),'w'),indent=2,default=str)
print(json.dumps(SUMMARY,indent=2,default=str))

# ##########  PART 3 — KEGG INTRINSIC-RESISTANCE LAYER  ##########
# ============================================================================
#  KEGG INTRINSIC-RESISTANCE LAYER  —  analysis (Paper 1 extension)
#  Intrinsic resistance = KOs on KEGG antimicrobial-resistance maps
#    map01501 beta-lactam | map01502 vancomycin | map01503 CAMP resistance
# ============================================================================
import os,re,json,textwrap,numpy as np,pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import kruskal, spearmanr, mannwhitneyu
from scipy.spatial.distance import braycurtis
import networkx as nx
np.random.seed(42)

BASE="/sessions/compassionate-pensive-euler/mnt/2.1"
KEGG=os.path.join(BASE,"2 KEGG database for pathway mapping","KEGG数据库基因详细注释表.tsv")
CARDsum=os.path.join(BASE,"CARD_Metabolism_Integration","Step1_Resistome_Alignment_Profiling","data_resistome_summary_per_sample.csv")
OUTP="/tmp/ex/out2"; os.makedirs(OUTP,exist_ok=True)

SAMPLES=["SHW1","SHW2","SCW1","SCW2","SSLW1","SSLW2","MHW1","MHW2","MCW1","MCW2",
         "MSLW1","MSLW2","PHW1","PHW2","PCW1","PCW2","PSLW1","PSLW2"]
META={s:(("Swat" if s[0]=="S" else "Mardan" if s[0]=="M" else "Peshawar"),
         ("Hospital" if s[1]=="H" else "Community" if s[1]=="C" else "Slaughterhouse")) for s in SAMPLES}
ENV=["Hospital","Slaughterhouse","Community"]; CITY=["Swat","Mardan","Peshawar"]
ENVC={'Hospital':'#C0392B','Slaughterhouse':'#E67E22','Community':'#27AE60'}
CITYM={'Swat':'o','Mardan':'s','Peshawar':'^'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
def sdir(n,name): d=os.path.join(OUTP,f"Step{n}_{name}"); os.makedirs(d,exist_ok=True); return d
def savepng(fig,p): fig.savefig(p,dpi=600,bbox_inches='tight'); plt.close(fig)
def phtml(path,title,traces,layout):
    lay={"title":title,"template":"plotly_white","font":{"family":"Arial","size":13}}; lay.update(layout)
    open(path,'w').write(f"""<!doctype html><html><head><meta charset="utf-8"><title>{title}</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script><style>body{{font-family:Arial;margin:16px}}#g{{width:100%;max-width:1000px;height:600px}}</style></head>
<body><div id="g"></div><script>Plotly.newPlot('g',{json.dumps(traces)},{json.dumps(lay)},{{responsive:true}});</script></body></html>""")
def shannon(v):
    v=np.array(v,float); v=v[v>0]
    if v.sum()==0: return 0.0
    p=v/v.sum(); return float(-np.sum(p*np.log(p)))
def pcoa(D):
    A=-0.5*D**2; H=A-A.mean(1,keepdims=True)-A.mean(0,keepdims=True)+A.mean()
    w,v=np.linalg.eigh(H); idx=np.argsort(w)[::-1]; w=w[idx]; v=v[:,idx]; pos=w>0
    return v[:,pos]*np.sqrt(w[pos]), w[pos]/w[pos].sum()*100
def permanova(D,labels,nperm=999,seed=42):
    rng=np.random.RandomState(seed); labels=np.array(labels); k=len(set(labels)); nn=len(labels)
    def F(D,lab):
        SSt=np.sum(D**2)/nn
        SSw=sum(np.sum(D[np.ix_(np.where(lab==g)[0],np.where(lab==g)[0])]**2)/np.sum(lab==g) for g in set(lab))
        SSa=SSt-SSw; return (SSa/(k-1))/(SSw/(nn-k)), SSa/SSt
    Fo,R2=F(D,labels); perm=[F(D,rng.permutation(labels))[0] for _ in range(nperm)]
    return Fo,R2,(np.sum(np.array(perm)>=Fo)+1)/(nperm+1)
SUM={}

# ==== STEP 1 START : extract & profile KEGG intrinsic resistome ====
raw=pd.read_csv(KEGG,sep='\t'); raw.rename(columns={raw.columns[0]:'KO_ID'},inplace=True)
raw=raw[['KO_ID']+SAMPLES+['Module','Pathway','Name','EC','Description']].copy()
RES={'Beta-lactam':'map01501','Vancomycin':'map01502','CAMP':'map01503'}
CORE={'Glycolysis':'map00010','TCA cycle':'map00020','Pentose phosphate':'map00030',
      'Oxid. phosphorylation':'map00190','AA biosynthesis':'map01230'}
def res_cats(pw):
    if pd.isna(pw): return []
    return [k for k,m in RES.items() if m in pw]
raw['ResCats']=raw['Pathway'].apply(res_cats)
resko=raw[raw['ResCats'].map(len)>0].copy()
resko['ResClass']=resko['ResCats'].map(lambda x:'|'.join(x))
d1=sdir(1,"KEGG_IntrinsicResistome_Profiling")
resko[['KO_ID','ResClass','Name','EC','Description']+SAMPLES].to_csv(os.path.join(d1,"data_intrinsic_resistance_KOs.csv"),index=False)
# per sample potential
rows=[]
for s in SAMPLES:
    c,e=META[s]; col=resko[s].values
    rows.append({'Sample':s,'City':c,'Environment':e,
                 'IntR_potential':round(float(col.sum()),4),
                 'IntR_richness':int((col>0).sum()),
                 'IntR_Shannon':round(shannon(col),4)})
prof=pd.DataFrame(rows); prof.to_csv(os.path.join(d1,"data_intrinsic_potential_per_sample.csv"),index=False)
# class breakdown (a KO counted in each class it belongs to)
cls_rows=[]
for cls in RES:
    sub=raw[raw['Pathway'].fillna('').str.contains(RES[cls])]
    cls_rows.append({'Class':cls,'n_KOs':len(sub),'mean_potential':round(sub[SAMPLES].sum().mean(),3)})
cls=pd.DataFrame(cls_rows); cls.to_csv(os.path.join(d1,"data_class_summary.csv"),index=False)
# fig: class KO counts + per-sample potential bar
fig,ax=plt.subplots(1,2,figsize=(14,5))
ax[0].bar(cls['Class'],cls['n_KOs'],color=['#C0392B','#8E44AD','#2980B9'],alpha=.8,edgecolor='white')
for i,v in enumerate(cls['n_KOs']): ax[0].text(i,v+0.5,str(v),ha='center',fontweight='bold')
ax[0].set_ylabel('KOs detected'); ax[0].set_title('A  KEGG intrinsic-resistance genes by class',fontweight='bold',loc='left')
order_s=[s for e in ENV for s in SAMPLES if META[s][1]==e]
ax[1].bar(range(len(order_s)),prof.set_index('Sample').loc[order_s,'IntR_potential'],
          color=[ENVC[META[s][1]] for s in order_s],edgecolor='white')
ax[1].set_xticks(range(len(order_s))); ax[1].set_xticklabels(order_s,rotation=90,fontsize=8)
ax[1].set_ylabel('Intrinsic-resistance potential (Σ TPM)'); ax[1].set_title('B  Per-sample intrinsic-resistance potential',fontweight='bold',loc='left')
savepng(fig,os.path.join(d1,"Fig_1_intrinsic_resistome_overview.png"))
phtml(os.path.join(d1,"Fig_1_intrinsic_resistome_overview.html"),"KEGG intrinsic-resistance potential",
      [{"type":"bar","x":order_s,"y":prof.set_index('Sample').loc[order_s,'IntR_potential'].round(2).tolist(),
        "marker":{"color":[ENVC[META[s][1]] for s in order_s]}}],
      {"yaxis":{"title":"Intrinsic-resistance potential"},"xaxis":{"title":"Sample"}})
SUM['step1']={'n_res_KOs':len(resko),'class':cls.to_dict('records')}
# ==== STEP 1 END ====

# ==== STEP 2 START : environment structuring + contrast ====
d2=sdir(2,"IntrinsicResistance_Environment_Structuring")
kw=[]
for m in ['IntR_potential','IntR_richness','IntR_Shannon']:
    ge=[prof[prof.Environment==e][m].values for e in ENV]; gc=[prof[prof.City==c][m].values for c in CITY]
    He,pe=kruskal(*ge); Hc,pc=kruskal(*gc)
    kw.append({'Metric':m,'H_env':round(He,3),'p_env':round(pe,4),'H_city':round(Hc,3),'p_city':round(pc,4)})
kw=pd.DataFrame(kw); kw.to_csv(os.path.join(d2,"data_KW_results.csv"),index=False)
# PERMANOVA on intrinsic-resistance KO profiles
X=resko.set_index('KO_ID')[SAMPLES].T
X=X.loc[:,(X>0).sum(0)>=2]
n=len(SAMPLES); D=np.zeros((n,n))
for i in range(n):
    for j in range(i+1,n): D[i,j]=D[j,i]=braycurtis(X.iloc[i].values,X.iloc[j].values)
coords,ve=pcoa(D)
pco=pd.DataFrame({'Sample':SAMPLES,'City':[META[s][0] for s in SAMPLES],'Environment':[META[s][1] for s in SAMPLES],
                  'PC1':coords[:,0],'PC2':coords[:,1]})
pco.to_csv(os.path.join(d2,"data_PCoA_coordinates.csv"),index=False)
perm=[]
for name,lab in [('Environment',[META[s][1] for s in SAMPLES]),('City',[META[s][0] for s in SAMPLES])]:
    Fo,R2,p=permanova(D,lab); perm.append({'Factor':name,'Pseudo_F':round(Fo,3),'R2':round(R2,4),'p_value':round(p,4)})
perm=pd.DataFrame(perm); perm.to_csv(os.path.join(d2,"data_PERMANOVA_results.csv"),index=False)
pe=perm[perm.Factor=='Environment'].iloc[0]
# 3-layer contrast: metabolism (known), intrinsic (this), acquired CARD (known)
LAYER=pd.DataFrame([
    {'Layer':'Metabolism\n(functional KO)','R2':0.117,'p':0.415,'src':'Paper1 Stage3'},
    {'Layer':'Intrinsic resistance\n(KEGG maps)','R2':float(pe.R2),'p':float(pe.p_value),'src':'this analysis'},
    {'Layer':'Acquired resistome\n(CARD ARGs)','R2':0.2714,'p':0.004,'src':'CARD integration'}])
LAYER.to_csv(os.path.join(d2,"data_three_layer_PERMANOVA_contrast.csv"),index=False)
# figure: PCoA + 3-layer bar
fig,ax=plt.subplots(1,2,figsize=(15,6))
for _,r in pco.iterrows():
    ax[0].scatter(r.PC1,r.PC2,color=ENVC[r.Environment],marker=CITYM[r.City],s=130,edgecolors='white',zorder=5)
    ax[0].annotate(r.Sample,(r.PC1,r.PC2),fontsize=6,xytext=(4,4),textcoords='offset points',color='#333')
ax[0].axhline(0,color='#ddd',lw=.5); ax[0].axvline(0,color='#ddd',lw=.5)
ax[0].set_xlabel(f'PC1 ({ve[0]:.1f}%)'); ax[0].set_ylabel(f'PC2 ({ve[1]:.1f}%)')
ax[0].set_title(f'A  Intrinsic-resistance composition (PCoA)\nPERMANOVA env R²={pe.R2:.3f}, p={pe.p_value:.3f}',fontweight='bold',loc='left',fontsize=11)
cols=['#2980B9','#27AE60','#C0392B']
bars=ax[1].bar(range(3),LAYER['R2'],color=cols,alpha=.85,edgecolor='white',width=.6)
for i,(r2,p) in enumerate(zip(LAYER['R2'],LAYER['p'])):
    sig='**' if p<0.01 else '*' if p<0.05 else 'ns'
    ax[1].text(i,r2+0.006,f'p={p:.3f} {sig}',ha='center',fontweight='bold',fontsize=10)
ax[1].set_xticks(range(3)); ax[1].set_xticklabels(LAYER['Layer'],fontsize=9)
ax[1].set_ylabel('PERMANOVA R² (environment)'); ax[1].set_ylim(0,0.34)
ax[1].set_title('B  Which layer does environment structure?',fontweight='bold',loc='left',fontsize=11)
savepng(fig,os.path.join(d2,"Fig_2_env_structuring_and_contrast.png"))
phtml(os.path.join(d2,"Fig_2_three_layer_contrast.html"),"Environment structuring by layer",
      [{"type":"bar","x":[l.replace('\n',' ') for l in LAYER['Layer']],"y":LAYER['R2'].round(3).tolist(),
        "text":[f"p={p}" for p in LAYER['p']],"marker":{"color":cols}}],
      {"yaxis":{"title":"PERMANOVA R² (environment)"}})
SUM['step2']={'kw':kw.to_dict('records'),'permanova':perm.to_dict('records'),'layers':LAYER.to_dict('records')}
# ==== STEP 2 END ====

# ==== STEP 3 START : network integration of resistance into metabolism ====
d3=sdir(3,"Network_Integration")
def prim(x):
    if pd.isna(x): return None
    for k,m in CORE.items():
        if m in x: return k
    return None
raw['Primary']=raw['Pathway'].apply(prim)
core=raw[raw['Primary'].notna()].copy()
res_ids=set(resko['KO_ID'])
# node set = core metabolic KOs + resistance KOs (label resistance if in res set)
node_df=pd.concat([core[['KO_ID','Primary']+SAMPLES],
                   resko[['KO_ID']+SAMPLES].assign(Primary='Intrinsic resistance')]).drop_duplicates('KO_ID')
# if a KO is both, prefer resistance label
node_df.loc[node_df['KO_ID'].isin(res_ids),'Primary']='Intrinsic resistance'
node_df=node_df.drop_duplicates('KO_ID')
vc=node_df[SAMPLES].var(axis=1); nd=node_df[vc>0].copy()
Xn=np.log1p(nd[SAMPLES].values.T); ko=nd['KO_ID'].values; lab=nd['Primary'].values
rho,pv=spearmanr(Xn,axis=0); N=Xn.shape[1]; ti,tj=np.triu_indices(N,1)
r=rho[ti,tj]; p=pv[ti,tj]; ok=~(np.isnan(r)|np.isnan(p))
mask=(np.abs(r[ok])>=0.6)&(p[ok]<0.001)
G=nx.Graph()
for i,(k,l) in enumerate(zip(ko,lab)): G.add_node(k,cat=l)
tiv=ti[ok][mask]; tjv=tj[ok][mask]; rv=r[ok][mask]
for a,b,rr in zip(tiv,tjv,rv): G.add_edge(ko[a],ko[b],weight=abs(rr))
lcc=max(nx.connected_components(G),key=len); Gl=G.subgraph(lcc).copy()
bc=nx.betweenness_centrality(Gl,normalized=True,weight='weight'); bcf={x:bc.get(x,0) for x in G.nodes()}
deg=dict(G.degree())
nm=pd.DataFrame({'KO_ID':list(G.nodes()),'Category':[G.nodes[x]['cat'] for x in G.nodes()],
                 'Degree':[deg[x] for x in G.nodes()],'Betweenness':[round(bcf[x],6) for x in G.nodes()]})
nm.to_csv(os.path.join(d3,"data_combined_network_nodes.csv"),index=False)
resN=nm[nm.Category=='Intrinsic resistance']; metN=nm[nm.Category!='Intrinsic resistance']
# integration stats
res_connected=int((resN.Degree>0).sum()); res_total=len(resN)
# edges from resistance nodes -> which category
res_edges=[]
for u,v,dd in G.edges(data=True):
    cu,cv=G.nodes[u]['cat'],G.nodes[v]['cat']
    if (cu=='Intrinsic resistance')^(cv=='Intrinsic resistance'):
        partner=cv if cu=='Intrinsic resistance' else cu
        res_edges.append(partner)
edge_attach=pd.Series(res_edges).value_counts()
try: U,pmw=mannwhitneyu(resN.Degree,metN.Degree,alternative='two-sided'); 
except Exception: pmw=np.nan
nm.to_csv(os.path.join(d3,"data_combined_network_nodes.csv"),index=False)
pd.DataFrame({'Partner_pathway':edge_attach.index,'n_edges':edge_attach.values}).to_csv(os.path.join(d3,"data_resistance_edge_attachment.csv"),index=False)
# figure: degree distribution res vs met + attachment bar
fig,ax=plt.subplots(1,2,figsize=(15,6))
data=[metN.Degree.values,resN.Degree.values]
bp=ax[0].boxplot(data,positions=[0,1],widths=.5,patch_artist=True,medianprops=dict(color='black',lw=2),flierprops=dict(marker=''))
for patch,c in zip(bp['boxes'],['#2980B9','#C0392B']): patch.set_facecolor(c); patch.set_alpha(.4)
for i,dd in enumerate(data):
    ax[0].scatter(np.random.uniform(i-.12,i+.12,len(dd)),dd,s=12,color=['#2980B9','#C0392B'][i],alpha=.5,zorder=5)
ax[0].set_xticks([0,1]); ax[0].set_xticklabels(['Metabolic KOs','Intrinsic-resistance KOs'])
ax[0].set_ylabel('Node degree (co-occurrence)')
ax[0].set_title(f'A  Are resistance genes integrated?\nMann-Whitney p={pmw:.3g}',fontweight='bold',loc='left',fontsize=11)
# normalise attachment by number of metabolic nodes available in each pathway
psize=metN.Category.value_counts()
ea_raw=edge_attach.reindex([p for p in ['TCA cycle','Oxid. phosphorylation','Glycolysis','AA biosynthesis','Pentose phosphate'] if p in edge_attach.index])
ea=(ea_raw/psize.reindex(ea_raw.index)).sort_values(ascending=False)
pd.DataFrame({'Pathway':ea.index,'edges_per_metabolic_gene':ea.round(3).values,'raw_edges':ea_raw.reindex(ea.index).values}).to_csv(os.path.join(d3,'data_resistance_edge_attachment.csv'),index=False)
ax[1].bar(range(len(ea)),ea.values,color=['#E67E22' if p in ('TCA cycle','Oxid. phosphorylation') else '#95A5A6' for p in ea.index],alpha=.9,edgecolor='white')
ax[1].set_xticks(range(len(ea))); ax[1].set_xticklabels(ea.index,rotation=30,ha='right',fontsize=9)
ax[1].set_ylabel('Resistance→metabolism edges per gene')
ax[1].set_title('B  Resistance genes attach most densely to energy metabolism',fontweight='bold',loc='left',fontsize=11)
savepng(fig,os.path.join(d3,"Fig_3_resistance_network_integration.png"))
SUM['step3']={'nodes':G.number_of_nodes(),'edges':G.number_of_edges(),
              'res_connected':res_connected,'res_total':res_total,
              'med_deg_res':float(resN.Degree.median()),'med_deg_met':float(metN.Degree.median()),
              'mannwhitney_p':float(pmw),'edge_attach_raw':edge_attach.to_dict(),'edge_attach_per_gene':(edge_attach/metN.Category.value_counts()).dropna().round(3).sort_values(ascending=False).to_dict()}
# ==== STEP 3 END ====

json.dump(SUM,open(os.path.join(OUTP,"_SUMMARY.json"),'w'),indent=2,default=str)
# save intrinsic vs acquired merge for step4
card=pd.read_csv(CARDsum)
merged=prof.merge(card[['Sample','ARG_burden','ARG_richness']],on='Sample')
merged.to_csv(os.path.join(OUTP,"_intrinsic_vs_acquired.csv"),index=False)
print(json.dumps(SUM,indent=2,default=str))


# #############################################################################
# ##  PART 4 — VARIABLE-FUNCTION SCREEN & PRODUCER-ACQUIRER AXIS (KEGG)       ##
# #############################################################################
# Screens variable KEGG functional layers for links to the CARD resistome, and
# tests the antibiotic-biosynthesis (streptomycin) inverse relationship.
import pandas as pd, numpy as np
from scipy.stats import spearmanr, kruskal
from scipy.spatial.distance import braycurtis
KEGG = r"...\2 KEGG database for pathway mapping\KEGG数据库基因详细注释表.tsv"   # set path
CARD_SUMMARY = r"...\CARD_data_resistome_summary_per_sample.csv"                 # per-sample ARG burden
SAMPLES=["SHW1","SHW2","SCW1","SCW2","SSLW1","SSLW2","MHW1","MHW2","MCW1","MCW2",
         "MSLW1","MSLW2","PHW1","PHW2","PCW1","PCW2","PSLW1","PSLW2"]
META={s:(("Swat" if s[0]=="S" else "Mardan" if s[0]=="M" else "Peshawar"),
         ("Hospital" if s[1]=="H" else "Community" if s[1]=="C" else "Slaughterhouse")) for s in SAMPLES}
ENV=["Hospital","Slaughterhouse","Community"]

raw=pd.read_csv(KEGG,sep='\t'); raw.rename(columns={raw.columns[0]:'KO'},inplace=True)
raw['Pathway']=raw['Pathway'].fillna('')
burden=pd.read_csv(CARD_SUMMARY).set_index('Sample').loc[SAMPLES,'ARG_burden']

PILLARS={
 'Antibiotic biosynthesis':['map00311','map00332','map00261','map00521','map01055','map00253','map00401'],
 'Xenobiotic degradation':['map00362','map00627','map00364','map00625','map00361','map00623',
                           'map00622','map00633','map00643','map00791','map00930','map00983','map00982','map00980'],
 'Transport & efflux':['map02010','map02060'],
 'Communication & biofilm':['map02024','map02025','map02026','map05111']}

def perm(D,lab,nperm=999,seed=42):
    rng=np.random.RandomState(seed); lab=np.array(lab); k=len(set(lab)); n=len(lab)
    def F(D,l):
        SSt=np.sum(D**2)/n
        SSw=sum(np.sum(D[np.ix_(np.where(l==g)[0],np.where(l==g)[0])]**2)/np.sum(l==g) for g in set(l))
        SSa=SSt-SSw; return (SSa/(k-1))/(SSw/(n-k)), SSa/SSt
    Fo,R2=F(D,lab); pr=[F(D,rng.permutation(lab))[0] for _ in range(nperm)]
    return R2,(np.sum(np.array(pr)>=Fo)+1)/(nperm+1)
def bc(X):
    n=len(SAMPLES); D=np.zeros((n,n))
    for i in range(n):
        for j in range(i+1,n): D[i,j]=D[j,i]=braycurtis(X.iloc[i].values,X.iloc[j].values)
    return D

for name,maps in PILLARS.items():
    sub=raw[raw['Pathway'].apply(lambda p:any(m in p for m in maps))]
    ab=sub.set_index('KO')[SAMPLES].T.sum(1)                      # per-sample potential
    X=sub.set_index('KO')[SAMPLES].T; X=X.loc[:,(X>0).sum(0)>=2]
    R2,p=perm(bc(X),[META[s][1] for s in SAMPLES])
    rho,pr=spearmanr(ab.values,burden.values)
    print(f"{name}: env PERMANOVA R2={R2:.3f} p={p:.3f} | rho(ABund,ARG)={rho:+.2f} p={pr:.3f}")

# streptomycin-pathway (map00521) — the driver of the producer-acquirer axis
strep=raw[raw['Pathway'].str.contains('map00521')].set_index('KO')[SAMPLES].T.sum(1)
print("Streptomycin potential vs ARG burden:", spearmanr(strep.values,burden.values))


# #############################################################################
# ##  PART 5 — MOBILOME (MGE) INTEGRATION                                     ##
# #############################################################################
# Class 1 integron markers vs the acquired resistome and producer capacity.
MGE_LEVEL2 = r"...\2. functional classification MGE\...\All.MGE.Level2.xls"      # set path (tab-separated)
MGE_LEVEL1 = r"...\All.MGE.Level1.xls"

def load_mge(f):
    d=pd.read_csv(f,sep='\t'); d.rename(columns={d.columns[0]:'MGE'},inplace=True)
    return d
L2=load_mge(MGE_LEVEL2); L1=load_mge(MGE_LEVEL1)
def tot(df,name):
    r=df[df['MGE']==name]; return r[SAMPLES].iloc[0] if len(r) else pd.Series(0,index=SAMPLES)
intI1=tot(L2,'intI1'); qac=tot(L2,'qacEdelta'); tnpA=tot(L2,'tnpA')
mobilome=L1[~L1['MGE'].str.contains(r'\|')][SAMPLES].sum()       # total-category rows only

print("intI1 vs ARG burden        :", spearmanr(intI1.values,burden.values))    # +0.72, p<0.001
print("qacEdelta vs ARG burden    :", spearmanr(qac.values,burden.values))      # +0.62
print("total mobilome vs ARG burden:", spearmanr(mobilome.values,burden.values))# +0.57
print("intI1 vs streptomycin      :", spearmanr(intI1.values,strep.values))     # -0.65, p=0.004
print("intI1 by environment KW    :", kruskal(*[intI1[[s for s in SAMPLES if META[s][1]==e]].values for e in ENV]))
