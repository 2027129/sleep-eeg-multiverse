import numpy as np, pandas as pd, mne, warnings
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
warnings.filterwarnings('ignore')
plt.rcParams['font.family'] = 'Liberation Serif'; plt.rcParams['font.size'] = 8

ROI = {'frontal': ['Fp1','Fpz','Fp2','AF7','AF3','AF4','AF8','F7','F5','F3','F1','Fz','F2','F4','F6','F8'],
       'centro_temporal': ['FT7','FC5','FC3','FC1','FC2','FC4','FC6','FT8','T7','C5','C3','C1','Cz','C2','C4','C6','T8',
                           'TP7','CP5','CP3','CP1','CPz','CP2','CP4','CP6','TP8','TP9','TP10'],
       'parieto_occipital': ['P7','P5','P3','P1','Pz','P2','P4','P6','P8','PO7','PO3','POz','PO4','PO8','O1','Oz','O2']}
B = pd.read_csv('results/coreB_epoch2_specification_curve.csv')
t = B[(B.band == 'theta') & (B.outcome == 'mean')]
med = t.groupby('roi').dz.median().to_dict(); sig = (100 * t.groupby('roi').sig.mean()).round(0).to_dict()

m = mne.channels.make_standard_montage('standard_1020')
pos3 = m.get_positions()['ch_pos']
def proj(p):
    x, y, z = p; r = np.sqrt(x*x + y*y + z*z); x, y, z = x/r, y/r, z/r
    theta = np.arccos(z); phi = np.arctan2(y, x)
    return theta * np.cos(phi), theta * np.sin(phi)
xy = {ch: proj(pos3[ch]) for r in ROI.values() for ch in r}
rmax = max(np.hypot(*v) for v in xy.values())
cmap = LinearSegmentedColormap.from_list('theta', ['#F4F7F8', '#9CBFC8', '#4E7C8A', '#2E4F5A'])
norm = Normalize(vmin=0.0, vmax=0.5)

fig, ax = plt.subplots(figsize=(3.35, 3.0), dpi=300)
ax.set_aspect('equal'); ax.axis('off')
R = rmax * 1.08
ax.add_patch(Circle((0, 0), R, fc='white', ec='black', lw=0.9, zorder=1))
ax.add_patch(Polygon([(-0.13*R, R*0.985), (0, R*1.12), (0.13*R, R*0.985)], closed=True, fc='white', ec='black', lw=0.9, zorder=0))
for sgn in (-1, 1):
    ax.add_patch(Polygon([(sgn*R*0.99, 0.18*R), (sgn*R*1.09, 0.12*R), (sgn*R*1.09, -0.12*R), (sgn*R*0.99, -0.18*R)], closed=True, fc='white', ec='black', lw=0.9, zorder=0))
for roi, chs in ROI.items():
    c = cmap(norm(med[roi]))
    for ch in chs:
        x, y = xy[ch]
        ax.add_patch(Circle((x, y), 0.055, fc=c, ec='black', lw=0.4, zorder=3))
lab = {'frontal': ('Frontal', (1.55*R, 0.72*R), (0.55*R, 0.72*R)),
       'centro_temporal': ('Centro-temporal', (1.55*R, 0.0), (1.02*R, 0.0)),
       'parieto_occipital': ('Parieto-occipital', (1.55*R, -0.72*R), (0.55*R, -0.72*R))}
for roi, (name, (tx, ty), (lx, ly)) in lab.items():
    ax.plot([lx, tx - 0.05*R], [ly, ty], color='grey', lw=0.6, zorder=2)
    ax.text(tx, ty, f"{name}\n$d_z$ = {med[roi]:.2f}\n{sig[roi]:.0f}% significant", ha='left', va='center', fontsize=7.2, zorder=4)
ax.set_xlim(-1.2*R, 2.75*R); ax.set_ylim(-1.2*R, 1.25*R)
ax.set_title('B. Theta effect by region, fixed sample', fontsize=7.4, fontweight='bold', loc='left', x=0.0)
sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
cb = fig.colorbar(sm, ax=ax, orientation='horizontal', fraction=0.045, pad=0.02, shrink=0.5, anchor=(0.0, 1.0))
cb.set_label("Median Cohen's $d_z$", fontsize=7.5); cb.ax.tick_params(labelsize=7)
fig.savefig('figures/headmap.png', bbox_inches='tight', facecolor='white')
print('saved')
