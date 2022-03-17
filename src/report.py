import sys

path = sys.argv[1]

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from stylia import create_figure
from stylia.colors.colors import Palette

df = pd.read_csv(os.path.join(path, "rdkit-descriptors.csv"))

fig, axs = create_figure(nrows=4, ncols=4, support="paper", area_proportion=1.5)
properties = ["MolWt", "MolLogP", "TPSA", "FractionCSP3", "qed", "NumRotatableBonds", "NumHAcceptors", "NumHDonors"]

fig.suptitle(path.split("/")[-1].upper())

palette = Palette(shuffle=False)

limits = {
    "MolWt": (100, 400),
    "MolLogP": (-3, 5),
    "TPSA": (0, 100),
    "FractionCSP3": (0,1),
    "qed": (0,1),
    "NumRotatableBonds": (0,5),
    "NumHAcceptors": (0,5),
    "NumHDonors": (0,5)
}

import collections
def barplot(ax, vals, xlim, color):
    counts = collections.defaultdict(int)
    for v in vals:
        counts[int(v)] += 1
    print(counts)
    x = []
    y = []
    for k,v in counts.items():
        x += [k]
        y += [v]
    ax.scatter(x, y, color=color)
    for x_, y_ in zip(x,y):
        ax.plot([x_, x_], [0, y_], color=color)

for p in properties:
    ax = axs.next()
    xlim = limits[p]
    
    vals = list(df[p])
    if "Num" in p:
        barplot(ax, vals, xlim, palette.next())
    else:
        vals = np.clip(vals, xlim[0], xlim[1])
        ax.hist(vals, 20, cumulative=False, color=palette.next())
    p50 = np.median(vals)
    p25 = np.percentile(vals, 25)
    p75 = np.percentile(vals, 75)
    ax.set_title("Med: {0} IQR: ({1}, {2})".format(np.round(p50, 1), np.round(p25, 1), np.round(p75, 1)))
    ax.set_xlabel(p)
    ax.set_ylabel("Counts")
    margin = xlim[1] - xlim[0]
    ax.set_xlim(xlim[0]-margin*0.05, xlim[1]+margin*0.05)


embedding = "grover-embedding"

def embeddings_plots(embedding, axs):

    with open(os.path.join(path, embedding+"_Xt.npy"), "rb") as f:
        Xt = np.load(f)

    with open(os.path.join(path, embedding+"_Xp.npy"), "rb") as f:
        Xp = np.load(f)

    with open(os.path.join(path, embedding+"_distances.npy"), "rb") as f:
        distances = np.load(f)

    palette = Palette(shuffle=False)

    ax = axs.next()
    ax.hist(distances[:,1], 100, histtype="step", color=palette.next(), cumulative=True, label="Neigh 1")
    ax.hist(distances[:,2], 100, histtype="step", color=palette.next(), cumulative=True, label="Neigh 2")
    ax.hist(distances[:,3], 100, histtype="step", color=palette.next(), cumulative=True, label="Neigh 3")
    if "grover" in embedding:
        ax.set_xlim(0, 0.2)
    else:
        ax.set_xlim(0, 0.6)
    ax.set_xlabel("Cosine distance")
    ax.set_ylabel("Counts")
    ax.set_title(embedding)

    ax = axs.next()
    ax.scatter(Xt[:,0], Xt[:,1], color=palette.next(), s=3)
    ax.set_title("UMAP")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")

    ax = axs.next()
    ax.scatter(Xp[:,0], Xp[:,1], color=palette.next(), s=3)
    ax.set_title("PCA")
    ax.set_xlabel("PC 1")
    ax.set_ylabel("PC 2")

    ax = axs.next()
    ax.scatter(Xp[:,2], Xp[:,3], color=palette.next(), s=3)
    ax.set_title("PCA")
    ax.set_xlabel("PC 3")
    ax.set_ylabel("PC 4")

embeddings_plots("grover-embedding", axs)
embeddings_plots("cc-signaturizer", axs)

plt.tight_layout()
plt.savefig(os.path.join(path, "plot.png"), dpi=600)


from rdkit import Chem
from rdkit.Chem import Draw
import os

mols = [Chem.MolFromSmiles(x) for x in list(df["input"])]
import random
mols = random.sample(mols, 100)
core = Chem.MolFromSmiles( 'c1ncc2nc[nH]c2n1' )
img = Draw.MolsToGridImage( mols, molsPerRow=10, useSVG=False)
img.save(os.path.join(path, "molecules.png"))