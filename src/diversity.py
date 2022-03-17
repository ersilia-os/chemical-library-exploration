import h5py
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
from umap import UMAP
import os
import numpy as np
import sys

path = sys.argv[1]

def distances_and_map(embedding):

    with h5py.File(os.path.join(path, "{0}.h5".format(embedding)), "r") as f:
        X = f["Values"][:]

    nn = NearestNeighbors(n_neighbors=4, metric="cosine")
    nn.fit(X)
    distances, _ = nn.kneighbors(X)

    reducer = UMAP(densmap=True, n_neighbors=15,
                metric="euclidean", low_memory=False,
                repulsion_strength=3, negative_sample_rate=50,
                verbose=True)

    reducer.fit(X)
    Xt = reducer.transform(X)

    pca = PCA(n_components=4, whiten=True)
    pca.fit(X)
    Xp = pca.transform(X)
    
    with open(os.path.join(path, embedding+"_distances.npy"), "wb") as f:
        np.save(f, distances)

    with open(os.path.join(path, embedding+"_Xt.npy"), "wb") as f:
        np.save(f, Xt)

    with open(os.path.join(path, embedding+"_Xp.npy"), "wb") as f:
        np.save(f, Xp)

distances_and_map("cc-signaturizer")
distances_and_map("grover-embedding")