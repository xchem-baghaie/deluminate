from typing import List, Tuple, Union

from configs.visualizations import EMBEDDINGS, DIMENSIONALITY_REDUCTION_METHODS, TMAP_L
from transforms.embeddings import parallel_map_get_embeddings


def calculate_nodes_and_edges(
    smiles_list: List[str],
    embedding: str,
    dimensionality_reduction: str,
    intermediate_pca_reduction: bool = False,
    intermediate_pca_n_dimensions: int = 50,
) -> Tuple[List, List, Union[List, None], Union[List, None]]:
    """
    Args:
    -smiles_list (List of SMILES strings)
    -embedding (the type of embedding to use for dimensionality reduction). One of:
        -ECFP2
        -ECFP4
        -ECFP6
        -MACCS
        -Daylight
        -Morgan
        -PubChem
        -MHFP6
        -MAP4
        -MXFP
        -MQN
    -dimensionality_reduction (the dimensionality reduction method). One of:
        -TMAP
        -UMAP
        -t-SNE
        -PCA
        -SOM
    -intermediate_pca_reduction: bool, whether to first reduce the dimensionality of the fingerprint embeddings using PCA
    =intermediate_pca_n_dimensions: int, number of dimensions to reduce to, if intermediate_pca_reduction = True

    Returns:
    -List of x coordinates
    -List of y coordinates
    -List of start nodes (None if dimensionality reduction is not TMAP. Otherwise, list of type: tmap.VectorUInt)
    -List of to nodes (None if dimensionality reduction is not TMAP. Otherwise, list of type: tmap.VectorUInt)
    """
    assert embedding in EMBEDDINGS
    assert dimensionality_reduction in DIMENSIONALITY_REDUCTION_METHODS

    fps = parallel_map_get_embeddings(smiles_list, embedding)

    if intermediate_pca_reduction:
        from sklearn.decomposition import PCA

        pca = PCA(n_components=intermediate_pca_n_dimensions)
        fps = pca.fit_transform(fps)

    if dimensionality_reduction == "TMAP":
        import tmap as tm

        d = len(fps[0])  # number of hash functions
        l = TMAP_L[embedding]  # number of prefix trees
        print(f"embedding={embedding}, d={d}, l={l}")
        tmap_fps = [tm.VectorUint(x) for x in fps]
        lf = tm.LSHForest(d, l)
        lf.batch_add(tmap_fps)
        lf.index()
        cfg = tm.LayoutConfiguration()
        cfg.k = 100
        cfg.sl_repeats = 2
        cfg.mmm_repeats = 2
        cfg.node_size = 0.015
        cfg.fme_randomize = False
        x, y, s, t, _ = tm.layout_from_lsh_forest(lf, config=cfg)
        return list(x), list(y), s, t
    elif dimensionality_reduction == "UMAP":
        import umap

        reducer = umap.UMAP()
        coords = reducer.fit_transform(fps)
        return coords[:, 0], coords[:, 1], None, None
    elif dimensionality_reduction == "t-SNE":

        from sklearn.manifold import TSNE

        reducer = TSNE(n_components=2)
        coords = reducer.fit_transform(fps)

        return coords[:, 0], coords[:, 1], None, None
    elif dimensionality_reduction == "PCA":
        from sklearn.decomposition import PCA

        reducer = PCA(n_components=2)
        coords = reducer.fit_transform(fps)
        return coords[:, 0], coords[:, 1], None, None
    elif dimensionality_reduction == "SOM":
        from minisom import MiniSom

        reducer = MiniSom(
            x=10,
            y=10,
            input_len=len(fps[0]),
            sigma=0.3,
            learning_rate=0.5,
            random_seed=42,
        )
        reducer.train_random(fps, num_iteration=20000)

        x = []
        y = []
        for fp in fps:
            x_val, y_val = reducer.winner(fp)
            x.append(x_val)
            y.append(y_val)

        return x, y, None, None
    else:
        raise Exception(
            f"Unexpected dimensionality reduction method, {dimensionality_reduction}, was used"
        )


if __name__ == "__main__":
    import pandas as pd
    from time import time

    IN_FILE = '/mnt/p/discovery chemistry/people/ryan/shuttle between computers/trxb2_combined_preds_top_10k.csv'
    EMBEDDING ='ECFP2'
    DIMENSIONALITY_REDUCTION = 'UMAP'
    OUT_FILE = f'/mnt/p/discovery chemistry/people/ryan/shuttle between computers/trxb2_combined_preds_top_10k_{EMBEDDING}_{DIMENSIONALITY_REDUCTION}.csv'
    
    df = pd.read_csv(IN_FILE).fillna('')

    n = len(df)
    smiles_ls = [x.replace('"', "").replace("'", "") for x in list(df["SMILES"])]

    tick = time()
    x, y, s, t = calculate_nodes_and_edges(
        smiles_ls,
        embedding=EMBEDDING,
        dimensionality_reduction=DIMENSIONALITY_REDUCTION,
    )
    tock = time()
    took = int((tock - tick) / 60)
    print(f'Took {took} minutes for {n} compounds')
    df[f'{EMBEDDING}_X'] = list(x)
    df[f'{EMBEDDING}_Y'] = list(y)
    if EMBEDDING == 'TMAP':
        df[f'{EMBEDDING}_S'] = [str(x) for x in s]
        df[f'{EMBEDDING}_T'] = [str(x) for x in t]
    df.to_csv(OUT_FILE, index=False)

