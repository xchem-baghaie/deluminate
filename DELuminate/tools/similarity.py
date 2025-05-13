from itertools import chain
from typing import Optional, Union
import math

import numpy as np
from mhfp.encoder import MHFPEncoder
from rdkit import DataStructs

def adaptive_batch_size_matrix(feat_array_size: int) -> int:
    """
    The batch size is computed as follows:
    - If array size is less than 1000, the batch size is 100 or half the size of the
        feature array, whichever is lower.
    - If array size is greater than 1000, the batch size starts at 256 and increases
        exponentially to a maximum of 1024 at array size 1e4.
    """
    if feat_array_size < int(1e3):
        return min(feat_array_size // 2, 100)
    else:
        return min(int(2 ** (2 * (math.log10(feat_array_size) + 1))), 1024)


def adaptive_batch_size_triu(feat_array_size: int) -> int:
    if feat_array_size < int(1e3):
        return min(feat_array_size // 2, 100)
    elif int(1e3) <= feat_array_size < int(1e4):
        return 256
    else:
        return 512


def tanimoto_complementary(
    fingerprint_x: Union[DataStructs.cDataStructs.ExplicitBitVect, np.ndarray],
    fingerprint_y: Union[DataStructs.cDataStructs.ExplicitBitVect, np.ndarray],
) -> float:
    """
    Compute the Tanimoto complementary between two fingerprints. Suitable for "ECFP2", "ECFP4",
    "ECFP6" and "Morgan" fingerprints.
    """
    if isinstance(fingerprint_x, np.ndarray):  # convert to bitvect
        fingerprint_x = DataStructs.cDataStructs.CreateFromBitString(
            "".join(fingerprint_x.astype(str))
        )
    if isinstance(fingerprint_y, np.ndarray):  # convert to bitvect
        fingerprint_y = DataStructs.cDataStructs.CreateFromBitString(
            "".join(fingerprint_y.astype(str))
        )
    return 1 - DataStructs.TanimotoSimilarity(fingerprint_x, fingerprint_y)


def mhfp_distance(
    fingerprint_x: Union[DataStructs.cDataStructs.ExplicitBitVect, np.ndarray],
    fingerprint_y: Union[DataStructs.cDataStructs.ExplicitBitVect, np.ndarray],
) -> float:
    """
    Compute the MHFP distance between two fingerprints. Suitable for "MHFP6" and "MAP4" fingerprints.
    """
    return MHFPEncoder.distance(fingerprint_x, fingerprint_y)


def parallel_mhfp_complementary_triu(
    fingerprint_array: np.ndarray, batch_size: int = None
) -> np.ndarray:
    """
    Optimized implementation, leveraging parallelism and symmetry.
    Returns the flattened upper triangle of a mhfp similarity matrix.
    """
    from metaflow import parallel_map

    total_len = len(fingerprint_array)
    triu_len = int(total_len * (total_len - 1) / 2)
    batch_size = adaptive_batch_size_triu(triu_len)

    # Contains N(N-1)/2 off-diagonal mcs ratios when flattened
    triu_feat = [
        [fingerprint_array[i], fingerprint_array[j]]
        for i, j in zip(*np.triu_indices(total_len, k=1))
    ]
    mhfp_complementary_triu_nested = parallel_map(
        lambda rolling_idx: [
            1.0 - mhfp_distance(fingerprint_x=feat_A, fingerprint_y=feat_B)
            for feat_A, feat_B in triu_feat[
                rolling_idx : min(rolling_idx + batch_size, triu_len)
            ]
        ],
        range(0, triu_len, batch_size),
    )

    return np.array(
        list(chain.from_iterable(mhfp_complementary_triu_nested)), dtype=float
    )


def parallel_mhfp_complementary_matrix(
    seed_fingerprint_array: np.ndarray,
    fingerprint_array: np.ndarray,
    batch_size: int = None,
) -> np.ndarray:
    """
    Optimized implementation, assuming M (# of rows) << (N - M) (# of cols) to leverage parallelism.
    Suitable for "MHFP6" and "MAP4" fingerprints represented as integer numpy arrays
    with dim (*, 2048).
    """
    from metaflow import parallel_map

    num_cols = len(fingerprint_array)
    batch_size = adaptive_batch_size_matrix(num_cols)

    # Mx(N-M) matrix of MHFP complementary similarities, built in batches of columns (N-M)
    mhfp_complementary_nested = parallel_map(
        lambda rolling_idx: [
            [
                1.0
                - mhfp_distance(
                    fingerprint_x=seed_finger_print, fingerprint_y=finger_print
                )
                for seed_finger_print in seed_fingerprint_array
            ]
            for finger_print in fingerprint_array[
                rolling_idx : min(rolling_idx + batch_size, num_cols)
            ]
        ],
        range(0, num_cols, batch_size),
    )
    return np.array(list(chain.from_iterable(mhfp_complementary_nested)), dtype=float).T


def parallel_sklearn_pairwise_distances(
    feature_array: np.ndarray,
    metric_type: str,
    y_array: Optional[np.ndarray] = None,
) -> np.ndarray:
    """

    Parallel map for sklearn metrics. y_array is an optional argumen, if not passed metrics compared in feature array only

    Args:
        feature_array (np.ndarray): Feature array used
        metric_type (str): Metric type, only manhattan and euclidean supported
        y_array (Optional[np.ndarray], optional): second array to compare to. Defaults to None.
        axis (int) = 1: axis to compare on
        metric_kwargs = None: kwargs to pass to metric

    Raises:
        ValueError: If metric type not supported

    Returns:
        np.array: Similarity matrix
    """
    from tools.parallel_map import parallel_map_numpy_array

    metric_type = metric_type.lower()

    if metric_type == "manhattan":
        from sklearn.metrics.pairwise import manhattan_distances

        call_fn = manhattan_distances
    elif metric_type == "euclidean":
        from sklearn.metrics.pairwise import euclidean_distances

        call_fn = euclidean_distances
    else:
        raise ValueError(
            f"Metric type {metric_type} not supported, only euclidean and manhattan "
        )

    pairwise_nested = parallel_map_numpy_array(
        array=feature_array,
        call_fn=call_fn,
        Y=y_array,  # Y allows for other arrays to be compared
    )

    return np.concatenate(list(chain.from_iterable(pairwise_nested)), axis=0)

if __name__ == '__main__':
    csv_path = '/mnt/p/Discovery Chemistry/Library Synthesis/Enumeration/Library Enumeration/XPL0195/presynthesis/output/diversity_test.csv'
    embedding = 'ECFP2'
    n_to_pick = 2000
    
    import pandas as pd
    
    from configs.similarity import MaxMinSampler
    from transforms.embeddings import parallel_map_get_embeddings
    
    df = pd.read_csv(csv_path)
    fps = parallel_map_get_embeddings(list(df['SMILES']), embedding=embedding)

    if isinstance(fps[0], np.ndarray):  # convert to bitvect
        fps = [
            DataStructs.cDataStructs.CreateFromBitString("".join(fps[i].astype(str)))
            for i in range(len(fps))
        ]

    sampler = MaxMinSampler(
        fps = fps, 
        embedding=embedding,
        n_to_pick = n_to_pick
    )

    df['Keep'] = sampler.lazypick()
    df.to_csv(csv_path, index=False)


