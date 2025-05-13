from typing import List, Union, Optional

import numpy as np
from rdkit.SimDivFilters import MaxMinPicker
from rdkit import DataStructs

from configs import Config
from tools.similarity import tanimoto_complementary, mhfp_distance


FINGERPRINT_SIMILARITY_METRIC_DICT = {
    'ECFP2':'tanimoto',
    'ECFP4':'tanimoto',
    'ECFP6':'tanimoto',
    'MACCS':'tanimoto',
    'Daylight':'tanimoto',
    'RDKit2D':'tanimoto',
    'Morgan':'tanimoto',
    'PubChem':'tanimoto',
    'MHFP6':'mhfp-distance',
    'MAP4':'mhfp-distance',
}

class MaxMinSampler(Config):
    """
    Performs a diversity selection using the MaxMin algorithm

    Attributes:
    ___________

    -fps, Union[List[DataStructs.cDataStructs.ExplicitBitVect], np.ndarray]
    -embedding, str
    -n_to_pick, int
    -random_seed, Optional[int]

    Methods:
    _______

    lazypick():
        Args: None
        Returns a list of booleans, corresponding to whether the compounds were selected or not
    
    _distance_metric()
        Args:
            -i, int
            -j, int
        Returns computed distance as a float
    """

    def __init__(
        self,
        fps: Union[List[DataStructs.cDataStructs.ExplicitBitVect], np.ndarray],
        embedding:str,
        n_to_pick: int,
        random_seed: Optional[int] = 42,
    ):
        self.fps = fps
        self.embedding = embedding
        self.n_to_pick = n_to_pick
        self.random_seed = random_seed
        self._maxmin_picker = MaxMinPicker()

    def lazypick(self) -> List[bool]:
        """
        Returns a list of bools, corresponding to whether the compound was selected in the diversity pick or not.

        MaxMinPicker.LazyPick() takes arguments:
            - distFunc: a function that should take two indices and return the distance between those two points.
            - poolSize: number of items in the pool
            - pickSize: number of items to pick from the pool
            - seed: random seed
        """
        poolSize = len(self.fps)

        if self.n_to_pick > poolSize:
            raise Exception(f'Cannot choose {self.n_to_pick} compounds from {poolSize}')

        selection = self._maxmin_picker.LazyPick(
            distFunc=self._distance_metric,
            poolSize=poolSize,
            pickSize=self.n_to_pick,
            seed=self.random_seed,
        )

        return [i in selection for i in range(poolSize)]

    def _distance_metric(self, i: int, j: int) -> float:
        """
        Compute the distance between two fingerprints.
        """
        try:
            FINGERPRINT_SIMILARITY_METRIC_DICT[self.embedding]
        except KeyError:
            raise ValueError(f'{self.embedding} is not supported for similarity comparisons')
        if FINGERPRINT_SIMILARITY_METRIC_DICT[self.embedding]  == "tanimoto":
            return tanimoto_complementary(self.fps[i], self.fps[j])
        elif FINGERPRINT_SIMILARITY_METRIC_DICT[self.embedding] == "mhfp-distance":
            return mhfp_distance(self.fps[i], self.fps[j])
        else:
            raise ValueError(
                f"Unknown metric: {FINGERPRINT_SIMILARITY_METRIC_DICT[self.embedding]}"
            )