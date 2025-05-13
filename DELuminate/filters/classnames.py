from typing import List

import pandas as pd

from configs.classnames import UNDESIRABLE_CLASSNAMES


def filter_undesirable_classnames(
    df: pd.DataFrame, undesirable_classnames: List[str] = UNDESIRABLE_CLASSNAMES
) -> pd.DataFrame:
    """
    Labels rows of `df` which have a classname in `undesirable_classnames`. Returns the labeled `df`
    """

    df["Undesirable Classname"] = [""] * len(df)
    for i in range(len(df)):
        if df["CLASSNAME"][i] in UNDESIRABLE_CLASSNAMES:
            df.loc[i, ["Include"]] = "No"
            df.loc[i, ["Undesirable Classname"]] = "Yes"
    return df
