import pandas as pd

from configs.property_profiles import MAX_ALLOWED_MONOSYNTHON_MW_NON_BRO5


def filter_mw(
    df: pd.DataFrame, mw: int = MAX_ALLOWED_MONOSYNTHON_MW_NON_BRO5
) -> pd.DataFrame:
    """
    Labels rows of `df` which have a monosynthon MW > `mw`. Returns the labeled `df`
    """
    for i in range(len(df)):
        if df["Monosynthon MW"][i] > mw:
            df.loc[i, ["Include"]] = "No"
    return df
