import pandas as pd

def filter_virtual_bbs(df: pd.DataFrame) -> pd.DataFrame:
    """
    Labels rows of `df` in which building blocks are virtual.
    Returns the labeled `df`
    """

    for i in range(len(df)):
        if df["IS_VIRTUAL_BB"][i] == "Y":
            df.loc[i, ["Include"]] = "No"
    return df