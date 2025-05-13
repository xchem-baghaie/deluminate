import pandas as pd


def filter_tbbs(df: pd.DataFrame) -> pd.DataFrame:
    """
    Labels rows of `df` in which the XBBID has a "TBB" prefix.
    Returns the labeled `df`
    """

    df["TBB"] = [""] * len(df)
    for i in range(len(df)):
        if df["XBBID"][i].startswith("TBB"):
            df.loc[i, ["Include"]] = "No"
            df.loc[i, ["TBB"]] = "Yes"
    return df
