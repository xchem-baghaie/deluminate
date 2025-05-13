import pandas as pd


def filter_df_for_invalid_xsmiles(df: pd.DataFrame) -> pd.DataFrame:
    """
    Labels rows of `df` which do not have valid XSMILES. Returns the labeled `df`
    """

    for i in range(len(df)):
        if not df["XSMILES Valid"][i]:
            df.loc[i, ["Include"]] = "No"
    return df
