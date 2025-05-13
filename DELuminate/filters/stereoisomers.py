import pandas as pd


def filter_stereoisomers(df: pd.DataFrame, stereoisomer_threshold: int) -> pd.DataFrame:
    """
    Labels rows of `df` which have more stereoisomers in the monosynthon than `stereoisomer_threshold`. Returns the labeled `df`
    """

    df[f"More than {str(stereoisomer_threshold)} Stereoisomers in Monosynthon"] = [
        ""
    ] * len(df)
    for i in range(len(df)):
        if df["Number of Stereoisomers in Monosynthon"][i] > stereoisomer_threshold:
            df.loc[i, ["Include"]] = "No"
            df.loc[
                i,
                [
                    f"More than {str(stereoisomer_threshold)} Stereoisomers in Monosynthon"
                ],
            ] = "Yes"
    return df
