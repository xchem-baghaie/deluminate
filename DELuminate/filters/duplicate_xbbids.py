import pandas as pd


def filter_duplicate_xbbids(df: pd.DataFrame) -> pd.DataFrame:
    """
    Labels rows of `df` which have duplicate XBBIDs, which should only occur if the BB can be used in more than one reaction sequence. Returns the labeled `df`
    """

    df["XBB Used In Another Reaction Sequence"] = [""] * len(df)

    # Sort the dataframe in ascending order by the prevalence of each reaction branch, so if duplicates are found, we remove them from the less exotic reactions
    df["Sort Col"] = df["Reaction Sequence ID"].map(
        df["Reaction Sequence ID"].value_counts()
    )
    df = (
        df.sort_values(by="Sort Col", ascending=True)
        .drop("Sort Col", axis=1)
        .reset_index(drop=True)
    )

    ids_seen = []
    for i in range(len(df)):
        id = df["XBBID"][i]
        if id in ids_seen:
            df.loc[i, ["Include"]] = "No"
            df.loc[i, ["XBB Used In Another Reaction Sequence"]] = (
                "Yes: Sequence" + str(df["Reaction Sequence ID"][i])
            )
        else:
            ids_seen.append(id)
    return df
