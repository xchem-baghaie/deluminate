import pandas as pd


def filter_deck(df: pd.DataFrame, deck: str) -> pd.DataFrame:
    """
    Labels rows of `df` which are not assigned to `deck`. Returns the labeled `df`
    """

    df["Deck Alert"] = [""] * len(df)
    for i in range(len(df)):
        if df["LIBRARY_DECK_NAME"][i] != deck:
            df.loc[i, ["Include"]] = "No"
            df.loc[i, ["Deck Alert"]] = f"BB not assigned to {deck}"
    return df
