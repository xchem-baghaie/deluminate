import pandas as pd


def filter_validation_for_virtual_libraries(df: pd.DataFrame) -> pd.DataFrame:
    """
    Labels rows of `df` in which there is experimental failing validation data.
    Returns the labeled `df`
    """

    for i in range(len(df)):
        if (
            df["Validation Outcome"][i] == "Fail"
            and df["Validation Outcome Source"][i] == "Experimental"
        ):
            df.loc[i, ["Include"]] = "No"
    return df
