from typing import List

import pandas as pd
from rdkit import Chem

from configs.inputs import ReactionConfig, FunctionalGroupConfig

from tools.substructure_searching import (
    parallel_map_generate_fg_match_dicts,
    parallel_map_count_substructure_matches,
)
from xviewer_queries.functional_groups import generate_dict_of_functional_groups


def filter_reactivity(
    df: pd.DataFrame,
    bbs: List[Chem.rdchem.Mol],
    reactions: List[ReactionConfig],
    override_multiple_substructure_matches: bool = False,
) -> pd.DataFrame:
    """
    Labels rows of `df` which have either zero or more than 1 reactive functional groups.
    Returns the labeled `df`
    """

    df["Reactive FG"] = [""] * len(df)
    df["Reactive FG Issue"] = [""] * len(df)
    for reaction in reactions:
        for fg in reaction.allowed_fgs_reactant2:
            substructure_match_counts = parallel_map_count_substructure_matches(
                bbs, fg.smarts
            )
            correct_reaction_mask = [
                df["Reaction 1"][i] == reaction.rxn_name for i in range(len(df))
            ]
            for i in range(len(df)):
                if correct_reaction_mask[i]:
                    if substructure_match_counts[i] == 1:
                        if df["Reactive FG"][i] == "":
                            df.loc[i, ["Reactive FG"]] = fg.fg_name
                        else:
                            df.loc[i, ["Reactive FG"]] = "Multiple"
                    elif substructure_match_counts[i] > 1:
                        if override_multiple_substructure_matches:
                            df.loc[i, ["Reactive FG"]] = fg.fg_name
                        else:
                            df.loc[i, ["Reactive FG"]] = "Multiple"

    # incompatible FGS
    for i in range(len(df)):
        if df["Reactive FG"][i] == "":
            df.loc[i, ["Include"]] = "No"
            df.loc[i, ["Reactive FG Issue"]] = "No reactive functional group"
        elif df["Reactive FG"][i] == "Multiple":
            df.loc[i, ["Include"]] = "No"
            df.loc[i, ["Reactive FG Issue"]] = "Multiple reactive functional groups"
    return df


def filter_undesirable_substructures(
    df: pd.DataFrame,
    monosynthons: List[Chem.rdchem.Mol],
    fg_list: List[FunctionalGroupConfig],
) -> pd.DataFrame:
    """
    Labels rows of `df` which match a Warn or Alert FG, and marks Alerts for exclusion.
    Returns the labeled `df`
    """

    if "Include" not in list(df.columns):
        df["Include"] = [""] * len(df)
    alertfgs = [fg for fg in fg_list if fg.compound_proposal_rule == "Alert"]
    warningfgs = [fg for fg in fg_list if fg.compound_proposal_rule == "Warn"]
    list_of_alert_matches = parallel_map_generate_fg_match_dicts(alertfgs, monosynthons)
    list_of_warning_matches = parallel_map_generate_fg_match_dicts(
        warningfgs, monosynthons
    )
    alert_labels = [""] * len(df)
    warning_labels = [""] * len(df)
    for j in range(len(list_of_alert_matches)):
        mask = list(list_of_alert_matches[j].values())[0]
        if any(mask):
            for i in range(len(df)):
                if mask[i]:
                    alert_labels[i] = (
                        alert_labels[i]
                        + "~ "
                        + list(list_of_alert_matches[j].keys())[0]
                        + " ~"
                    )
    for j in range(len(list_of_warning_matches)):
        mask = list(list_of_warning_matches[j].values())[0]
        if any(mask):
            for i in range(len(df)):
                if mask[i]:
                    warning_labels[i] = (
                        warning_labels[i]
                        + "~ "
                        + list(list_of_warning_matches[j].keys())[0]
                        + " ~"
                    )

    df["Alert FGs"] = alert_labels
    df["Warning FGs"] = warning_labels

    df["Exclude Due to Alert Substructure"] = [""] * len(df)
    for i in range(len(df)):
        if df["Alert FGs"][i] != "":
            df.loc[i, ["Include"]] = "No"
            df.loc[i, ["Exclude Due to Alert Substructure"]] = "Yes"
    return df


if __name__ == "__main__":
    from time import time

    smiles_list = [
        "[2H]CN(C(C1=CC=CC=C1)=O)C",
        "[2H]CN(C(C1=CC=C(O)C=C1)=O)C",
        "[2H]CN(C(C1=CC=C(O)C=C1)=O)CCSC",
        "[2H]CN(C(C1=CC=C(O)C=C1)=O)CCS",
        "[2H]CN(C(C1=CC=CC=C1)=O)CCS",
        "[2H]CN(C(C1=CC=CC(CSSC)=C1)=O)CCS",
    ]
    monosynthons = [Chem.MolFromSmiles(x) for x in smiles_list] * 5

    df = pd.DataFrame({"Idx": [1, 2, 3, 4, 5, 6] * 5, "SMILES": smiles_list * 5})

    fgs = generate_dict_of_functional_groups().values()
    tick = time()
    df = filter_undesirable_substructures(df, monosynthons, fgs)
    tock = time()
    df.to_csv("Testing_undesirable_substructure_filtering.csv", index=False)
    took = round((tock - tick) / 60, 2)
    print(
        f"filter_undesirable_substructures took {took} minutes for {len(df)} molecules and {len(fgs)} functional groups"
    )
    print(df)
    # filter_undesirable_substructures took 1.8 minutes for 30 molecules and 611 functional groups (20Apr2023)

    #         Idx                             SMILES Include               Alert FGs              Warning FGs Exclude Due to Alert Substructure
    # 0     1          [2H]CN(C(C1=CC=CC=C1)=O)C
    # 1     2       [2H]CN(C(C1=CC=C(O)C=C1)=O)C                                               ~ Phenol ~
    # 2     3    [2H]CN(C(C1=CC=C(O)C=C1)=O)CCSC                                  ~ Phenol ~~ Thioether ~
    # 3     4     [2H]CN(C(C1=CC=C(O)C=C1)=O)CCS      No               ~ Thiol ~               ~ Phenol ~                               Yes
    # 4     5        [2H]CN(C(C1=CC=CC=C1)=O)CCS      No               ~ Thiol ~                                                        Yes
    # 5     6  [2H]CN(C(C1=CC=CC(CSSC)=C1)=O)CCS      No  ~ Disulfide ~~ Thiol ~                                                        Yes
    # 6     1          [2H]CN(C(C1=CC=CC=C1)=O)C
    # 7     2       [2H]CN(C(C1=CC=C(O)C=C1)=O)C                                               ~ Phenol ~
    # 8     3    [2H]CN(C(C1=CC=C(O)C=C1)=O)CCSC                                  ~ Phenol ~~ Thioether ~
    # 9     4     [2H]CN(C(C1=CC=C(O)C=C1)=O)CCS      No               ~ Thiol ~               ~ Phenol ~                               Yes
    # 10    5        [2H]CN(C(C1=CC=CC=C1)=O)CCS      No               ~ Thiol ~                                                        Yes
    # 11    6  [2H]CN(C(C1=CC=CC(CSSC)=C1)=O)CCS      No  ~ Disulfide ~~ Thiol ~                                                        Yes
    # 12    1          [2H]CN(C(C1=CC=CC=C1)=O)C
    # 13    2       [2H]CN(C(C1=CC=C(O)C=C1)=O)C                                               ~ Phenol ~
    # 14    3    [2H]CN(C(C1=CC=C(O)C=C1)=O)CCSC                                  ~ Phenol ~~ Thioether ~
    # 15    4     [2H]CN(C(C1=CC=C(O)C=C1)=O)CCS      No               ~ Thiol ~               ~ Phenol ~                               Yes
    # 16    5        [2H]CN(C(C1=CC=CC=C1)=O)CCS      No               ~ Thiol ~                                                        Yes
    # 17    6  [2H]CN(C(C1=CC=CC(CSSC)=C1)=O)CCS      No  ~ Disulfide ~~ Thiol ~                                                        Yes
    # 18    1          [2H]CN(C(C1=CC=CC=C1)=O)C
    # 19    2       [2H]CN(C(C1=CC=C(O)C=C1)=O)C                                               ~ Phenol ~
    # 20    3    [2H]CN(C(C1=CC=C(O)C=C1)=O)CCSC                                  ~ Phenol ~~ Thioether ~
    # 21    4     [2H]CN(C(C1=CC=C(O)C=C1)=O)CCS      No               ~ Thiol ~               ~ Phenol ~                               Yes
    # 22    5        [2H]CN(C(C1=CC=CC=C1)=O)CCS      No               ~ Thiol ~                                                        Yes
    # 23    6  [2H]CN(C(C1=CC=CC(CSSC)=C1)=O)CCS      No  ~ Disulfide ~~ Thiol ~                                                        Yes
    # 24    1          [2H]CN(C(C1=CC=CC=C1)=O)C
    # 25    2       [2H]CN(C(C1=CC=C(O)C=C1)=O)C                                               ~ Phenol ~
    # 26    3    [2H]CN(C(C1=CC=C(O)C=C1)=O)CCSC                                  ~ Phenol ~~ Thioether ~
    # 27    4     [2H]CN(C(C1=CC=C(O)C=C1)=O)CCS      No               ~ Thiol ~               ~ Phenol ~                               Yes
    # 28    5        [2H]CN(C(C1=CC=CC=C1)=O)CCS      No               ~ Thiol ~                                                        Yes
    # 29    6  [2H]CN(C(C1=CC=CC(CSSC)=C1)=O)CCS      No  ~ Disulfide ~~ Thiol ~                                                        Yes
