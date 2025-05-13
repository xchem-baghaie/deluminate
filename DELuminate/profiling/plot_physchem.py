from typing import Dict, List
import os
from os.path import exists
import ast

import pandas as pd
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
from PIL import Image
import ternary
from ternary.helpers import simplex_iterator

from configs.property_profiles import (
    RDKIT_PROPERTIES,
    VALID_PROPERTIES,
    PMI_XY_COLUMN_NAMES,
    PMI_RDS_COLUMN_NAMES,
    PROPERTIES_FOR_GRID,
    BIN_DCT,
    COLOR_DCT,
    ROUNDING_DCT,
)


def plot_generic(
    data: List, prop: str, property_profile_directory: str, bins=None, kde: bool = False
) -> None:
    fig, ax = plt.subplots()
    if prop in ["pKa", "pKb", "Lundbeck CNS MPO Score", "Pfizer CNS MPO Score"]:
        data = [
            float(x)
            for x in data
            if x != "No value calculated"
            and x != "Calculation not licensed"
            and x != "Unable to calculate"
        ]

    if prop in ["Rod-Likeness", "Disc-Likeness", "Sphere-Likeness"]:
        expanded_data = []
        for i in range(len(data)):
            if isinstance(data[i], List):
                expanded_data.extend(
                    [np.float64(x) for x in data[i] if not np.isnan(x)]
                )
            elif isinstance(data[i], float):
                if not np.isnan(data[i]):
                    expanded_data.extend([data[i]])
                else:
                    pass
            else:
                raise Exception(
                    f"Unexpected data type encountered, {data[i]}, on row {i}"
                )
        data = expanded_data

    if bins is not None:
        # Adjust the domain of the plot, if more than 1% of it falls outside.
        # First we move to the left if necessary, then we move to the right if necessary.
        min = bins[0]
        max = bins[-1]
        delta = bins[1] - bins[0]
        while np.count_nonzero(np.array(data) < min) / len(data) > 0.01:
            # Shift the domain to the left, by delta
            bins -= delta
            # Recalculate min and max
            min = bins[0]
            max = bins[-1]
        while np.count_nonzero(np.array(data) > max) / len(data) > 0.01:
            # Shift the domain to the right, by delta
            bins += delta
            # Recalculate min and max
            min = bins[0]
            max = bins[-1]
    ax = sns.histplot(data, bins=bins, kde=kde, color=COLOR_DCT[prop])
    ax.set_xticks(bins)
    plt.title(prop + " Distribution", fontsize=18)
    if prop in ["pKa", "pKb"]:
        plt.title(prop + " Distribution (Ionizable only)", fontsize=15)
    plt.yticks([])
    ax.get_yaxis().set_visible(False)
    fig.savefig(property_profile_directory + prop + ".png", bbox_inches="tight")
    return None


def plot_lipinski_compliance(
    data: list, prop: str, property_profile_directory: str
) -> None:
    fig, _ = plt.subplots()
    unique_data = sorted(list(set(data)), reverse=True)
    ls = ["PASS", "EXCEPTIONS: 1", "EXCEPTIONS: 2", "EXCEPTIONS: 3", "EXCEPTIONS: 4"]
    unique_data = sorted(list(set(unique_data).intersection(ls)))
    del unique_data[-1]
    unique_data = ["PASS"] + unique_data
    percentages = {}
    for outcome in unique_data:
        percentages[outcome] = data.count(outcome) / len(data)
    labels = percentages.keys()
    pc = percentages.values()
    colors = [
        "#2ecc71",  # Green for no violations
        "#ffc3a0",  # Tan for one violation
        "#f08080",  # Red for two violations
        "#808080",  # Gray for three violations
        "#000000",  # Black for four violations
    ]
    plt.pie(
        pc, labels=labels, colors=colors
    )  # , autopct='%.0f%%') #Uncomment to display the percentages on the plot
    plt.title("Lipinski Ro5 Compliance", fontsize=18)
    fig.savefig(property_profile_directory + prop + ".png", bbox_inches="tight")
    return None


def plot_cns_mpo(data: list, prop: str, property_profile_directory: str):
    data = [
        float(x)
        for x in data
        if x != "No value calculated"
        and x != "Calculation not licensed"
        and x != "Unable to calculate"
    ]

    fig, _ = plt.subplots()
    thresholded_data = []
    for d in data:
        if d >= 4:
            thresholded_data.append("> 4")
        else:
            thresholded_data.append("< 4")
    percentages = {}
    for outcome in ["> 4", "< 4"]:
        percentages[outcome] = thresholded_data.count(outcome) / len(thresholded_data)
    labels = percentages.keys()
    pc = percentages.values()
    colors = [
        "#2ecc71",  # Green for > 4
        "#008FFF",  # Light Blue for < 4
    ]
    plt.pie(
        pc, labels=labels, colors=colors
    )  # , autopct='%.0f%%') #Uncomment to display the percentages on the plot
    plt.title(prop, fontsize=18)
    fig.savefig(property_profile_directory + prop + "_Pie.png", bbox_inches="tight")


def plot_pmi(
    npr1: List,
    npr2: List,
    property_profile_directory: str,
    plot_scatter: bool = True,
    plot_binned: bool = True,
) -> None:
    npr1 = [
        x
        for x in npr1
        if ((isinstance(x, float) and not np.isnan(x)) or isinstance(x, List))
    ]
    npr2 = [
        x
        for x in npr2
        if ((isinstance(x, float) and not np.isnan(x)) or isinstance(x, List))
    ]
    try:
        assert len(npr1) == len(npr2)
    except:
        raise AssertionError(
            f"{len(npr1)} elements were passed in for npr1, but {len(npr2)} elements were passed in for npr2"
        )

    npr1_list = []
    npr2_list = []
    if isinstance(npr1[0], float) and isinstance(npr2[0], float):
        # If the first row is both floats, we assume the rest of them are also floats.
        npr1_list = npr1
        npr2_list = npr2
    elif isinstance(npr1[0], List) and isinstance(npr2[0], List):
        # If the first row contains both lists, we assume the rest of them are also lists
        for i in range(len(npr1)):
            npr1_new_values = [np.float64(x) for x in npr1[i] if not np.isnan(x)]
            npr2_new_values = [np.float64(x) for x in npr2[i] if not np.isnan(x)]
            assert len(npr1_new_values) == len(npr2_new_values)
            for j in range(len(npr1_new_values)):
                npr1_list.append(npr1_new_values[j])
                npr2_list.append(npr2_new_values[j])
    else:
        raise Exception(
            "Unexpected list types were passed in to plot_pmi. Expecting both to be lists of floats, or both to be lists of lists"
        )

    rod_likeness = [npr2_list[i] - npr1_list[i] for i in range(len(npr1_list))]
    disc_likeness = [2 - (2 * npr2_list[i]) for i in range(len(npr1_list))]
    sphere_likeness = [npr1_list[i] + npr2_list[i] - 1 for i in range(len(npr1_list))]

    pmi_x = np.asarray(npr1_list)
    pmi_y = np.asarray(npr2_list)

    if plot_scatter:
        fig, ax = plt.subplots()
        xy = np.vstack((pmi_x, pmi_y))
        z = gaussian_kde(xy)(xy)
        ax = sns.scatterplot(
            x=pmi_x, y=pmi_y, c=z, cmap="jet", s=5, linewidth=0, alpha=1
        )
        x1, y1 = [np.float64(0.5), np.float64(0.0)], [np.float64(0.5), np.float64(1.0)]
        x2, y2 = [np.float64(0.5), np.float64(1.0)], [np.float64(0.5), np.float64(1.0)]
        x3, y3 = [np.float64(0.0), np.float64(1.0)], [np.float64(1.0), np.float64(1.0)]
        plt.plot(x1, y1, x2, y2, x3, y3, c="gray", ls="--", lw=2)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.text(
            0,
            1.02,
            s="Rod",
            fontsize=16,
            horizontalalignment="center",
            verticalalignment="center",
            fontweight="bold",
        )
        plt.text(
            1,
            1.02,
            s="Sphere",
            fontsize=16,
            horizontalalignment="center",
            verticalalignment="center",
            fontweight="bold",
        )
        plt.text(
            0.56,
            0.49,
            s="Disc",
            fontsize=16,
            horizontalalignment="center",
            verticalalignment="center",
            fontweight="bold",
        )

        plt.tick_params("both", width=2, labelsize=14)
        plt.axis("off")
        plt.title("Principal Moments of Inertia", fontsize=18)
        plt.tight_layout()
        fig.savefig(property_profile_directory + "PMI_Scatter.png", bbox_inches="tight")
    if plot_binned:
        fig, ax = plt.subplots()
        temp_ternary_list = list(zip(sphere_likeness, disc_likeness, rod_likeness))

        scale = 25
        bins = np.linspace(0, 1, scale + 1)[1:]

        # digitize returns the index of the bin to which a value was assigned
        x_inds = list(np.digitize(np.array([i[0] for i in temp_ternary_list]), bins))
        y_inds = list(np.digitize(np.array([i[1] for i in temp_ternary_list]), bins))
        z_inds = list(np.digitize(np.array([i[2] for i in temp_ternary_list]), bins))
        ternary_list_binned = list(zip(x_inds, y_inds, z_inds))

        # Populate ternary_dict with {(i,j,k):frequency}
        ternary_dict = dict()

        # Initiate all possible (i,j,k) vertices (only keep i and j as k is implied by these)
        for i, j, k in simplex_iterator(scale):
            ternary_dict[(i, j)] = 0
        # Count number of occurences of each (i,j,k) in binned data
        for i, j, k in ternary_list_binned:
            ternary_dict[(i, j)] += 1

        tax = ternary.TernaryAxesSubplot(ax=ax, scale=scale)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        plt.tick_params("both", width=2, labelsize=14)
        plt.axis("off")
        plt.title("Principal Moments of Inertia", fontsize=18)
        plt.tight_layout()

        plt.text(
            0,
            -1,
            s="Rod",
            fontsize=16,
            horizontalalignment="center",
            verticalalignment="center",
            fontweight="bold",
        )
        plt.text(
            25,
            -1,
            s="Sphere",
            fontsize=16,
            horizontalalignment="center",
            verticalalignment="center",
            fontweight="bold",
        )
        plt.text(
            14,
            22,
            s="Disc",
            fontsize=16,
            horizontalalignment="center",
            verticalalignment="center",
            fontweight="bold",
        )

        tax.set_title("Principal Moments of Inertia", fontsize=18)
        tax.boundary(linewidth=2)

        # Get colormap with white below vmin (eps)
        cmap_whiteZero = matplotlib.colormaps[
            "jet"
        ]  # If reverse palette is desired, add "_r" to the end of the colormap name
        cmap_whiteZero.set_under("w")
        eps = np.spacing(
            0.0
        )  # set to zero so we never show the white. If this is increased, we would start to see the white background.

        # Plot ternary heatmap
        tax.heatmap(
            ternary_dict, style="t", cmap=cmap_whiteZero, vmin=eps, colorbar=False
        )  # style can be t, d, or h (see documentation)
        fig.gca().invert_yaxis()  # turn plot upside down

        plt.savefig(property_profile_directory + "PMI_Binned.png", bbox_inches="tight")
    return None


def plot_property_profiles(
    df: pd.DataFrame,
    property_profile_directory: str,
    properties: Dict = VALID_PROPERTIES,
    keep_autogenerated_property_columns: bool = False,
) -> pd.DataFrame:
    """
    Exports plots of physchem property profiles.
    Returns the `df` with property columns
    """

    # Create the output path if it doesn't exist
    if not exists(property_profile_directory):
        os.makedirs(property_profile_directory)
    df = df.fillna("")

    # Determine which properties to use (if the dataframe contains the property or its alias as a column, the property is used)
    props_to_use = [
        prop
        for prop in list(set(list(properties.keys())))
        if (prop in list(df.columns) or properties[prop] in list(df.columns))
    ]

    # Convert to list in case we had imported a csv which contained string literals
    for col_name in PMI_XY_COLUMN_NAMES + PMI_RDS_COLUMN_NAMES:
        first_row = list(df[col_name])[0]
        if isinstance(first_row, str) and len(first_row) > 0:
            if first_row[0] == "[":
                df[col_name] = [ast.literal_eval(x) for x in list(df[col_name])]

    # Create an empty dictionary to export statistics
    stats_for_export_dict = {
        "Property": [],
        "Mean": [],
        "Stdev": [],
        "Median": [],
        "n": [],
    }

    # Generate plots for each property.
    for prop in props_to_use:
        data = df[prop].tolist()
        prop = properties[prop]  # Use its alias
        if prop in PMI_XY_COLUMN_NAMES:
            if prop == PMI_XY_COLUMN_NAMES[0]:
                # Regardless of what "prop" is, this will look for the "PMI_X_All" and "PMI_Y_All" columns only.
                if "PMI_X_All" not in list(df.columns) or "PMI_Y_All" not in list(
                    df.columns
                ):
                    raise Exception(
                        'PMI columns were named something different than "PMI_X_All" and "PMI_Y_All"'
                    )
                plot_pmi(
                    npr1=list(df["PMI_X_All"]),
                    npr2=list(df["PMI_Y_All"]),
                    property_profile_directory=property_profile_directory,
                )
            elif prop in ["Rod-Likeness", "Disc-Likeness", "Sphere-Likeness"]:
                plot_generic(
                    data=data,
                    prop=prop,
                    property_profile_directory=property_profile_directory,
                    bins=BIN_DCT[prop],
                )
            else:
                continue  # We don't need to repeat the same PMI plotting for the other PMI X/Y columns.
        elif prop == "Lipinski Compliance":
            plot_lipinski_compliance(data, prop, property_profile_directory)
        elif prop in ["Lundbeck CNS MPO Score", "Pfizer CNS MPO Score"]:
            plot_cns_mpo(data, prop, property_profile_directory)
            plot_generic(
                data=data,
                prop=prop,
                property_profile_directory=property_profile_directory,
                bins=BIN_DCT[prop],
            )
        else:
            plot_generic(
                data=data,
                prop=prop,
                property_profile_directory=property_profile_directory,
                bins=BIN_DCT[prop],
            )
            filtered_data = [
                x for x in data if isinstance(x, float) or isinstance(x, int)
            ]
            if len(filtered_data) > 0:
                stats_for_export_dict["Property"].append(prop)
                stats_for_export_dict["Mean"].append(
                    round(np.mean(filtered_data), ROUNDING_DCT[prop])
                )
                stats_for_export_dict["Stdev"].append(
                    round(np.std(filtered_data), ROUNDING_DCT[prop])
                )
                stats_for_export_dict["Median"].append(
                    round(np.median(filtered_data), ROUNDING_DCT[prop])
                )
                stats_for_export_dict["n"].append(len(filtered_data))

    # Write property statistics .csv file
    stats_for_export_df = pd.DataFrame(stats_for_export_dict)
    stats_for_export_df["Property"] = pd.Categorical(
        stats_for_export_df["Property"], categories=list(ROUNDING_DCT.keys())
    )
    stats_for_export_df.sort_values(by=["Property"], inplace=True)
    stats_for_export_df.to_csv(
        property_profile_directory + "property_statistics.csv", index=False
    )

    if all([(x in list(df.columns)) for x in PROPERTIES_FOR_GRID]):
        mw = Image.open(property_profile_directory + properties["CA_MW"] + ".png")
        clogp = Image.open(property_profile_directory + properties["CA_cLogP"] + ".png")
        clogd = Image.open(
            property_profile_directory + properties["CA_cLogD_7.4"] + ".png"
        )
        fsp3 = Image.open(property_profile_directory + properties["CA_FSP3"] + ".png")
        hbd = Image.open(property_profile_directory + properties["CA_HBD"] + ".png")
        hba = Image.open(property_profile_directory + properties["CA_HBA"] + ".png")
        tpsa = Image.open(property_profile_directory + properties["CA_tPSA"] + ".png")
        rotb = Image.open(property_profile_directory + properties["CA_ROT"] + ".png")

        x = mw.size[0]
        y = mw.size[1]
        grid = Image.new("RGB", (4 * x, 2 * y), (250, 250, 250))
        grid.paste(mw, (0, 0))
        grid.paste(clogp, (x, 0))
        grid.paste(clogd, (2 * x, 0))
        grid.paste(fsp3, (3 * x, 0))
        grid.paste(hbd, (0, y))
        grid.paste(hba, (x, y))
        grid.paste(tpsa, (2 * x, y))
        grid.paste(rotb, (3 * x, y))
        grid.save(property_profile_directory + "grid.png")

    plt.close("all")
    if keep_autogenerated_property_columns:
        cols_to_drop = ["Mols"]
        for property in RDKIT_PROPERTIES:
            df.rename(columns={property: VALID_PROPERTIES[property]}, inplace=True)
    else:
        cols_to_drop = ["Mols"] + RDKIT_PROPERTIES
    for col in cols_to_drop:
        if col in list(df.columns):
            df = df.drop(columns=[col])
    return df
