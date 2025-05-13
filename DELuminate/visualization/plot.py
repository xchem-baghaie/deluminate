from typing import List, Tuple, Union
from copy import deepcopy
from itertools import compress

import pandas as pd
import seaborn as sns
from faerun import Faerun
from matplotlib import pyplot as plt
import numpy as np

from configs.visualizations import (
    SeabornConfig,
    FaerunConfig,
    ColumnToColorConfig,
    MAX_N_VALUES_TO_PROHIBIT_PLOTTING_ALL_SEPARATE_LAYERS,
    SEABORN_PLOT_DPI,
    RANDOM_STATE,
)

np.random.seed(RANDOM_STATE)


def randomize_point_order(
    x_coords: List, y_coords: List, color_data: List
) -> Tuple[List]:
    assert len(x_coords) == len(y_coords)
    if len(x_coords) > 0:
        df = pd.DataFrame({"X": x_coords, "Y": y_coords, "C": color_data})
        df = df.sample(frac=1, random_state=RANDOM_STATE).reset_index(drop=True)

        return list(df["X"]), list(df["Y"]), list(df["C"])
    else:
        return [], [], []


def randomize_point_order_readidels_in_front(
    x_coords: List, y_coords: List, color_data: List
) -> Tuple[List]:
    assert len(x_coords) == len(y_coords)
    if len(x_coords) > 0:
        df = pd.DataFrame({"X": x_coords, "Y": y_coords, "C": color_data})
        df = df.sample(frac=1, random_state=RANDOM_STATE).reset_index(drop=True)

        df_readidels = df[df["C"] == "ReadiDELs"].reset_index(drop=True)
        df_not_readidels = df[df["C"] != "ReadiDELs"].reset_index(drop=True)
        df = pd.concat([df_not_readidels, df_readidels])

        return list(df["X"]), list(df["Y"]), list(df["C"])
    else:
        return [], [], []


def filter_column_values(
    x_coords: List,
    y_coords: List,
    column: ColumnToColorConfig,
    cfg: Union[SeabornConfig, FaerunConfig],
):

    if column.background_values is not None or column.foreground_values is not None:
        filtered_x_coords = deepcopy(x_coords)
        filtered_y_coords = deepcopy(y_coords)
        column_copy = deepcopy(column)
        layer_cfg = deepcopy(cfg)
        color_data = column_copy.data
        if hasattr(layer_cfg, "labels"):
            labels = layer_cfg.labels
        else:
            labels = None
        if column_copy.not_background_and_only_foreground_mask is not None:
            if (
                column_copy.not_background_and_only_foreground_mask.count(True) > 0
                and column_copy.not_background_and_only_foreground_mask.count(False) > 0
            ):
                filtered_x_coords = list(
                    compress(
                        filtered_x_coords,
                        column_copy.not_background_and_only_foreground_mask,
                    )
                )
                filtered_y_coords = list(
                    compress(
                        filtered_y_coords,
                        column_copy.not_background_and_only_foreground_mask,
                    )
                )
                color_data = list(
                    compress(
                        color_data, column_copy.not_background_and_only_foreground_mask
                    )
                )
                if labels is not None:
                    labels = list(
                        compress(
                            labels, column_copy.not_background_and_only_foreground_mask
                        )
                    )
            elif all(column_copy.not_background_and_only_foreground_mask):
                print(
                    f'All datapoints will be plotted for the "{column_copy.column_name}" column.'
                )
            else:
                print(
                    f'No datapoints will be plotted for the "{column_copy.column_name}" column.'
                )
                filtered_x_coords = []
                filtered_y_coords = []
                color_data = []
                if labels is not None:
                    labels = []
        elif column_copy.not_background_value_mask is not None:
            if (
                column_copy.not_background_value_mask.count(True) > 0
                and column_copy.not_background_value_mask.count(False) > 0
            ):
                filtered_x_coords = list(
                    compress(filtered_x_coords, column_copy.not_background_value_mask)
                )
                filtered_y_coords = list(
                    compress(filtered_y_coords, column_copy.not_background_value_mask)
                )
                color_data = list(
                    compress(color_data, column_copy.not_background_value_mask)
                )
                if labels is not None:
                    labels = list(
                        compress(labels, column_copy.not_background_value_mask)
                    )
            elif all(column_copy.not_background_value_mask):
                print(
                    f'All datapoints will be plotted for the "{column_copy.column_name}" column.'
                )
            else:
                print(
                    f'No datapoints will be plotted for the "{column_copy.column_name}" column.'
                )
                filtered_x_coords = []
                filtered_y_coords = []
                color_data = []
                if labels is not None:
                    labels = []
        elif column_copy.only_foreground_value_mask is not None:
            if (
                column_copy.only_foreground_value_mask.count(True) > 0
                and column_copy.only_foreground_value_mask.count(False) > 0
            ):
                filtered_x_coords = list(
                    compress(filtered_x_coords, column_copy.only_foreground_value_mask)
                )
                filtered_y_coords = list(
                    compress(filtered_y_coords, column_copy.only_foreground_value_mask)
                )
                color_data = list(
                    compress(color_data, column_copy.only_foreground_value_mask)
                )
                if labels is not None:
                    labels = list(
                        compress(labels, column_copy.only_foreground_value_mask)
                    )
            elif all(column_copy.only_foreground_value_mask):
                print(
                    f'All datapoints will be plotted for the "{column_copy.column_name}" column.'
                )
            else:
                print(
                    f'No datapoints will be plotted for the "{column_copy.column_name}" column.'
                )
                filtered_x_coords = []
                filtered_y_coords = []
                color_data = []
                if labels is not None:
                    labels = []
        else:
            print(
                f'All datapoints will be plotted for the "{column_copy.column_name}" column.'
            )

        # Normalize the values in color_data to floats in [0,1]
        unique_color_data_values = list(set(color_data))
        unique_color_data_values.sort()

        # Find the unused color values
        if column_copy.background_values is not None:
            unused_color_data_values = [
                x
                for x in list(set(column_copy.data))
                if x in column_copy.background_values
            ]
        else:
            unused_color_data_values = []
        if column_copy.foreground_values is not None:
            unused_color_data_values = list(
                set(
                    unused_color_data_values
                    + [
                        x
                        for x in list(set(column_copy.data))
                        if x not in unique_color_data_values
                        and x not in column_copy.foreground_values
                    ]
                )
            )

        unused_color_data_values.sort()

        setattr(column_copy, "unique_values", unique_color_data_values)
        int_color_data = [i for i in range(len(column_copy.unique_values))]

        if len(int_color_data) > 0:
            max_color_data = max(int_color_data)
            if max_color_data > 0:
                color_idxs = [float(i) / max_color_data for i in int_color_data]
            elif len(int_color_data) == 1:
                color_idxs = [0]
            else:
                raise Exception("Error normalizing color data")

            # Set palette based on the remaining column values
            column_copy.set_palette()

            # Use the normalized values to get positions in column.palette
            palette = [column_copy.palette(x) for x in color_idxs] + [
                cfg.inner_color
            ] * len(unused_color_data_values)

            # Update the set of unique values to include the unused color data values
            setattr(
                column_copy,
                "unique_values",
                column_copy.unique_values + unused_color_data_values,
            )

            # Update palette and hue_order attributes for the layer configuration
            setattr(layer_cfg, "palette", palette)
            setattr(layer_cfg, "hue_order", column_copy.unique_values)
            if hasattr(layer_cfg, "labels"):
                setattr(layer_cfg, "labels", labels)
            setattr(column_copy, "data", color_data)
            setattr(column_copy, "palette", palette)
        else:
            setattr(column_copy, "unique_values", unused_color_data_values)
            palette = [cfg.inner_color] * len(unused_color_data_values)
            # Update palette and hue_order attributes for the layer configuration
            setattr(layer_cfg, "palette", palette)
            setattr(layer_cfg, "hue_order", column_copy.unique_values)
            if hasattr(layer_cfg, "labels"):
                setattr(layer_cfg, "labels", labels)
            setattr(column_copy, "data", color_data)
            setattr(column_copy, "palette", palette)
        return filtered_x_coords, filtered_y_coords, column_copy, layer_cfg


def filter_column_values_complement(
    x_coords: List,
    y_coords: List,
    column: ColumnToColorConfig,
    cfg: Union[SeabornConfig, FaerunConfig],
):

    if column.background_values is not None or column.foreground_values is not None:
        filtered_x_coords = deepcopy(x_coords)
        filtered_y_coords = deepcopy(y_coords)
        column_copy = deepcopy(column)
        layer_cfg = deepcopy(cfg)
        color_data = column_copy.data
        if hasattr(layer_cfg, "labels"):
            labels = layer_cfg.labels
        else:
            labels = None
        if column_copy.not_background_and_only_foreground_mask is not None:
            if (
                column_copy.not_background_and_only_foreground_mask.count(True) > 0
                and column_copy.not_background_and_only_foreground_mask.count(False) > 0
            ):
                filtered_x_coords = list(
                    compress(
                        filtered_x_coords,
                        [
                            not x
                            for x in column_copy.not_background_and_only_foreground_mask
                        ],
                    )
                )
                filtered_y_coords = list(
                    compress(
                        filtered_y_coords,
                        [
                            not x
                            for x in column_copy.not_background_and_only_foreground_mask
                        ],
                    )
                )
                color_data = list(
                    compress(
                        color_data,
                        [
                            not x
                            for x in column_copy.not_background_and_only_foreground_mask
                        ],
                    )
                )
                if labels is not None:
                    labels = list(
                        compress(
                            labels,
                            [
                                not x
                                for x in column_copy.not_background_and_only_foreground_mask
                            ],
                        )
                    )
            elif all(column_copy.not_background_and_only_foreground_mask):
                filtered_x_coords = []
                filtered_y_coords = []
                color_data = []
                if labels is not None:
                    labels = []
            else:
                pass
        elif column_copy.not_background_value_mask is not None:
            if (
                column_copy.not_background_value_mask.count(True) > 0
                and column_copy.not_background_value_mask.count(False) > 0
            ):
                filtered_x_coords = list(
                    compress(
                        filtered_x_coords,
                        [not x for x in column_copy.not_background_value_mask],
                    )
                )
                filtered_y_coords = list(
                    compress(
                        filtered_y_coords,
                        [not x for x in column_copy.not_background_value_mask],
                    )
                )
                color_data = list(
                    compress(
                        color_data,
                        [not x for x in column_copy.not_background_value_mask],
                    )
                )
                if labels is not None:
                    labels = list(
                        compress(
                            labels,
                            [not x for x in column_copy.not_background_value_mask],
                        )
                    )
            elif all(column_copy.not_background_value_mask):
                filtered_x_coords = []
                filtered_y_coords = []
                color_data = []
                if labels is not None:
                    labels = []
            else:
                pass
        elif column_copy.only_foreground_value_mask is not None:
            if (
                column_copy.only_foreground_value_mask.count(True) > 0
                and column_copy.only_foreground_value_mask.count(False) > 0
            ):
                filtered_x_coords = list(
                    compress(
                        filtered_x_coords,
                        [not x for x in column_copy.only_foreground_value_mask],
                    )
                )
                filtered_y_coords = list(
                    compress(
                        filtered_y_coords,
                        [not x for x in column_copy.only_foreground_value_mask],
                    )
                )
                color_data = list(
                    compress(
                        color_data,
                        [not x for x in column_copy.only_foreground_value_mask],
                    )
                )
                if labels is not None:
                    labels = list(
                        compress(
                            labels,
                            [not x for x in column_copy.only_foreground_value_mask],
                        )
                    )
            elif all(column_copy.only_foreground_value_mask):
                filtered_x_coords = []
                filtered_y_coords = []
                color_data = []
                if labels is not None:
                    labels = []
            else:
                pass
        else:
            filtered_x_coords = []
            filtered_y_coords = []
            color_data = []
            if labels is not None:
                labels = []

        # Normalize the values in color_data to floats in [0,1]
        unique_color_data_values = list(set(color_data))
        unique_color_data_values.sort()

        # Find the unused color values
        if column_copy.background_values is not None:
            unused_color_data_values = [
                x
                for x in list(set(column_copy.data))
                if x in column_copy.background_values
            ]
        else:
            unused_color_data_values = []
        if column_copy.foreground_values is not None:
            unused_color_data_values = list(
                set(
                    unused_color_data_values
                    + [
                        x
                        for x in list(set(column_copy.data))
                        if x not in unique_color_data_values
                        and x not in column_copy.foreground_values
                    ]
                )
            )

        unused_color_data_values.sort()

        setattr(column_copy, "unique_values", unique_color_data_values)
        int_color_data = [i for i in range(len(column_copy.unique_values))]

        if len(int_color_data) > 0:
            max_color_data = max(int_color_data)
            if max_color_data > 0:
                color_idxs = [float(i) / max_color_data for i in int_color_data]
            elif len(int_color_data) == 1:
                color_idxs = [0]
            else:
                raise Exception("Error normalizing color data")

            # Set palette based on the remaining column values
            column_copy.set_palette()

            # Use the normalized values to get positions in column.palette
            palette = [column_copy.palette(x) for x in color_idxs] + [
                cfg.inner_color
            ] * len(unused_color_data_values)

            # Update the set of unique values to include the unused color data values
            setattr(
                column_copy,
                "unique_values",
                column_copy.unique_values + unused_color_data_values,
            )

            # Update palette and hue_order attributes for the layer configuration
            setattr(layer_cfg, "palette", palette)
            setattr(layer_cfg, "hue_order", column_copy.unique_values)
            if hasattr(layer_cfg, "labels"):
                setattr(layer_cfg, "labels", labels)
            setattr(column_copy, "data", color_data)
            setattr(column_copy, "palette", palette)
        else:
            setattr(column_copy, "unique_values", unused_color_data_values)
            palette = [cfg.inner_color] * len(unused_color_data_values)
            # Update palette and hue_order attributes for the layer configuration
            setattr(layer_cfg, "palette", palette)
            setattr(layer_cfg, "hue_order", column_copy.unique_values)
            if hasattr(layer_cfg, "labels"):
                setattr(layer_cfg, "labels", labels)
            setattr(column_copy, "data", color_data)
            setattr(column_copy, "palette", palette)

        return filtered_x_coords, filtered_y_coords, column_copy, layer_cfg


def generate_legend(palette, unique_color_data_values):
    fig_legend, ax_legend = plt.subplots()

    f = lambda m, c: plt.plot([], [], marker=m, color=c, ls="none")[0]
    handles = [f("s", palette[i]) for i in range(len(palette))]
    legend = plt.legend(
        handles,
        unique_color_data_values,
        loc="center",
        ncol=max(1, int(len(unique_color_data_values) / 10)),
    )
    ax_legend.set(xlabel="", ylabel="", xticklabels=[], yticklabels=[])
    sns.despine(left=True, bottom=True)
    ax_legend.tick_params(left=False, bottom=False)
    fig_legend = legend.figure
    fig_legend.canvas.draw()
    bbox = legend.get_window_extent()
    bbox = bbox.from_extents(*(bbox.extents + np.array([-5, -5, 5, 5])))
    bbox = bbox.transformed(fig_legend.dpi_scale_trans.inverted())
    return fig_legend, bbox


def format_seaborn_plot(fig, ax) -> None:
    ax.set_box_aspect(1)
    ax.set(xlabel="", ylabel="", xticklabels=[], yticklabels=[])
    sns.despine(left=True, bottom=True)
    ax.tick_params(left=False, bottom=False)
    ax.legend([], [], frameon=False)
    fig.tight_layout()

    return None


def plot_base_seaborn(
    x_coords: List,
    y_coords: List,
    cfg: SeabornConfig,
):

    if cfg.plot_borders_base:
        ax = sns.scatterplot(
            x=x_coords,
            y=y_coords,
            hue=[""] * len(x_coords),
            palette=[cfg.outer_color],
            s=cfg.outer_size_base,
            linewidth=cfg.linewidth,
            alpha=cfg.alpha,
        )
        ax = sns.scatterplot(
            x=x_coords,
            y=y_coords,
            hue=[""] * len(x_coords),
            palette=[cfg.inner_color],
            s=cfg.inner_size_base,
            linewidth=cfg.linewidth,
            alpha=cfg.alpha,
        )
    else:
        ax = sns.scatterplot(
            x=x_coords,
            y=y_coords,
            hue=[""] * len(x_coords),
            palette=[cfg.gray],
            s=cfg.outer_size_base,
            linewidth=cfg.linewidth,
            alpha=cfg.alpha,
        )
    return ax


def plot_layer_seaborn(
    x_coords: List,
    y_coords: List,
    color_data: List = None,
    cfg: SeabornConfig = None,
):

    if cfg.plot_borders_layer:
        ax = sns.scatterplot(
            x=x_coords,
            y=y_coords,
            hue=[""] * len(x_coords),
            palette=[cfg.outer_color],
            s=cfg.outer_size_layer,
            linewidth=cfg.linewidth,
            alpha=cfg.alpha * cfg.alpha_outer_scale_factor_for_foreground,
        )
        ax = sns.scatterplot(
            x=x_coords,
            y=y_coords,
            hue=color_data,
            palette=cfg.palette,
            s=cfg.inner_size_layer,
            linewidth=cfg.linewidth,
            alpha=cfg.alpha * cfg.alpha_inner_scale_factor_for_foreground,
            hue_order=cfg.hue_order,
        )
    else:
        ax = sns.scatterplot(
            x=x_coords,
            y=y_coords,
            hue=color_data,
            palette=cfg.palette,
            s=cfg.outer_size_layer,
            linewidth=cfg.linewidth,
            alpha=cfg.alpha * cfg.alpha_inner_scale_factor_for_foreground,
            hue_order=cfg.hue_order,
        )

    return ax


def plot_seaborn(
    x_coords: List,
    y_coords: List,
    cfg: SeabornConfig = SeabornConfig(),
) -> None:

    fig, ax = plt.subplots()
    ax = plot_base_seaborn(x_coords, y_coords, cfg)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    format_seaborn_plot(fig, ax)

    if cfg.columns_to_color is not None:
        if 'png' in cfg.formats:
            plt.savefig(
                f"{cfg.base_layer_all_points_path}{cfg.name}_base.png",
                dpi=SEABORN_PLOT_DPI,
                bbox_inches="tight",
            )
        if 'eps' in cfg.formats:
            plt.savefig(
                f"{cfg.base_layer_all_points_path}{cfg.name}_base.eps",
                dpi=SEABORN_PLOT_DPI,
                bbox_inches="tight",
                format='eps',
            )
        for column in cfg.columns_to_color:
            # Get only the background points
            fig_background, ax_background = plt.subplots()
            background_x_coords, background_y_coords, _, _ = (
                filter_column_values_complement(
                    x_coords=x_coords,
                    y_coords=y_coords,
                    column=column,
                    cfg=cfg,
                )
            )
            ax_background = plot_base_seaborn(
                background_x_coords, background_y_coords, cfg
            )
            ax_background.set_xlim(xlim)
            ax_background.set_ylim(ylim)
        
            fig_layer = deepcopy(fig_background)
            ax_layer = deepcopy(ax_background)

            filtered_x_coords, filtered_y_coords, column, layer_cfg = (
                filter_column_values(
                    x_coords=x_coords,
                    y_coords=y_coords,
                    column=column,
                    cfg=cfg,
                )
            )

            if column.unique_values == ["ChemBL", "Enamine REAL", "ReadiDELs"]:
                filtered_x_coords, filtered_y_coords, color_data = (
                    randomize_point_order_readidels_in_front(
                        filtered_x_coords, filtered_y_coords, column.data
                    )
                )
            elif column.unique_values == ["Synthesized", "Virtual"]:
                color_data = column.data
            else:
                # Randomize the order in which points are drawn, so certain classes are not always plotted in front
                filtered_x_coords, filtered_y_coords, color_data = (
                    randomize_point_order(
                        filtered_x_coords, filtered_y_coords, column.data
                    )
                )

            ax_layer = plot_layer_seaborn(
                filtered_x_coords, filtered_y_coords, color_data, layer_cfg
            )

            format_seaborn_plot(fig_layer, ax_layer)

            if 'png' in cfg.formats:
                plt.savefig(
                    f"{cfg.output_path}{cfg.name}_{column.column_name}.png",
                    dpi=SEABORN_PLOT_DPI,
                    bbox_inches="tight",
                )
            if 'eps' in cfg.formats:
                plt.savefig(
                    f"{cfg.output_path}{cfg.name}_{column.column_name}.eps",
                    dpi=SEABORN_PLOT_DPI,
                    bbox_inches="tight",
                    format='eps',
                )                

            # Save the legend separately
            fig_legend, legend_bbox = generate_legend(
                column.palette, column.unique_values
            )

            if 'png' in cfg.formats:
                fig_legend.savefig(
                    f"{cfg.legend_path}{cfg.name}_{column.column_name}_Legend.png",
                    dpi=SEABORN_PLOT_DPI,
                    bbox_inches=legend_bbox,
                )
            if 'eps' in cfg.formats:
                fig_legend.savefig(
                    f"{cfg.legend_path}{cfg.name}_{column.column_name}_Legend.eps",
                    dpi=SEABORN_PLOT_DPI,
                    bbox_inches=legend_bbox,
                    format = 'eps',
                )                

            # Get only the background points
            fig_background, ax_background = plt.subplots()
            background_x_coords, background_y_coords, _, _ = (
                filter_column_values_complement(
                    x_coords=x_coords,
                    y_coords=y_coords,
                    column=column,
                    cfg=cfg,
                )
            )
            ax_background = plot_base_seaborn(
                background_x_coords, background_y_coords, cfg
            )
            ax_background.set_xlim(xlim)
            ax_background.set_ylim(ylim)
            format_seaborn_plot(fig_background, ax_background)
            if 'png' in cfg.formats:
                plt.savefig(
                    f"{cfg.base_layer_only_background_points_path}{cfg.name}_{column.column_name}_base.png",
                    dpi=SEABORN_PLOT_DPI,
                    bbox_inches="tight",
                )
            if 'eps' in cfg.formats:
                plt.savefig(
                    f"{cfg.base_layer_only_background_points_path}{cfg.name}_{column.column_name}_base.eps",
                    dpi=SEABORN_PLOT_DPI,
                    bbox_inches="tight",
                    format = 'eps',
                )                

            if column.plot_each_value_as_its_own_layer:
                if (
                    len(column.unique_values)
                    > MAX_N_VALUES_TO_PROHIBIT_PLOTTING_ALL_SEPARATE_LAYERS
                ):
                    print(
                        f'There are {len(column.unique_values)} unique values in the "{column.column_name}" column, which is more than the maximum allowed for plotting each separately ({MAX_N_VALUES_TO_PROHIBIT_PLOTTING_ALL_SEPARATE_LAYERS}). The separate layers will not be plotted.'
                    )
                else:
                    for i in range(len(column.unique_values)):
                        fig_individual_layer = deepcopy(fig_background)
                        ax_individual_layer = deepcopy(ax_background)
                        individual_layer_cfg = deepcopy(layer_cfg)

                        value = column.unique_values[i]
                        setattr(
                            individual_layer_cfg,
                            "palette",
                            [column.palette[i]] * len(column.unique_values),
                        )
                        individual_layer_x_coords = []
                        individual_layer_y_coords = []
                        individual_layer_color_data = []

                        for j in range(len(filtered_x_coords)):
                            if color_data[j] == value:
                                individual_layer_x_coords.append(filtered_x_coords[j])
                                individual_layer_y_coords.append(filtered_y_coords[j])
                                individual_layer_color_data.append(color_data[j])

                        ax_individual_layer = plot_layer_seaborn(
                            individual_layer_x_coords,
                            individual_layer_y_coords,
                            individual_layer_color_data,
                            individual_layer_cfg,
                        )

                        format_seaborn_plot(fig_individual_layer, ax_individual_layer)

                        if 'png' in cfg.formats:
                            plt.savefig(
                                f"{cfg.output_path}{cfg.name}_{column.column_name}_{value}.png",
                                dpi=SEABORN_PLOT_DPI,
                                bbox_inches="tight",
                            )
                        if 'eps' in cfg.formats:
                            plt.savefig(
                                f"{cfg.output_path}{cfg.name}_{column.column_name}_{value}.eps",
                                dpi=SEABORN_PLOT_DPI,
                                bbox_inches="tight",
                                format='eps',
                            )                            

                        individual_layer_legend, individual_layer_bbox = (
                            generate_legend(
                                palette=[column.palette[i]],
                                unique_color_data_values=[value],
                            )
                        )
                        if 'png' in cfg.formats:
                            individual_layer_legend.savefig(
                                f"{cfg.legend_path}{cfg.name}_{column.column_name}_{value}_Legend.png",
                                dpi=SEABORN_PLOT_DPI,
                                bbox_inches=individual_layer_bbox,
                            )
                        if 'eps' in cfg.formats:
                            individual_layer_legend.savefig(
                                f"{cfg.legend_path}{cfg.name}_{column.column_name}_{value}_Legend.eps",
                                dpi=SEABORN_PLOT_DPI,
                                bbox_inches=individual_layer_bbox,
                                format='eps',
                            )                            
    else:
        if 'png' in cfg.formats:
            plt.savefig(
                f"{cfg.output_path}{cfg.name}.png",
                dpi=SEABORN_PLOT_DPI,
                bbox_inches="tight",
            )
        if 'eps' in cfg.formats:
            plt.savefig(
                f"{cfg.output_path}{cfg.name}.eps",
                dpi=SEABORN_PLOT_DPI,
                bbox_inches="tight",
                format='eps',
            )
    return None


def plot_base_faerun(
    x_coords: List,
    y_coords: List,
    cfg: FaerunConfig,
):
    f = Faerun(
        title=cfg.title,
        clear_color=cfg.clear_color,
        coords=cfg.coords,
        coords_color=cfg.coords_color,
        coords_box=cfg.coords_box,
        coords_ticks=cfg.coords_ticks,
        coords_grid=cfg.coords_grid,
        coords_tick_count=cfg.coords_tick_count,
        coords_tick_length=cfg.coords_tick_length,
        coords_offset=cfg.coords_offset,
        x_title=cfg.x_title,
        y_title=cfg.y_title,
        show_legend=cfg.show_legend,
        legend_title=cfg.legend_title,
        legend_orientation=cfg.legend_orientation,
        legend_number_format=cfg.legend_number_format,
        view=cfg.view,
        scale=cfg.scale,
        alpha_blending=cfg.alpha_blending,
        anti_aliasing=cfg.anti_aliasing,
        style=cfg.style,
        impress=cfg.impress,
        thumbnail_width=cfg.thumbnail_width,
        thumbnail_fixed=cfg.thumbnail_fixed,
    )
    if cfg.plot_borders_base:
        f.add_scatter(
            name="Outer_Base",
            data={
                "x": x_coords,
                "y": y_coords,
                "c": [0] + [1] * (len(x_coords) - 1),
                "labels": cfg.labels,
            },
            colormap=[cfg.outer_color] * len(x_coords),
            shader=cfg.shader,
            point_scale=cfg.outer_size_base,
            has_legend=False,
        )
        f.add_scatter(
            name="Inner_Base",
            data={
                "x": x_coords,
                "y": y_coords,
                "c": [0] + [1] * (len(x_coords) - 1),
                "labels": cfg.labels,
            },
            colormap=[cfg.inner_color] * len(x_coords),
            shader=cfg.shader,
            point_scale=cfg.inner_size_base,
            has_legend=False,
        )
    else:
        f.add_scatter(
            name="Base",
            data={
                "x": x_coords,
                "y": y_coords,
                "c": [0] + [1] * (len(x_coords) - 1),
                "labels": cfg.labels,
            },
            colormap=[cfg.gray] * len(x_coords),
            shader=cfg.shader,
            point_scale=cfg.inner_size_base,
            has_legend=False,
        )
    f.add_scatter(
        name=cfg.name,
        data={
            "x": x_coords,
            "y": y_coords,
            "c": [column.faerun_data for column in cfg.columns_to_color],
            "labels": cfg.labels,
        },
        mapping=cfg.mapping,
        colormap=[column.palette for column in cfg.columns_to_color],
        shader=cfg.shader,
        point_scale=cfg.inner_size_base,
        fog_intensity=cfg.fog_intensity,
        saturation_limit=cfg.saturation_limit,
        categorical=cfg.categorical,
        interactive=cfg.interactive,
        has_legend=cfg.has_legend,
        legend_title=cfg.legend_title,
        legend_labels=[column.legend_label for column in cfg.columns_to_color],
        max_legend_label=[column.max_legend_label for column in cfg.columns_to_color],
        min_legend_label=[column.min_legend_label for column in cfg.columns_to_color],
        series_title=[column.series_title for column in cfg.columns_to_color],
        ondblclick=cfg.ondblclick,
        selected_labels=cfg.selected_labels,
        label_index=cfg.label_index,
        title_index=cfg.title_index,
    )
    return f


def plot_layer_faerun(
    x_coords: List,
    y_coords: List,
    cfg: FaerunConfig,
    f,
):

    f.add_scatter(
        name=cfg.name,
        data={
            "x": x_coords,
            "y": y_coords,
            "c": [column.faerun_data for column in cfg.columns_to_color],
            "labels": cfg.labels,
        },
        mapping=cfg.mapping,
        colormap=[column.palette for column in cfg.columns_to_color],
        shader=cfg.shader,
        point_scale=cfg.inner_size_base,
        fog_intensity=cfg.fog_intensity,
        saturation_limit=cfg.saturation_limit,
        categorical=cfg.categorical,
        interactive=cfg.interactive,
        has_legend=cfg.has_legend,
        legend_title=cfg.legend_title,
        legend_labels=[column.legend_label for column in cfg.columns_to_color],
        max_legend_label=[column.max_legend_label for column in cfg.columns_to_color],
        min_legend_label=[column.min_legend_label for column in cfg.columns_to_color],
        series_title=[column.series_title for column in cfg.columns_to_color],
        ondblclick=cfg.ondblclick,
        selected_labels=cfg.selected_labels,
        label_index=cfg.label_index,
        title_index=cfg.title_index,
    )
    return f


def add_tree_faerun(point_helper_name: str, start_nodes, end_nodes, f):
    if start_nodes is None or end_nodes is None:
        pass
    else:
        import tmap as tm

        if not isinstance(start_nodes, tm.VectorUint):
            start_nodes = tm.VectorUint(start_nodes)
        if not isinstance(end_nodes, tm.VectorUint):
            end_nodes = tm.VectorUint(end_nodes)
        f.add_tree(
            "Tree",
            {"from": start_nodes, "to": end_nodes},
            point_helper=point_helper_name,
        )

    return f


def plot_faerun(
    x_coords: List,
    y_coords: List,
    start_nodes=None,
    end_nodes=None,
    cfg: FaerunConfig = FaerunConfig(),
) -> None:
    if cfg.columns_to_color is not None:
        for column in cfg.columns_to_color:
            column.set_palette()
            column.update_faerun_data()
        cfg.prune_problematic_columns()
    cfg.generate_labels()

    f = plot_base_faerun(
        x_coords,
        y_coords,
        cfg,
    )

    # TODO plotting each layer individually

    f = add_tree_faerun(
        point_helper_name=cfg.name, start_nodes=start_nodes, end_nodes=end_nodes, f=f
    )
    f.plot(cfg.name, path=cfg.output_path, template="smiles")

    return None


def plot(
    x_coords: List[Union[int, float]],
    y_coords: List[Union[int, float]],
    start_nodes=None,  # List[Union[int, tmap.VectorUint]]
    end_nodes=None,  # List[Union[int, tmap.VectorUint]]
    cfg: Union[SeabornConfig, FaerunConfig] = SeabornConfig(),
) -> None:

    try:
        assert len(x_coords) == len(y_coords)
    except AssertionError:
        raise Exception(
            f"Mismatched number of coordinates ({len(x_coords)} X coordinates, {len(y_coords)} Y coordinates)"
        )

    if isinstance(cfg, SeabornConfig):
        plot_seaborn(
            x_coords=x_coords,
            y_coords=y_coords,
            cfg=cfg,
        )

    elif isinstance(cfg, FaerunConfig):
        plot_faerun(
            x_coords=x_coords,
            y_coords=y_coords,
            start_nodes=start_nodes,
            end_nodes=end_nodes,
            cfg=cfg,
        )
    else:
        raise Exception(f"Unexpected plotting configuration, {type(cfg)}, was used")

    return None
