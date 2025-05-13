import json
import os

import pandas as pd
from metaflow import FlowSpec, Parameter, step

from configs.visualizations import (
    SeabornConfig,
    FaerunConfig,
    ColumnToColorConfig,
    CUSTOM_COLORMAPS,
    COLUMN_COLORING_CONFIG_FILEPATH
)
from visualization.calculate_coordinates import calculate_nodes_and_edges
from visualization.plot import plot

class Visualize_Chemical_Space(FlowSpec):
    filename = Parameter(
        "filename",
        help='The name of the .csv file containing the data. ".csv" should be included.',
        default=None,
    )

    parent_directory = Parameter(
        "parent_directory",
        help="The path to the parent directory",
        default="/mnt/p/Discovery Chemistry/Chemical Space/",
    )

    generate_fresh_coordinates = Parameter(
        "generate_fresh_coordinates",
        help="If the X and Y coordinates are found in the input file, do not use them for the plot. Instead, calculate fresh coordinates, and use those.",
        default=False,
    )

    overwrite_input_file_to_append_coordinates = Parameter(
        "overwrite_input_file_to_append_coordinates",
        help="Whether to overwrite the input file with the same fields + coordinate fields. If the coordinate fields did not exist in the input file, they will be appended, and if they existed in the input file, they will be overwritten.",
        default=True,
    )

    dimensionality_reduction = Parameter(
        "dimensionality_reduction",
        help="The dimensionality reduction method",
        default="UMAP",
    )

    embedding = Parameter(
        "embedding",
        help="The embedding used for dimensionality reduction (usually a type of fingerprint)",
        default="ECFP4",
    )

    plot_config_string = Parameter(
        "plot_config_string",
        help='The plotting configuration. Currently, only "seaborn" and "faerun" are supported.',
        default="seaborn",
    )

    change_inner_points_to_white_for_tmaps = Parameter(
        "change_inner_points_to_white_for_tmaps",
        help="Whether to change the color of inner points to white for tmap plots",
        default=True,
    )

    save_eps_files = Parameter(
        "save_eps_files",
        help="Whether to save .eps files in addition to the .png files that are exported",
        default=False,
    )

    apply_transparency_to_background_and_foreground = Parameter(
        "apply_transparency_to_background_and_foreground",
        help="Whether to apply equal transparency to all background and foreground points",
        default = False,
    )

    apply_additional_transparency_to_foreground = Parameter(
        "apply_additional_transparency_to_layered_points",
        help="Whether to apply additional transparency to foreground points, independent of whether `apply_transparency_to_background_and_foreground` = True",
        default = False,
    )

    @step
    def start(self):
        with open(COLUMN_COLORING_CONFIG_FILEPATH, "r") as f:
            # Config lists are of the form: [
            # name:str,
            # List[background_values:str],
            # List[foreground_values:str],
            # plot_each_value_as_its_own_layer:bool
            # ]
            self.color_configs = json.load(f)[self.filename]

        if self.plot_config_string == "seaborn":
            if self.apply_transparency_to_background_and_foreground:
                alpha = 0.2
            else:
                alpha = 1.0
            if self.apply_additional_transparency_to_foreground:
                alpha_outer_scale_factor_for_foreground = 0.02
                alpha_inner_scale_factor_for_foreground = 0.2
            else:
                alpha_outer_scale_factor_for_foreground = 1.0
                alpha_inner_scale_factor_for_foreground = 1.0             
            self.plot_config = SeabornConfig(
                outer_size_base=6,
                inner_size_base=2,
                outer_size_layer=12,
                inner_size_layer=4,
                alpha=alpha,
                alpha_outer_scale_factor_for_foreground=alpha_outer_scale_factor_for_foreground,
                alpha_inner_scale_factor_for_foreground=alpha_inner_scale_factor_for_foreground,
            )
        elif self.plot_config_string == "faerun":
            self.plot_config = FaerunConfig(
                outer_size_base=8,
                inner_size_base=4,
                outer_size_layer=20,
                inner_size_layer=10,
            )
        else:
            raise Exception(
                f"Unexpected plotting config, {self.plot_config_string} was passed."
            )

        print(
            f"Processing {self.filename}, {self.dimensionality_reduction} dimensionality reduction, {self.embedding} embedding, {self.plot_config_string} plotting"
        )
        self.filepath = self.parent_directory + self.filename
        self.df = pd.read_csv(self.filepath).fillna("")

        self.next(self.calculate_coordinates)

    @step
    def calculate_coordinates(self):
        self.smiles_list = list(self.df["SMILES"])

        if (
            self.generate_fresh_coordinates
            or f"{self.dimensionality_reduction}_{self.embedding}_X"
            not in list(self.df.columns)
            or f"{self.dimensionality_reduction}_{self.embedding}_Y"
            not in list(self.df.columns)
        ):
            self.x, self.y, s, t = calculate_nodes_and_edges(
                smiles_list=self.smiles_list,
                embedding=self.embedding,
                dimensionality_reduction=self.dimensionality_reduction,
            )
            if s is not None:
                self.s = list(s)
            else:
                self.s = s
            if t is not None:
                self.t = list(t)
            else:
                self.t = t
            self.df[f"{self.dimensionality_reduction}_{self.embedding}_X"] = self.x
            self.df[f"{self.dimensionality_reduction}_{self.embedding}_Y"] = self.y
            if self.dimensionality_reduction.upper() == "TMAP":
                self.df[f"{self.dimensionality_reduction}_{self.embedding}_S"] = (
                    self.s + [""]
                )  # number of edges = number of nodes - 1
                self.df[f"{self.dimensionality_reduction}_{self.embedding}_T"] = (
                    self.t + [""]
                )  # number of edges = number of nodes - 1
            if self.overwrite_input_file_to_append_coordinates:
                self.df.to_csv(self.filepath, index=False)
        else:
            self.x = list(
                self.df[f"{self.dimensionality_reduction}_{self.embedding}_X"]
            )
            self.y = list(
                self.df[f"{self.dimensionality_reduction}_{self.embedding}_Y"]
            )
            if (
                self.dimensionality_reduction.upper() == "TMAP"
                and f"{self.dimensionality_reduction}_{self.embedding}_S"
                in list(self.df.columns)
                and f"{self.dimensionality_reduction}_{self.embedding}_T"
                in list(self.df.columns)
            ):
                self.s = [
                    int(x)
                    for x in list(
                        self.df[f"{self.dimensionality_reduction}_{self.embedding}_S"]
                    )
                    if x != ""
                ]
                self.t = [
                    int(x)
                    for x in list(
                        self.df[f"{self.dimensionality_reduction}_{self.embedding}_T"]
                    )
                    if x != ""
                ]
            else:
                self.s = None
                self.t = None
        self.next(self.plot)

    @step
    def plot(self):

        columns_to_color = [
            ColumnToColorConfig(
                column_name=config[0],
                data=list(self.df[config[0]]),
                background_values=config[1],
                foreground_values=config[2],
                plot_each_value_as_its_own_layer=config[3],
            )
            for config in self.color_configs
        ]

        name = self.filename.split(".")[0]  # Remove file extensions from the name

        path = self.parent_directory + name + "/"
        if isinstance(self.plot_config, FaerunConfig):
            path = path + "Interactive Plots/"
        else:
            path = path + "Static Plots/"
        output_path = path + self.dimensionality_reduction + "/" + self.embedding + "/"
        legend_path = output_path + "Legends/"
        base_layer_path = output_path + "Base Layers/"
        base_layer_all_points_path = base_layer_path + "All Points/"
        base_layer_only_background_points_path = (
            base_layer_path + "Only Background Points/"
        )
        if not isinstance(self.plot_config, FaerunConfig):
            for dir in [
                path,
                output_path,
                legend_path,
                base_layer_path,
                base_layer_all_points_path,
                base_layer_only_background_points_path,
            ]:
                if not os.path.exists(dir):
                    os.makedirs(dir)
        else:
            for dir in [path, output_path]:
                if not os.path.exists(dir):
                    os.makedirs(dir)
        setattr(self.plot_config, "name", name)
        if isinstance(self.plot_config, SeabornConfig) and self.save_eps_files:
            setattr(self.plot_config, "formats", ["png", "eps"])
        if isinstance(self.plot_config, FaerunConfig):
            setattr(
                self.plot_config, "title", ""
            )  # Title is overlaid on the plot and used to populate the tab title in web browser...setting it to '' removes the overlay, and tab title defaults to name.
        setattr(self.plot_config, "columns_to_color", columns_to_color)
        setattr(self.plot_config, "smiles_list", self.smiles_list)
        setattr(self.plot_config, "output_path", output_path)
        setattr(self.plot_config, "legend_path", legend_path)
        setattr(self.plot_config, "base_layer_path", base_layer_path)
        setattr(
            self.plot_config, "base_layer_all_points_path", base_layer_all_points_path
        )
        setattr(
            self.plot_config,
            "base_layer_only_background_points_path",
            base_layer_only_background_points_path,
        )

        if (
            self.change_inner_points_to_white_for_tmaps
            and self.dimensionality_reduction.upper() == "TMAP"
        ):
            if isinstance(self.plot_config, SeabornConfig):
                setattr(
                    self.plot_config, "inner_color", "#FFFFFF"
                )  # Set inner color to white
            elif isinstance(self.plot_config, FaerunConfig):
                setattr(
                    self.plot_config, "inner_color", CUSTOM_COLORMAPS["white"]
                )  # Set inner color to white
            else:
                pass

        plot(
            x_coords=self.x,
            y_coords=self.y,
            start_nodes=self.s,
            end_nodes=self.t,
            cfg=self.plot_config,
        )

        self.next(self.end)

    @step
    def end(self):
        pass


if __name__ == "__main__":
    Visualize_Chemical_Space()
