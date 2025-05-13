from typing import Dict, List, Any, Union, Optional

from faerun import Faerun
import numpy as np
from matplotlib.cm import get_cmap
from matplotlib.colors import Colormap, LinearSegmentedColormap, ListedColormap
from scipy import stats as ss

from configs import Config
from configs.property_profiles import VALID_PROPERTIES

RANDOM_STATE = 42

TDC_FEATURES = [
    "ECFP2",
    "ECFP4",
    "ECFP6",
    "MACCS",
    "Daylight",
    "PubChem",
]

EMBEDDINGS = [
    "ECFP2",
    "ECFP4",
    "ECFP6",
    "MACCS",
    "Daylight",
    "Morgan",
    "PubChem",
    "MHFP6",  # MinHash FingerPrint, up to 6 bonds
    "MAP4",  # Minhash Atom pair FingerPrint with a radius of 2
    "MXFP",  # Macromolecule eXtended FingerPrint
    "MQN",  # Molecular Quantum Numbers
    #'ISIDA', #TODO
]

DIMENSIONALITY_REDUCTION_METHODS = [
    "TMAP",  # Tree-map
    "UMAP",  # Uniform Manifold Approximation & Projection
    "t-SNE",  # t-distributed Stochastic Neighbor Embedding
    "PCA",  # Principal Component Analysis
    "SOM",  # Self-Organizing Map
    # GTM, #Generative Topographic Mapping #TODO
    # SPE, #Stochastic Proximity Embedding #TODO
]

TMAP_L = {
    "ECFP2": 8,
    "ECFP4": 4,
    "ECFP6": 2,
    "MACCS": 8,
    "Daylight": 8,
    "Morgan": 8,
    "PubChem": 8,
    "MHFP6": 128,
    "MAP4": 128,
    "MXFP": 8,
    "MQN": 8,
}

PLOTTING_LIBRARIES = [
    "faerun",
    "seaborn",
    #'bokeh', #TODO
]

CUSTOM_COLORMAPS = {
    "categorical_cmap": LinearSegmentedColormap.from_list(
        "categorical_cmap",
        [
            "#993300",  # Brown
            "#FF002A",  # Red
            "#FF6C00",  # Orange
            "#00E600",  # Light Green
            "#00B050",  # Dark Green
            "#008FFF",  # Light Blue
            "#0008FB",  # Dark Blue
            "#AC00FF",  # Purple
            "#FF00BE",  # Pink
        ],
    ),
    "unary_cmap": ListedColormap(["#FF00BE"], name="unary_cmap"),  # Pink
    "binary_cmap": LinearSegmentedColormap.from_list(
        "binary_cmap",
        [
            "#FF6C00",  # Orange
            "#008FFF",  # Light Blue
        ],
    ),
    "ternary_cmap": LinearSegmentedColormap.from_list(
        "ternary_cmap",
        [
            "#00A8F1",  # X-Chem Blue
            "#2FA74E",  # X-Chem Green
            "#F24C80",  # X-Chem Red
        ],
    ),
    "continuous_cmap": get_cmap("jet"),
    "black": ListedColormap(["#000000"], name="black"),
    "white": ListedColormap(["#FFFFFF"], name="white"),
    "gray": ListedColormap(["#BFBFBF"], name="gray"),
}

MAX_N_VALUES_TO_BE_CONSIDERED_CATEGORICAL = 30

MAX_N_VALUES_TO_PROHIBIT_PLOTTING_ALL_SEPARATE_LAYERS = 100

SEABORN_PLOT_DPI = 300

COLUMN_COLORING_CONFIG_FILEPATH = "./DELuminate/configs/visualization_column_coloring.json"

class ColumnToColorConfig(Config):
    """
    Attributes:

    column_name: str
    data: List
    unique_values: List. If None, deduced from the data attribute, and sorted in ascending order.
    plot_each_value_as_its_own_layer:bool. If True, each unique value will be plotted as its own layer over the base
    background_values: List[Any]. If not None, rows with values in this list will not be plotted
    foreground_values: List[Any]. If not None, only rows with values in this list will be plotted
    idxs_not_background_values: List[int]. Indexes of rows which are not considered background values
    idxs_only_foreground_values: List[int]. Indexes of rows which are considered foreground values
    not_background_value_mask: A mask with True for rows which are not considered background values
    only_foreground_value_mask: A mask with True for rows which are considered foreground values
    not_background_and_only_foreground_mask: A mask with True for rows wich are not considered background values, and are considered foreground values

    max_equals_min: bool. Whether the maximum value of the column equals the minimum value of the column

    If background_values is not None and/or foreground_values is not None, idxs_not_background_values, idxs_only_foreground_values,
    not_background_value_mask, only_foreground_value_mask and not_background_and_only_foreground_value_mask are deduced from those attributes
    as part of __init__().

    Methods:
    set_palette(self): Sets the palette from self.unique_values, based on whether the column is binary, categorical, or continuous
    update_faerun_data(self): Formats the data for plotting using the Faerun library

    """

    def __init__(
        self,
        column_name: str,
        data: List,
        faerun_data: List[int] = None,
        unique_values: List = None,
        plot_each_value_as_its_own_layer: bool = False,
        background_values: Optional[List[Any]] = None,
        foreground_values: Optional[List[Any]] = None,
        palette: Union[Colormap, LinearSegmentedColormap] = None,
        legend_label: str = None,
        max_legend_label: str = None,
        min_legend_label: str = None,
        series_title: str = None,
        idxs_not_background_values: Optional[List[int]] = None,
        idxs_only_foreground_values: Optional[List[int]] = None,
        not_background_value_mask: Optional[List[bool]] = None,
        only_foreground_value_mask: Optional[List[bool]] = None,
        not_background_and_only_foreground_mask: Optional[List[bool]] = None,
        max_equals_min: bool = None,
    ):
        self.column_name = column_name
        self.data = data
        self.faerun_data = faerun_data
        self.unique_values = unique_values
        self.plot_each_value_as_its_own_layer = plot_each_value_as_its_own_layer
        self.background_values = background_values
        self.foreground_values = foreground_values
        self.palette = palette

        self.legend_label = legend_label
        self.max_legend_label = max_legend_label
        self.min_legend_label = min_legend_label
        self.series_title = series_title

        self.idxs_not_background_values = idxs_not_background_values
        self.idxs_only_foreground_values = idxs_only_foreground_values
        self.not_background_value_mask = not_background_value_mask
        self.only_foreground_value_mask = only_foreground_value_mask
        self.not_background_and_only_foreground_mask = (
            not_background_and_only_foreground_mask
        )

        self.max_equals_min = max_equals_min

        # get unique values
        if self.unique_values is None:
            unique_values = list(set(data))
            unique_values.sort()
            self.unique_values = unique_values

        # get masks for filtering values
        if (
            self.idxs_not_background_values is None
            and self.background_values is not None
        ):
            if len(self.background_values) > 0:
                self.idxs_not_background_values = [
                    i
                    for i in range(len(self.data))
                    if not self.data[i] in self.background_values
                ]
        if (
            self.idxs_only_foreground_values is None
            and self.foreground_values is not None
        ):
            if len(self.foreground_values) > 0:
                self.idxs_only_foreground_values = [
                    i
                    for i in range(len(self.data))
                    if self.data[i] in self.foreground_values
                ]

        if (
            self.not_background_value_mask is None
            and self.idxs_not_background_values is not None
        ):
            self.not_background_value_mask = [
                i in self.idxs_not_background_values for i in range(len(self.data))
            ]
        if (
            self.only_foreground_value_mask is None
            and self.idxs_only_foreground_values is not None
        ):
            self.only_foreground_value_mask = [
                i in self.idxs_only_foreground_values for i in range(len(self.data))
            ]

        if (
            self.not_background_and_only_foreground_mask is None
            and self.not_background_value_mask is not None
            and self.only_foreground_value_mask is not None
        ):
            self.not_background_and_only_foreground_mask = [
                (
                    self.not_background_value_mask[i]
                    and self.only_foreground_value_mask[i]
                )
                for i in range(len(self.data))
            ]

    def set_palette(self) -> None:
        if len(self.unique_values) == 1:
            self.palette = CUSTOM_COLORMAPS["unary_cmap"]
        elif len(self.unique_values) == 2:
            self.palette = CUSTOM_COLORMAPS["binary_cmap"]
        elif len(self.unique_values) == 3:
            self.palette = CUSTOM_COLORMAPS["ternary_cmap"]
        elif len(self.unique_values) <= MAX_N_VALUES_TO_BE_CONSIDERED_CATEGORICAL:
            self.palette = CUSTOM_COLORMAPS["categorical_cmap"]
        else:
            self.palette = CUSTOM_COLORMAPS["continuous_cmap"]
        return None

    def update_faerun_data(self):
        if isinstance(self.data[0], str):
            max_legend_label = sorted(self.data)[-1]
            min_legend_label = sorted(self.data)[0]
            labels_data, faerun_data = Faerun.create_categories(self.data)
            # consider the data to be categorical if it contains fewer than MAX_N_VALUES_TO_BE_CONSIDERED_CATEGORICAL
            if len(set(self.data)) > MAX_N_VALUES_TO_BE_CONSIDERED_CATEGORICAL:
                # call Faerun.create_categories again, this time with int data (from the previous call)
                faerun_data = ss.rankdata(
                    np.array(faerun_data) / np.max(faerun_data)
                ) / len(faerun_data)
                labels_data = None
        else:
            faerun_data = ss.rankdata(np.array(self.data) / np.max(self.data)) / len(
                self.data
            )
            labels_data = None
            maximum_value = np.max(self.data)
            minimum_value = np.min(self.data)
            if isinstance(maximum_value, int) or (
                isinstance(maximum_value, float) and maximum_value.is_integer()
            ):
                max_legend_label = str(int(maximum_value))
                min_legend_label = str(int(minimum_value))
            else:
                max_legend_label = str(round(maximum_value, 2))
                min_legend_label = str(round(minimum_value, 2))
        if self.column_name in list(VALID_PROPERTIES.keys()):
            series_title = VALID_PROPERTIES[self.column_name]
        else:
            series_title = self.column_name

        # Determine whether a divide by zero exception will be thrown when generating the scatter plot
        self.max_equals_min = max_legend_label == min_legend_label

        self.faerun_data = faerun_data
        self.legend_label = labels_data
        self.series_title = series_title
        self.max_legend_label = max_legend_label
        self.min_legend_label = min_legend_label


class SeabornConfig(Config):
    """
    Attributes:

    columns_to_color: List[ColumnToColorConfig]
    palette: color palette to be used for plotting
    plot_borders: bool, whether to plot each point once (False), or twice, once with the outer_color border and outer_size, and again with the inner_color_border and inner_size (True)
    outer_size: int, the size of the points to be plotted. If plot_borders = True, this is the border color.
    inner_size: int, only used if plot_borders = True. The size of the inner points
    linewidth: int, linewidth (alternative border strategy to the double plotting)
    alpha: float, transparency of each point
    alpha_outer_scale_factor_for_foreground: float, how much the transparency of the outer points of the layers will be scaled down from `alpha`
    alpha_inner_scale_factor_for_foreground: float, how much the transparency of the inner points of the layers will be scaled down from `alpha`
    outer_color: str, the color of the outer points
    inner_color: str, the color of the inner points
    gray: str, the desired gray hue
    output_path: str, the path to which the output files will be saved
    legend_path: str, the path to which the legend files will be saved
    base_layer_path: str, the path to which the base layers will be saved
    base_layer_all_points_path: str, the path to which base layers where all points are plotted will be saved
    base_layer_only_background_points_path: str, the path to which base layers where only the background points are plotted will be saved
    smiles_list: str, the list of SMILES strings used to generate this plot
    formats: List[str], any combination of "png" and "eps". Exported files will be saved in these formats

    """

    def __init__(
        self,
        name: str = "All_Compounds",
        columns_to_color: List[ColumnToColorConfig] = None,
        palette: Union[Colormap, LinearSegmentedColormap] = CUSTOM_COLORMAPS[
            "categorical_cmap"
        ],
        plot_borders_base: bool = True,
        plot_borders_layer: bool = True,
        outer_size_base: int = 3,
        inner_size_base: int = 1,
        outer_size_layer: int = 10,
        inner_size_layer: int = 4,
        linewidth: int = 0,
        alpha: float = 1.0,
        alpha_outer_scale_factor_for_foreground: float = 1.0,
        alpha_inner_scale_factor_for_foreground: float = 1.0,
        outer_color: str = "#000000",  # Black
        inner_color: str = "#D3D3D3",  # Gray
        gray: str = "#D3D3D3",  # Gray
        output_path: str = "./",
        legend_path: str = "./",
        base_layer_path: str = "./",
        base_layer_all_points_path: str = "./",
        base_layer_only_background_points_path: str = "./",
        smiles_list: List[str] = None,
        formats:List[str] = ["png"],
        #### Kwargs from sns.scatterplot() ####
        data=None,
        *,
        x=None,
        y=None,
        hue=None,
        size=5,
        # palette = None,
        style=None,
        hue_order=None,
        hue_norm=None,
        sizes=None,
        size_order=None,
        size_norm=None,
        markers=True,
        style_order=None,
        legend=False,  # "auto",
        ax=None,
        **kwargs,
        #### End Kwargs from sns.scatterplot() ####
    ):
        self.name = name
        self.columns_to_color = columns_to_color
        self.palette = palette
        self.plot_borders_base = plot_borders_base
        self.plot_borders_layer = plot_borders_layer
        self.outer_size_base = outer_size_base
        self.inner_size_base = inner_size_base
        self.outer_size_layer = outer_size_layer
        self.inner_size_layer = inner_size_layer
        self.linewidth = linewidth
        self.alpha = alpha
        self.alpha_outer_scale_factor_for_foreground = alpha_outer_scale_factor_for_foreground
        self.alpha_inner_scale_factor_for_foreground = alpha_inner_scale_factor_for_foreground
        self.outer_color = outer_color
        self.inner_color = inner_color
        self.gray = gray
        self.output_path = output_path
        self.legend_path = legend_path
        self.base_layer_path = base_layer_path
        self.base_layer_all_points_path = base_layer_all_points_path
        self.base_layer_only_background_points_path = (
            base_layer_only_background_points_path
        )
        self.smiles_list = smiles_list
        self.formats = formats

        #### From sns.scatterplot() ####
        self.data = data
        self.x = x
        self.y = y
        self.hue = hue
        self.size = size
        self.style = style
        self.hue_order = hue_order
        self.hue_norm = hue_norm
        self.sizes = sizes
        self.size_order = size_order
        self.size_norm = size_norm
        self.markers = markers
        self.style_order = style_order
        self.legend = legend
        self.ax = ax
        self.kwargs = kwargs
        #### End From sns.scatterplot() ####


class FaerunConfig(Config):
    """
    Attributes:

    Kwargs from Faerun base class:
    title (:obj:`str`, optional): The plot title
    clear_color (:obj:`str`, optional): The background color of the plot
    coords (:obj:`bool`, optional): Show the coordinate axes in the plot
    coords_color (:obj:`str`, optional): The color of the coordinate axes
    coords_box (:obj:`bool`, optional): Show a box around the coordinate axes
    coords_tick (:obj:`bool`, optional): Show ticks on coordinate axes
    coords_grid (:obj:`bool`, optional): Extend ticks to create a grid
    coords_tick_count (:obj:`int`, optional): The number of ticks to display per axis
    coords_tick_length (:obj:`float`, optional): The length of the coordinate ticks
    coords_offset (:obj:`float`, optional): An offset added to the coordinate axes
    x_title (:obj:`str`, optional): The title of the x-axis
    y_title (:obj:`str`, optional): The title of the y-axis
    show_legend (:obj:`bool`, optional): Whether or not to show the legend
    legend_title (:obj:`str`, optional): The legend title
    legend_orientation (:obj:`str`, optional): The orientation of the legend ('vertical' or 'horizontal')
    legend_number_format (:obj:`str`, optional): A format string applied to the numbers displayed in the legend
    view (:obj:`str`, optional): The view (front, back, top, bottom, left, right, free)
    scale (:obj:`float`, optional): To what size to scale the coordinates (which are normalized)
    alpha_blending (:obj:`bool`, optional): Whether to activate alpha blending (required for smoothCircle shader)
    anti_aliasing (:obj:`bool`, optional): Whether to activate anti-aliasing. Might improve quality at the cost of (substantial) rendering performance
    style (:obj:`Dict[str, Dict[str, Any]]`, optional): The css styles to apply to the HTML elements
    impress (:obj:`str`, optional): A short message that is shown on the HTML page
    thumbnail_width (:obj:`int`, optional): The width of the thumbnail images. Defaults to 250.
    thumbnail_fixed (:obj:`bool`, optional): Whather to show the thumbnail on the top instead as next to the mouse Mainly used for reactions. Defaults to False.

    Kwargs from Faerun.plot_scatter():
    mapping (:obj:`dict`, optional): The keys which contain the data in the input dict or the column names in the pandas :obj:`DataFrame`
    colormap (:obj:`str`, :obj:`Colormap`, :obj:`List[str]`, or :obj:`List[Colormap]` optional): The name of the colormap (can also be a matplotlib Colormap object). A list when visualizing multiple series
    shader (:obj:`str`, optional): The name of the shader to use for the data point visualization
    point_scale (:obj:`float`, optional): The relative size of the data points
    max_point_size (:obj:`int`, optional): The maximum size of the data points when zooming in
    fog_intensity (:obj:`float`, optional): The intensity of the distance fog
    saturation_limit (:obj:`float` or :obj:`List[float]`, optional): The minimum saturation to avoid "gray soup". A list when visualizing multiple series
    categorical (:obj:`bool` or :obj:`List[bool]`, optional): Whether this scatter layer is categorical. A list when visualizing multiple series
    interactive (:obj:`bool`, optional): Whether this scatter layer is interactive
    has_legend (:obj:`bool`, optional): Whether or not to draw a legend
    legend_title (:obj:`str` or :obj:`List[str]`, optional): The title of the legend. A list when visualizing multiple series
    legend_labels (:obj:`Dict` or :obj:`List[Dict]`, optional): A dict mapping values to legend labels. A list when visualizing multiple series
    min_legend_label (:obj:`str`, :obj:`float`, :obj:`List[str]` or :obj:`List[float]`, optional): The label used for the miminum value in a ranged (non-categorical) legend. A list when visualizing multiple series
    max_legend_label (:obj:`str`, :obj:`float`, :obj:`List[str]` or :obj:`List[float]`, optional): The label used for the maximum value in a ranged (non-categorical) legend. A list when visualizing multiple series
    series_title (:obj:`str` or :obj:`List[str]`, optional): The name of the series (used when multiple properites supplied). A list when visualizing multiple series
    ondblclick (:obj:`str` or :obj:`List[str]`, optional): A JavaScript snippet that is executed on double-clicking on a data point. A list when visualizing multiple series
    selected_labels: (:obj:`Dict` or :obj:`List[Dict]`, optional): A list of label values to show in the selected box. A list when visualizing multiple series
    label_index: (:obj:`int` or :obj:`List[int]`, optional): The index of the label value to use as the actual label (when __ is used to specify multiple values). A list when visualizing multiple series
    title_index: (:obj:`int` or :obj:`List[int]`, optional): The index of the label value to use as the selected title (when __ is used to specify multiple values). A list when visualizing multiple series

    columns_to_color: Dict, of the form {column_name:column_values}
    output_path: str, the path to which the output files will be saved
    legend_path: str, the path to which the legend files will be saved
    base_layer_path: str, the path to which the base layers will be saved
    base_layer_all_points_path: str, the path to which base layers where all points are plotted will be saved
    base_layer_only_background_points_path: str, the path to which base layers where only the background points are plotted will be saved
    palette: color palette to be used for plotting
    plot_borders: bool, whether to plot each point once (False), or twice, once with the outer_color border and outer_size, and again with the inner_color_border and inner_size (True)
    outer_size: int, the size of the points to be plotted. If plot_borders = True, this is the border color.
    inner_size: int, only used if plot_borders = True. The size of the inner points
    outer_color: str, the color of the outer points
    inner_color: str, the color of the inner points
    gray: str, the desired gray hue
    column_color_configs: Dict, the configurations for the column colorings in the tmap
    labels: str, the labels to be used for the scatter plot
    smiles_list: List[str], used to construct the labels for the scatter plot
    """

    def __init__(
        self,
        #### From Faerun Base class ####
        title: str = "",
        clear_color: str = "#FFFFFF",
        coords: bool = False,
        coords_color: str = "#888888",
        coords_box: bool = False,
        coords_ticks: bool = True,
        coords_grid: bool = False,
        coords_tick_count: int = 10,
        coords_tick_length: float = 2.0,
        coords_offset: float = 5.0,
        x_title: str = "",
        y_title: str = "",
        show_legend: bool = True,
        legend_title: str = "Legend",
        legend_orientation: str = "vertical",  # this is an attribute in the base Faerun class, but isn't actually used in that class!
        legend_number_format: str = "{:.0f}",
        view: str = "front",
        scale: float = 750.0,
        alpha_blending=True,
        anti_aliasing=True,
        style: Dict[str, Dict[str, Any]] = {},
        impress: str = None,
        thumbnail_width: int = 250,
        thumbnail_fixed: bool = False,
        #### End From Faerun Base Class ####
        #### From Faerun.add_scatter() ####
        name: str = "All_Compounds",
        mapping: Dict = {
            "x": "x",
            "y": "y",
            "z": "z",
            "c": "c",
            "cs": "cs",
            "s": "s",
            "labels": "labels",
            "knn": "knn",
        },
        colormap: Union[str, Colormap, List[str], List[Colormap]] = "plasma",
        shader: str = "smoothCircle",
        # point_scale: float = 2.,
        max_point_size: float = 100.0,
        fog_intensity: float = 0.0,
        saturation_limit: Union[float, List[float]] = 0.2,
        categorical: Union[bool, List[bool]] = False,
        interactive: bool = True,
        has_legend: bool = True,
        legend_labels: Union[Dict, List[Dict]] = None,
        min_legend_label: Union[str, float, List[str], List[float]] = None,
        max_legend_label: Union[str, float, List[str], List[float]] = None,
        series_title: Union[str, List[str]] = None,
        ondblclick: Union[str, List[str]] = None,
        selected_labels: Union[List, List[List]] = None,
        label_index: Union[int, List[int]] = 0,
        title_index: Union[int, List[int]] = 2,
        #### End From Faerun.add_scatter() ####
        columns_to_color: List[ColumnToColorConfig] = None,
        output_path: str = "./",
        legend_path: str = "./",
        base_layer_path: str = "./",
        base_layer_all_points_path: str = "./",
        base_layer_only_background_points_path: str = "./",
        palette: Union[Colormap, LinearSegmentedColormap] = CUSTOM_COLORMAPS[
            "categorical_cmap"
        ],
        plot_borders_base: bool = True,
        plot_borders_layer: bool = True,
        outer_size_base: int = 4,
        inner_size_base: int = 2,
        outer_size_layer: int = 10,
        inner_size_layer: int = 5,
        outer_color: Union[
            Colormap, LinearSegmentedColormap, ListedColormap
        ] = CUSTOM_COLORMAPS["black"],
        inner_color: Union[
            Colormap, LinearSegmentedColormap, ListedColormap
        ] = CUSTOM_COLORMAPS["gray"],
        gray: Union[
            Colormap, LinearSegmentedColormap, ListedColormap
        ] = CUSTOM_COLORMAPS["gray"],
        column_color_configs: Optional[Dict] = None,
        labels=None,
        smiles_list: List[str] = None,
    ):

        #### From Faerun Base Class ####
        self.title = title
        self.clear_color = clear_color
        self.coords = coords
        self.coords_color = coords_color
        self.coords_box = coords_box
        self.coords_ticks = coords_ticks
        self.coords_grid = coords_grid
        self.coords_tick_count = coords_tick_count
        self.coords_tick_length = coords_tick_length
        self.coords_offset = coords_offset
        self.x_title = x_title
        self.y_title = y_title
        self.show_legend = show_legend
        self.legend_title = legend_title
        self.legend_orientation = legend_orientation
        self.legend_number_format = legend_number_format
        self.view = view
        self.scale = scale
        self.alpha_blending = alpha_blending
        self.anti_aliasing = anti_aliasing
        self.style = style
        self.impress = impress
        self.thumbnail_width = thumbnail_width
        self.thumbnail_fixed = thumbnail_fixed
        #### End From Faerun Base Class ####

        #### From Faerun.add_scatter() ####
        self.name = name
        self.mapping = mapping
        self.colormap = colormap
        self.shader = shader
        # self.point_scale = point_scale
        self.max_point_size = max_point_size
        self.fog_intensity = fog_intensity
        self.saturation_limit = saturation_limit
        self.categorical = categorical
        self.interactive = interactive
        self.has_legend = has_legend
        self.legend_labels = legend_labels
        self.min_legend_label = min_legend_label
        self.max_legend_label = max_legend_label
        self.series_title = series_title
        self.ondblclick = ondblclick
        self.selected_labels = selected_labels
        self.label_index = label_index
        self.title_index = title_index
        #### End From Faerun.add_scatter() ####

        self.columns_to_color = columns_to_color
        self.output_path = output_path
        self.legend_path = legend_path
        self.base_layer_path = base_layer_path
        self.base_layer_all_points_path = base_layer_all_points_path
        self.base_layer_only_background_points_path = (
            base_layer_only_background_points_path
        )
        self.palette = palette
        self.plot_borders_base = plot_borders_base
        self.plot_borders_layer = plot_borders_layer
        self.outer_size_base = outer_size_base
        self.inner_size_base = inner_size_base
        self.outer_size_layer = outer_size_layer
        self.inner_size_layer = inner_size_layer
        self.outer_color = outer_color
        self.inner_color = inner_color
        self.gray = gray
        self.column_color_configs = column_color_configs
        self.labels = labels
        self.smiles_list = smiles_list

        if (
            self.labels is None
            and self.columns_to_color is not None
            and self.smiles_list is not None
        ):
            # Start label text with only the compound SMILES
            self.labels = [self.smiles_list[i] for i in range(len(self.smiles_list))]
            # Append to the label text values for each field in columns_to_color
            for column in self.columns_to_color:
                assert len(column.data) == len(self.smiles_list)
                self.labels = [
                    f"{self.labels[i]}__</a>__{column.column_name}: {str(column.data[i])}"
                    for i in range(len(column.data))
                ]
        else:
            pass  # Don't overwrite the labels that were passed in

        # Set the palette for each column, and update the faerun metadata
        if self.columns_to_color is not None:
            for column in self.columns_to_color:
                column.set_palette()
                column.update_faerun_data()

    def prune_problematic_columns(self):
        updated_columns_to_color = []
        for column in self.columns_to_color:
            if column.max_equals_min:
                print(
                    f"{column.column_name} column maximum value equals the minimum value, so will not be included in the plot"
                )
            else:
                updated_columns_to_color.append(column)
        self.columns_to_color = updated_columns_to_color

    def generate_labels(self):
        self.labels = [self.smiles_list[i] for i in range(len(self.smiles_list))]
        if self.columns_to_color is not None:
            for column_to_color in self.columns_to_color:
                assert len(column_to_color.data) == len(self.smiles_list)
                self.labels = [
                    f"{self.labels[i]}__</a>__{column_to_color.column_name}: {str(column_to_color.data[i])}"
                    for i in range(len(column_to_color.data))
                ]
