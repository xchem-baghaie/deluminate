# DELuminate 🧬✨
The repository contains workflows for the following use cases:
- [Presynthesis enumeration](#Presynthesis-enumeration)
- [Postsynthesis enumeration](#Postsynthesis-enumeration)
- [Instance exemplar enumeration](#Instance-exemplar-enumeration)
- [Physicochemical property profile generation](#Physicochemical-property-profile-generation)
- [Chemical space visualization](#Chemical-space-visualization)
- [Exhaustive disynthon enumeration](#Exhaustive-disynthon-enumeration)
- [Exhaustive instance enumeration](#Exhaustive-instance-enumeration)

## Local installation
This installation was tested on Windows Subsystem for Linux (WSL) version 1. WSL is only required because of the tmap-viz dependency, though it should presumably work on Linux as well.
Windows installation is not currently supported, but if this repository were to be forked and only modules which do not rely on the tmap-viz dependency were retained in the forked version, a Windows installation of that version should be feasible.

Ensure the latest version of poetry is installed, using `pip install poetry`. Then run the following commands:

```
git clone https://github.com/RyanTWalsh/DELuminate.git
cd DELuminate
poetry install
```

Then, add the following line to the `DELuminate/.venv/bin/activate` file, substituting "YOUR_USERNAME" for your own username.
This will ensure that relative imports work within the cloned repository

`export PYTHONPATH="${PYTHONPATH}:/home/YOUR_USERNAME/DELuminate/DELuminate/"`

## File requirements
- Libraries are defined using a populated form of `P:\Discovery Chemistry\Library Synthesis\Enumeration\Library Enumeration\Enumeration Information Template.xlsx`, which is stored in:
    - `P:\Discovery Chemistry\Library Synthesis\Enumeration\Library Enumeration\{library_name}\presyntheis\input\{library_name}_presynthesis.xlsx` for presynthesis
    - `P:\Discovery Chemistry\Library Synthesis\Enumeration\Library Enumeration\{library_name}\postsyntheis\input\{library_name}_postsynthesis.xlsx` for postsynthesis
- In `./DELuminate/xviewer_queries/`:
    - BB_XSMILES.csv is the file output by the "All XSMILES V2" department query in X-Viewer
    - BB_XSMILES_Placeholders.csv is the file output by the "All XSMILES Placeholders" department query in X-Viewer
    - XBBs_{DATE}.csv is an export of all building blocks in X-Viewer with X-Chem permissions
        - When this file is updated, the `BB_FILEPATH` variable in `./DELuminate/configs/library_definition/` should be updated to the new filename
        - If a different BB set is required (e.g. partner BBs, or commercial BBs), that file can be deposited in this folder, and the `BB_FILEPATH` variable may be updated (temporarily!) to use that file.
    - Functional_Groups.csv is a list of functional groups (a superset of the output of the "All Functional Groups" department query in X-Viewer), also containing some additional information not yet captured in X-Viewer. It is maintained manually.
    - Reactions.csv is a list of reactions, with necessary information about tolerated/untolerated functional groups, etc. It is generated and maintained manually, and is not yet stored in X-Viewer.

## Usage
Any of the workflows in `/DELuminate/DELuminate/workflows/` can be run using the following commands, assuming the user is in the root `./DELuminate` directory and has activated the virtual environment using `source .venv/bin/activate`. "WORKFLOW.py" should be subsituted with the filename of the appropriate workflow (e.g. "presynthesis.py", without quotation marks):

`python ./DELuminate/workflows/WORKFLOW.py run --PARAM1 PARAM1_VALUE --PARAM2 PARAM2_VALUE` (etc. for other parameters. If a parameter is not included in the command, its default setting will be used.)

Input parameters and output files for individual workflows are described below.

### Presynthesis enumeration
|Parameter|Type|Default|Description|
|---------|----|-------|-----------|
|library_name|str|""|The XLIB or XPL number of the library. For example, XLIB0001 or XPL0002|
|num_instances|int|3000|The number of library instances to enumerate|
|generate_physchem_property_profiles|bool|True|Whether to calculate and generate physicochemical property profiles for the enumerated instances|
|enumerate_qc_set|bool|True|Whether to enumerate a QC set, in which every building block is used in at least one instance. This is run in addition to the main instance enumeration, which will enumerate `num_instances`|
|compute_xsmiles|bool|True|Whether to generate and QC building block XSMILES strings|
|xsmiles_directory|str|"/mnt/d/XSMILES/BB_XSMILES/"|The directory in which the building block XSMILES files for further review, and ultimately, upload to X-Viewer will be written|
|write_xsmiles|bool|False|Whether to write XSMILES lists to the XSMILES directory (overwriting existing files for the library, if they exist)|

![presynthesis_workflow](https://github.com/RyanTWalsh/DELuminate/assets/93160524/09d87698-9383-401e-91c3-ec501ac47849)

`./DELuminate/workflows/presynthesis.py` will generate and QC building block lists _de novo_, given the library definition provided in `P:\Discovery Chemistry\Library Synthesis\Enumeration\Library Enumeration\{library_name}\presyntheis\input\{library_name}_presynthesis.xlsx` and the full suite of available filters, and will output a log file, an output .xlsx file containing filtered BB lists for each cycle and instance enumerations, and (optionally) a populated `property_profiles` subdirectory to `P:\Discovery Chemistry\Library Synthesis\Enumeration\Library Enumeration\{library_name}\presyntheis\output\`


### Postsynthesis enumeration
|Parameter|Type|Default|Description|
|---------|----|-------|-----------|
|library_name|str|""|The XLIB or XPL number of the library. For example, XLIB0001 or XPL0002|
|num_instances|int|3000|The number of library instances to enumerate|
|generate_physchem_property_profiles|bool|True|Whether to calculate and generate physicochemical property profiles for the enumerated instances|
|enumerate_qc_set|bool|True|Whether to enumerate a QC set, in which every building block is used in at least one instance. This is run in addition to the main instance enumeration, which will enumerate `num_instances`|
|compute_xsmiles|bool|True|Whether to generate and QC building block XSMILES strings|
|xsmiles_directory|str|"/mnt/d/XSMILES/BB_XSMILES/"|The directory in which the building block XSMILES files for further review, and ultimately, upload to X-Viewer will be written|
|write_xsmiles|bool|False|Whether to write XSMILES lists to the XSMILES directory (overwriting existing files for the library, if they exist)|

![postsynthesis_workflow](https://github.com/RyanTWalsh/DELuminate/assets/93160524/800b1585-a64b-4f8a-8675-5ccb5f8f90b9)

`./DELuminate/workflows/postsynthesis.py` will QC user-supplied building block lists using a relevant subset of the available filters, given the library definition provided in `P:\Discovery Chemistry\Library Synthesis\Enumeration\Library Enumeration\{library_name}\postsyntheis\input\{library_name}_postsynthesis.xlsx`, and will output a log file, an output .xlsx file containing filtered BB lists for each cycle and instance enumerations, and (optionally) a populated `property_profiles` subdirectory to `P:\Discovery Chemistry\Library Synthesis\Enumeration\Library Enumeration\{library_name}\postsyntheis\output\`

### Instance exemplar enumeration
|Parameter|Type|Default|Description|
|---------|----|-------|-----------|
|library_name|str|""|The XLIB or XPL number of the library. For example, XLIB0001 or XPL0002|
|num_instances|int|3000|The number of library instances to enumerate|
|generate_physchem_property_profiles|bool|True|Whether to calculate and generate physicochemical property profiles for the enumerated instances|
|enumerate_qc_set|bool|True|Whether to enumerate a QC set, in which every building block is used in at least one instance. This is run in addition to the main instance enumeration, which will enumerate `num_instances`|

![enumerate_instances_workflow](https://github.com/RyanTWalsh/DELuminate/assets/93160524/5f97f97d-01e5-424e-8164-ee489989a625)

`./DELuminate/workflows/enumerate_instances.py` will behave the same as `postsynthesis.py`, but will not QC and filter the building block list before attempting to enumerate instances.

### Physicochemical property profile generation
|Parameter|Type|Default|Description|
|---------|----|-------|-----------|
|library_name|str|""|The XLIB or XPL number of the library. For example, XLIB0001 or XPL0002|
|mode|str|"postsynthesis"|"presynthesis" or "postsynthesis"|

![property_profiles_workflow](https://github.com/RyanTWalsh/DELuminate/assets/93160524/17585ed4-984e-42a1-b28e-a75d30f5e921)

`./DELuminate/workflows/property_profiles.py` will generate a list of instances using the "Instance Enumeration" tab of `P:\Discovery Chemistry\Library Synthesis\Enumeration\Library Enumeration\{library_name}\postsyntheis\output\{library_name}_output.xlsx`. It then calculates any physchem properties which are not present in that dataset, and outputs plots of all properties to `P:\Discovery Chemistry\Library Synthesis\Enumeration\Library Enumeration\{library_name}\postsyntheis\output\property_profiles\` as .png files, overwriting any existing files with the same names if they are present.

### Chemical space visualization
|Parameter|Type|Default|Description|
|---------|----|-------|-----------|
|filename|str|None|The name of the .csv file containing the data. ".csv" should be included.|
|parent_directory|str|"/mnt/p/Discovery Chemistry/Chemical Space/"|The path to the parent directory|
|generate_fresh_coordinates|bool|False|If the X and Y coordinates are found in the input file, do not use them for the plot. Instead, calculate fresh coordinates, and use those.|
|overwrite_input_file_to_append_coordinates|bool|True|Whether to overwrite the input file with the same fields + coordinate fields. If the coordinate fields did not exist in the input file, they will be appended, and if they existed in the input file, they will be overwritten.|
|dimensionality_reduction|str|"UMAP"|The dimensionality reduction method|
|embedding|str|"ECFP4"|The embedding used for dimensionality reduction (usually a type of fingerprint)|
|plot_config_string|str|"seaborn"|The plotting configuration. Currently, only "seaborn" and "faerun" are supported.|
|change_inner_points_to_white_for_tmaps|bool|False|Whether to change the color of the inner points to white for tmap plots|

![visualize_chemical_space_workflow](https://github.com/RyanTWalsh/DELuminate/assets/93160524/05809457-af1e-4e96-8f2b-55d478eab9a5)

`./DELuminate/workflows/visualize_chemical_space.py` will read the file at `{parent_directory}/{filename}` and, assuming a column header of "SMILES" is present, will calculate embeddings for all of the molecules. It will then use those embeddings to calculate coordinates for each molecule. If the X and Y coordinates were already present in the input file and `generate_fresh_coordinates` = False, this step is skipped and the existing coordinates are used. If `overwrite_input_file_to_append_coordinates` = True, the input file is overwritten at this point.

In the plotting step, the plotting configuration and the configuration for each column to color is built. Currently, the column color configurations are provided to the script by manually editing the `color_config_dict` dictionary at the top of the script, but in a future update, a mechanism to pass in that information as a parameter may be implemented. Plots are then output to `{parent_directory}/{filename_without_the_.csv_suffix}/`, filed by their interactivity, dimensionality reduction method, and embedding.

### Exhaustive disynthon enumeration
|Parameter|Type|Default|Description|
|---------|----|-------|-----------|
|library_name|str|""|The XLIB or XPL number of the library. For example, XLIB0001 or XPL0002|
|enum_destination|str|"X-Viewer"|"X-Viewer" or "S3"|

![xsmiles_disynthon_enumeration_workflow](https://github.com/RyanTWalsh/DELuminate/assets/93160524/191e7494-af9e-4e1f-8fce-ab8f2be4f5cc)

`./DELuminate/workflows/xsmiles_disynthon_enumeration.py` will define the library scheme using `P:\Discovery Chemistry\Library Synthesis\Enumeration\Library Enumeration\{library_name}\postsyntheis\output\{library_name}_output.xlsx`, lookup the associated BB XSMILES and placeholder XSMILES from `./DELuminate/xviewer_queries/BB_XSMILES.csv` and `./DELuminate/xviewer_queries/BB_XSMILES_Placeholders.csv`, respectively, then enumerate all disynthon combinations which do not use any BB with an instruction of "FAIL<>" via concatenation of the BB and placeholder XSMILES strings. 

If `enum_destination` = "X-Viewer", it will output each disynthon file to `/mnt/d/XSMILES/Disynthon_Enumeration/{library_name}/{library_name}_{disynthon_type}.csv`, which should then be passed into the KNIME workflow located at `/mnt/d/knime-workspace/XSMILES_V2_Disynthon_Enumeration_QC` for SMILES canonicalization and QC (visual inspection, ensuring the expected number of disynthons were enumerated, and that all disynthons are valid molecules), prior to uploading the disynthons to X-Viewer. 

If `enum_destination` = "S3", it will output each disynthon file to `/mnt/d/Virtual_Libraries/DELs/{library_deck}/{n_cycles}-cycle/{library_name}/Disynthons/{disynthon_type}.parquet`, which is a directory containing individual .parquet files of approximately 10K rows each. As these files are not subjected to the same QC as the .csv files destined for X-Viewer upload, it is recommended that they are only generated after confirming that all disynthons in the .csv file to be uploaded to X-Viewer passed QC.

### Exhaustive instance enumeration
|Parameter|Type|Default|Description|
|---------|----|-------|-----------|
|library_name|str|""|The XLIB or XPL number of the library. For example, XLIB0001 or XPL0002|

![xsmiles_instance_enumeration_workflow](https://github.com/RyanTWalsh/DELuminate/assets/93160524/b86576d9-f06f-4b84-a958-d624541db14e)

`./DELuminate/workflows/xsmiles_disynthon_enumeration.py` will define the library scheme using `P:\Discovery Chemistry\Library Synthesis\Enumeration\Library Enumeration\{library_name}\postsyntheis\output\{library_name}_output.xlsx`, lookup the associated BB XSMILES from `./DELuminate/xviewer_queries/BB_XSMILES.csv`, and enumerate all instance combinations which do not use any BB with an instruction of "FAIL<>" via concatenation of the BB XSMILES strings in chunks of about 100M rows each. Each chunk is written to `/mnt/d/Virtual_Libraries/DELs/{library_deck}/{n_cycles}-cycle/{library_name}/Instances/{cycle_strings}.parquet`, which is a directory containing individual .parquet files of approximately 10K rows each. As these files are not subjected to any QC beyond confirming during workflow execution that the total number of rows across all .parquet files matches the expected library size, it is recommended that the full instance files are only generated once all disynthon files for the library have passed QC.
