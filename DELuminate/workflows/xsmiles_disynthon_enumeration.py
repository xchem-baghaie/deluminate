import inspect
from itertools import combinations
import os

import pandas as pd
import pyarrow as pa
from pyarrow import csv
from pyarrow import dataset as pda
from metaflow import FlowSpec, Parameter, step

from tools.directory_handling import create_path
from transforms.rdkit_mols import (
    parallel_map_canonicalize_smiles_in_chunks
)
from xviewer_queries.library_definition import define_library

import logging


class XSMILES_Disynthon_Enumeration(FlowSpec):
    library_name = Parameter(
        "library_name",
        help="The XLIB or XPL number of the library. For example, XLIB0001 or XPL0002",
        default="",
    )

    enum_destination = Parameter(
        "enum_destination", help="X-Viewer or S3", default="X-Viewer"
    )

    canonicalize_smiles = Parameter(
        "canonicalize_smiles",
        help="Whether to canonicalize the SMILES strings",
        default=True,
    )

    use_xviewer_connection = Parameter(
        "use_xviewer_connection",
        help = "Whether to use the XViewer connection. If False, a static export will be used.",
        default = True,
    )

    def generate_logfile_path(self, enum_destination: str):
        if enum_destination == "X-Viewer":
            return create_path(self.output_directory) + self.library_name + "_log.txt"
        elif enum_destination == "S3":
            return (
                create_path(self.output_directory + "logs/")
                + self.library_name
                + "_disynthon_log.txt"
            )

    def log(self, message: str, level: str):
        acceptable_levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
        try:
            assert level.upper() in acceptable_levels
        except AssertionError:
            raise AssertionError(
                f"Provided level, {level}, is invalid. Expecting a value in {acceptable_levels}"
            )

        logfile = self.generate_logfile_path(self.enum_destination)

        logger = logging.getLogger("Filter_BB_Lists")
        logger.setLevel(logging.INFO)

        ch = logging.StreamHandler()
        ch.setLevel(logging.WARNING)

        fh = logging.FileHandler(logfile)
        fh.setLevel(logging.INFO)
        fh.filemode = "a"

        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )

        ch.setFormatter(formatter)
        fh.setFormatter(formatter)

        logger.addHandler(ch)
        logger.addHandler(fh)
        if level.upper() == "DEBUG":
            logger.debug(message)
        if level.upper() == "INFO":
            logger.info(message)
        if level.upper() == "WARNING":
            logger.warning(message)
        if level.upper() == "ERROR":
            logger.error(message)
        if level.upper() == "CRITICAL":
            logger.critical(message)
        return None

    @step
    def start(self):
        print(self.library_name)
        print(f"Enumeration destination: {self.enum_destination}")
        self.input_directory = create_path(
            f"/mnt/p/Discovery Chemistry/Library Synthesis/Enumeration/Library Enumeration/{self.library_name}/postsynthesis/input/"
        )
        if self.enum_destination == "X-Viewer":
            self.output_directory = create_path(
                f"/mnt/d/XSMILES/Disynthon_Enumeration/{self.library_name}/"
            )
        elif self.enum_destination == "S3":
            self.output_directory = create_path(f"/mnt/d/Virtual_Libraries/DELs/")
        else:
            raise Exception(
                f"enum_destination must be either X-Viewer or S3, but {self.enum_destination} was input."
            )
        self.input_filepath = (
            self.input_directory + self.library_name + "_postsynthesis.xlsx"
        )
        logfile = self.generate_logfile_path(enum_destination=self.enum_destination)
        with open(logfile, "w") as _:
            pass  # creates file, or overwrites it if it exists already

        self.library = define_library(
            self.library_name,
            self.input_filepath,
            compute_xsmiles=True,
            is_presynthesis=False,
            use_xviewer_connection=self.use_xviewer_connection,
        )

        def generate_library_initialization_report() -> str:
            string = f"""
            {self.library_name} was initialized successfully!
            Enumeration destination: {self.enum_destination}
            Deck: {self.library.deck}
            Cycles: {str(tuple([cycle.cycle_name for cycle in self.library.cycles]))}
            Using XViewer connection: {self.use_xviewer_connection}
            Synonym: {self.library.synonym}

            """

            return string

        library_initialization_report = generate_library_initialization_report()
        self.log(library_initialization_report, "info")
        self.cycles = self.library.cycles
        self.next(self.initialize_linker_cycle)

    @step
    def initialize_linker_cycle(self):
        self.library.linker_cycle.initialize_bb_list()
        self.log(
            f"Finished step {inspect.stack()[0][3]}.",
            "info",
        )
        self.next(self.linker_monosynthon_enumeration)

    @step
    def linker_monosynthon_enumeration(self):
        cycle_names = tuple([cycle.cycle_name for cycle in self.library.cycles])
        placeholder_reaction_branches = tuple(
            [cycle.reaction_branches[0] for cycle in self.library.cycles]
        )
        placeholder_bb_exemplars = tuple(
            [branch.monosynthon_exemplar for branch in placeholder_reaction_branches]
        )
        placeholder_xsmiles = tuple(
            [cycle.placeholder_xsmiles for cycle in self.library.cycles]
        )
        self.library.linker_cycle.monosynthon_enumeration(
            cycle_names = cycle_names,
            placeholder_reaction_branches = placeholder_reaction_branches,
            placeholder_bb_exemplars = placeholder_bb_exemplars,
            compute_xsmiles = True,
            placeholder_xsmiles = placeholder_xsmiles,
        )
        self.log(
            f"Finished step {inspect.stack()[0][3]}.",
            "info",
        )
        self.next(self.initialize_bb_list, foreach="cycles")

    @step
    def initialize_bb_list(self):
        self.cycle = self.input
        if self.cycle.forced_inclusion_bb_dict is None:
            raise ValueError(
                f"ERROR: BB list was not provided for cycle {self.cycle.cycle_name}"
            )
        self.cycle.initialize_bb_list()
        self.cycle.working_df["Include"] = ["Yes"] * len(self.cycle.working_df)
        self.log(
            f"Finished step {inspect.stack()[0][3]} for cycle {self.cycle.cycle_name}. Current split size  = {self.cycle.split_size}",
            "info",
        )
        self.next(self.get_bb_and_placeholder_xsmiles)

    @step
    def get_bb_and_placeholder_xsmiles(self):
        if self.enum_destination == "X-Viewer":
            use_int_bbids = False
        elif self.enum_destination == "S3":
            use_int_bbids = True
        else:
            raise Exception(
                f"enum_destination must be either X-Viewer or S3, but {self.enum_destination} was input."
            )

        if (
            self.library.linker_cycle.encoding_tag_set is not None
            and self.cycle.cycle_name
            == [cycle.cycle_name for cycle in self.library.cycles][0]
        ):
            prefix = self.library.linker_cycle.placeholder_xsmiles
        else:
            prefix = ""

        self.cycle.lookup_bb_xsmiles(use_int_bbids=use_int_bbids, prefix=prefix)
        self.cycle.lookup_placeholder_xsmiles()
        self.log(
            f"Finished step {inspect.stack()[0][3]} for cycle {self.cycle.cycle_name}. Current split size  = {self.cycle.split_size}. Placeholder XSMILES = {self.cycle.placeholder_xsmiles}",
            "info",
        )
        self.next(self.join_1)

    @step
    def join_1(self, inputs):
        self.merge_artifacts(
            inputs, include=["library", "input_filepath", "output_directory"]
        )
        cycle_list = [input.cycle for input in inputs]
        cycle_list.sort(key=lambda x: x.cycle_name)
        cycle_names = sorted([cycle.cycle_name for cycle in cycle_list])
        self.disynthon_types = [e for e in list(combinations(cycle_names, 2))]
        print("Disynthon types:", self.disynthon_types)
        self.library.cycles = tuple(cycle_list)
        self.log(f"Finished step {inspect.stack()[0][3]}", "info")
        self.next(self.move_to_next_step)

    @step
    def move_to_next_step(self):
        pass  # Workaround because the join step above cannot also be a split step
        self.next(self.enumerate, foreach="disynthon_types")

    @step
    def enumerate(self):
        self.disynthon_type = self.input  # tuple of the form (cycle1, cycle2)
        self.disynthon_string = (
            self.disynthon_type[0] + self.disynthon_type[1]
        )  # string
        expected_size = 1  # init
        for cycle in self.library.cycles:
            if cycle.cycle_name in self.disynthon_type:
                expected_size *= cycle.split_size
        self.disynthon_tab = self.library.enumerate_disynthons(
            self.disynthon_type, self.enum_destination
        )
        if self.enum_destination == "X-Viewer":
            self.output_filepath = (
                self.output_directory
                + self.library_name
                + "_"
                + self.disynthon_string
                + ".csv"
            )
        elif self.enum_destination == "S3":
            self.output_filepath = f"{self.output_directory}{self.library.deck}/{len(self.library.cycles)}-cycle/{self.library_name}/Disynthons/{self.disynthon_string}.parquet"
        self.log(
            f"Finished step {inspect.stack()[0][3]}. \n\nDisynthon type: {self.disynthon_string} \nExpected size = {expected_size:,} \nActual size = {len(self.disynthon_tab):,} \nOutput filepath: {self.output_filepath}\n",
            "info",
        )
        self.next(self.canonicalize)

    @step
    def canonicalize(self):
        failure_out_path = (
            self.output_directory
            + "INVALID_"
            + self.library_name
            + "_"
            + self.disynthon_string
            + ".csv"
        )
        if self.canonicalize_smiles:
            canonical_smiles = parallel_map_canonicalize_smiles_in_chunks([str(x) for x in self.disynthon_tab["SMILES"]])
            
            invalid_mask = [s == "" for s in canonical_smiles]

            if any(invalid_mask):
                # Delete any "valid" output file, if it is present
                try:
                    os.remove(self.output_filepath)
                except OSError:
                    pass

                # Write out the invalid rows, and raise an exception
                failures = self.disynthon_tab.filter(invalid_mask).to_pandas()
                failures.to_csv(failure_out_path, index=False)
                raise Exception(
                    f"{len(failures)} {self.disynthon_string} disynthons were invalid. They were written to {failure_out_path} for inspection."
                )
            else:
                # Delete any "invalid" output file, if it is present
                try:
                    os.remove(failure_out_path)
                except OSError:
                    pass

                # Overwrite the existing SMILES column with the canonical SMILES
                self.disynthon_tab = self.disynthon_tab.set_column(
                    self.disynthon_tab.schema.get_field_index("SMILES"),
                    "SMILES",
                    pa.array(canonical_smiles),
                )
        self.next(self.write)

    @step
    def write(self):
        if self.enum_destination == "X-Viewer":
            # Write subset for visual inspection
            with_instruction_df = self.disynthon_tab.filter(
                pa.compute.field("Contains Instruction") == "Y"
            ).to_pandas()
            no_instruction_df = self.disynthon_tab.filter(
                pa.compute.field("Contains Instruction") != "Y"
            ).to_pandas()
            if len(with_instruction_df) > 10:
                with_instruction_df = with_instruction_df.sample(10, random_state=42)
            if len(no_instruction_df) > 10:
                no_instruction_df = no_instruction_df.sample(10, random_state=42)
            visual_inspection_df = pd.concat([with_instruction_df, no_instruction_df])
            visual_inspection_df.to_csv(
                self.output_directory
                + self.library_name
                + "_"
                + self.disynthon_string
                + "_visual_inspection.csv",
                index=False,
            )

            # Write full table
            self.disynthon_tab = self.disynthon_tab.select(
                [f"BB{self.disynthon_type[0]}", f"BB{self.disynthon_type[1]}", "SMILES"]
            )
            csv.write_csv(self.disynthon_tab, self.output_filepath)
        elif self.enum_destination == "S3":
            pda.write_dataset(
                self.disynthon_tab,
                self.output_filepath,
                format="parquet",
                max_rows_per_file=10000,
                max_rows_per_group=10000,
                existing_data_behavior="delete_matching",
            )
        else:
            raise Exception(
                f"enum_destination must be either X-Viewer or S3, but {self.enum_destination} was input."
            )
        self.log(f"Finished step {inspect.stack()[0][3]}", "info")
        self.next(self.join_2)

    @step
    def join_2(self, inputs):
        self.merge_artifacts(
            inputs, include=["library", "input_filepath", "output_directory"]
        )
        self.log(f"Finished step {inspect.stack()[0][3]}", "info")
        self.next(self.end)

    @step
    def end(self):
        self.log("Finished!", "info")


if __name__ == "__main__":
    XSMILES_Disynthon_Enumeration()
