import inspect

from metaflow import FlowSpec, Parameter, step

from tools.directory_handling import create_path
from xviewer_queries.library_definition import define_library

import logging


class XSMILES_Instance_Enumeration(FlowSpec):
    library_name = Parameter(
        "library_name",
        help="The XLIB or XPL number of the library. For example, XLIB0001 or XPL0002",
        default="",
    )

    canonicalize_smiles = Parameter(
        "canonicalize_smiles",
        help="Whether to canonicalize the SMILES strings",
        default=False,
    )

    use_xviewer_connection = Parameter(
        "use_xviewer_connection",
        help = "Whether to use the XViewer connection. If False, a static export will be used.",
        default = True,
    )

    def generate_logfile_path(self):
        return (
            create_path(self.output_directory + "logs/")
            + self.library_name
            + "_instance_log.txt"
        )

    def log(self, message: str, level: str):
        acceptable_levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
        try:
            assert level.upper() in acceptable_levels
        except AssertionError:
            raise AssertionError(
                f"Provided level, {level}, is invalid. Expecting a value in {acceptable_levels}"
            )

        logfile = self.generate_logfile_path()

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
        self.input_directory = create_path(
            f"/mnt/p/Discovery Chemistry/Library Synthesis/Enumeration/Library Enumeration/{self.library_name}/postsynthesis/input/"
        )
        self.output_directory = create_path(f"/mnt/d/Virtual_Libraries/DELs/")
        self.input_filepath = (
            self.input_directory + self.library_name + "_postsynthesis.xlsx"
        )
        logfile = self.generate_logfile_path()
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
        self.next(self.get_bb_xsmiles)

    @step
    def get_bb_xsmiles(self):
        if (
            self.library.linker_cycle.encoding_tag_set is not None
            and self.cycle.cycle_name
            == [cycle.cycle_name for cycle in self.library.cycles][0]
        ):
            prefix = self.library.linker_cycle.placeholder_xsmiles
        else:
            prefix = ""

        self.cycle.lookup_bb_xsmiles(use_int_bbids=True, prefix=prefix)
        self.log(
            f"Finished step {inspect.stack()[0][3]} for cycle {self.cycle.cycle_name}. Current split size  = {self.cycle.split_size}",
            "info",
        )
        self.next(self.join)

    @step
    def join(self, inputs):
        self.merge_artifacts(
            inputs, include=["library", "input_filepath", "output_directory"]
        )
        cycle_list = [input.cycle for input in inputs]
        cycle_list.sort(key=lambda x: x.cycle_name)
        self.library.cycles = tuple(cycle_list)
        cycle_string = ""  # init
        expected_size = 1  # init
        for cycle in self.library.cycles:
            cycle_string = cycle_string + cycle.cycle_name
            expected_size *= cycle.split_size
        self.output_filepath = f"{self.output_directory}{self.library.deck}/{len(self.library.cycles)}-cycle/{self.library_name}/Instances/{cycle_string}.parquet"
        self.log(
            f"Finished step {inspect.stack()[0][3]}. \n\nCycles: {cycle_string} \nExpected size: {expected_size:,} \nOutput filepath: {self.output_filepath}\n",
            "info",
        )
        self.next(self.enumerate_and_write)

    @step
    def enumerate_and_write(self):
        self.library.enumerate_and_write_instances(
            self.output_filepath, canonicalize_smiles=self.canonicalize_smiles
        )
        self.log(f"Finished step {inspect.stack()[0][3]}", "info")
        self.next(self.end)

    @step
    def end(self):
        self.log("Finished!", "info")


if __name__ == "__main__":
    XSMILES_Instance_Enumeration()
