import inspect
import logging

from metaflow import FlowSpec, Parameter, step

from configs.inputs import ReactionConfig, ReactionBranchConfig, ValidationStepConfig, CycleConfig, LibraryConfig
from tools.directory_handling import create_path
from xviewer_queries.library_definition import define_library


class Enumerate_Instance_Exemplars(FlowSpec):
    library_name = Parameter(
        "library_name",
        help="The XLIB or XPL number of the library. For example, XLIB0001 or XPL0002",
        default="",
    )

    num_instances = Parameter(
        "num_instances",
        help="The number of library instances to enumerate",
        default=3000,
    )

    generate_physchem_property_profiles = Parameter(
        "generate_physchem_property_profiles",
        help="Whether to calculate and generate physicochemical property profiles for the enumerated instances",
        default=True,
    )

    enumerate_qc_set = Parameter(
        "enumerate_qc_set",
        help="Whether to enumerate a QC set, in which every building block is used in at least one instance. This is run in addition to the main instance enumeration, which will enumerate num_instances",
        default=True,
    )

    use_xviewer_connection = Parameter(
        "use_xviewer_connection",
        help = "Whether to use the XViewer connection. If False, a static export will be used.",
        default = True,
    )

    add_substructure_warnings_to_instance_enumerations = Parameter(
        "add_substructure_warnings_to_instance_enumerations",
        help = "Whether to add substructure warnings to the instance enumerations",
        default = True,
    )

    def generate_logfile_path(self):
        return self.output_directory + self.library_name + "_log.txt"

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
        self.output_directory = create_path(
            f"/mnt/p/Discovery Chemistry/Library Synthesis/Enumeration/Library Enumeration/{self.library_name}/postsynthesis/output/"
        )
        self.input_filepath = (
            self.input_directory + self.library_name + "_postsynthesis.xlsx"
        )
        logfile = self.generate_logfile_path()
        with open(logfile, "w") as _:
            pass  # creates file, or overwrites it if it exists already
        self.library = define_library(
            self.library_name,
            self.input_filepath,
            compute_xsmiles=False,
            is_presynthesis=False,
            use_xviewer_connection=self.use_xviewer_connection,
        )

        def construct_reaction_str(reaction: ReactionConfig) -> str:
            if reaction.allowed_fgs_reactant1 is not None:
                allowed_fgs_reactant1 = [
                    fg.fg_name for fg in reaction.allowed_fgs_reactant1
                ]
            else:
                allowed_fgs_reactant1 = None
            if reaction.forbidden_fgs_reactant1 is not None:
                forbidden_fgs_reactant1 = [
                    fg.fg_name for fg in reaction.forbidden_fgs_reactant1
                ]
            else:
                forbidden_fgs_reactant1 = None
            if reaction.allowed_fgs_reactant2 is not None:
                allowed_fgs_reactant2 = [
                    fg.fg_name for fg in reaction.allowed_fgs_reactant2
                ]
            else:
                allowed_fgs_reactant2 = None
            if reaction.forbidden_fgs_reactant2 is not None:
                forbidden_fgs_reactant2 = [
                    fg.fg_name for fg in reaction.forbidden_fgs_reactant2
                ]
            else:
                forbidden_fgs_reactant2 = None
            if reaction.manually_added_fgs_reactant2 is not None:
                manually_added_fgs_reactant2 = [
                    fg.fg_name for fg in reaction.manually_added_fgs_reactant2
                ]
            else:
                manually_added_fgs_reactant2 = None
            if reaction.manually_excluded_fgs_reactant2 is not None:
                manually_excluded_fgs_reactant2 = [
                    fg.fg_name for fg in reaction.manually_excluded_fgs_reactant2
                ]
            else:
                manually_excluded_fgs_reactant2 = None
            if reaction.product_fgs is not None:
                product_fgs = [fg.fg_name for fg in reaction.product_fgs]
            else:
                product_fgs = None

            return f"""
            ----- {reaction.rxn_name} -----

            Branch {reaction.branch_step[0]}, step {reaction.branch_step[1]}
            SMARTS: {reaction.smarts}
            Number of reactants: {reaction.num_reactants}
            Intermediate isolated: {reaction.intermediate_isolated}
            Example Reactant 1: {reaction.example_reactant1}
            Example Reactant 2: {reaction.example_reactant2}
            Example Product: {reaction.example_product}
            Tolerated FGs in Reactant 1: {allowed_fgs_reactant1}
            Untolerated FGs in Reactant 1: {forbidden_fgs_reactant1}
            Tolerated FGs in Reactant 2: {allowed_fgs_reactant2}
            Untolerated FGs in Reactant 2: {forbidden_fgs_reactant2}
            Product FGs: {product_fgs}
            Manually Added FGs in Reactant 2: {manually_added_fgs_reactant2}
            Manually Excluded FGs from Reactant 2: {manually_excluded_fgs_reactant2}
            
            """

        def construct_validation_step_str(validation_step: ValidationStepConfig) -> str:
            return f"""
            _____ Validation Step {validation_step.validation_step_idx} _____

            Assays: {validation_step.assay_names}
            Assay to use when results are equal: {validation_step.assay_to_use_when_equal}
            Cutoff: {validation_step.validation_cutoff}
            Include unvalidated BBs: {validation_step.include_unvalidated_bbs}

            """

        def construct_reaction_branch_str(reaction_branch: ReactionBranchConfig) -> str:
            string = f"""
            @@@@@ Reaction Sequence ID {reaction_branch.reaction_branch_idx} @@@@@
            
            Monosynthon exemplar: {reaction_branch.monosynthon_exemplar}
            
            """
            for reaction in reaction_branch.reaction_sequence:
                string = string + construct_reaction_str(reaction)
            if reaction_branch.validation_sequence is not None:
                for validation_step in reaction_branch.validation_sequence:
                    string = string + construct_validation_step_str(validation_step)
            return string

        def construct_cycle_str(cycle: CycleConfig) -> str:
            contains_included_bbs = cycle.forced_inclusion_bb_dict is not None
            contains_excluded_bbs = cycle.forced_exclusion_bbs is not None

            string = f"""
            ##### Cycle {cycle.cycle_name} #####
            
            Operating off of a user-defined BB list: {contains_included_bbs}
            Filtering user-defined BBs out of the BB list for this cycle: {contains_excluded_bbs}
            
            """
            for reaction_branch in cycle.reaction_branches:
                string = string + construct_reaction_branch_str(reaction_branch)
            return string

        def generate_library_initialization_report(library: LibraryConfig) -> str:
            string = f"""
            All building blocks will be subjected to filtering by monosynthon enumeration ONLY.
            
            {self.library_name} was initialized successfully!
            Using XViewer connection: {self.use_xviewer_connection}
            Synonym: {self.library.synonym}
            Deck: {self.library.deck}
            Linker Encoding Tag Set: {self.library.linker_cycle.encoding_tag_set}
            Linker XBBIDs and SMILES: {self.library.linker_cycle.xbbid_smiles_dict}
            Linker Default SMILES: {self.library.linker_cycle.default_smiles}
            Cycles: {str(tuple([cycle.cycle_name for cycle in self.library.cycles]))}

            """
            for cycle in library.cycles:
                string = string + construct_cycle_str(cycle)

            return string

        library_initialization_report = generate_library_initialization_report(
            self.library
        )
        self.log(library_initialization_report, "info")
        self.cycles = self.library.cycles
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
        self.next(self.join)

    @step
    def join(self, inputs):
        self.merge_artifacts(
            inputs, include=["library", "input_filepath", "output_directory"]
        )
        cycle_list = [input.cycle for input in inputs]
        cycle_list.sort(key=lambda x: x.cycle_name)
        self.library.cycles = tuple(cycle_list)
        self.log(f"Finished step {inspect.stack()[0][3]}", "info")
        self.next(self.instance_enumeration)

    @step
    def instance_enumeration(self):
        self.library.enumerate_instances(self.num_instances, add_substructure_warnings=self.add_substructure_warnings_to_instance_enumerations)
        if self.enumerate_qc_set:
            self.library.enumerate_qc_set()
        self.log(f"Finished step {inspect.stack()[0][3]}", "info")
        self.next(self.calculate_physchem_properties)

    @step
    def calculate_physchem_properties(self):
        if self.generate_physchem_property_profiles:
            self.library.calculate_physchem_properties()
            self.log(f"Finished step {inspect.stack()[0][3]}", "info")
        self.next(self.plot_physchem_property_profiles)

    @step
    def plot_physchem_property_profiles(self):
        if self.generate_physchem_property_profiles:
            property_profile_directory = create_path(
                f"{self.output_directory}/property_profiles/"
            )
            self.library.plot_physchem_property_profiles(property_profile_directory)
            self.log(f"Finished step {inspect.stack()[0][3]}", "info")
        self.next(self.write_lists)

    @step
    def write_lists(self):
        self.library.write_excel_file(
            self.input_filepath,
            self.output_directory
            + self.library.lib_name.replace("_postsynthesis", "").replace(
                "_presynthesis", ""
            )
            + "_output.xlsx",
            color_list_tabs=True,
            is_presynthesis=False,
        )
        self.log(f"Finished step {inspect.stack()[0][3]}\n", "info")
        self.next(self.end)

    @step
    def end(self):
        self.log("Finished!", "info")


if __name__ == "__main__":
    Enumerate_Instance_Exemplars()
