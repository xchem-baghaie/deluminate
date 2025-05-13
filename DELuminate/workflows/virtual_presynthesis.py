import inspect

from metaflow import FlowSpec, Parameter, step
import pandas as pd

from configs.inputs import ReactionConfig, ReactionBranchConfig, ValidationStepConfig, CycleConfig, LibraryConfig
from configs.property_profiles import (
    MAX_ALLOWED_MONOSYNTHON_MW_NON_BRO5,
    MAX_ALLOWED_MONOSYNTHON_MW_BRO5,
)
from configs.library_definition import (
    BRO5_LIBRARIES,
    SYMMETRICAL_MONO_BOC_TRIAMINE_LIBRARY_IDS,
)
from tools.directory_handling import create_path
from transforms.rdkit_mols import parallel_map_compute_mol_from_smiles
from xviewer_queries.library_definition import define_library
from xviewer_queries.functional_groups import (
    generate_dict_of_functional_groups,
    query_functional_groups_by_first_reaction_in_reaction_branch,
)

import logging


class Virtual_Library_Enumeration(FlowSpec):
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

    compute_xsmiles = Parameter(
        "compute_xsmiles",
        help="Whether to generate and QC building block XSMILES strings",
        default=True,
    )

    xsmiles_directory = Parameter(
        "xsmiles_directory",
        help="The directory in which the building block XSMILES files for further review, and ultimately, upload to X-Viewer will be written",
        default="/mnt/d/XSMILES/BB_XSMILES/Virtual/",
    )

    write_xsmiles = Parameter(
        "write_xsmiles",
        help="Whether to write XSMILES lists to the XSMILES directory (overwriting existing files for the library, if they exist)",
        default=True,
    )

    bb_list_by_reaction_dir = Parameter(
        "bb_list_by_reaction_dir",
        help="Where to find the list of compatible building blocks by reaction",
        default="/mnt/p/discovery chemistry/people/ryan/tech dev - cheminformatics/enumeration/compatible_bbs_by_reaction/predictions/",
    )

    use_xbbs_only = Parameter(
        "use_xbbs_only",
        help="Whether to only consider XBBs as candidate BBs",
        default=False,
    )

    apply_mw_filter = Parameter(
        "apply_mw_filter",
        help='Whether to apply a molecular weight filter on the monosynthons',
        default = True,
    )

    prohibit_exact_overlap_in_instance_enum = Parameter(
        "prohibit_exact_overlap_in_instance_enum",
        help = "Whether to prohibit exact overlap in instance enumerations",
        default = True,
    )

    use_xviewer_connection = Parameter(
        "use_xviewer_connection",
        help = "Whether to use the XViewer connection. If False, a static export will be used.",
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
            f"/mnt/p/Discovery Chemistry/Library Synthesis/Enumeration/Library Enumeration/{self.library_name}/virtual/input/"
        )
        if self.use_xbbs_only:
            self.output_directory = create_path(
                f"/mnt/p/Discovery Chemistry/Library Synthesis/Enumeration/Library Enumeration/{self.library_name}/virtual/output/XBBs_only"
            )
        else:
            self.output_directory = create_path(
                f"/mnt/p/Discovery Chemistry/Library Synthesis/Enumeration/Library Enumeration/{self.library_name}/virtual/output/"
            )            
        self.input_filepath = self.input_directory + self.library_name + "_virtual.xlsx"
        logfile = self.generate_logfile_path()
        with open(logfile, "w") as _:
            pass  # creates file, or overwrites it if it exists already
        self.library = define_library(
            self.library_name,
            self.input_filepath,
            self.compute_xsmiles,
            is_presynthesis=True,
            is_virtual=True,
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
            
            Filtering for deck in this cycle: False
            Perform a diversity selection for this cycle (manually as post-processing -- not done automatically): {cycle.perform_diversity_selection}
            Target number of BBs: {cycle.target_number_of_bbs}            
            Operating off of a user-defined BB list: {contains_included_bbs}
            Filtering user-defined BBs out of the BB list for this cycle: {contains_excluded_bbs}
            Placeholder XSMILES Precursor: {cycle.placeholder_xsmiles_precursor}
            Placeholder XSMILES: {cycle.placeholder_xsmiles}

            """
            for reaction_branch in cycle.reaction_branches:
                string = string + construct_reaction_branch_str(reaction_branch)
            return string

        def generate_library_initialization_report(library: LibraryConfig) -> str:
            string = f"""
            All building blocks will be subjected to the full suite of filters, including those based on monosynthon enumeration.
            
            {self.library_name} was initialized successfully!
            Using XViewer connection: {self.use_xviewer_connection}
            Synonym: {self.library.synonym}
            Deck: {self.library.deck}
            Linker Encoding Tag Set: {self.library.linker_cycle.encoding_tag_set}
            Linker XBBIDs and SMILES: {self.library.linker_cycle.xbbid_smiles_dict}
            Linker Default SMILES: {self.library.linker_cycle.default_smiles}
            Cycles: {str(tuple([cycle.cycle_name for cycle in self.library.cycles]))}

            {str(self.num_instances)} instances will be enumerated.
            Generate physiochemical property profiles: {str(self.generate_physchem_property_profiles)}
            
            """
            for cycle in library.cycles:
                string = string + construct_cycle_str(cycle)

            return string

        library_initialization_report = generate_library_initialization_report(
            self.library
        )
        self.log(library_initialization_report, "info")
        self.next(self.initialize_linker_cycle)
        

    @step
    def initialize_linker_cycle(self):
        self.library.linker_cycle.initialize_bb_list()
        self.log(
            f"Finished step {inspect.stack()[0][3]}.",
            "info",
        )
        self.cycles = self.library.cycles
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
            compute_xsmiles = self.compute_xsmiles,
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

        assert self.cycle.forced_inclusion_bb_dict is None

        max_num_reactions = max(
            [len(branch.reaction_sequence) for branch in self.cycle.reaction_branches]
        )
        df_list = []
        if self.use_xbbs_only:
            counter = 1
            synthesized_library = define_library(
                self.library_name,
                filename=f"/mnt/p/Discovery Chemistry/Library Synthesis/Enumeration/Library Enumeration/{self.library_name}/postsynthesis/input/{self.library_name}_postsynthesis.xlsx",
                compute_xsmiles=False,
                is_presynthesis=False,
                is_virtual=False,
            )
            synthesized_library_rbs = [cycle for cycle in synthesized_library.cycles if cycle.cycle_name == self.cycle.cycle_name][0].reaction_branches
            synthesized_library_reaction_sequences = [rb.reaction_sequence for rb in synthesized_library_rbs]
            synthesized_library_reaction_name_lists = [[rxn.rxn_name for rxn in rxn_seq]for rxn_seq in synthesized_library_reaction_sequences]

            rbs = list(self.cycle.reaction_branches)
            rbs_to_keep = []

            for rb in rbs:
                virtual_library_reaction_name_list = [rxn.rxn_name for rxn in rb.reaction_sequence]
                if "Primary Amines" in virtual_library_reaction_name_list[0
                    ] or "Secondary Amines" in virtual_library_reaction_name_list[0
                    ] or "Suzuki" in virtual_library_reaction_name_list[0
                    ] or virtual_library_reaction_name_list in synthesized_library_reaction_name_lists:
                    setattr(rb, 'reaction_branch_idx', counter)
                    rbs_to_keep.append(rb)
                    counter+=1
            assert len(rbs_to_keep) > 0
            setattr(self.cycle, "reaction_branches", tuple(rbs_to_keep))

        for rb in self.cycle.reaction_branches:
            rb_df = pd.read_csv(
                self.bb_list_by_reaction_dir + rb.reaction_sequence[0].rxn_name + ".csv"
            ).fillna("")
            if self.use_xbbs_only:
                rb_df = rb_df[rb_df['Set']=='XBBs']
                rb_df = rb_df.reset_index(drop=True)
            # ID column name handling
            if "ID" in list(rb_df.columns):
                rb_df.rename(columns={"ID": "XBBID"}, inplace=True)

            # Populate info about the reaction branch
            rb_df["Reaction Sequence ID"] = [rb.reaction_branch_idx] * len(rb_df)

            for i in range(max_num_reactions):
                rb_df["Reaction " + str(i + 1)] = [""] * len(rb_df)
            for j in range(len(rb.reaction_sequence)):
                rb_df["Reaction " + str(j + 1)] = [
                    rb.reaction_sequence[j].rxn_name
                ] * len(rb_df)

            bbs = parallel_map_compute_mol_from_smiles(list(rb_df["SMILES"]))
            rb_df = query_functional_groups_by_first_reaction_in_reaction_branch(
                rb,
                rb_df,
                bbs,
                overwrite_existing_reaction_assignments=False,
                forced_inclusion_bb_dict=self.cycle.forced_inclusion_bb_dict,
            )
            assert len(rb_df) > 0
            df_list.append(rb_df)

        self.cycle.working_df = pd.concat(df_list, axis=0)
        self.cycle.working_df = self.cycle.working_df.reset_index(drop=True)
        self.cycle.split_size = len(self.cycle.working_df)
        self.cycle.bbs = parallel_map_compute_mol_from_smiles(
            list(self.cycle.working_df["SMILES"])
        )

        # Now exclude any BBs which the user wanted to force exclusion for by their IDs
        if self.cycle.forced_exclusion_bbs is not None:
            for i in range(len(self.cycle.working_df)):
                if self.cycle.working_df["XBBID"][i] in self.cycle.forced_exclusion_bbs:
                    self.cycle.working_df.loc[i, "Excluded By User (XBBID)"] = "Yes"
                    self.cycle.working_df.loc[i, "Include"] = "No"
        self.cycle.filtered_OUT_df = None
        self.cycle.update_working_df()

        self.log(
            f"Finished step {inspect.stack()[0][3]} for cycle {self.cycle.cycle_name}. Current split size  = {self.cycle.split_size}",
            "info",
        )
        self.next(self.filter_validation)

    @step
    def filter_validation(self):
        self.cycle.filter_validation_for_virtual_libraries()
        self.log(
            f"Finished step {inspect.stack()[0][3]} for cycle {self.cycle.cycle_name}. Current split size  = {self.cycle.split_size}",
            "info",
        )
        self.next(self.convert_smiles_to_mols)

    @step
    def convert_smiles_to_mols(self):
        self.cycle.convert_smiles_to_mols()
        self.log(
            f"Finished step {inspect.stack()[0][3]} for cycle {self.cycle.cycle_name}. Current split size  = {self.cycle.split_size}",
            "info",
        )
        self.next(self.filter_reactivity)

    @step
    def filter_reactivity(self):
        override_multiple_substructure_matches = (
            self.library_name in SYMMETRICAL_MONO_BOC_TRIAMINE_LIBRARY_IDS
            and self.cycle.cycle_name == "A"
        )  # bool
        self.cycle.filter_reactivity(
            override_multiple_substructure_matches=override_multiple_substructure_matches
        )
        self.log(
            f"Finished step {inspect.stack()[0][3]} for cycle {self.cycle.cycle_name}. Current split size  = {self.cycle.split_size}",
            "info",
        )
        self.next(self.monosynthon_enumeration)

    @step
    def monosynthon_enumeration(self):
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
        self.cycle.monosynthon_enumeration(
            self.library.linker_cycle.default_smiles,
            cycle_names,
            placeholder_reaction_branches,
            placeholder_bb_exemplars,
            self.compute_xsmiles,
            placeholder_xsmiles,
            linker_cycle_encoded = self.library.linker_cycle.encoding_tag_set is not None,
            linker_placeholder_xsmiles = self.library.linker_cycle.placeholder_xsmiles
        )
        self.log(
            f"Finished step {inspect.stack()[0][3]} for cycle {self.cycle.cycle_name}. Current split size  = {self.cycle.split_size}",
            "info",
        )
        self.next(self.filter_duplicate_enums)

    @step
    def filter_duplicate_enums(self):
        self.cycle.filter_duplicate_enums()
        self.log(
            f"Finished step {inspect.stack()[0][3]} for cycle {self.cycle.cycle_name}. Current split size  = {self.cycle.split_size}",
            "info",
        )
        self.next(self.filter_undesirable_substructures)

    @step
    def filter_undesirable_substructures(self):
        fgs_to_search = [
            fg
            for fg in generate_dict_of_functional_groups().values()
            if fg.compound_proposal_rule in ["Warn", "Alert"]
        ]
        self.cycle.filter_undesirable_substructures(fgs_to_search)
        self.log(
            f"Finished step {inspect.stack()[0][3]} for cycle {self.cycle.cycle_name}. Current split size  = {self.cycle.split_size}",
            "info",
        )
        self.next(self.filter_duplicate_ids)

    @step
    def filter_duplicate_ids(self):
        self.cycle.filter_duplicate_ids()
        self.log(
            f"Finished step {inspect.stack()[0][3]} for cycle {self.cycle.cycle_name}. Current split size  = {self.cycle.split_size}",
            "info",
        )
        self.next(self.filter_invalid_xsmiles)

    @step
    def filter_invalid_xsmiles(self):
        self.cycle.filter_invalid_xsmiles()
        self.log(
            f"Finished step {inspect.stack()[0][3]} for cycle {self.cycle.cycle_name}. Current split size  = {self.cycle.split_size}",
            "info",
        )
        self.next(self.filter_molecular_weight)

    @step
    def filter_molecular_weight(self):
        if self.apply_mw_filter:
            if self.library_name in BRO5_LIBRARIES:
                self.cycle.filter_mw(MAX_ALLOWED_MONOSYNTHON_MW_BRO5)
            else:
                self.cycle.filter_mw(MAX_ALLOWED_MONOSYNTHON_MW_NON_BRO5)
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
        if self.prohibit_exact_overlap_in_instance_enum:
            self.library.enumerate_virtual_instances(self.num_instances)
        else:
            self.library.enumerate_instances(self.num_instances)
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
            self.library.plot_physchem_property_profiles(
                property_profile_directory, keep_autogenerated_property_columns=True
            )
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
            + "_virtual_output.xlsx",
            self.output_directory
            + self.library.lib_name.replace("_postsynthesis", "").replace(
                "_presynthesis", ""
            )
            + "_virtual_QC_set.csv",            
            color_list_tabs=False,
            is_presynthesis=True,
        )
        self.log(f"Finished step {inspect.stack()[0][3]}\n", "info")
        self.next(self.write_xsmiles_lists)

    @step
    def write_xsmiles_lists(self):
        if self.compute_xsmiles and self.write_xsmiles:
            xsmiles_directory = create_path(
                f"{self.xsmiles_directory}{self.library_name}/"
            )
            self.library.write_xsmiles_lists(xsmiles_directory)
            self.log(f"Finished step {inspect.stack()[0][3]}\n", "info")
        self.next(self.end)

    @step
    def end(self):
        split_size_ls = [cycle.split_size for cycle in self.library.cycles]
        total_size = 1  # init
        for split_size in split_size_ls:
            total_size *= split_size
        self.log(f"Finished! Max library size = product of {split_size_ls} = {total_size:,}","info")


if __name__ == "__main__":
    Virtual_Library_Enumeration()
