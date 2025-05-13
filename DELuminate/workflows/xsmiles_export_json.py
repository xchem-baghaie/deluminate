from copy import deepcopy
import json

from metaflow import FlowSpec, Parameter, step

from configs.xsmiles import CONNECTIVITY_HANDLES
from tools.directory_handling import create_path
from transforms.rdkit_mols import compute_mol_from_smiles, compute_smiles_from_mol
from xviewer_queries.library_definition import define_library
from xviewer_queries.xsmiles import query_xsmiles


class XSMILES_Export_JSON(FlowSpec):
    library_name = Parameter(
        "library_name",
        help="The XLIB or XPL number of the library. For example, XLIB0001 or XPL0002",
        default="",
    )

    output_directory = Parameter(
        "output_directory",
        help="The output directory in which the JSON file will be exported to",
        default='../Thompson_Sampling_DELs/data/X-Chem/',
    )

    use_xviewer_connection = Parameter(
        "use_xviewer_connection",
        help = "Whether to use the XViewer connection. If False, a static export will be used.",
        default = True,
    )

    @step
    def start(self):
        self.input_directory = create_path(
            f"/mnt/p/Discovery Chemistry/Library Synthesis/Enumeration/Library Enumeration/{self.library_name}/postsynthesis/input/"
        )
        self.input_filepath = (
            self.input_directory + self.library_name + "_postsynthesis.xlsx"
        )

        self.library = define_library(
            self.library_name,
            self.input_filepath,
            compute_xsmiles=True,
            is_presynthesis=False,
            use_xviewer_connection=self.use_xviewer_connection,
        )

        self.cycles = self.library.cycles
        self.library.linker_cycle.initialize_bb_list()
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
        self.next(self.process_cycles, foreach="cycles")

    @step
    def process_cycles(self):
        self.cycle = self.input
        if self.cycle.forced_inclusion_bb_dict is None:
            raise ValueError(
                f"ERROR: BB list was not provided for cycle {self.cycle.cycle_name}"
            )
        self.cycle.initialize_bb_list()
        self.cycle.working_df["Include"] = ["Yes"] * len(self.cycle.working_df)
        self.cycle.lookup_bb_xsmiles(use_int_bbids=False, prefix="")
        self.next(self.join)

    @step
    def join(self, inputs):
        self.merge_artifacts(
            inputs, include=["library"]
        )
        cycle_list = [input.cycle for input in inputs]
        cycle_list.sort(key=lambda x: x.cycle_name)
        self.library.cycles = tuple(cycle_list)

        xsmiles_closure_element_pattern_dict = {value:key for key, value in CONNECTIVITY_HANDLES.items()}

        xsmiles_dict = {self.library_name: {}} #init
        
        #Determine whether there are linker XSMILES which will not be captured just by looping through the other cycles
        if "A" not in [cycle.cycle_name for cycle in self.library.cycles] and self.library.linker_cycle.encoding_tag_set is not None:
            xsmiles_dict[self.library_name]["A"] = {}
            linker_xbbid_smiles_dict = {self.library.linker_cycle.df['XBBID'][i]: self.library.linker_cycle.df['SMILES'][i] for i in range(len(self.library.linker_cycle.df))}
            cyc_a_xsmiles_df = query_xsmiles(self.library_name, "A", self.library.xviewer_connection)
            for i in range(len(cyc_a_xsmiles_df)):
                xbbid = cyc_a_xsmiles_df['XBBID'][i]
                smiles = compute_smiles_from_mol(compute_mol_from_smiles(self.library.linker_cycle.xbbid_smiles_dict[xbbid]))
                linker_xsmiles= cyc_a_xsmiles_df['XSMILES_ORIG'][i]
                linker_xsmiles_instruction = cyc_a_xsmiles_df['INSTRUCTION'][i]
                status = cyc_a_xsmiles_df['XSMILES_STATUS'][i]
                linker_xsmiles_precursor = deepcopy(linker_xsmiles)
                for pattern, elements in xsmiles_closure_element_pattern_dict.items():
                    if linker_xsmiles_precursor.count(pattern) == 1:
                        # If there are multiple of the pattern in the XSMILES, the connection is already handled, so should not be replaced.
                        linker_xsmiles_precursor = linker_xsmiles_precursor.replace(f"-{pattern}", f"(-{elements})")
                        linker_xsmiles_precursor = linker_xsmiles_precursor.replace(f":{pattern}", f"(:{elements})")
                        linker_xsmiles_precursor = linker_xsmiles_precursor.replace(f"={pattern}", f"(={elements})")
                        linker_xsmiles_precursor = linker_xsmiles_precursor.replace(f"#{pattern}", f"(#{elements})")
                        linker_xsmiles_precursor = linker_xsmiles_precursor.replace(f"/{pattern}", f"(/{elements})")
                        linker_xsmiles_precursor = linker_xsmiles_precursor.replace(f"\{pattern}", f"(\{elements})")
                        linker_xsmiles_precursor = linker_xsmiles_precursor.replace(f"\\{pattern}", f"(\\{elements})")
                        linker_xsmiles_precursor = linker_xsmiles_precursor.replace(f"\\\{pattern}", f"(\\\{elements})")
                        linker_xsmiles_precursor = linker_xsmiles_precursor.replace(f"\\\\{pattern}", f"(\\\\{elements})")
                        linker_xsmiles_precursor = linker_xsmiles_precursor.replace(pattern, f'({elements})')
                if len(linker_xsmiles_precursor) > 0:
                    #Replace the trailing period if it exists
                    if linker_xsmiles_precursor[-1] == '.':
                        linker_xsmiles_precursor = linker_xsmiles_precursor[:-1]
                linker_xsmiles_precursor = compute_smiles_from_mol(compute_mol_from_smiles(linker_xsmiles_precursor))
                if linker_xsmiles != '' and linker_xsmiles_precursor == '':
                    raise Exception(f'XSMILES precursor generation failed for {self.library_name} cycle A {xbbid} with XSMILES "{linker_xsmiles}", instruction "{linker_xsmiles_instruction}" and SMILES "{smiles}"')
                if linker_xsmiles_instruction == "FAIL<>" or xbbid == "" or status != "VALID":
                    continue
                xsmiles_dict[self.library_name]["A"][xbbid] = {
                    "SMILES": smiles,
                    "XVIEWER_SMILES_CANONICALIZED": compute_smiles_from_mol(compute_mol_from_smiles(linker_xbbid_smiles_dict[xbbid])),
                    "XSMILES": linker_xsmiles,
                    "XSMILES_PRECURSOR": linker_xsmiles_precursor,
                    "XSMILES_INSTRUCTION": linker_xsmiles_instruction,
                    "REACTION_SEQUENCE_ID": 1,
                    "REACTION_SEQUENCE": []
                }
        
        for cycle in self.library.cycles:
            xsmiles_dict[self.library_name][cycle.cycle_name] = {}
            xbbid_smiles_dict = {cycle.working_df['XBBID'][i]: cycle.working_df['SMILES'][i] for i in range(len(cycle.working_df))}
            xbbid_xsmiles_dict = cycle.xsmiles_dict_with_instruction | cycle.xsmiles_dict_no_instruction #Merge dictionaries
            #use xbbid_xsmiles_dict as the keys to iterate through because BBs with an instruction of FAIL<> are not included, but will be if xbbid_smiles_dict is used.
            reaction_information_dict = {}
            for reaction_branch in cycle.reaction_branches:
                reaction_information_dict[reaction_branch.reaction_branch_idx] = []
                for reaction in reaction_branch.reaction_sequence:
                    reaction_information_dict[reaction_branch.reaction_branch_idx].append(
                        {
                            "STEP": reaction.branch_step[1],
                            "REACTION_NAME": reaction.rxn_name,
                            "REACTION_SMARTS": reaction.smarts,
                        }
                    )
            for xbbid in xbbid_xsmiles_dict.keys():
                row_index = cycle.working_df.index[cycle.working_df['XBBID'] == xbbid].tolist()[0]
                reaction_sequence_id = cycle.working_df['Reaction Sequence ID'][row_index]
                xsmiles = xbbid_xsmiles_dict[xbbid][0]
                xsmiles_instruction = xbbid_xsmiles_dict[xbbid][1]
                xsmiles_precursor = deepcopy(xsmiles)
                for pattern, elements in xsmiles_closure_element_pattern_dict.items():
                    if xsmiles_precursor.count(pattern) == 1:
                        # If there are multiple of the pattern in the XSMILES, the connection is already handled, so should not be replaced.
                        xsmiles_precursor = xsmiles_precursor.replace(f"-{pattern}", f"(-{elements})")
                        xsmiles_precursor = xsmiles_precursor.replace(f":{pattern}", f"(:{elements})")
                        xsmiles_precursor = xsmiles_precursor.replace(f"={pattern}", f"(={elements})")
                        xsmiles_precursor = xsmiles_precursor.replace(f"#{pattern}", f"(#{elements})")
                        xsmiles_precursor = xsmiles_precursor.replace(f"/{pattern}", f"(/{elements})")
                        xsmiles_precursor = xsmiles_precursor.replace(f"\{pattern}", f"(\{elements})")
                        xsmiles_precursor = xsmiles_precursor.replace(f"\\{pattern}", f"(\\{elements})")
                        xsmiles_precursor = xsmiles_precursor.replace(f"\\\{pattern}", f"(\\\{elements})")
                        xsmiles_precursor = xsmiles_precursor.replace(f"\\\\{pattern}", f"(\\\\{elements})")
                        xsmiles_precursor = xsmiles_precursor.replace(pattern, f'({elements})')
                if len(xsmiles_precursor) > 0:
                    #Replace the trailing period if it exists
                    if xsmiles_precursor[-1] == '.':
                        xsmiles_precursor = xsmiles_precursor[:-1]
                xsmiles_precursor = compute_smiles_from_mol(compute_mol_from_smiles(xsmiles_precursor))
                if xsmiles !='' and xsmiles_precursor == '':
                    raise Exception(f'XSMILES precursor generation failed for {self.library_name} cycle {cycle.cycle_name} {xbbid} with XSMILES "{xsmiles}", instruction "{xsmiles_instruction}" and SMILES "{smiles}"')
                xsmiles_dict[self.library_name][cycle.cycle_name][xbbid] = {
                    "SMILES": compute_smiles_from_mol(compute_mol_from_smiles(xbbid_smiles_dict[xbbid])),
                    "XSMILES": xsmiles,
                    "XSMILES_PRECURSOR": xsmiles_precursor,
                    "XSMILES_INSTRUCTION": xsmiles_instruction,
                    "REACTION_SEQUENCE_ID": reaction_sequence_id,
                    "REACTION_SEQUENCE": reaction_information_dict[reaction_sequence_id]
                }
        with open(f'{self.output_directory}{self.library_name}.json', 'w') as f:
            json.dump(xsmiles_dict, f, indent=4)
        self.next(self.end)

    @step
    def end(self):
        pass

if __name__ == "__main__":
    XSMILES_Export_JSON()
