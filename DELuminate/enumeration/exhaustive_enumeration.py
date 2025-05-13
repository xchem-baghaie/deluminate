from itertools import chain
from typing import List, Tuple

from tools.parallel_map import parallel_map_generic
from transforms.xsmiles import parse_xsmiles_instructions


def concatenate_xsmiles_no_instruction(xsmiles_strs: List[str]) -> str:
    """
    Args:
    _____
    -strs: A list of strings, which will be joined by concatenation
    """
    return "".join(xsmiles_strs)


def concatenate_xsmiles_with_instruction(
    xsmiles: List[List[str]], cycle_name_list: List[str]
) -> str:
    """
    Args:
    _____
    -xsmiles: a 2-tuple, first element is a list of XSMILES strings, second element is a list of instruction strings
    -cycle_name_list: a list of strings corresponding to the names of the cycles used in this library.

    next_cycle below is used for the implementation of CONCATENATE_AND_SKIP_TO_CYCLE and SKIP_TO_CYCLE.
    Start with next_cycle = the first cycle
    Only allow the enumeration to proceed if the current cycle matches the next cycle from the previous iteration.
    If the CONCATENATE_AND_SKIP_TO_CYCLE or SKIP_TO_CYCLE instructions are encountered, next_cycle is set to that instruction's argument,
    to prohibit the intermediate cycles' enumerations from occurring. It also prevents the application of any instructions following
    CONCATENATE_AND_SKIP_TO_CYCLE or SKIP_TO_CYCLE to curr_smiles.
    """
    curr_xsmiles = ""
    xsmiles_strs = xsmiles[0]
    instructions = xsmiles[1]
    next_cycle = cycle_name_list[0]
    for i in range(len(cycle_name_list)):
        if instructions[i] == "FAIL<>":
            return ""
        else:
            curr_cycle = cycle_name_list[i]
            if curr_cycle == next_cycle:
                parsed_instructions, parsed_args = parse_xsmiles_instructions(
                    instructions[i]
                )
                if parsed_instructions is not None:
                    # Resolve the string replacement on the current XSMILES before concatenating
                    if (
                        parsed_instructions[0] == "REPLACE_PREV_END"
                    ):  # First instruction string
                        args = parsed_args[0]  # First set of arguments
                        n_chars = len(args[0])
                        if curr_xsmiles[-n_chars:] == args[0]:
                            curr_xsmiles = curr_xsmiles[:-n_chars] + args[1]
                    else:
                        pass

                    # Concatenate
                    curr_xsmiles = curr_xsmiles + xsmiles_strs[i]

                    # Resolve remaining instructions
                    for j in range(len(parsed_instructions)):
                        instruction = parsed_instructions[j]  # Instruction string
                        args = parsed_args[j]  # Instruction arguments

                        if instruction == "REPLACE":
                            curr_xsmiles = curr_xsmiles.replace(args[0], args[1])
                        elif instruction == "REPLACE_PREV_END":
                            if j == 0:
                                continue  # Already resolved
                            else:
                                raise Exception(
                                    f"REPLACE_PREV_END is only supported if it is the first instruction, but it was used in position {j+1}"
                                )
                        elif instruction == "CONCATENATE_AND_STOP":
                            return curr_xsmiles
                        elif instruction in [
                            "CONCATENATE_AND_SKIP_TO_CYCLE",
                            "SKIP_TO_CYCLE",
                        ]:
                            next_cycle = args[0]
                        else:
                            raise Exception(
                                f"{instruction} is not implemented as an instruction"
                            )
                else:
                    # Concatenate only
                    curr_xsmiles = curr_xsmiles + xsmiles_strs[i]
                if curr_cycle == next_cycle and curr_cycle != cycle_name_list[-1]:
                    next_cycle = cycle_name_list[i + 1]
                else:
                    pass
            else:
                continue

    return curr_xsmiles


def parallel_map_concatenate_xsmiles(
    combos: List[Tuple[Tuple]],
    cycle_name_list: List[str],
    use_instruction: bool = False,
):
    """
    Args:
    _____
    -combos: List of tuple of tuples.
        The outer list has length equal to the number of combinations
        The outer tuple has length equal to the number of cycles. Each element of this tuple is one instance.
        The inner tuple has length equal to 3, and is of the form (XBBID: str, XSMILES: str, Instruction: str = '')
    -use_instruction: bool
    """
    str_list = []
    for combo in combos:
        str_list.append([elm[1] for elm in combo])
    if use_instruction:
        instructions = []
        for combo in combos:
            instructions.append([elm[2] for elm in combo])
        xsmiles = [(str_list[i], instructions[i]) for i in range(len(str_list))]
        return list(
            chain.from_iterable(
                parallel_map_generic(
                    xsmiles,
                    concatenate_xsmiles_with_instruction,
                    cycle_name_list=cycle_name_list,
                )
            )
        )
    else:
        return list(
            chain.from_iterable(
                parallel_map_generic(str_list, concatenate_xsmiles_no_instruction)
            )
        )


if __name__ == "__main__":
    xs_list = ["A.", r"B=%16%17%25%31%34=%37%40%45%46.", "C.", "D.", "E"]

    expected_results = {
        "": r"A.B=%16%17%25%31%34=%37%40%45%46.C.D.E",
        "CONCATENATE_AND_STOP<>": r"A.B=%16%17%25%31%34=%37%40%45%46.",
        "REPLACE<%34|(C#C)>CONCATENATE_AND_STOP<>": r"A.B=%16%17%25%31(C#C)=%37%40%45%46.",
        "REPLACE<%37|>CONCATENATE_AND_STOP<>": r"A.B=%16%17%25%31%34=%40%45%46.",
        "REPLACE<%37|>": r"A.B=%16%17%25%31%34=%40%45%46.C.D.E",
        "REPLACE<%37|%45>": r"A.B=%16%17%25%31%34=%45%40%45%46.C.D.E",
        "FAIL<>": r"",
        "CONCATENATE_AND_SKIP_TO_CYCLE<E>": r"A.B=%16%17%25%31%34=%37%40%45%46.E",
        "REPLACE<%45|%51>": r"A.B=%16%17%25%31%34=%37%40%51%46.C.D.E",
        "REPLACE_PREV_END<.|>REPLACE<%40|(=O)>CONCATENATE_AND_STOP<>": r"AB=%16%17%25%31%34=%37(=O)%45%46.",
        "REPLACE<%37|%51>SKIP_TO_CYCLE<E>": r"A.B=%16%17%25%31%34=%51%40%45%46.E",
        "REPLACE<%45|(=O)>CONCATENATE_AND_STOP<>": r"A.B=%16%17%25%31%34=%37%40(=O)%46.",
        "REPLACE_PREV_END<.|>REPLACE<%16|(=O)>REPLACE<%17|(=O)>CONCATENATE_AND_STOP<>": r"AB=(=O)(=O)%25%31%34=%37%40%45%46.",
        "REPLACE_PREV_END<.|>REPLACE<%45|(=O)>": r"AB=%16%17%25%31%34=%37%40(=O)%46.C.D.E",
        "REPLACE_PREV_END<.|>REPLACE<%40|>": r"AB=%16%17%25%31%34=%37%45%46.C.D.E",
        "REPLACE_PREV_END<.|>REPLACE<%16|>": r"AB=%17%25%31%34=%37%40%45%46.C.D.E",
        "REPLACE_PREV_END<.|>REPLACE<=%16|>REPLACE<%37|>": r"AB%17%25%31%34=%40%45%46.C.D.E",
        "REPLACE<%25|%37>": r"A.B=%16%17%37%31%34=%37%40%45%46.C.D.E",
        "REPLACE_PREV_END<.|>REPLACE<%45|>": r"AB=%16%17%25%31%34=%37%40%46.C.D.E",
        "REPLACE<%31|>CONCATENATE_AND_STOP<>": r"A.B=%16%17%25%34=%37%40%45%46.",
        "REPLACE<%31|>": r"A.B=%16%17%25%34=%37%40%45%46.C.D.E",
        "REPLACE_PREV_END<.|>REPLACE<=%16|>REPLACE<%17|>": r"AB%25%31%34=%37%40%45%46.C.D.E",
        "REPLACE<=%37|>REPLACE<%16|>": r"A.B=%17%25%31%34%40%45%46.C.D.E",
        "REPLACE_PREV_END<.|>REPLACE<%40|(=O)>": r"AB=%16%17%25%31%34=%37(=O)%45%46.C.D.E",
        "REPLACE<%37|(=O)>": r"A.B=%16%17%25%31%34=(=O)%40%45%46.C.D.E",
        "REPLACE_PREV_END<.|>REPLACE<%37|>CONCATENATE_AND_STOP<>": r"AB=%16%17%25%31%34=%40%45%46.",
        "REPLACE_PREV_END<.|>REPLACE<%45|(=[N+]=[N-])>REPLACE<%46|>": r"AB=%16%17%25%31%34=%37%40(=[N+]=[N-]).C.D.E",
        "REPLACE<=%16|(=O)>REPLACE<%37|>": r"A.B(=O)%17%25%31%34=%40%45%46.C.D.E",
        "REPLACE_PREV_END<.|>REPLACE<%45|(=O)>REPLACE<%46|(=O)>": r"AB=%16%17%25%31%34=%37%40(=O)(=O).C.D.E",
        "REPLACE_PREV_END<.|>REPLACE<%37|(=O)>": r"AB=%16%17%25%31%34=(=O)%40%45%46.C.D.E",
    }

    for instruction in list(expected_results.keys()):
        xs_enum = concatenate_xsmiles_with_instruction(
            (xs_list, ["", instruction, "", "", ""]), ["A", "B", "C", "D", "E"]
        )
        try:
            assert xs_enum == expected_results[instruction]
        except:
            raise Exception(
                f"Instruction {list(expected_results.keys()).index(instruction)+1} of {len(expected_results)}, {instruction}, produced unexpected behavior. Expected and actual results are displayed below\n{expected_results[instruction]}\n{xs_enum}"
            )
    print("All instructions produced the expected behavior")

    from time import time
    from statistics import mean, stdev

    for n in [
        # 1,
        # 10,
        # 100,
        # 1000,
        # 10000,
        # 100000,
        # 1000000,
        # 10000000,
        # 100000000
    ]:
        combos = [
            (
                ("XBB001324", r"C%13.N%37%13C(c1nc(C)cs1)C(F)(F)F.", ""),
                ("Plc", r"c1%45ccc(C(=O)%37)cc1.", ""),
                ("XBB011740", r"c1%45c(F)cccc1F", ""),
            )
        ] * n
        t_list = []
        for _ in range(5):
            t = time()
            rs = parallel_map_concatenate_xsmiles(combos, use_instruction=False)
            t_list.append((time() - t) / 60)
        print(
            f"n = {n}, avg runtime (5 runs) = {mean(t_list)} min, stdev runtime (5 runs) = {stdev(t_list)} min"
        )
    # (1) n = 1, avg runtime (5 runs) = 8.344650268554687e-08 min, stdev runtime (5 runs) = 4.955118462777108e-08 min
    # (10) n = 10, avg runtime (5 runs) = 1.3748804728190104e-07 min, stdev runtime (5 runs) = 5.837334292238189e-08 min
    # (100) n = 100, avg runtime (5 runs) = 0.034953460693359376 min, stdev runtime (5 runs) = 0.000927458840413892 min
    # (1K) n = 1000, avg runtime (5 runs) = 0.030045989354451498 min, stdev runtime (5 runs) = 0.0006193290956932797 min
    # (10K) n = 10000, avg runtime (5 runs) = 0.029808342456817627 min, stdev runtime (5 runs) = 0.0006168006254449203 min
    # (100K) n = 100000, avg runtime (5 runs) = 0.033507176240285236 min, stdev runtime (5 runs) = 0.0009239435202416883 min
    # (1M) n = 1000000, avg runtime (5 runs) = 0.06002254327138265 min, stdev runtime (5 runs) = 0.00237064446099151 min
    # (10M) n = 10000000, avg runtime (5 runs) = 0.31629184166590374 min, stdev runtime (5 runs) = 0.029654960153288768 min
    # (100M) n = 100000000, avg runtime (5 runs) = 3.3262389659881593 min, stdev runtime (5 runs) = 0.5728247384536419 min
