from typing import List, Tuple, Union
from itertools import product


def generate_unique_combinations_with_at_least_one_instruction(
    with_instruction_rows: Union[List[List[Tuple]], List[List]],
    no_instruction_rows: Union[List[List[Tuple]], List[List]],
) -> List[Tuple[Tuple]]:
    """
    Args:
    _____
    -with_instruction_rows: Outer list has length equal to the number of cycles
        Inner list has length equal to the number of building blocks with an instruction. If there are no such building blocks, this will be an empty list (so will not contain an inner tuple)
        Inner Tuple has length equal to 3, and is of the form (XBBID: Union[str,int], XSMILES: str, Instruction: str = '')
    -no_instruction_rows: Outer list has length equal to the number of cycles
        Inner list has length equal to the number of building blocks with no instruction. If there are no such building blocks, this will be an empty list (so will not contain an inner tuple)
        Inner Tuple has length equal to 3, and is of the form (XBBID: Union[str,int], XSMILES: str, Instruction: str = '')

    Returns:
    ________
    List of tuple of tuples.
        The outer list has length equal to the number of combinations
        The outer tuple has length equal to the number of cycles. Each element of this tuple is one instance.
        The inner tuple has length equal to 3, and is of the form (XBBID: Union[str,int], XSMILES: str, Instruction: str = '')
    If no combinations are generated, None is returned.

    The function operates under the assumption that no instruction is present in any row of the no_instruction_rows argument, and instructions are present in all rows of the with_instruction_rows argument.

    It creates boolean masks corresponding to the different possible combinations of whether an instance contains a BB with an instruction or not, stored in the "cases" variable.
    Cases will have length = 2^(n_cycles) - 1
    Each case element will have length = n_cycles, and will be a unique combination of True or False, with the all False case excluded.

    For example:
    If n_cycles = 2, cases = [(True, True), (True, False), (False, True)]
    If n_cycles = 3, cases = [(True, True, True), (True, True, False), (True, False, True), (True, False, False), (False, True, True), (False, True, False), (False, False, True)]
    etc.

    For each case, if the value at a position is True, the corresponding list from with_instruction_rows will be used. Otherwise, the corresponding list from no_instruction_rows will be used.

    This ensures that only unique combos are output.

    """
    try:
        assert len(with_instruction_rows) == len(no_instruction_rows)
    except:
        raise Exception(
            f"Mismatched number of cycles: {len(with_instruction_rows)} with instruction, {len(no_instruction_rows)} without instruction."
        )

    n_cycles = len(with_instruction_rows)

    cases = [c for c in list(product([True, False], repeat=n_cycles)) if any(c)]
    assert len(cases) == (2**n_cycles) - 1

    combos = []
    for case in cases:
        list_args = []
        for i in range(n_cycles):
            if case[i]:
                list_args.append(with_instruction_rows[i])
            else:
                list_args.append(no_instruction_rows[i])
        combos.extend(list(product(*list_args)))
    if len(combos) == 0:
        return None
    else:
        return combos


def generate_combinations_with_no_instruction(
    no_instruction_rows: Union[List[List[Tuple]]],
) -> List[Tuple[Tuple]]:
    """
    Args:
    _____
    -no_instruction_rows: Outer list has length equal to the number of cycles
        Inner list has length equal to the number of building blocks with no instruction. If there are no such building blocks, this will be an empty list (so will not contain an inner tuple)
        Inner Tuple has length equal to 3, and is of the form (XBBID: Union[str,int], XSMILES: str, Instruction: str = '')

    Returns:
    ________
    List of tuple of tuples.
        The outer list has length equal to the number of combinations
        The outer tuple has length equal to the number of cycles. Each element of this tuple is one instance.
        The inner tuple has length equal to 3, and is of the form (XBBID: Union[str,int], XSMILES: str, Instruction: str = '')
    If no combinations are generated, None is returned.

    The function operates under the assumption that no instruction is present in any row of the no_instruction_rows argument.
    """
    combos = list(product(*no_instruction_rows))
    if len(combos) == 0:
        return None
    else:
        return combos


if __name__ == "__main__":
    import numpy as np

    # (2,1,3)
    with_instruction_rows_1 = [
        [("BID1", "BXS1", "BI1"), ("BID2", "BXS2", "BI2")],
        [("CID1", "CXS1", "CI1")],
        [("DID1", "DXS1", "DI1"), ("DID2", "DXS2", "DI2"), ("DID3", "DXS3", "DI3")],
    ]
    # (0,1,3)
    with_instruction_rows_2 = [
        [],
        [("CID1", "CXS1", "CI1")],
        [("DID1", "DXS1", "DI1"), ("DID2", "DXS2", "DI2"), ("DID3", "DXS3", "DI3")],
    ]
    # (0,1,0)
    with_instruction_rows_3 = [[], [("CID1", "CXS1", "CI1")], []]
    # (0,0,0)
    with_instruction_rows_4 = [[], [], []]
    # (1,3,2)
    no_instruction_rows_1 = [
        [("BID3", "BXS3", "")],
        [("CID2", "CXS2", ""), ("CID3", "CXS3", ""), ("CID4", "CXS4", "")],
        [("DID4", "DXS4", ""), ("DID5", "DXS5", "")],
    ]
    # (0,3,2)
    no_instruction_rows_2 = [
        [],
        [("CID2", "CXS2", ""), ("CID3", "CXS3", ""), ("CID4", "CXS4", "")],
        [("DID4", "DXS4", ""), ("DID5", "DXS5", "")],
    ]
    # (1,1,2)
    no_instruction_rows_3 = [
        [("BID4", "BXS4", "")],
        [("Plc", "Plc_SMILES", "")],
        [("DID6", "DXS6", ""), ("DID7", "DXS7", "")],
    ]
    # (3,2)
    no_instruction_rows_4 = [
        [("BID5", "BXS5", ""), ("BID6", "BXS6", ""), ("BID7", "BXS7", "")],
        [("CID5", "CXS5", ""), ("CID6", "CXS6", "")],
    ]

    with_instruction_rows = with_instruction_rows_1
    no_instruction_rows = no_instruction_rows_1

    count_no_instruction = [
        len(no_instruction_rows[i]) for i in range(len(with_instruction_rows))
    ]
    combined_list = [
        len(with_instruction_rows[i]) + len(no_instruction_rows[i])
        for i in range(len(with_instruction_rows))
    ]
    expected_total_no_instruction = np.prod(count_no_instruction)
    expected_total = np.prod(combined_list)
    expected_total_with_instruction = expected_total - expected_total_no_instruction
    combos1 = generate_unique_combinations_with_at_least_one_instruction(
        with_instruction_rows, no_instruction_rows
    )
    combos2 = generate_combinations_with_no_instruction(no_instruction_rows)
    combos = combos1 + combos2
    print(
        f"At least one instruction: {len(combos1)}, expecting {expected_total_with_instruction}.\nNo instruction: {len(combos2)}, expecting {expected_total_no_instruction}. \nTotal: {len(combos1) + len(combos2)}, expecting {expected_total}."
    )
