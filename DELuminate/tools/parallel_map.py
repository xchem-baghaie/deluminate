import os
from typing import Optional, Callable, List, Union

import numpy as np
from metaflow import parallel_map

BATCH_SIZE = os.cpu_count() * 1


def select_adaptive_batch_size(total_len: int, batch_size: Optional[int] = None):
    if batch_size is None:
        batch_size = total_len // BATCH_SIZE
    else:
        batch_size = min(total_len, batch_size)
    return batch_size


def parallel_map_generic(
    array: Union[List, np.ndarray],
    call_fn: Callable,
    batch_size: Optional[int] = None,
    **kwargs
) -> Union[List, np.ndarray]:
    """
    Args:
        -array (Union[List, np.ndarray]): Numpy array or list
        -call_fn (Callable): Function being used
        -batch_size (Optional[int]): Batch size used for parallel map. If not based, defaults to CPU count
        -kwargs (Any): Used in function call

    Returns:
        List: Combined list of function called
    """
    total_len = len(array)
    batch_size = select_adaptive_batch_size(total_len, batch_size=batch_size)
    if batch_size > 0:
        return parallel_map(
            lambda roll_idx: [
                call_fn(x, **kwargs)
                for x in array[roll_idx : min(roll_idx + batch_size, total_len)]
            ],
            range(0, total_len, batch_size),
        )
    else:
        return [[call_fn(x, **kwargs) for x in array]]

def parallel_map_numpy_array(
    array: Union[List, np.ndarray],
    call_fn: Callable,
    batch_size: Optional[int] = None,
    **kwargs
) -> Union[List, np.ndarray]:
    """
    Parallel map used when we want a numpy array input, i.e. 2D array, and not a list input parallel_map_generic

    Args:
        array (Union[List, np.ndarray]): Numpy array or list
        call_fn (Callable): Function being used
        batch_size (Optional[int]): Batch size used for parallel map. If not based, defaults to CPU count
        kwargs (Any): Used in function call

    Returns:
        List: Combined list of function called
    """
    from metaflow import parallel_map
    import numpy as np

    total_len = len(array)
    batch_size = select_adaptive_batch_size(total_len, batch_size=batch_size)
    return parallel_map(
        lambda roll_idx: [
            call_fn(np.reshape(x, (1, -1)), **kwargs)
            for x in array[roll_idx : min(roll_idx + batch_size, total_len)]
        ],
        range(0, total_len, batch_size),
    )
