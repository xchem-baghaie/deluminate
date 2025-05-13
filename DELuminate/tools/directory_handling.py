import os


def create_path(directory: str) -> str:
    if directory[-1] not in ["/", "\\"]:
        directory = directory + "/"
    if not os.path.exists(directory):
        os.makedirs(directory)
    return directory
