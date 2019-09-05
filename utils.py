import os
import pathlib
import requests


def download_file(
    url: str,
    filename: str,
    base: str = ".",
    dir: str = "data",
    overwrite: bool = False,
):
    """Method for downloading data

    Arguments:
        filename {str} -- File access of the ENCODE data file

    Keyword Arguments:
        base {str} -- Base directory (default: {"."})
        dir {str} -- Download directory (default: {"data"})
        overwrite {bool} -- If {True} existing files with be overwritten (default: {False})

    Returns:
        {str} -- Returns a pointer to `filename`.
    """
    filepath = os.path.join(base, dir, filename)

    if pathlib.Path(filepath).is_file() and not overwrite:
        print("File already exist. To overwrite pass `overwrite=True`")
        return filepath

    chunkSize = 1024
    name, _ = os.path.splitext(filename)
    r = requests.get(url, stream=True)

    print("Download {}...".format(filename), end='')

    with open(filepath, "wb") as f:
        for chunk in r.iter_content(chunk_size=chunkSize):
            if chunk:  # filter out keep-alive new chunks
                f.write(chunk)

    print(" done!")

    return filepath
