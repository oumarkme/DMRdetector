"""Bundled BED files (data/assets/), read through importlib.resources."""

from importlib.resources import files


def asset(name: str) -> str:
    """Returns the path of a bundled file in data/assets/.

    Args:
        name: File name, for example par_hg38.bed or icr_gene_list.bed.

    Returns:
        The absolute path, which works from a checkout and from an installed
        copy.
    """
    return str(files(__package__) / "assets" / name)
