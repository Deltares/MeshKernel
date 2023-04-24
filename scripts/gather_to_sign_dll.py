#!/usr/bin/python3
"""
Script to harvest and restore dlls for signing purposes.
"""

__author__ = "Maarten Tegelaers"
__copyright__ = "Copyright (C) Stichting Deltares, 2020"
__version__ = "0.1.1"
__maintainer__ = "Maarten Tegelaers"
__email__ = "Maarten.Tegelaers@deltares.nl"
__status__ = "Development"


from pathlib import Path
from typing import Generator
import shutil
import argparse
import json


def determine_to_sign_dll_names(build_dir: Path) -> Generator[str, None, None]:
    """
    Gather all the dll names which need to be signed.

    We assume that the only dll's produced by the D-HYDRO solution are the
    ones created by their corresponding .vcxproj. Any other dll is the result of
    a NuGet package, which should already be signed to begin with.

    Args:
        build_dir (Path): Path to the build dir

    Returns:
        Generator[str, None, None]: Generator containing the dll names to sign.
    """
    libs_vcxprojs = (build_dir / Path("libs")).glob("**/*.vcxproj")
    return (x.with_suffix(".dll").name for x in libs_vcxprojs)


def get_corresponding_path(build_dir: Path, dll_name: str) -> Path:
    """
    Get the dll path corresponding with the dll_name from the provided build_dir.

    Args:
        build_dir (Path): Path to the build dir
        dll_name (str): Name of the dll
    """
    # this hack is required to keep the capitalisation
    return next(
        (build_dir / Path("libs")).glob("**/Release/{}".format(dll_name))
    ).parent / Path(dll_name)


def has_corresponding_path(build_dir: Path, dll_name: str) -> bool:
    """
    Verify whether the provided dll has a corresponding path.
    Args:
        build_dir (Path): Path to the build dir
        dll_name (Generator[str, None, None]): Generator containing the dll names to sign.

    Returns:
        bool: True if a corresponding path can be found; False otherwise.
    """
    return any((build_dir / Path("libs")).glob(f"**/Release/{dll_name}"))


def find_all_dll_paths(
    build_dir: Path, dll_names: Generator[str, None, None]
) -> Generator[Path, None, None]:
    """
    Find the dll paths of the dll names provided in dll_names.

    Args:
        build_dir (Path): Path to the build dir
        dll_names (Generator[str, None, None]): Generator containing the dll names to sign.

    Returns:
        Generator[Path, None, None]: Generator containing the dll paths to sign.
    """
    return (
        get_corresponding_path(build_dir, x)
        for x in dll_names
        if has_corresponding_path(build_dir, x)
    )


def generate_json_mapping_content(
    build_dir: Path, dll_paths: Generator[Path, None, None]
) -> dict:
    """
    Generate the json mapping file, with which the signed dlls can be restored to their original paths.

    Args:
        build_dir (Path): Path to the build dir
        dll_paths (Generator[Path, None, None]): Generator containing the dll names to sign.

    Returns:
        dict: Dict describing the mapping
    """
    return {
        "mapping": list(
            ({"dll": p.name, "path": str(p.relative_to(build_dir))} for p in dll_paths)
        )
    }


DLL_MAPPING_FILE_NAME = "dll_mapping.json"


def write_json_mapping_file(sign_dir: Path, content: dict) -> None:
    """
    Write the specified content to the mapping file within sign_dir.

    Args:
        sign_dir (Path): Path to the directory containing the dlls to sign.
        content (dict): dictionary describing the dll mapping.
    """
    json_mapping_file = sign_dir / Path(DLL_MAPPING_FILE_NAME)
    with json_mapping_file.open("w") as f:
        json.dump(content, f, indent=4)


def read_json_mapping_file(sign_dir: Path) -> dict:
    """
    Load the dll mapping file from the sign_dir.
    Args:
        sign_dir (Path): Path to the directory containing the signed dlls.
    Returns:
        dict: A dictionary containing the mapping of the dlls.
    """
    json_mapping_file = sign_dir / Path(DLL_MAPPING_FILE_NAME)
    with json_mapping_file.open("r") as f:
        return json.load(f)


def harvest_to_sign_dlls(
    dll_paths: Generator[Path, None, None], sign_dir: Path
) -> None:
    """
    Move the dlls specified by dll_paths to the sign_dir.

    Args:
        dll_paths: The dlls to copy
        sign_dir (Path): Path to the directory containing the dlls to sign.
    """
    for p in dll_paths:
        shutil.move(str(p), str(sign_dir))


def restore_signed_dlls(build_dir: Path, dll_mapping: dict, sign_dir) -> None:
    """
    Restore the signed dlls to their original locations.
    Args:
        build_dir (Path): The path to the build dir.
        dll_mapping (dict): Dictionary describing the mapping of dlls
        sign_dir (Path): Path to the directory containing the signed dlls.
    """
    for dll_info in dll_mapping["mapping"]:
        shutil.move(
            str(sign_dir / dll_info["dll"]), str(build_dir / Path(dll_info["path"]))
        )


TO_SIGN_DIR_NAME = "to_sign"


def harvest(build_dir: Path) -> None:
    """
    Create a to_sign folder containing all the dlls to sign, and their corresponding mapping

    Args:
        build_dir: Path to the build dir.
    """
    dll_names = determine_to_sign_dll_names(build_dir)
    dll_paths = list(find_all_dll_paths(build_dir, dll_names))
    for x in dll_paths:
        print("path: ", x)

    sign_dir = build_dir / Path(TO_SIGN_DIR_NAME)
    if not (sign_dir.exists() and sign_dir.is_dir()):
        sign_dir.mkdir(parents=True)

    dll_mapping = generate_json_mapping_content(build_dir, dll_paths)
    write_json_mapping_file(sign_dir, dll_mapping)

    harvest_to_sign_dlls(dll_paths, sign_dir)


def restore(build_dir: Path) -> None:
    """
    Restore the signed dlls to their original locations.
    Args:
        build_dir: Path to the build dir.
    """
    sign_dir = build_dir / Path(TO_SIGN_DIR_NAME)
    dll_mapping = read_json_mapping_file(sign_dir)

    restore_signed_dlls(build_dir, dll_mapping, sign_dir)


def parse_arguments():
    """
    Parse the arguments with which this script was called through
    argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("build_dir_path", help="Path to the root of the repository.")
    action_group = parser.add_mutually_exclusive_group(required=True)
    action_group.add_argument("--harvest", action="store_true")
    action_group.add_argument("--restore", action="store_false")

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()
    build_dir = Path(args.build_dir_path)

    if args.harvest:
        harvest(build_dir)
    else:
        restore(build_dir)
