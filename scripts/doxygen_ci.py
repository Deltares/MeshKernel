import subprocess
import sys
from pathlib import Path
import argparse


def configure_doxyfile(
    doxyfile_in: Path,
    doxyfile_path: Path,
    doxygen_log_path: Path,
    doxygen_input_list: str,
    output_dir: Path,
    docs_dir: Path,
    fail_with_warnings: bool,
):
    with open(doxyfile_in, "r") as file:
        doxyfile_data = file.read()

    doxyfile_data = doxyfile_data.replace("@DOXYGEN_INPUT_LIST@", doxygen_input_list)
    doxyfile_data = doxyfile_data.replace("@DOXYGEN_OUTPUT_DIR@", str(output_dir))
    doxyfile_data = doxyfile_data.replace("@CMAKE_CURRENT_SOURCE_DIR@", str(docs_dir))

    doxyfile_data = doxyfile_data.replace(
        "@DOXYGEN_WARN_LOG_FILE@", str(doxygen_log_path)
    )

    if fail_with_warnings:
        doxyfile_data = doxyfile_data.replace("@DOXYGEN_EXTRACT_PRIVATE@", "TRUE")
    else:
        doxyfile_data = doxyfile_data.replace("@DOXYGEN_EXTRACT_PRIVATE@", "FALSE")

    with open(doxyfile_path, "w") as file:
        file.write(doxyfile_data)


def file_empty(file: Path):
    return file.stat().st_size == 0


def print_file(file: Path):
    with open(file, "r") as f:
        print(f.read())


# Set dirs
root_dir = Path(__file__).parent.parent
meshkernel_include_dir = root_dir / "include" / "MeshKernel"
meshkernelapi_include_dir = root_dir / "include" / "MeshKernelApi"
output_dir = root_dir / "build" / "docs"
docs_dir = root_dir / "docs"

# Set paths
main_page_path = docs_dir / "main_page.md"
doxyfile_in_path = docs_dir / "Doxyfile.in"
doxyfile_path = output_dir / "Doxyfile"
doxygen_log_path = output_dir / "Doxygen_log.txt"
doxygen_input_list = (
    f"{meshkernel_include_dir} {meshkernelapi_include_dir} {main_page_path}"
)

# The dir tree is not created automatically
output_dir.parent.mkdir(exist_ok=True)
output_dir.mkdir(exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument(
    "--fail-with-warnings",
    action="store_true",
    help="let the script exit with code 1 with doxygen warnings",
)
args = parser.parse_args()

configure_doxyfile(
    doxyfile_in_path,
    doxyfile_path,
    doxygen_log_path,
    doxygen_input_list,
    output_dir,
    docs_dir,
    args.fail_with_warnings,
)

# Call doxygen
subprocess.call(f"doxygen {doxyfile_path}", shell=True)


if args.fail_with_warnings:
    # Check if file is empty
    # If it is not, there were warnings and we exit with 1
    if not file_empty(doxygen_log_path):
        print("There were warnings")
        print_file(doxygen_log_path)
        sys.exit(1)
