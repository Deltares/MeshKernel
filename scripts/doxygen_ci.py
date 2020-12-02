import subprocess
import sys
from pathlib import Path
import argparse


def configure_doxyfile(
    doxyfile_in_path: Path,
    doxyfile_path: Path,
    doxygen_log_path: Path,
    input_dir: Path,
    output_dir: Path,
    docs_dir: Path,
):
    with open(doxyfile_in_path, "r") as file:
        doxyfile_data = file.read()

    doxyfile_data = doxyfile_data.replace("@DOXYGEN_INPUT_DIR@", str(input_dir))
    doxyfile_data = doxyfile_data.replace("@DOXYGEN_OUTPUT_DIR@", str(output_dir))
    doxyfile_data = doxyfile_data.replace("@CMAKE_CURRENT_SOURCE_DIR@", str(docs_dir))

    doxyfile_data = doxyfile_data.replace(
        "@DOXYGEN_WARN_LOG_FILE@", str(doxygen_log_path)
    )

    with open(doxyfile_path, "w") as file:
        file.write(doxyfile_data)


def file_empty(file: Path):
    return file.stat().st_size == 0


def print_file(file: Path):
    with open(file, "r") as f:
        print(f.read())


# Set paths
root_dir = Path(__file__).parent.parent
input_dir = root_dir / "include" / "MeshKernel"
output_dir = root_dir / "build" / "docs"
docs_dir = root_dir / "docs"
doxyfile_in_path = docs_dir / "Doxyfile.in"
doxyfile_path = output_dir / "Doxyfile"
doxygen_log_path = output_dir / "Doxygen_log.txt"

# The dir tree is not created automatically
output_dir.parent.mkdir(exist_ok=True)
output_dir.mkdir(exist_ok=True)


configure_doxyfile(
    doxyfile_in_path, doxyfile_path, doxygen_log_path, input_dir, output_dir, docs_dir
)

# Call doxygen
subprocess.call(f"doxygen {doxyfile_path}", shell=True)

parser = argparse.ArgumentParser()
parser.add_argument(
    "--fail-with-warnings",
    action="store_true",
    help="let the script exit with code 1 with doxygen warnings",
)
args = parser.parse_args()

if args.fail_with_warnings:
    # Check if file is empty
    # If it is not, there were warnings and we exit with 1
    if not file_empty(doxygen_log_path):
        print("There were warnings")
        print_file(doxygen_log_path)
        sys.exit(1)
