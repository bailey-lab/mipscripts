import csv
from mipscripts import utils
import pytest

def test_to_snake_case():
    assert utils.to_snake_case("Sample.set") == "sample_set"
    assert utils.to_snake_case("Sample_set") == "sample_set"
    assert utils.to_snake_case("SampleSet") == "sample_set"
    assert utils.to_snake_case("sampleSet") == "sample_set"
    assert utils.to_snake_case("sample Set") == "sample_set"
    assert utils.to_snake_case("sample\nSet") == "sample_set"

def test_header_to_snake_case(tmp_path):
    tmp_dir = tmp_path / "header"
    tmp_dir.mkdir()
    tmp_file = tmp_dir / "sample_list.tsv"

    # Create initial header
    init_header = [
        "sample_name",
        "sampleSet",
        "probe.set",
        "replicate",
        "fw",
        "rev",
        "owner",
        "Capture plate-name",
        "Capture.Plate Location",
        "library_prep",
    ]

    # Write header to file
    with open(tmp_file, mode="w") as file:
        file_content = csv.writer(file, delimiter="\t")
        file_content.writerow(init_header)

    # Convert header to snake case
    utils.header_to_snake_case(tmp_file, overwrite=True)

    # Read result
    with open(tmp_file) as file:
        file_content = csv.reader(file, delimiter="\t")
        res_header = next(file_content)

    # Define expected header
    expect_header = [
        "sample_name",
        "sample_set",
        "probe_set",
        "replicate",
        "fw",
        "rev",
        "owner",
        "capture_plate_name",
        "capture_plate_location",
        "library_prep",
    ]

    # Check that headers match
    assert res_header == expect_header
