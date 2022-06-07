from distutils import dir_util
from mipscripts.seqrun_stats import seqrun_stats
from snapshottest.file import FileSnapshot
import pytest
import os


def test_incorrect_args(snapshot, tmpdir):
    # Copy test data into a temporary directory and find sample list path
    dir_util.copy_tree("tests/test_data", str(tmpdir))
    path = str(tmpdir.join("sample_list.tsv"))

    # Error if groups are wrong
    with pytest.raises(SystemExit):
        seqrun_stats([path], "incorrect", None)
    with pytest.raises(SystemExit):
        seqrun_stats([path], "sample_set", "incorrect")
    with pytest.raises(SystemExit):
        seqrun_stats([path], "incorrect", "incorrect")


def test_seqrun_stats(snapshot, tmpdir):
    # Copy test data into a temporary directory and find sample list path
    dir_util.copy_tree("tests/test_data", str(tmpdir))
    path = str(tmpdir.join("sample_list.tsv"))

    # Run seqrun_stats
    seqrun_stats([path], "sample_set", None)

    # Check snapshot
    snapshot.assert_match(FileSnapshot("_readcnt".join(os.path.splitext(path))))
