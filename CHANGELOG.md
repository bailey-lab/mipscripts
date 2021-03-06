# Changelog

## mipscripts 0.3.4

- Improve searching for files. Previously, users would need to provide the
  complete path to input files or file searching would fail. This is no longer
  the case and relative paths may now be input.
- Add short options to the CLI.
- Warn user when converting to snake case.
- Add a function to create a directory.
- Add tests for the subcommands.
- Improve documentation and docstrings.
- Use `argparse` to handle subcommands.

## mipscripts 0.3.3

- Remove unused packages and modules.
- Print out per sample set summary information.
- Throw error if sample sheet column names are duplicated
  ([#4](https://github.com/bailey-lab/mipscripts/issues/4)).
- Fix check to see if FASTQ names are the duplicated. Previously, similar names
  would occasionally trigger an error. Only identical names should trigger an
  error ([#2](https://github.com/bailey-lab/mipscripts/issues/2)).
- Ensure standard column names by converting names to snake case
  ([#3](https://github.com/bailey-lab/mipscripts/issues/3)).

## mipscripts 0.3.2

- Install dependencies when installing package.
- Add a progress bar when counting the number of FASTQ reads.
- Improve error messages.

## mipscripts 0.3.0

- Initial public release.