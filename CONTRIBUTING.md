# Contributing to xFischer

Thank you for helping improve xFischer. Please open an issue before a substantial change so its scientific scope and compatibility can be discussed.

## Development setup

```bash
conda env create -f environment.yml
conda activate xfischer
python -m unittest discover -s tests -v
```

To verify package construction:

```bash
python -m pip install -e ".[dev]"
python -m build
```

## Pull requests

- Keep changes focused and document user-visible behavior in `CHANGELOG.md`.
- Do not change equations, filter semantics, defaults, or statistical definitions without a clear derivation and comparison against an independently calculated reference case.
- Update both the user and scientific documentation when behavior changes.
- Ensure the test suite passes on a clean environment.

By contributing, you agree that your contribution is licensed under the repository's MIT License.

