# Validation notes

This package was cleaned and checked at the package level.

## Checks performed during assembly

- Python helper files were syntax-checked with `python -m py_compile`.
- Shell entry points were syntax-checked with `bash -n`.
- Synthetic helper tests were included under `tests/`.

## Limitation of the assembly environment

The helper scripts depend on `astropy`, and the assembly environment used to prepare this tarball did not have `astropy` installed. Because of that, the synthetic tests were discoverable and executed by `unittest`, but they were skipped rather than run to completion.

After installing `requirements/python.txt`, run:

```bash
python -m unittest discover -s tests -v
```

for the helper tests, and then validate the full pipeline against SAS on the real observation tree.
