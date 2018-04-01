# Performance notes
At the time, things done to ensure performance include:
* Ensuring that any transformation of variants is only performed once per variant
* Removing the need for lookup dictionaries by establishing the order of comparison within each line of the VCF file right on the header line

## Additional possibilities:
* Refactor `DefaultComparer.update_comparison` so that is does not require a set. This should reduce some overhead and possibly speedup things a bit more
* Use `itertools.islice` to avoid creating a copy of the list of fields when we extract the metadata and data fields for each line, in `Handler.run` (see the line containing the code `data_fields = fields[FIXED_FIELDS_LENGTH:]`)
    * Note that this will require some refactoring regarding the `process_variant` method, which rewrites the value in the list
