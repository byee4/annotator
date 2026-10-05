# Small Ensembl GTF example

This seven-feature GTF uses Ensembl's `gene_biotype` and
`transcript_biotype` attributes. Neither transcript has a separate `gene`
feature, which exercises [issue #3](https://github.com/byee4/annotator/issues/3).
The coding transcript is illustrative. The second transcript uses the ID,
coordinates, and exon structure returned by the
[Ensembl GRCh37 lookup for ENST00000606034](https://grch37.rest.ensembl.org/lookup/id/ENST00000606034?expand=1),
the transcript named in the issue; this is a curated subset, not the complete
release 87 GTF.
The committed `mini_ensembl.gtf.db` is a gffutils SQLite database built from
`mini_ensembl.gtf`.

From the repository root, rebuild the database with:

```sh
build_gffutils_db \
  --annotation_file examples/ensembl_gtf/mini_ensembl.gtf \
  --db_file examples/ensembl_gtf/mini_ensembl.gtf.db \
  --force --disable_infer_genes --disable_infer_transcripts
```

Disabling inference keeps the database free of a synthetic `gene` feature.
Run the example with:

```sh
annotator \
  --input examples/ensembl_gtf/queries.bed \
  --output /tmp/mini_ensembl.annotated.bed \
  --gtfdb examples/ensembl_gtf/mini_ensembl.gtf.db \
  --format ensembl
```

The `Genic Region` column (column 9) should read `CDS`, `five_prime_utr`,
`intergenic`, and `noncoding_exon` for the four queries. The BED file uses
zero-based starts; the GTF uses one-based starts. The source GTF and database
are test fixtures, not a biological reference.

To verify the committed database and annotations:

```sh
python -m pytest -q unit_tests/test_ensembl_example.py
```
