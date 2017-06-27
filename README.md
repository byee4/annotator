# annotator

example usage: 

```
python annotate.py \
--input BEDFILE \
--output ANNOTATEDFILE \
--gtfdb gencode.v19.annotation.gtf.db \
--limit-chroms-to chr7
```

(```--limit-chroms-to``` will limit the dictionary build to only include these chromosomes for faster processing and less memory footprint. Leave blank to hash all chromosomes in the db file)

outputs bedfile + first priority annotation + all overlapping annotations

Let me know if you have issues/questions: bay001@ucsd.edu
