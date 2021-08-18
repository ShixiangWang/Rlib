# Open R Function Library

Hi, this is a set of R function for handling different topics in basic data analysis and bioinformatics.

## Usage

R package [`{modules}`](https://github.com/klmr/modules) is recommended to be used here importing different R files as different modules like Python.

```R
install.packages("modules")
```

For example, use the TCGA functions:

```R
> tcga <- modules::use("tcga.R")
> tcga$filterReplicates(tsb = c("TCGA-55-7913-01B-11D-2237-01", "TCGA-55-7913-01B-11X-2237-01", "TCGA-55-7913-01B-11D-2237-01"))
ooo Filter barcodes successfully!
[1] "TCGA-55-7913-01B-11D-2237-01"
```

## LICENSE

[GPL-3](LICENSE)

Copyright (C) Shixiang Wang <w_shixiang@163.com>
