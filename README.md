# Open R Function Library

Hi, this is a set of R function for handling different topics in basic data analysis and bioinformatics.

## Usage

R package [`{modules}`](https://github.com/klmr/modules) is recommended to be used here importing different R files as different modules like Python.

```R
install.packages("modules")
```

For example, use the TCGA functions:

```R
> tcga <- modules::use("https://biosisyphus.github.io/Rlib/tcga.R")
> tcga$filterReplicates(tsb = c("TCGA-55-7913-01B-11D-2237-01", "TCGA-55-7913-01B-11X-2237-01", "TCGA-55-7913-01B-11D-2237-01"))
ooo Filter barcodes successfully!
[1] "TCGA-55-7913-01B-11D-2237-01"
```

You can also clone the repo to your local machine and then load the specific module with `modules::use("<the_path_to_module_file>")`.

## Modules

Modules are saved as single files in the repo. All dependencies are well controlled or would be installed after starting the function
you want to run.

Use the following way to check public functions available in a module:

```R
> install
clone:
function(url, local_path, gitee = FALSE, reset_remote = FALSE, ...)


download:
function(repo, destdir, release = NULL, gitee = FALSE, ...)


install:
function(pkg, gitee = FALSE, ...)


save:
function()
```

You may need to check the source file in the GitHub/Gitee repo to read the detail documentation.

**Please file an issue if you find any bugs.**

Module list:

- `tcga.R` -  contains functions useful for processing TCGA data.
- `install.R` - contains functions useful for installing packages.

## LICENSE

[GPL-3](LICENSE)

Copyright (C) Shixiang Wang <w_shixiang@163.com>
