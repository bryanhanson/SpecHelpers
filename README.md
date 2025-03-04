
## How to install SpecHelpers

### To install from CRAN:

````r
install.packages("SpecHelpers")
library("SpecHelpers")
````

### To install from Github using R:

````r
install.packages("remotes")
library("remotes")
install_github(repo = "bryanhanson/SpecHelpers@main")
library("SpecHelpers")
````

If you use `@some_other_branch` you can download other branches that might be available.  They may or may not pass CRAN checks and thus may not install automatically using the method above.  Check the NEWS file to see what's up.

Questions?  hanson@depauw.edu
