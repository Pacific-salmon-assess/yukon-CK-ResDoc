# For Brendan to run once, then delete this code until line 9
# This script will replace "load.R"
# renv will create an renv folder, which will be tracked. Many files within will be untracked (gitignored). Leave defaults for now
install.packages("renv")
library(renv)
init()


## Recreate the R environment that analyses ran in.
## We use the renv package, which tracks the exact packages, package versions, and R version we used.
## Using this will not change the environment of any of your other projects.

install.packages("renv")
library(renv)
renv::restore()
