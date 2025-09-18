## Recreate the R environment that analyses ran in.
## We use the renv package, which tracks the exact packages, package versions, and R version we used.
## Using this will not change the environment of any of your other projects.

install.packages("renv")
renv::restore()

# It is a good idea to restart your R session at this point.
