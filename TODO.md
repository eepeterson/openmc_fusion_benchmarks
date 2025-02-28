# to do
- the to_hdf method should be revised: one should be able to provide the "when" and "where" and a default hdf_file name should be added, think about the xaxis thing as well (?)
- the to_hdf file should be more resilient to errors
- move the tally_h5 functiono out of the OpenmcResults class for more versatility (made necessary by the fng_duct heating tally)

## tests

## notebooks
- show how the Benchmark class and factory pattern works (the different ways to call a benchmark), a bit of focus on run_option
- Show how the model objects in src/ofb/benchmarks/ are structured and can be used
- the irdff function

## clean postprocessing
- being able to extract tally from statepoint files
- pick tally from database
- mean, std_dev, rel_std_dev, ratio_of_means, over_lethargy, heating
