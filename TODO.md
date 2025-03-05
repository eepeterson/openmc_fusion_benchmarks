# to do
### download geometry
- better (more resilient) way to name the downloaded file?
- add the possibility to chose the path to download it in
- implement the rtt-to-h5m workflow (`rtt_to_pydagmc` and `pydagmc_to_h5m` scripts) 


- the to_hdf method should be revised: one should be able to provide the "when" and "where" and a default hdf_file name should be added, think about the xaxis thing as well (?)
- the to_hdf file should be more resilient to errors
- move the tally_h5 functiono out of the OpenmcResults class for more versatility (made necessary by the fng_duct heating tally)

## tests
- add tests to everything

## notebooks
- show how the Benchmark class and factory pattern works (the different ways to call a benchmark), a bit of focus on run_option
- Show how the model objects in src/ofb/benchmarks/ are structured and can be used
- the irdff function
- download geometry

## clean postprocessing
- rework the whole statepoint file thing
