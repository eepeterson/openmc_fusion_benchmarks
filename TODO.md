# To do
### Manage CAD geometry
- Better (more resilient) way to name the downloaded file?
- Add the possibility to chose the path to download it in
- Implement the `rtt_to_h5m` workflow
- implement choice of rtt workflow or directly h5m in Benchmark class
- `is_file_in_folder()` on h5m/rtt to avoid download everytime the file
- use the above function in `get_model()` of the `Benchmark` class if files are not in `cwd`, download them (in `cwd`)

### Benchmark modules
- Figure out how to import mesh for the DAGMC geometry (download h5m and json for tallies? download rtt and go `rtt_to_h5m` for tallies? Or download h5m and save tallies directly in the python api model?)
- implement UM tallies choice
- check if wwinp file is there, if not, run analog
- print the python api model input file in cwd (?)

## Tests
- Add tests to everything

## Notebooks
- Show how the `Benchmark` class and factory pattern works (the different ways to call a benchmark), a bit of focus on run_option
- Show how the model objects in src/ofb/benchmarks/ are structured and can be used
- The irdff function
- Download geometry

## clean postprocessing
- Rework the whole statepoint file thing

## Miscellanea
- The to_hdf method should be revised: one should be able to provide the "when" and "where" and a default hdf_file name should be added, think about the xaxis thing as well (?)
- The to_hdf file should be more resilient to errors
- Move the tally_h5 functiono out of the OpenmcResults class for more versatility (made necessary by the fng_duct heating tally)
- Add a `run_on_hpc` method to the `Benchmark` class that parses a `run_model.sh` file to to run the model as a slurm job


## Questions
- rtt/h5m files should be automatically downloaded during `Benchmark.get_model()` or only when necessary (e.g. `Benchmark.get_model.run()`) or requested (e.g. `Benchmark.download_h5m_file()`)?
- Should we implement the `rtt_to_h5m` workflow? 
- Should we provide both rtt and h5m meshes or just one of the two would be sufficient? Because w/ the `rtt_to_h5m` workflow we could provide just one as long as models and scripts don't get modified
- Split `rtt_to_h5m` to branch out the rtt parser?