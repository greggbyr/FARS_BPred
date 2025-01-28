# FARS_BPred
Contains code and benchmark tests for Fully Adiabatic, Reversible, and Superscalar Branch Predictor implementations in the SimpleScalar simulator.

To replicate results, run the following command in simulator/ after cloning:

  source <(cat bench_runs/*); source <(cat extract_data/*)

Results can be found in simulator/compiled_results/*.csv
Note: there are some benchmarks that are not yet fully enabled to run in reverse mode. These are either completely missing bpred_* data points or have all rev data points at 0. 
Known good benchmarks are gcc, go, ijpeg, li, and perl.
