These files provide synthetic test inputs for the exported OptiDesign functions.

Files for prep_famoptg()
- prep_famoptg_input.csv
- prep_famoptg_id_map.csv
- prep_famoptg_GRM.csv
- prep_famoptg_A.csv
- prep_famoptg_K.csv

Files for alpha_rc_stream()
- alpha_input.csv
- alpha_id_map.csv
- alpha_GRM.csv
- alpha_A.csv
- alpha_K.csv

Matrix files use the first column as row names (LineID).
In R, read them with for example:

mat <- as.matrix(read.csv(system.file("extdata", "prep_famoptg_GRM.csv", package = "OptiDesign"), row.names = 1, check.names = FALSE))
