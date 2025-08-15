Nanopore_m6A_Finder can be used to identify m6A sites and m6A levels 
in both coding region and poly(A) region of Nanopore DRS reads.

We trained a XGBoost model and a Random Forest model using in vitro
transcribed 5mers with As or with m6As. The trained models were saved
in google drive.

The training data were saved as 81_1_1500_top.sam and 82_1_1500.sam.

Steps to identify m6A sites in coding region need to firstly run "get_pos_CDs.py" 
to extract the coding region from reads using a aligned sam file.

Then, you may use the output of "get_pos_CDs.py" and DRS fast5 files to run 
"Extract_feature_per_5mer.py".

Then you may run "basecall_feature.py" to encode the feature directory which 
were as outputs of "Extract_feature_per_5mer.py".

The encoded features of the 5mers from reads were input of " xgb_model_pred.py"

The output of "xgb_model_pred.py" was the input of "get_chr_pos.py"

After "get_chr_pos.py", we run "get_interval_and_read_count.py" to get the m6A sites of each dataset.

Then we run "calculate_coverage.sh" to calculate the coverage of each potential m6A site

After get the coverage, we run "merge_to_get_coverage.py" to calculate each sites' m6A level, and 
"merge_WT_KO.py" to combine the WT and m6A writer deletion dataset together.

Then we run "bootstrap_new.py" with the existence of "config_bootstrap.py" to get the reliable 
m6A sites whose WT m6A levels significantly higher than the KO.

We screen the output sites with the column significant. We output the reliable m6A sites and levels.

For polyA region identification. We follow the similar process.

We first run "get_pos_polyA.py" with two aligned sam files.

Then we run the sam feature extraction and encoding steps.

Lastly, we run "RF_combine_xgb_predict.py" to assemble the 
trained XGBoost model and random forest model, getting the final prediction results.
