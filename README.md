# Testing pseudo bayes grouping algorithm and patterns

Run the code:

       cmsRun L1Trigger/DTTriggerPhase2/test/test_primitivesPhase2Prod_DEBUG_Grouping_cfg.py > dump.txt

Create plots and log:

       cd L1Trigger/DTTriggerPhase2/test/macros/
       python compare_groupings_StdToBayes.py
       python plotGroupings_StdToBayes.py
       python compare_number_of_fits.py
       cd -



