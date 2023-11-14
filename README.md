# Common Tools

## Setting the environment

```
cmsrel CMSSW_13_0_13
cd CMSSW_13_0_13/src
cmsenv
git clone https://github.com/Ma128-bit/xgboost_to_ULs/ .
scram b -j20
```

## Run the script
Example: 
```
python3 plots_for_tau3mu.py --config ../../../xgboost/CMSSW_13_0_13/src/xgboost/results/BDT/20231111-232226/Run2022_config.json --categories "Cat_A Cat_B Cat_C" --index 1
```
```
python3 BDT_optimal_cut.py --config ../../../xgboost/CMSSW_13_0_13/src/xgboost/results/BDT/20231111-232226/Run2022_config.json
```
