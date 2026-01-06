# CGOSS code and related dataset
This all of relevant Codes and used Dataset for CGOSS. 

Dataset 'climate_change_impact_on_agriculture.csv' is a public synthetic dataset was collected from kaggle with URL:https://www.kaggle.com/datasets/waqi786/climate-change-impact-on-agriculture.

'energydata_complete.csv' is a real dataset which collected from UCI public dataset with DOI: https://doi.org/10.24432/C5VC8G.

People can find three R codes do the comparison between CGOSS GOSS OSS GIBOSS and IBOSS and two cpp codes which is the Performance-critical components in the CGOSS code file. The 'simulation.R' is do the comparison in the simulation dataset that mentioned in the paper. 'publicdataset.R' and 'realdataset.R' used same main code apply to two different datasets, 'climate_change_impact_on_agriculture.csv' and'energydata_complete.csv'.

All of those hyper-parameters are not theoreticaly optimized for all cases which means for different problem people should use different configure of hyper-parameters.
