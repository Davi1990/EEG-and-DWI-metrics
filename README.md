# EEG-and-DWI-metrics

This repository includes the code required to reproduce the results in: "Network-level Macroscale Structural Connectivity Predicts Propagation of Transcranial Magnetic Stimulation" Momi D., Ozdemir R., Tadayon E., Boucher P., Shafi M., Pascual-Leone A., Santarnecchi E.. NeuroImage (2020).

<p align="center">
    <img src="https://github.com/Davi1990/EEG-and-DWI-metrics/blob/main/Figure_1.png" width="1000"/>
</p>


# TMS-EEG 
All EEG data pre-processing was performed offline using EEGLAB 14.1 (https://github.com/sccn/eeglab) (Delorme and Makeig, 2004) and customized script running in Matlab R2017b (Math-Works Inc., USA). All TMS-evoked EEG source reconstruction and analysis was performed using Brainstorm (https://github.com/brainstorm-tools/brainstorm3) (Tadel et al., 2019) and customized script running in Matlab R2017b (Math-Works Inc., USA). 

For EEG data preprocessing please use:
```matlab
    'TMS_EEG_preprocessing.m' 
```

For EEG data analysis please use 
```matlab
    'EEG_metric_extraction.m ' 
```

# DWI
All DWI metrics were extracted using a customized script running in Matlab R2017b (Math-Works Inc., USA). 
For DWI data analysis please use 
```matlab
    'DWI_metrics_extraction.m' 
```

There is an updated version of the code at the following link: https://github.com/Davi1990/DissNet

# TMS Targeting 
For TMS targets please use: 
```matlab
    'yeo2subject.sh subj nNet' 
```

