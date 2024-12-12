# MT3
## Exploring Rainfall Propagation Off Sumatran Coast
This project aims to analyze the propagation of rainfall off the coast of Sumatra by mapping relevant variables such as precipitation, winds, temperature, and cloud cover, with the goal of better understanding and perhaps better predicting precipitation and its movements in an area that faces particularly intense rainfall.  
This project makes use of a 10-year ERA5 Reanalysis dataset, spanning 2010-2019, which is accessible via the EMCWF, and a precipitation dataset derived from that dataset via code by Jingyi Hu with the Pennsylvania State University.
### Introductory Files
- 1hr_pcp.py  
This file plots one hour of precipitation data.
- 10m_wind.py  
This file plots one hour of 10 meter winds.
- pcpwind_anysrc_anyhr_cbar.py  
This file plots one hour of precipitation and one hour of 10 meter winds together, with a unique colorbar. It works with multiple sources of precipitation (total, MCS, deep convection, and congestus).

### Files that Cover Entire Dataset
- 10yravgpcp.py  
This file averages the precipitation.
- 10yrUVavg.py  
This file averages the U and V winds.
- plot_10yrpcp.py  
This file plots the averaged precipitation.
- plot_10yrwind.py  
This file plots the averaged winds. 
- ccavg.py  
This file averages the cloud cover.
- plot_CCavg.py  
This file averages the cloud cover.
- tempavg.py  
This file averages the temperature.
- plot_Tempavg.py  
This file plots the averaged temperature.
- MAM_plot.py  
This file plots exclusively the averaged precipitation for March, April, and May of the 10 years.
- DJF_plot.py
This file plots exclusively the averaged precipitation for December, January, and February of the 10 years.
