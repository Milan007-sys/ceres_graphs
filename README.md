# A Python script for visualization of the global radiation at the top of the atmosphere from CERES

[alt The final image](https://github.com/Milan007-sys/ceres_graphs/blob/main/CERES_radiation_fluxes_globe_ENG.png?raw=true)

The final picture shows the 12-month running mean of the CERES radiation components
along with the GISTEMP global tempareture dataset.

The GISTEMP data can be obtained here:
https://data.giss.nasa.gov/gistemp/graphs_v4/graph_data/Monthly_Mean_Global_Surface_Temperature/graph.txt
The last data are up to 2024-07 and are edited to match the dates of the CERES data. 

You have to download the CERES data yourself as they are too large to store at 
Github. This script works with the dataset CERES_EBAF-TOA_Ed4.2_Subset_200003-202407.nc
The source of the CERES data: 
https://ceres.larc.nasa.gov/data/#ebaftoa-level-3

This script is tested on WSL/Ubuntu but should run on any system able to host
Anaconda/Miniconda environments

In Miniconda, you have to install some packages, namely pyplot, Matplotlib, xarray, dask, netCDF4 and bottleneck. 

You can use and improve the script in any way, if you find some errors or shortcoming, let me know.

Enjoy!

Milan Salek 
