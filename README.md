# Quantifying the impact of black carbon pollution on snow melt rates using an energy balance model in the Mont Blanc massif

## Project Description

This program estimates snowmelt over the Mont-Blanc region using publicly available meteorological data. It computes the snow surface energy balance to quantify melt rates under realistic atmospheric conditions. By comparing clean and polluted snow albedo, it evaluates the impact of surface impurities like black carbon on melting. Results are visualized in 3D, highlighting areas where pollution accelerates snowmelt and showing spatial melt patterns over time.

### Input files

The main input file used by the program, mont_blanc_weather_2025_day.csv, is too large to be uploaded to GitHub. For this reason, it is provided separately by email in order to allow the code to be executed. Nevertheless, the [data/](data/) directory contains all the scripts that were used to generate this dataset. Each script is prefixed with a numerical index indicating the order in which it must be executed to reproduce the full data processing pipeline.

The first step consists in defining the spatial domain of the study area (centered on Mont Blanc) with the script [1_coordinates_creation.py](data/Preprocessing/1_coordinates_creation.py). It outputs a CSV file containing all coordinate pairs of the grid. This procedure resulted in a total of 2160 grid points, which was considered an appropriate spatial resolution for the objectives of this project.

In the second step, the script [2_meteo_data.py](data/Preprocessing/2_meteo_data.py) is used to retrieve meteorological data for each grid point. This script relies on the Open-Meteo Historical Forecast API to download the atmospheric variables required to compute the snow energy balance. The study period was chosen as the first week of May 2025, a time when snowmelt typically occurs in alpine environments. The resulting file, mont_blanc_weather_2025.csv, contains 211,680 rows and includes the following variables: air temperature at 2 meters, relative humidity, mean sea-level pressure, surface pressure, wind speed at 10 meters, incoming shortwave radiation, rainfall, snowfall, snow depth, and a binary indicator specifying whether the time step corresponds to daytime or nighttime.

Finally, the script [3_Filter_table.m](data/Preprocessing/3_Filter_table.m) is used to post-process the meteorological dataset. This step removes all nighttime data by keeping only rows where the is_day variable equals 1. This simplification is justified because several terms of the snow energy balance depend directly on solar radiation, and snowmelt during the night can be considered negligible due to low temperatures. In addition, the relative humidity values are converted from percentages to fractions by dividing by 100, ensuring consistency with the physical formulas used later in the model.

After completing these steps, the final file mont_blanc_weather_2025_day.csv is obtained. This dataset contains all the necessary input variables, properly filtered and formatted, to compute the snow energy balance and analyze the impact of black carbon pollution on snowmelt rates over the Mont-Blanc area.


### Output files

The program generates a file named results_melt.csv, which contains several columns with calculated values from the energy balance model applied over the Mont-Blanc area (this file is also provided separately by email, as it is too large to be uploaded to GitHub). 

This data will be used to create two 3D mountain plots, one showing pure snow melt and the other showing polluted snow melt. These plots are interactive, allowing users to slide through time, with each step corresponding to an hour of the first week of May 2025 (when there is sunlight). Below the 3D plots, a time series graph will be displayed showing the melt rates (polluted vs clean over the course of the week). All these visualizations are saved in the [Plot_result.fig](results/Plot_result.fig) file. 

Additionally, the program generates a [melt_volume_results.tex](results/melt_volume_results.tex) file, which summarizes the total snowmelt volume for both clean and polluted snow and highlights the extra melt caused by pollution. The document also compares melt volumes, quantifies the additional melt in Olympic-size swimming pools, and reports the overall percentage increase under polluted conditions, illustrating the impact of black carbon on snowmelt rates.


### Report

Templates for .pdf, .tex and .odt formats are provided in [docs/](docs/). The final report should be placed in "docs/" as "report.pdf".

## Running the program

### Dependencies

Data preprocessing and acquisition are performed using Python (version ≥ 3.8). The Python scripts require the following external packages: requests (for accessing the Open-Meteo historical weather API), csv (standard Python library, no installation required), time and os (standard Python libraries). These packages can be installed using pip if not already available in the environment.

The numerical core of the project is written in C and compiled using GCC (GNU Compiler Collection, version ≥ 9.0). The C code relies only on standard libraries (stdio.h, stdlib.h, math.h, string.h) and therefore does not require any external scientific libraries. Compilation is performed using the -lm flag to link the standard math library.

Data filtering and visualization are performed using MATLAB (version ≥ R2021a). MATLAB is required for: filtering and preprocessing large CSV datasets, running the 3D visualization and time-dependent melt plots and generating interactive figures with sliders. No additional MATLAB toolboxes beyond the base MATLAB installation are required. However, a recent MATLAB version is recommended to ensure full compatibility with the datetime, uicontrol, and 3D plotting functionalities used in the scripts.


### Build

The snow energy balance model is implemented in C [snow_energy_balance.c](src/snow_energy_balance.c). To compile it, use the GNU Compiler Collection (gcc) from the root of the project (`Project-CMT`):

```bash
gcc src/snow_energy_balance.c -o bin/snow_energy_balance.out -lm
```

### Execute

The complete workflow, including compilation of the C program, computation of snow melt, and plotting of results, is automated in the MATLAB script [plot_results.m](src/plot_results.m).

To run it, open MATLAB and press **Run**. All output files (CSV, LaTeX, figures) will be generated in the `results/` folder.

Alternatively, from a terminal in the project root, you can execute the workflow in batch mode with:

```bash
matlab -batch "src/plot_results"
```
## Contributors

Adrien Mathieu & Nelson Almeida Faxiolo

## Acknowledgments

### Data sources

The meteorological input data used in this project are obtained from the [Open-Meteo Historical Forecast API](https://open-meteo.com/en/docs/historical-weather-api), an open-access weather data service providing free and publicly available meteorological datasets. Open-Meteo aggregates numerical weather prediction (NWP) outputs from major forecasting models and reanalysis products, and exposes them through a standardized REST API. In this project, the historical forecast endpoint is used to retrieve hourly meteorological variables for a user-defined spatial grid and time period.

### Code

The core scientific methodology, physical assumptions, and overall program architecture were designed and implemented by the authors of this project. All numerical models, parameter choices, and scientific interpretations are the result of the authors’ own work.

Some portions of the code and documentation benefited from assistance provided by a large language model (LLM), specifically OpenAI’s ChatGPT. All LLM-generated suggestions were reviewed, adapted, and validated by the authors, and responsibility for the correctness of the final implementation remains entirely with them.

No external proprietary codebases were directly copied into the project. 

