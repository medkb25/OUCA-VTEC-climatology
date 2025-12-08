# OUCA VTEC Climatology – Processing and Plotting Code

This repository contains the code used in:

> Kaab, M., et al. (2025). Quiet-time GNSS-based VTEC climatology over Oukaimeden (Morocco), 2015–2025 (QSL–GIM20) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.17144341

It provides the processing pipeline and plotting scripts for the quiet-time GNSS-based VTEC climatology over Oukaimeden Observatory (OUCA), Morocco, for the period October 2015 to September 2025.

## Contents

Main notebooks:

- "00_ZENODO_Generate_Master_files.ipynb"  
  Builds the master dataframes and intermediate products from raw inputs:
  - RINEX observation files → Parquet epoch–satellite tables
  - Daily / monthly DCB master tables (satellite and receiver)
  - Navigation master table from broadcast ephemerides
  - Daily receiver bias estimates (WLS and MS)
  - Co-located GIM VTEC time series and daily offsets (ΔVTEC)
  - Solar and geomagnetic indices (F10.7, Kp, Dst, etc.)
  - GFZ Q/D labels and derived quiet flags (QSL–GIM20)

- "00_ZENODO_TEC_Plots.ipynb"  
  Generates all figures and key diagnostics used in the JGR: Space Physics manuscript and in the Zenodo dataset description, including:
  - Diurnal, monthly, and seasonal climatologies
  - Slope diagnostics and EIA-crest relation
  - Histograms of daily maximum time and amplitude
  - Boxplots, ECDFs, exceedance curves, and yearly anomalies

Core parsers:

- "rinex2_obs_lean.py"  
  Minimal RINEX-2 GPS observation parser for VTEC work. Returns epoch–satellite "P2 - C1" time series (in meters) as a pandas DataFrame.
- "nav_rinex2_to_df.py"  
  RINEX-2 navigation parser. Reads broadcast ephemerides into a DataFrame (satellite orbital parameters, clocks).
- "latest_bsx_dsb_dcb_sat_c2w-c1c_parser_corr.py"  
  SINEX-BIAS / DCB parser. Extracts satellite differential code biases (C2W–C1C) for use in the P2–C1 observable model.

## Requirements

This code targets Python 3.10–3.11 and relies on standard scientific Python packages.

Recommended installation is via Conda using the provided "environment.yml":

"""bash
conda env create -f TEC_environment.yml
conda activate ouca_vtec

The environment includes:Python, NumPy, pandas, Matplotlib, Seaborn, pytz, pyarrow, Jupyter / IPython

## Usage
- Prepare raw inputs
	Organize the raw data locally (paths must match those configured inside the notebooks):
	+ RINEX-2 observation files (*.YYo or *.YYO)
	+ RINEX-2 navigation files (*.YYn or *.YYN)
	+ SINEX-BIAS / DCB products (daily and/or monthly)
	+ Global Ionospheric Maps (IONEX)
	+ Daily solar and geomagnetic indices (F10.7, Kp, Dst, etc.)
	+ GFZ/ISGI quiet/disturbed day lists, where used
- Generate master files
	Open and run:
	+ 00_ZENODO_Generate_Master_files.ipynb
	This notebook:
	+ Parses RINEX, navigation, and DCB products
	+ Estimates daily receiver biases (WLS and MS)
	+ Applies mapping to VTEC
	+ Builds:
		* 30-minute VTEC master file with columns: time, VTEC_median, VTEC_mean, VTEC_q25, VTEC_q75, N, vtec_gim, date_utc
		* Daily VTEC master file with columns: date_utc, VTEC_median, VTEC_mean, VTEC_min, VTEC_max, VTEC_q25, VTEC_q75, VTEC_std, N, VTEC_max_from_		30min, max_bin_utc, max_time_utc_str, max_hour_utc, max_ts_utc, gim_offset_tecu, offset_flag_abs_ge_20, mitigation_used, f107_obs, f107_adj, 		kp_daily_mean, kp_daily_max, Ap, SN, Dst_min, Dst_mean, Dst_n, date, D_rx_m_MS, N_epochs, Nsat_mean, solar_label, gfz_label, gfz_rank, 			gfz_flag, gfz_code, geomag_label_gfz_QDNQ
- Reproduce figures and diagnostics
	Open and run:
	+ 00_ZENODO_TEC_Plots.ipynb
	This notebook reads the master files and produces the figures and numerical diagnostics used in the manuscript and dataset description, including:
	+ Diurnal climatologies (monthly, seasonal)
	+ Boxplots by local-time bin
	+ Daily maximum-time and amplitude distributions
	+ Yearly anomalies relative to the quiet-time baseline
	+ ECDFs and exceedance probabilities
- Path configuration
	Both notebooks assume specific directory layouts for raw and processed data (e.g., local paths to RINEX, DCB, IONEX, and output folders). Before running, update these paths in the notebooks to match your local file system.

## Related data
	The processed VTEC datasets produced by this code are archived at Zenodo:

	Kaab, M., et al. (2025). Quiet-time GNSS-based VTEC climatology over Oukaimeden (Morocco), 2015–2025 (QSL–GIM20) [Data set]. Zenodo. 	https://doi.org/10.5281/zenodo.17144341

	The Zenodo record contains:
	+ 30-minute VTEC master file for OUCA
	+ Daily VTEC master file with solar/geomagnetic indices and flags
	+ QSL–GIM20 day list and supporting metadata

## Citation

If you use this code or the corresponding data products, please cite:
- The Zenodo dataset:
Kaab, M., et al. (2025). Quiet-time GNSS-based VTEC climatology over Oukaimeden (Morocco), 2015–2025 (QSL–GIM20) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.17144341
- And, once published, the associated JGR: Space Physics article (full reference to be added).
## License
This repository is released under the MIT License. See the LICENSE file for details.
