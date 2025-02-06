# IAV_Sources_and_Sinks

Data and R code used to run analyses for the manuscript "Sources and sinks of influenza A virus genomic diversity in swine from 2009 to 2022 in the United States". This manuscript is in prep.

## Description

This manuscript summarizes 13 years of influenza A virus evolution in the swine host, focusing on the interplay between ecological concepts and swine agricultural practices, and offers tools for modeling the spread of novel virus strains across the US.

## Getting Started

### Dependencies

* These scripts were run using R version 4.2.0 (2022-04-22 ucrt) on Windows x86_64.

## Outline

### Script 1
This script brings together a lot of data, cleans it, formats it, and produces a few data frames that the following scripts will use. This script must be run once to build the `IAV_Sources_and_Sinks.RData` file, which all following scripts will load into memory before progressing. Some summary analyses are also run in this script. While not strictly necessary, I recommend clearing the environment before moving on to Script 2 or Script 3, either by restarting RStudio or by using `rm(list=ls())`, as it will clear resources.

### Script 2
This script runs the Zeta diversity analyses and produces the associated plots. Core functionality is provided with the [`zetadiv` package](https://cran.r-project.org/web/packages/zetadiv/index.html). You must load `IAV_Sources_and_Sinks.RData` from Script 1 to run this script. This script saves nothing to the .RData file.

### Script 3
This script establishes a time series of detections and runs analyses and looks at temporal and geographic patterns of detections. You must load `IAV_Sources_and_Sinks.RData` from Script 1 to run this script.

### Script 4
This script uses Markov Chains to model state-to-state IAV transition probabilities, provides code for testing real and hypothetical IAV strain detection orders, and statistics to evaluate model prediction performance. You must load `IAV_Sources_and_Sinks.RData`, modified and saved in Script 3, to run this script.

### Microreact
The file `microreact_input_5y.csv` is written by script 1. The microreact figure in this manuscript can be made by following these steps:
 * Go to [the microreact website](https://microreact.org/upload).
 * Drag and drop the file into the page.
 * Select `continue`.
 * Select `continue`.
 * To make the legend, select `Legend` to pull up the Legend, download it by clicking on the three bars in that window and downloading it as PNG.
 * To make the timeline, select `Create new Timeline` from the three bars from the pencil icon at the top right. Change the Temporal Data Type to `One column: Formatted Values`, change the Temporal Data Column to `Date`, and click `close`. Download it by clicking on the three bars in that window and downloading it as PNG.
 * To format the map, select the sliders icon in that window and adjust as desired.
   
## Authors
Garrett M. Janzen
[@garrettjanzen](https://twitter.com/garrettjanzen)\
Blake T. Inderski\
Jennifer Chang\
Zebulun W. Arendsee\
Alicia Janas-Martindale\
Mia Kim Torchetti\
Amy L. Baker\
Tavis K. Anderson

## Acknowledgments
Thanks to co-authors Blake T. Inderski, Jennifer Chang, Zebulun W. Arendsee, Alicia Janas-Martindale, Mia Torchetti, Amy L. Baker, Tavis K. Anderson, and to all members of the Flu-Crew at the National Animal Disease Center, USDA-ARS.

