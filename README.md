# IAV_Sources_and_Sinks

Data and R code used to run analyses for the manuscript "Sources and sinks of influenza A virus genome diversity in swine in the USA from 2009 to 2022".

## Description

This manuscript summarizes 13 years of influenza A virus evolution in the swine host, focusing on the interplay between ecological concepts and swine agricultural practices, and offers tools for modeling the spread of novel virus strains across the US.

## Getting Started

### Dependencies

* These scripts were run using R version 4.2.0 (2022-04-22 ucrt) on Windows x86_64.


## Outline

### Script 1
This script brings together a lot of data, cleans it, formats it, and produces a few data frames that the following scripts will use. This script must be run once to build the .RData file, which all following scripts will load into memory before progressing. Some summary analyses are also run in this script.



### Script 2
This script runs the Zeta diversity analyses and produces the associated plots.

* How to run the program
* Step-by-step bullets
```
code blocks for commands
```

### Microreact
The file microreact_input_5y.csv is written by script 1. The microreact figure in this manuscript can be made by following these steps:
 * Go to [the microreact website](https://microreact.org/upload).
 * Drag and drop the file into the page.
 * Select `continue`.
 * Select `continue`.
 * To make the legend, select `Legend` to pull up the Legend, download it by clicking on the three bars in that window and downloading it as PNG.
 * To make the timeline, select `Create new Timeline` from the three bars from the pencil icon at the top right. Change the Temporal Data Type to `One column: Formatted Values`, change the Temporal Data Column to `Date`, and click `close`. Download it by clicking on the three bars in that window and downloading it as PNG.
 * To format the map, select the sliders icon in that window and adjust as desired.
   
## Authors

Contributors names and contact info

ex. Garrett Janzen
ex. [@garrettjanzen](https://twitter.com/garrettjanzen)

## Acknowledgments

Inspiration, code snippets, etc.
* [awesome-readme](https://github.com/matiassingers/awesome-readme)
* [PurpleBooth](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)
* [dbader](https://github.com/dbader/readme-template)
* [zenorocha](https://gist.github.com/zenorocha/4526327)
* [fvcproductions](https://gist.github.com/fvcproductions/1bfc2d4aecb01a834b46)
