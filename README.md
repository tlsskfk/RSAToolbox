# <center>RSA Toolbox

Introduced by [David Rothlein](david.rothlein@gmail.com)

Contributors: Jeremy Purcell, Sam Rosenberg, Stephen Shin

Maintainers: [Stephen Shin](shin.skfk@gmail.com), [David Rothlein](david.rothlein@gmail.com)
</center>

## Directory Structure

Top level folder is for any uses/scripts of the main functions introduced in the toolbox.  They live here for simplicity and for sharing.

`Functions` is the folder in which the base level matlab functions will live.

`fmriprep_table` is the folder in which a preliminary dataset should be created

## Usage
Create a matlab script that will use `./Functions` as a path (selected folders and subfolder).
Use script to process data.

## Suggested Workflow for NTR project
Use the `Full_Analysis_Word.m` to create a data table living in `fmriprep_table`
Use the created file for futher RSA analysis.  
