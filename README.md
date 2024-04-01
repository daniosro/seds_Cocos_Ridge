## Repository Architecture

The repository is split into four main directories. Please see each directory for information about each file.

### **`data`**
 All raw data collected from the analytical measurements in this work. 
 
 * **`raw`** 
 Includes the raw data, subdivided as:
 
 1. EA_CO2: blanks, linearity, and main data files for the measurements of TOC and C isotope composition of organic matter.
 
 2. EA_SO2_sulfate: 
 
 - blanks, linearity, and main data files for the measurements of S isotope composition of sulfate.
 - main data files for the measurements of O isotope composition of sulfate.
 
 3. ICP-MS: trace metal porewater compositions quantified by ICP-MS.
 
 4. quantifications: measurements of:
 -SR2113_IC: major ions quantified by IC.
 -S_intermediates: reduced sulfur species quantified by XEVO-TOF-UPLC with an MS and FLR detector.
 -POC fluxes: fluxes of POC measured on board of the SR2113 cruise.
  -Oxygen profiles: water column oxygen concentrations measured by CTD on the SR2113 cruise.
 -Elemental sulfur: solid-phase elemental sulfur concentrations measured with UV-HPLC.

 5. 16s_rRNA_seqs: 16S rRNA sequences for 4 sediment cores in the Eastern Pacific, one from the San Clemente Basin and 3 from Cocos Ridge.

 6. Extras: SR2113 sampling station coordinates and file used to construct the map.
 
* **`processed`** 

Includes the processed data, subdivided as:
 
 1. EA_CO2: measurements of TOC and C isotope composition of organic matter.
 
 2. EA_SO2_sulfate: measurements of S isotope composition of sulfate.

 3. ICP-MS: trace metal porewater compositions quantified by ICP-MS.

 4. 16s_rRNA_seqs: cleaned 16S rRNA sequences (primers removed with cutadapt)
 
### **`code`**
Files with the code to process the data generated in this work and reproduce the plots, subdivided as:

* **`analysis`**
Contains the R scripts used to analyze the sediment 16S rRNA data.

* **`processing`**
Contains the scripts used to generate the map with the SR2113 stations, and process the quantification and isotope data.

* **`figures`**
Contains the scripts used to generate the figures in the main text and the supplementary material.

### **`figures`**
Figures generated from the geochemical measurements in this work.

 

# License Information

All creative works (code, figures, etc) are licensed under the [Creative
Commons CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/) license. This work is published from: United States.
