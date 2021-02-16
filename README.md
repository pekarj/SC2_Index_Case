# Title: Timing the SARS-CoV-2 Index Case in Hubei Province using Coalescent Approaches
Authors: Jonathan Pekar<sup>1,2</sup>, Michael Worobey<sup>3,\*</sup>, Niema Moshiri<sup>4</sup>, Konrad Scheffler<sup>5</sup>, and Joel O. Wertheim<sup>6,\*</sup><br />

Affiliations:<br />
<sup>1</sup>Bioinformatics and Systems Biology Graduate Program, University of California San Diego, La Jolla, CA 92093, USA<br />
<sup>2</sup>Department of Biomedical Informatics, University of California San Diego, La Jolla, CA 92093, USA<br />
<sup>3</sup>Department of Ecology and Evolutionary Biology, University of Arizona, Tucson, AZ 85721, USA.<br />
<sup>4</sup>Department Computer Science & Engineering, University of California San Diego, La Jolla, CA 92093, USA<br />
<sup>5</sup>Illumina, Inc., San Diego, CA 92122, USA<br />
<sup>6</sup>Department of Medicine, University of California San Diego, La Jolla, CA 92093, USA<br />

<sup>\*</sup>Corresponding authors. Email: worobey@arizona.edu (MW); jwertheim@health.ucsd.edu (JOW)

The main and sensitivity analysis JSON files are used to run the initial FAVITES simulations. Then, for the established epidemics, the VirusTreeSimulator JSON file must be edited to include the proper paths to the contact and transmission networks in the FAVITES outputs from the initial simulations. "TF" stands for transmission factor (e.g., "0.33xTF" is 1/3 the transmission rate of the main analysis), "R" is ascertainment rate (e.g., "2xR" is 2-times the ascertainment rate of the main analysis), and "hosp" indicates modulated hospitalization rate (e.g., "hosp2x" is 2-times the hospitalization rate). There can be more than one of these modulations per CONFIG file, as indicated in the name of the file. Please refer to the supplementary tables in the manuscript for further details.  

The XML files are the BEAST inputs for the strict, relaxed clock, and Skygrid analyses. 
