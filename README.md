Code developed in SAS v9.4 <br>
Simulations used to estimate power for global rank composite versus decision rule, mimicking CITRIS-ALI.

**The global rank** <br>
Global Rank is a generic term encompassing a number of decision rules. The intention is to assign every patient a rank among the other patients that summarises their response.
The following equation describes the calculation.

The global rank composite is seen in the medical literature. For example, see:<br>
--Califf 1990:  "[Left Ventricular Ejection Fraction May Not Be Useful as an End Point of Thrombolytic Therapy Comparative Trials](https://pubmed.ncbi.nlm.nih.gov/2225381/)" 
--Felker & Maisel 2010: "[A Global Rank End Point for Clinical Trials in Acute Heart Failure](https://www.ahajournals.org/doi/full/10.1161/circheartfailure.109.926030)"

Input required for each outcome includes a base rate, effect due to intervention, and correlations among outcomes;
although from our experience, the correlations have an unimportant effect on power estimates.
We must also specify the order of outcomes for the hierarchy and cut-off's defining when to move to the 
next outcome in the hierarchy (i.e. cut-offs defining "failure").

Effect sizes and correlations can be found in simul_data.sas
Ordering of outcomes and cut-offs can be found in derive_GR.sas