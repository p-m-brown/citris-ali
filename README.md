Code developed in SAS v9.4 <br>
Simulations used to estimate power for global rank composite versus decision rule, mimicking CITRIS-ALI.

**The global rank** <br>
Global Rank is a generic term encompassing a number of composites. The intention is to assign every patient a rank among the other patients that summarises their response.
The following equation describes the calculation.



The global rank composite is promoted in the medical literature. For example, see:<br>
- Califf 1990:  "[Left Ventricular Ejection Fraction May Not Be Useful as an End Point of Thrombolytic Therapy Comparative Trials](https://pubmed.ncbi.nlm.nih.gov/2225381/)" <br>
- Felker & Maisel 2010: "[A Global Rank End Point for Clinical Trials in Acute Heart Failure](https://www.ahajournals.org/doi/full/10.1161/circheartfailure.109.926030)" <br>
- Packer 2016 "[Development and Evolution of a Hierarchical Clinical Composite End Point for the Evaluation of Drugs and Devices for Acute and Chronic Heart Failure](https://www.ahajournals.org/doi/10.1161/circulationaha.116.023538)"

Input required for each outcome includes a base rate, effect due to intervention, and correlations among outcomes;
although from our experience, correlations have an unimportant effect on power estimates.
We must also specify the order of outcomes for the hierarchy and cut-off's defining when to move to the 
next outcome in the hierarchy (i.e. cut-offs defining "failure").

- Effect sizes and correlations can be found in simul_data.sas
- Ordering of outcomes and cut-offs can be found in derive_GR.sas

For our simulations we used moderate effect sizes:
- Control: 
- Intervention:
- Initial working correlations: 

These guesstimates are based on:
- VITAMINS: "[Effect of Vitamin C, Hydrocortisone, and Thiamine vs Hydrocortisone Alone on Time Alive and Free of Vasopressor Support Among Patients With Septic Shock](https://jamanetwork.com/journals/jama/fullarticle/2759414)" <br>
- CITRIS-ALI: "[Effect of Vitamin C Infusion on Organ Failure and Biomarkers of Inflammation and Vascular Injury in Patients With Sepsis and Severe Acute Respiratory Failure](https://pubmed.ncbi.nlm.nih.gov/31573637/)"

These assumed effect estimates, cut-offs etc. can be modified to evaluate the effect this has on the power estimates.






















