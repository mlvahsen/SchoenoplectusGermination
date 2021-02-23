The metadata for the data files in this folder are as follows:

## Vahsen_etal_Pb210.csv
* **depth_cm** = middle depth (in cm) of the 2 cm soil layer				
* **core_year** = year that the core was extracted from Kirkpatrick Marsh				
* **core_month** = month that the core was extracted from Kirkpatrick Marsh				
* **sample** = core id (MX1, MX2, or MX4). MX stands for 'mixed' plots from a long-term elevated CO2 experiment (see Drake *et al.* (2014, doi:10.1111/gcb.12631) for details). Mixed plots contained a mix of two dominant marsh plant species *Schoenoplectus americanus* and *Spartina patens* at the start of the elevated CO2 experiment in 1987.
* **bulk_dens_210pb** =	bulk density g cm-3				
* **acc_rate** = calculated soil accretion rate using constant rate of supply model cm yr-1				
* **year** = calculated soil year (AD)				
* **se** = calculated soil year (AD) standard error	

## assay3_pretreatment.csv
* **ID** = unique ID for a seed
* **Assay_Num** = experimental assay. here these are all from assay = 3 which maps to assay 3 in the Supplementary Material in Vahsen et al.
* **Assay_Date** = year that the experimental assay was conducted. here these are all from 2010
* **Depth_Top** = depth (cm) from the marsh surface to the top of a slice of a soil core from which the seeds were recovered from
* **Core_Location** = (1-11) provenance of the core from which seeds were generated. core location numbers align with Figure 1 and Table 1 in Vahsen et al. (1 = Kirkpatrick Marsh, 2 = Corn Island, 3 = Hog Island, 3 = Virginia, 5 = Bay Bridge, 6 = Eastern Shore, 7 = Blackwater, 8 = Taylor Island, 9 = Delaware Bay, 10 = Sellman Creek, 11 = Greenhouse.
* **Temperature** = temperature regime of the experiment (here all are 30&deg;C)
* **HighTemp** = maximum temperature that seeds experience during the experiment (here all are 30&deg;C)
* **Light_Period** = number of daytime hours that the seeds experience (here all are 15)
* **Treatment** = one of six possible pre-treatments (Auxin, Bleach + Water, Ethylene, Gibberellic Acid, Water Bath, or Bleach) or no pre-treatment (None)
* **Germinated** = binary indicator for whether a seed germinated or not (0 = did not germinate, 1 = germinated)
