The metadata for the two data files in this folder are as follows:

## germination_data.csv
* Depth_Top = depth (cm) from the marsh surface to the top of a slice of a soil core from which the seeds were recovered from
* Assay_Num = (1-13) experimental assay that seeds were germinated in. assay numbers align with numbers in the Supplementary Materials in Vahsen et al. 
* Core_Location = (1-11) provenance of the core from which seeds were generated. core location numbers align with Figure 1 and Table 1 in Vahsen et al. (1 = Kirkpatrick Marsh, 2 = Corn Island, 3 = Hog Island, 3 = Virginia, 5 = Bay Bridge, 6 = Eastern Shore, 7 = Blackwater, 8 = Taylor Island, 9 = Delaware Bay, 10 = Sellman Creek, 11 = Greenhouse. Note: 11 is not on the map because seeds were collected from extant plants within a greenhouse.
* Assay_Date = year in which the experimental assay was conducted (2003-2019)
* Temperature = temperature regime that the seeds underwent (1 = constant 30&deg;C, 2 = fluctuating 20/15&deg;C, 3 = fluctuating 27/15&deg;C, 4 = constant 25&deg;C). Note: temperatures fluctuating in sync with photoperiod daytime/nighttime schedules.
* Media = medium on which seeds were germinated (1 = sand/soil mix, 2 = sand, 3 = growth media [Murashige and Skoog salt and vitamin, sucrose, and agar mixture]).
* Photoperiod = number of daytime and nighttime hours that seeds experiments in an experiment (1 = 15 daytime/9 nighttime, 2 = 12 daytime/12 nighttime, 3 = all dark)
* Treatment = indicator variable that indicates whether seeds experienced a pre-treatment (see Supplementary Materials for details) or not (0 = no, 1 = yes)
* Media1 = indicator variable to test for the effect of Media 2 vs Media 1 (reference level)
* Media2 = indicator variable to test for the effect of Media 3 vs Media 1 (reference level)
* Temp1 = indicator variable to test for the effect of Temp 1 vs Temp 3 (reference level)
* Temp2 = indicator variable to test for the effect of Temp 3 vs Temp 3 (reference level)
* Temp3 = indicator variable to test for the effect of Temp 4 vs Temp 3 (reference level)
* Photo1 = indicator variable to test for the effect of Photo 1 vs Photo 3 (reference level)
* Photo2 = indicator variable to test for the effect of Photo 2 vs Photo 3 (reference level)
* Germinated = the number of seeds that germinated for a given unique trial
* Seeds_Planted = the number of seeds that were planted for a given unique trial
