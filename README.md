Code repository "Confounding Contact: Bias Induced by Post Randomization Changes in Staff Behavior in Infection Prevention Trials"
----

The Python scripts in this repository require the 'os', 'numpy', 'stochpy' and 'random' libraries, all of which are available via pip. The R script requires the 'vioplot' library.

Contents:
----
* 'contact_mrsa.py': This script generates the primary results from the manuscript

* 'contact_mrsasensitivity.py': One of two scripts for the sensitivity analysis reported in the manuscript. This script generates the main results of the model allowing all parameters to vary +/-10% of their original values.

* 'contact_mrsasweep.py': The second of two scripts for the sensitivity analysis reported in the manuscript. This script varies the contact rate changes over a range from 0 to 50%.

* 'analysis.R': The R script for generating the statistical results and plots found in the manuscript.