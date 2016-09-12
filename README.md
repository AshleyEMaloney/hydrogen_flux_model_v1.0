# hydrogen_flux_model_v1.0
Stationary 6-box Matlab model of hydrogen and deuterium cycling in phytoplankton for testing the sensitivity of δ2Hlipid to metabolic parameters. 

Please cite the following paper if you use this code:

A.E. Maloney, A.L.C. Shinneman, K. Hemeon, and J.P. Sachs (2016). Exploring lipid 2H/1H fractionation mechanisms in response to salinity with continuous cultures of the diatom Thalassiosira pseudonana. Organic Geochemistry. http://dx.doi.org/10.1016/j.orggeochem.2016.08.015

There are 6 individual "Metabolic Parameter" scripts: use these for tinkering with with the model if you want to change alphas or flux porportionality assumptions since there is a mass balance check for each pool of hydrogen (if the fluxes into each pool do not equal the fluxes out of each pool than you have violated the steady-state assumption). A longer script is provided that was used for making the figures. Additionally there is a schematic to help illustrate what the model is doing and a "cheat sheet" to pair the numbers used in the script to each flux.

# Software
MatLab student version R2013a

# Acknowledgements
XXX vetted this code release. We also appreciate the assistance of Allison Smith.

This material is based upon work supported by the U.S. National Science Foundation under Grant No. OCE-1027079 (J.P.S.). The University of Washington Program on Climate Change Fellowship and IGERT Ocean Change Fellowship award #NSF1068839 provided partial support for A.E.M. We would like to thank P. MacCready for lending lab space, E. V. Armbrust and A. Ingalls for lending equipment, F. Ribolet, G. Hennon, L. Fisher, and O. Kawka for assisting with culture data interpretation, A. Chomos and N. Kato for assisting with culture maintenance, J.A. Gregersen for assisting with water isotope measurements, C. Paschall and R. D’Jay for assisting with Coulter Counter measurements, R.A. Cattolico and H. Hunsperger for assisting with flow cytometry measurements, D.B. Nelson and J.N. Richey for assisting with lipid analysis, and S.G. Warren and M.D. Wolhowe for writing edits. Thanks to K.A. Krogslund and A. Morello who analyzed nutrient and chlorophyll samples. Valuable discussions with S.R. Smith, J. Levering, A. Gruber, M.J. Behrenfeld, G. Rocap, E.V. Armbrust, P.D. Quay, L.A. Thompson, A. Ingalls, K.R. Heal, C. Deutch, N.S. Banas, A.E. Shao, C. Saenger, and L. Johnson aided model development. The comments of two anonymous reviewers and Associate Editor S. Schouten improved the manuscript.
