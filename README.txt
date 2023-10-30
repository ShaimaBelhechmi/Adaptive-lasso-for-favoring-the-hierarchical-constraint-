Source code for manuscript "Favoring the hierarchical constraint in penalized survival models for randomized trials in precision medicine" Shaima Belhechmi, Gwénaël Le Teuf, Riccardo De Bin, Federico Rotolo & Stefan Michiels. https://doi.org/10.1186/s12859-023-05162-x  

To replicate the results presented in the article, please execute the scripts located in the following directory: "./functions."

runsims.R: This script generates datasets based on the four scenarios outlined in the article.
The datas folder contains an example dataset for each of scenarios 1, 2, 3, and 4.

AL_LRT.R, AL_Wald.R: These scripts implement the Likelihood Ratio Test (LRT) and Single Wald weighting strategies for the Adaptive Lasso method.
Analyse_AL_LRT.R, Analyse_AL_Wald.R: These scripts apply the AL_LRT.R and AL_Wald.R functions to the simulated datasets.

SGL.R, GEL.R, cMCP.R, and ALRidge.R: These scripts implement different methods, namely SGL, GEL, cMCP, and ALRidge.
Analyse_SGL.R, Analyse_GEL.R, Analyse_cMCP.R, and Analyse_ALRidge.R: These scripts apply the SGL.R, GEL.R, cMCP.R, and ALRidge.R functions to the simulated datasets.

For questions, comments or remarks about the code please contact S. Michiels (stefan.michiels@gustaveroussy.fr) or S. Belhechmi (belhechmishaima@gmail.com).
