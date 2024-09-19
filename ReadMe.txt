This is a repository of the code for this paper: https://www.fortunejournals.com/articles/analysis-and-modeling-of-early.pdf.

We used the "Image_Processing" file to process the raw breast cancer cells, including the DAPI, Intraonic, and exonic channels. The outcome of this process is a Python dictionary saved as a pickle file. We used the pickle file from the "Image_Processing" to get the table of $Q_{i}$'s probability of activation for a single or a pair of alleles, compute mutual information, and perform further statistical analysis.

Next, we used the "Statistical_Analysis" Python file to compute $Q_{i}$ and mutual information. The outcome of this process is a data frame that you can save as a CSV or XLSX file for further visualization. 

Furthermore, for simulating our model for the population of cells, we used the "FV+E2_experiments_prediction" to simulate the model and the "E2_experiments_prediction" for the extended model where we have transcription before time 0. This simulation's outcome is also a Python dictionary saved as a pickle file. 

Then, we used the "FV+E2_Simulation_Analysis" and the "E2_Simulation_Analysis" to fit the result of the simulation to the probability of activation for a single allele extracted from the image processing. 

Finally, we used the "Independent_Model" file to simulate the possible solution for the "independent model" based on the independent equation (Equation 2), introduced in the paper. It requires only the CSV files of $Q_{i}$ to show that these two values do not have any overlapping. 
