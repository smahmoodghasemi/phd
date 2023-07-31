
First, we need to use the "Image Segmentation" file to process the raw breast cancer cell  images. The outcome of this process is a Python dictionary saved as a pickle file. 
To get the table of Qi's probability of activation for a single or a pair of alleles and mutual information, we used the pickle file from the "Image Segmentation." 
Next, we used the "Post Processing of Image Segmentation" Python file. The outcome of this process is a data frame that you can save as CSV or XLSX file.
Furthermore, for simulating our model for the population of cells, we used the "Simulation for Flavopiridol experiments" to simulate the model and the "Simulation for Non-Flavopirdol experiments" for the extended model.
This simulation's outcome is also a Python dictionary saved as a pickle file. Then, we used the "Post Processing of Flavopiridol Simulation" and the "Post Processing of Non-Flavopirdol Simulation" to fit the result of the simulation to the probability of activation for a single allele extracted from the image processing.
Finally, we used the "Simulation of Independent Model" file to simulate the possible solution for the independent model based on the independent equation, introduced in the paper.
It requires only the CSV files of Qi to show that these two values do not have any overlapping. 
