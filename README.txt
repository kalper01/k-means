Kenny Alperin
COMP 135
Assignment 3

README

All of the code for the assignment is located in kmc.java. In that file, each of
the functions is described and the code is extensively commented. The program is 
configured so that running it once does each experiment for each data set. To run the
program, in the top level directory enter in the command line: 
	java -jar kmc.jar.
Note that the names of the output text files are defined in the clustering function in 
kmc.java. The output files for 2-2 will end in '22.txt' and the output files for 2-3 will
end in '23.txt'. For 2-2, each line of the output file will include the run number, the 
CS value, and the NMI value. For 2-3, each line of the output file will include the
k value and the CS value. All of these output files will be generated in /output.
Assignment 3 Writeup.docx details the results from the experiments, and includes plots
generated from plotAssignment3.m. Note the data for the plots was manually entered into
the MATLAB script for section 2-2, but the data was parsed from the output text files
for section 2-3.
