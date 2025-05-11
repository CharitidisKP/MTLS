"""
Small script made to plot histograms of ORF lengths predicted by glimmer 
"""

# import matplotlib for making the histogram. 
import matplotlib.pyplot as plt

''' 
TODO:  Change the names of the input and output files here, so that the script uses your glimmer predictions as input,
and outputs a histogram named according to the input genome predictions.
'''
glimmer_pred = "<genome>.glimmer.predict"
out_file = "ORFlength_<genome>_hist.png"

# start an empty list that we can store genelengths in
orf_lengths = list()

# read file and loop over each line
with open(glimmer_pred, "r") as fil:
    # skip the first line to avoid the >name line
    next(fil)
    for line in fil:
        # split the line on whitespace. Note that glimmer returns a multispace seperated file not "\t" or ","
        # Note that python str.split() will split on and remove all whitespace from a string 
        Sline = line.split()

        # calculate the length of the predicted orf
        orf_len = float(Sline[2])-float(Sline[1])
        # Make the length absolute. Since some orfs are on the reverse strand, we do not want negative values.
        orf_len = abs(orf_len)
        # append the orf length to the list of lengths 
        orf_lengths.append(orf_len)

'''
TODO: use matplotlib to plot the orf_lengths in a histogram
'''

'''
TODO: Save the histogram to a png file named as the out_file variable
'''


