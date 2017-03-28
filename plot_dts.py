# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 12:27:29 2016

@author: ANI

This program compiles information from all tra-files into one file, and
then plots a filled contour plot from the information obtained.

If there are no tra-files available, it reads the data directly from
the text file, which should be located in the same directory.

It relies on a document dts_functions.py.

User defined parameters are:

*** general parameters ***
    language (str): language used in file names and on the plot.
        Choose 'en' for English or 'dk' for Danish.
    commas (int): defines decimal point in the file with the results.
        Choose 1 for comma and 0 for dot.
    langs (dict): dictionary with basic terms used for plotting and
        file names in the available languages.

*** contour plot ***
    minidno (int): defines the minimum ID number for which measurements
        should be included. All ID numbers below minidno are disregarded,
        minidno is still included on the graph.
        Change of this variable will affect how contour plot is defined.
    maxidno (int): defines the maximum ID number for which measurements
        should be included. All datapoints with ID number higher than maxidno
        are excluded from the contour plot (the value itself is included).
        If maxidno is higher than the total number of measurement points,
        all points are automatically included.
        Change of this variable will affect how contour plot is defined.
    timestep (int): defines how many labels should be skipped on the graph.
        timestep = n means that every nth label is included.
        For example, if n = 12, every 12th label is included corresponding
        to one y-tick for every two hours (with measurements taken every
        10 minutes).
    stepcb (float): the difference in Celsius degrees between
        individual ticks on the colorbar.
    clevels (float): number of contour levels, the higher the number,
        the smoother the transition between contours. clevels = 500
        gives a high quality plot.
    changelimits (list): defines whether user defined limits should be used
        for x-axis. If empty, limits on the x-axis are adjusted
        automatically. If not empty, the elements of the list are taken as
        xmin and xmax, respectively. The list should have either no elements,
        or 2 elements. Example: changelimits = [8, 1017].
        Change of this variable does not affect the contour plot itself,
        only which part of it should be visible on the graph.
    udef_xticks (list): list in which all elements are floats. User defined
        x-ticks. If empty, the ticks are applied automatically.
        Example: udef_xticks = [200, 400, 1000].
"""

# ### user defined input ####
# ## general
language = 'en'
commas = 0

# ## contour plot
minidno = 380
maxidno = 1002
timestep = 12
stepcb = 1
clevels = 10
changelimits = []
udef_xticks = []


# ## define labels in different languages
langs = dict()
# 0, 1, 2: used in file saving
# 3: label on a colorbar
# 4, 5: x and y axis
# 6, 7: plot title
langs['dk'] = ["gemt_den", "komma", "alle", "temperatur", "kabellængde (m)",
               "tid (hh:mm)", "start", "stop"]
langs['en'] = ["saved_on", "comma", "all", "temperature", "cabel length (m)",
               "time (hh:mm)", "start", "stop"]


##############
##############
##############

# ## import packages
import numpy as np
import matplotlib.pyplot as plt
import os
import time

import dts_functions as dfun

plt.ion()

tic = time.time()

# find all .tra-files from the current directory
# (the python script needs to be saved in the same folder as datafiles)
trafiles = [f for f in sorted(os.listdir(".")) if ".tra" in f]
txtfile = [f for f in sorted(os.listdir(".")) if ".txt" in f]


def get_info_txt(txtfile):
    return txtfile


if len(trafiles) > 0:
    print("read data from files")
    # initialise the files
    (eltime, timelabels, tempstack, tempstack_comma, lines,
     dataheader, savedate, initialtime) = dfun.get_data_initfile(trafiles[0])

    for trafile in trafiles[1:]:
        (eltime, timelabels, tempstack, tempstack_comma,
         dataheader) = dfun.get_data(trafile, dataheader, eltime,
                                     timelabels, initialtime,
                                     tempstack, tempstack_comma)

    print("save the results (text file)")

    idno = np.arange(0, lines+1)

    saveddata, savenamefig = dfun.savefile(lines, tempstack, tempstack_comma,
                                           commas, langs, language, savedate,
                                           dataheader, idno)

else:
    print("read data from the drive")

    (idno, eltime, tempstack, dataheader, timelabels,
     savenamefig) = dfun.get_data_from_drive(txtfile, langs, language)


print("plot and save the results (figure)")


fig, colorbar = dfun.plot_and_save(idno, eltime, tempstack,
                                   dataheader, timelabels, savenamefig,
                                   minidno, maxidno, timestep, stepcb,
                                   clevels, changelimits,  udef_xticks,
                                   langs, language)

toc = time.time()

print("***** end of script, running time: %.1f s *****" % (toc - tic))
