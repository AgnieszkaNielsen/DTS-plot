# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 17:31:47 2016

@author: ANI
"""

str1 = '.'
str2 = ','

import numpy as np
import sys
import time
import matplotlib.pyplot as plt

# ### functions used in plot_dts.py


def get_info(trafile, initialtime):
    """
    Returns basic information about the input file.
    Individual files (XXX.tra) correspond to measurements taken
    at a specific date and time.

    Args:
        trafile (str): data file to be investigated (usually a .tra-file).
        initialtime (str): time of the first measurement.
            Acceptable format: "0" for the first file,
            otherwise "YYYY-MM-DD HH:MM:SS".

    Returns:
        str: date of the measurement, "DD-MM-YYYY"
        str: time of the measurement, "HH:MM:SS"
        float: time in hours since initialtime (needed for plotting)
        int: number of data points in the file (0-indexing)
        str: date of the measurement, "YYYYMMDD" (used in filename)
    """

    timestamp_markers = ["Date.Year", "Date.Month", "Date.Day",
                         "Time.Hour", "Time.Minute", "Time.Second"]

    vals = dict()

    fh = open(trafile)
    lines = 0
    for line in fh:
        sp = line.split(";")
        if sp[0] not in timestamp_markers:
            try:
                lines = int(sp[0])
                # lines with measurements
            except:
                continue
                # initial lines with file information
        else:
            # lines with timestamp
            var = sp[1].strip()
            if len(var) == 1:
                var = "0" + var  # add a preceding zero for one-digit values

            # add to the dictionary
            vals[sp[0].split(".")[1]] = var

    fh.close()

    date = vals["Day"] + '-' + vals["Month"] + '-' + vals["Year"]
    time = vals["Hour"] + ':' + vals["Minute"] + ':' + vals["Second"]

    savedate = vals["Year"] + vals["Month"] + vals["Day"]

    date_eltime = vals["Year"] + '-' + vals["Month"] + '-' + vals["Day"]
    testtime = date_eltime + ' ' + time

    eltime = elapsedtime(testtime, initialtime)

    return date, time, eltime, lines, savedate


def append_info(trafile, dataheader, eltime,
                timelabels, initialtime):
    """
    Append information from the current file to the existing lists.

    Args:
        trafile (str): data file to be investigated (usually a .tra-file).
        dataheader (str): names of individual columns in the final file
            with the results.
            First column is "ID", the consequent columns are timestamps.
            Columns are separated by tabs:
            "ID\tDD-MM-YYYY HH:MM:SS\t...\tDD-MM-YYYY HH:MM:SS".
            For each file that is read, timestamp is appended at the end
            of the existing dataheader.
        eltime (list of floats): list with the time elapsed for individual
            measurements since initialtime (in hours). For each file
            that is read, elapsed time is appended at the end
            of the existing eltime-list. Used for plotting.
        timelabels (list of strings): list with timestamps of individual
            measurements ("HH:MM"), used as y-label on the plot.
            For each file that is read, time label is appended at the end
            of the existing timelabels-list.
        initialtime (str): time of the first measurement.
            Acceptable format: "0" for the first file,
            otherwise "YYYY-MM-DD HH:MM:SS".

    Returns:
        list of floats: same as eltime,
            but with new values appended from the input file.
        list of strings: same as timelabels,
            but with new values appended from the input file.
        list of floats: temperature measurements.
        list of strings: temperature measurements, but with comma
            as a decimal point.
        int: number of data points in the file (0-indexing)
        str: same as dataheader,
            but with new values appended from the input file.
        str: date of the measurement, "YYYYMMDD" (used in filename)
    """

    # print("open the file: %s" % trafile)

    (current_date, current_time, elapsed_time,
     lines, savedate) = get_info(trafile, initialtime)

    header_timestamp = current_date + ' ' + current_time
    dataheader += '\t' + header_timestamp

    temp = []
    temp_comma = []

    fh = open(trafile)
    for line in fh:
        try:
            int(line.split(";")[0])
        except:
            continue

        sp = line.split(";")

        # append data to lists
        temp.append(float(sp[2]))
        temp_comma.append(sp[2].replace(str1, str2))
    fh.close()

    eltime.append(elapsed_time)

    timelabels.append(current_time[:5])

    templist = []
    templistco = []

    templist.append(temp)
    templistco.append(temp_comma)

    temp = templist
    temp_comma = templistco

    return eltime, timelabels, temp, temp_comma, lines, dataheader, savedate


def elapsedtime(testtime, initialtime):
    """
    Return time in hours since the beginning of the test.
    First measurement is taken at t=0.

    Args:
        testtime (str): time of the current measurement.
            Acceptable format: "YYYY-MM-DD HH:MM:SS".
        initialtime (str): time of the first measurement.
            Acceptable format: "0" for the first file,
            otherwise "YYYY-MM-DD HH:MM:SS".

    Returns:
       float: time in hours since initialtime.
    """

    if initialtime == '0':
        timediff = 0
    else:
        time0 = initialtime
        time1 = testtime

        timediff = np.datetime64(time1) - np.datetime64(time0)
        timediff.astype('timedelta64[h]')
        timediff = timediff/np.timedelta64(1, 'h')

    return timediff


def get_data_initfile(trafile):
    """
    Initialises the list with the results, as well as
    defines some key parameters needed for plotting and final file
    with the results.

    Args:
        trafile (str): data file to be investigated (usually a .tra-file).

    Returns:
        list: one element list with the time elapsed since first
            measurement (in hours); here dt = 0.
        list: one element list with the time label of the first
            measurement (used for plotting).
        array: numpy array with temperature measurements, size n x 1,
            where n is number of measurement points. All measurements
            are saved as floats.
        array: numpy array with temperature measurements, size n x 1,
            where n is number of measurement points. All measurements
            are saved as strings with commas as decimal points.
        int: number of measurements per file.
        str: beginning of the header for the data file. 'ID\tTIMESTAMP',
            where TIMESTAMP is the date and time of the first measurement.
            Format: 'DD-MM-YYYY HH:MM:SS'.
        str: date of the first measurement, used in the file name.
            Format: 'YYYYMMDD'.
        str: date and time of the first measurement, used to calculate
            time difference between the consecutive measurements.
            Format: 'DD-MM-YYYY HH:MM:SS'.
    """

    (eltime, timelabels, temp, temp_comma, lines, dataheader,
     savedate) = append_info(trafile, "ID", [], [], '0')

    tempstack = np.column_stack((temp))
    tempstack_comma = np.column_stack((temp_comma))

    # savedate = dataheader.split("\t")[1].split()[0]
    initialtime = dataheader.split("\t")[1]
    initdate, inittime = initialtime.split()
    initday, initmonth, inityear = initdate.split("-")
    initialtime = inityear + '-' + initmonth + '-' + initday + ' ' + inittime

    return (eltime, timelabels, tempstack, tempstack_comma,
            lines, dataheader, savedate, initialtime)


def get_data(trafile, dataheader, eltime, timelabels,
             initialtime, tempstack, tempstack_comma):
    """
    Appends data to the lists initialised by get_data_initfile(...).

    Args:
        trafile (str): data file to be investigated (usually a .tra-file).
        dataheader (str): beginning of the header for the data file.
            'ID\tTIMESTAMP', where TIMESTAMP is the date and time
            of the first measurement. Format: 'DD-MM-YYYY HH:MM:SS'.
        eltime (list): one element list with the time elapsed since first
            measurement (in hours); here dt = 0.
        timelabels (list): one element list with the time label of the first
            measurement (used for plotting).
        initialtime (str): date and time of the first measurement,
            used to calculate time difference between the consecutive
            measurements. Format: 'DD-MM-YYYY HH:MM:SS'.
        tempstack (array): numpy array with temperature measurements,
            size n x 1, where n is number of measurement points.
            All measurements are saved as floats.
        tempstack_comma (array): numpy array with temperature measurements,
            size n x 1, where n is number of measurement points.
            All measurements are saved as strings
            with commas as decimal points.

    Returns:
        list: as eltime, but with new information appended
            (time elapsed since first measurement, in hours).
        list: as timelabels, but with new information appended.
            Time label of the first measurement (used for plotting).
            Format: 'HH:MM'.
        array: as tempstack, but with new information appended.
            Numpy array with temperature measurements, size n x m,
            where n is number of measurement points, and m is number
            of measurement times. All measurements are saved as floats.
        array: as tempstack_comma, but with new information appended.
            Numpy array with temperature measurements, size n x m,
            where n is number of measurement points, and m is number
            of measurement times. All measurements are saved as strings
            with commas as decimal points.
        str: as dataheader, but with new information appended.
            Header for the data file. 'ID\tT1\tT2\t...\tTm',
            where Ti's is the date and time of the measurement.
            Format: 'DD-MM-YYYY HH:MM:SS'.
    """

    (eltime, timelabels, temp, temp_comma, null1, dataheader,
     null2) = append_info(trafile,
                          dataheader, eltime, timelabels, initialtime)

    if np.shape(tempstack)[0] != np.shape(temp)[-1]:
        errormessage = "***** ERROR, files should have the same length "
        errormessage += "(and they do not). *****"
        sys.exit(errormessage)

    tempstack = np.column_stack((tempstack, np.array(temp).transpose()))
    tempstack_comma = np.column_stack((tempstack_comma,
                                       np.array(temp_comma).transpose()))

    return eltime, timelabels, tempstack, tempstack_comma, dataheader


def savefile(lines, tempstack, tempstack_comma, commas,
             langs, language, savedate, dataheader, idno):
    """
    Saves the final file with the results from all measurements.
    Additionally, it returns the filename used for saving figures.

    Args:
        lines (int): number of measurements per file.
        tempstack (array): numpy array with temperature measurements,
            size n x m, where n is number of measurement points,
            and m is number of measurement times.
            All measurements are saved as floats.
        tempstack_comma (array): numpy array with temperature measurements,
            size n x m, where n is number of measurement points,
            and m is number of measurement times.
            All measurements are saved as strings with commas as decimal point.
        commas (int): defines whether commas should be used as decimal point
            in the final file (1: yes, 0: no).
        langs (dict): dictionary with basic key words in different languages.
            Used for file name.
        language (str): defines which language should be used in filename
            ('en' for English and 'dk' for Danish).
        savedate (str): date of the first measurement, used in filename.
            Format: 'YYYYMMDD'.
        dataheader: header for the data file. 'ID\tT1\tT2\t...\tTm',
            where Ti's is the date and time of the measurement.
            Format: 'DD-MM-YYYY HH:MM:SS'.
        idno (array): numpy array with numbers from 0 to n, where n
            is variable lines.

    Returns:
        file: text file saved in the current directory, with n+1 columns
            and m+1 rows. n is number of measurement times, m is number
            of measurement points per measurement time.
            The first column is idno.
            The first row is dataheader.
        str: name used in figures as file name.
    """

    if commas == 0:
        results = np.column_stack((idno, tempstack))
    else:
        results = np.column_stack((idno, tempstack_comma))

    savename = "DTS_" + langs[language][2] + "_" + savedate + "_"
    savename += langs[language][0] + "_" + time.strftime("%Y%m%d")
    savenamefig = savename
    if commas == 1:
        savename += "_" + langs[language][1]

    # # if you want to change the directory, this should be added here,
    # # before savename
    # # example: mydata = np.savetxt("new_directory/" + savename + ....)
    saveddata = np.savetxt(savename + ".txt",
                           results, delimiter='\t', fmt="%s", comments='',
                           header=dataheader)

    return saveddata, savenamefig


def get_data_from_drive(txtfile, langs, language):
    """
    Extracts data from the text file with the results.
    The file has been previously created from the .tra-files.
    It is important, that ONLY ONE text file is saved in the same
    directory as the python scripts.

    Args:
        txtfile (str): data file, from which the data is to be extracted.
        langs (dict): dictionary with basic key words in different languages.
            Used for file name.
        language (str): defines which language should be used in filename
            ('en' for English and 'dk' for Danish).

    Returns:
        array: numpy array with numbers from 0 to n, where n+1
            is number of measurements (ID numbers).
        list of floats: list with the time elapsed for individual
            measurements since first measurement (in hours). Used for y-axis.
        array: numpy array with temperature measurements,
            size n x m, where n is number of measurement points,
            and m is number of measurement times.
            All measurements are saved as floats.
        str: names of individual columns in the final file
            with the results (data header).
            First column is "ID", the consequent columns are timestamps.
            Columns are separated by tabs:
            "ID\tDD-MM-YYYY HH:MM:SS\t...\tDD-MM-YYYY HH:MM:SS".
            Used to obtain timestamp for first and last measurement.
        list of strings: list with timestamps of individual
            measurements ("HH:MM"), used as y-label on the plot.
        str: string used in the figure name, includes basic
            information such as date and type.

    Raises:
       exception is raised if there are more than one text files in
       the same directory.
    """

    if len(txtfile) > 1:
        sys.exit("too many files, copy one file only to a separate folder")

    rows = 0

    idno = []
    eltime = []
    tempstack = []
    timelabels = []

    fh = open(txtfile[0])

    savedate = txtfile[0].split('_')[2]
    savenamefig = "DTS_" + langs[language][2] + "_" + savedate + "_"
    savenamefig += langs[language][0] + "_" + time.strftime("%Y%m%d")

    for line in fh:
        if rows == 0:
            dataheader = line
            rows += 1
            continue

        sp = line.split('\t')
        idno.append(int(sp[0].split('.')[0]))
        tempstack_onerow = []
        for temperature in sp[1:]:
            temperature = float(temperature.strip().replace(str2, str1))
            tempstack_onerow.append(temperature)

        tempstack.append(tempstack_onerow)
    fh.close()

    tempstack = np.array(tempstack)
    idno = np.array(idno)

    timestamps = dataheader.split('\t')[1:]
    item = 0
    for timestamp in timestamps:
        timestamp = timestamp.strip()
        dmy, hms = timestamp.split()
        d, m, y = dmy.split('-')
        s = '-'
        testtime = s.join([y, m, d]) + ' ' + hms

        if item == 0:
            eltime.append(0)
            initialtime = testtime
        else:
            eltime.append(elapsedtime(testtime, initialtime))

        timelabels.append(hms[:-3])
        item += 1

    return idno, eltime, tempstack, dataheader, timelabels, savenamefig


def plot_and_save(idno, eltime, tempstack,
                  dataheader, timelabels, savenamefig,
                  minidno, maxidno, timestep, stepcb, clevels,
                  changelimits,  udef_xticks,
                  langs, language):
    """
    Plots the results, saves the resulting figure on the drive.

    Args:
        idno (array): numpy array with numbers from 0 to n, where n+1
            is number of measurements.
        eltime (list of floats): list with the time elapsed for individual
            measurements since first measurement (in hours). Used for y-axis.
        tempstack (array): numpy array with temperature measurements,
            size n x m, where n is number of measurement points,
            and m is number of measurement times.
            All measurements are saved as floats.
        dataheader (str): names of individual columns in the final file
            with the results.
            First column is "ID", the consequent columns are timestamps.
            Columns are separated by tabs:
            "ID\tDD-MM-YYYY HH:MM:SS\t...\tDD-MM-YYYY HH:MM:SS".
            Used to obtain timestamp for first and last measurement.
        timelabels (list of strings): list with timestamps of individual
            measurements ("HH:MM"), used as y-label on the plot.
        savenamefig (str): string used in the figure name, includes basic
            information such as date and type.
        minidno (int): lower limit for measurement points on x-axis.
            All datapoints with ID number lower than minidno are excluded
            from the contour plot (the value itself is included).
        maxidno (int): upper limit for measurement points on x-axis.
            All datapoints with ID number higher than maxidno are excluded
            from the contour plot (the value itself is included).
            If maxidno is higher than the total number of measurement points,
            all points are automatically included.
        timestep (int): defines how many labels should be skipped on the graph.
            timestep = n means that every nth label is shown on the y-axis.
        stepcb (float): difference in Celsius degrees between consecutive
            y-ticks on the colorbar.
        clevels (float): number of contour levels.
        changelimits (list): defines whether user defined limits should be used
            for x-axis. If empty, limits on the x-axis are adjusted
            automatically. If not empty, the elements of the list are taken as
            xmin and xmax, respectively.
        udef_xticks (list): list in which all elements are floats. User defined
            x-ticks. If empty, the ticks are applied automatically.
        langs (dict): dictionary with basic key words in different languages.
            Used for file name.
        language (str): defines which language should be used in filename
            ('en' for English and 'dk' for Danish).

    Returns:
        figure: matplotlib.pyplot.subplot() with contour plots and colorbar.
        colorbar: colorbar used for the figure.
    """

    plt.close("all")

    if maxidno >= max(idno):
        xlist = np.linspace(minidno, max(idno), (max(idno)-minidno)+1)
    else:
        xlist = np.linspace(minidno, maxidno, (maxidno-minidno)+1)

    ylist = eltime

    X, Y = np.meshgrid(xlist, ylist)

    if maxidno >= max(idno):
        Z = tempstack.transpose()[:, minidno:]
    else:
        Z = tempstack.transpose()[:, minidno:maxidno+1]

    f, ax = plt.subplots()
    cp = ax.contourf(X, Y, Z, clevels, cmap=plt.cm.jet)

    # colorbar
    cbmax = np.ceil(tempstack.max())
    cbmin = np.floor(tempstack.min())
    # cbdiff = cbmax - (cbmin-1)
    cbtick = np.arange(cbmin, cbmax, stepcb)
    cb = plt.colorbar(cp, ax=ax,
                      label=langs[language][3] + " ($^\circ$C)", ticks=cbtick)

    ax.set_xlabel(langs[language][4])
    ax.set_ylabel(langs[language][5])

    starttime = dataheader.split("\t")[1][:-3]
    stoptime = dataheader.split("\t")[-1][:-3]
    ax.set_title("%s: %s, %s: %s" % (langs[language][6],
                                     starttime, langs[language][7], stoptime))

    ax.set_ylim(ylist[0], ylist[-1])
    # change x limits
    if changelimits == []:
        xmin = minidno
        xmax = xlist[-1]
    else:
        xmin = changelimits[0]
        xmax = changelimits[1]

    if udef_xticks == []:
        xticks = list(ax.get_xticks()) + [minidno]
    else:
        xticks = udef_xticks

    ax.set_xticks(xticks)

    ax.set_xlim(xmin, xmax)

    ylabels = []
    yticks = []

    for ii in range(int(len(ylist)/timestep)):
        ylabels.append(timelabels[ii*timestep][:5])
        yticks.append(ylist[ii*timestep])
    if np.mod(len(ylist), timestep) == 0:
        ylabels.append(timelabels[-1][:5])  # add the last one
        yticks.append(ylist[-1])
    else:
        ylabels.append(timelabels[(ii+1)*timestep][:5])  # add the last but one
        yticks.append(ylist[(ii+1)*timestep])

        # add the last one, if there is space
        if np.mod(len(ylist), timestep)-1 >= 4:
            ylabels.append(timelabels[-1][:5])  # add the last one
            yticks.append(ylist[-1])

    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)

    savenamefig += "_minidno_" + str(minidno)
    savenamefig += "_clevels_" + str(clevels) + "_x_" + str(xmin)
    savenamefig += "_" + str(xmax)
    savenamefig += "_stepcb_" + str(stepcb) + "_timestep_" + str(timestep)
    plt.savefig(savenamefig + ".png")

    return f, cb
