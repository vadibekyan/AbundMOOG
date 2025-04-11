#!/usr/bin/python
#Imports
import time
import numpy as np
import string
import os
import re
import math
import pandas as pd
from glob import glob
import sys

def mad(a, axis=None):
    """
    Compute *Median Absolute Deviation* of an array along given axis.
    """
    med = np.median(a, axis=axis)                # Median along given axis
    if axis is None:
        umed = med                              # med is a scalar
    else:
        umed = np.expand_dims(med, axis)         # Bring back the vatished axis
    mad = np.median(np.absolute(a - umed), axis=axis)  # MAD along given axis

    return mad

def getRoundedThresholdv1(a, MinClip):
    return round(float(a) / MinClip) * MinClip

def weighted_avg_and_std(values):
    """Get the weighted average and standard deviation.

    Input
    -----
    values : ndarray
      The list of values from which to calculate the weighted
      average and standard deviation

    Output
    ------
    average : float
      The weighted average
    std : float
      The weighted standard deviation
    """
    values = np.array(values)
    weights =  (np.abs(values-np.median(values))/(np.std(values)+1E-13)+0.25)
    weights_rounded = (np.round(weights/0.5+1E-13)*0.5)**(-1)
    average = round(np.average(values, weights=weights_rounded), 3)
    std = np.sqrt(np.average((values-average)**2, weights=weights_rounded))
    return average, std

def round_float(s):
    '''1. if s is float, round it to 3 decimals
       2. else return s as is
    '''
    import re
    m = re.match("(\d+\.\d+)",s.__str__())
    try:
        r = round(float(m.groups(0)[0]),3)
    except:
        r = s
    return r


def abfind_par(damping):
    '''Creates the abfind.par file.
    Damping = 0 uses Unsold aproximation, but if a factor is read from the line-list then,
        if factor greater than 10^-10, then multiply by the factor
        otherwise replace the Unsold value by the factor
    Damping = 1 uses the Barklem data, if no Barklem data, then goes to option 'damping = 0'
    Damping = 2 uses Unsold aproximation multiplied by a factor recommended by Blackwell group
    '''
    fout = "abfind\n"
    fout += "terminal    X11\n"
    fout += "atmosphere  1\n"
    fout += "molecules   1\n"
    fout += "lines       1\n"
    fout += "flux/int    0\n"
    fout += "plot        0\n"
    fout += "damping     %s\n" % damping
    fout += "units       0\n"
    fout += "standard_out  'standard_out'\n"
    fout += "summary_out   'outab.dat'\n"
    fout += "model_in      'out.atm'\n"
    fout += "lines_in      'lines.dat'\n"
    with open('abfind.par', 'w') as f:
        f.writelines(fout)

def blends_par(damping, elem):
    '''Creates the abfind.par file
    Damping = 0 uses Unsold aproximation, but if a factor is read from the line-list then,
        if factor greater than 10^-10, then multiply by the factor
        otherwise replace the Unsold value by the factor
    Damping = 1 uses the Barklem data, if no Barklem data, then goes to option 'damping = 0'
    Damping = 2 uses Unsold aproximation multiplied by a factor recommended by Blackwell group
    '''
    fout = "blends\n"
    fout += "terminal    X11\n"
    fout += "atmosphere  1\n"
    fout += "molecules   2\n"
    fout += "lines       1\n"
    fout += "freeform    0\n"
    fout += "flux/int    0\n"
    fout += "plot        0\n"
    fout += "damping     %s\n" % damping
    fout += "units       0\n"
    fout += "standard_out  'standard_out'\n"
    fout += "summary_out   'outab.dat'\n"
    fout += "model_in      'out.atm'\n"
    fout += "lines_in      'lines.dat'\n"

    if elem == 'ScI':
        fout += "blenlimits\n"
        fout += "   0.8  0.01   21\n"

    elif elem == 'ScII':
        fout += "blenlimits\n"
        fout += "   0.8  0.01   21\n"

    elif elem == 'VI':
        fout += "blenlimits\n"
        fout += "   0.8  0.01   23\n"

    elif elem == 'MnI':
        fout += "blenlimits\n"
        fout += "   0.8  0.01   25\n"

    elif elem == 'CoI':
        fout += "blenlimits\n"
        fout += "   0.8  0.01   27\n"

    elif elem == 'CuI':
        fout += "isotopes  2  1\n"
        fout += " 29.0630     1.4457\n"
        fout += " 29.0650     3.2436\n"
        fout += "blenlimits\n"
        fout += "   0.8  0.01   29\n"

    elif elem == 'BaII':
        fout += "isotopes  7  1\n"
        fout += " 56.1130   943.3962\n"
        fout += " 56.1132   990.0990\n"
        fout += " 56.1134    41.3736\n"
        fout += " 56.1135    15.1699\n"
        fout += " 56.1136    12.7324\n"
        fout += " 56.1137     8.9047\n"
        fout += " 56.1138     1.3947\n"
        fout += "blenlimits\n"
        fout += "   0.8  0.01   56\n"

    with open('blends.par', 'w') as f:
        f.writelines(fout)

def abfind(starname, param):
    """
    Derives abundances for a given star ('starname') and given 'parameter' (param).
    param = orig (the main abundances); teff (varying the teff) etc.
    """

    print (starname, param)
    os.system('mkdir moog_abundances')
    os.system('mkdir moog_abundances/moog_outab')
    os.system('mkdir moog_abundances/moog_std_out')
    os.system('mkdir moog_abundances/moog_outab_formatted')
    os.system('mkdir moog_abundances/moog_abstar')

    #create the abfind.par file
    abfind_par(1)

    #run MOOGSILENT2019
    os.system('cp atmospheres/%s_%s.atm out.atm' % (starname,param))
    os.system('cp ./ares_ew_corr/%s_aresout.dat lines.dat'  % starname)
    os.system('echo abfind.par | MOOGSILENT2019')


    #formatting the moog output
    ALPHA = string.ascii_letters
    with open('outab.dat', 'r') as fin, open('outab_formated.tmp', 'w') as fout:
        for line in fin:
            if line.lstrip():
                line = line.lstrip()
                if not line.startswith(tuple(ALPHA)):
                    new_line = re.sub(' +','\t', line)
                    fout.write(new_line)

    #open the formated moog output abudnances
    abund_star = pd.read_table("outab_formated.tmp", skiprows=1, usecols=(0,1,4,6,7), names=['wave', 'num', 'EWin',  'abund', "delavg"],
                          delimiter=r'\s+')
    abund_star.to_csv('./moog_abundances/moog_outab_formatted/%s_%s_outab_formated.dat' % (starname, param), sep = '\t', index = False)

    #List of elements that have lines in the moog output
    element_list_tmp = abund_star.groupby('num')['num'].count()
    # element_list has two columns: num - element's atomic number and nlines - number of lines
    element_list = pd.DataFrame([element_list_tmp.index, element_list_tmp.values], index = ['num', 'nlines']).transpose()


    element_info = pd.read_table("elements.rdb")
    element_abund = pd.merge(left=element_list,right=element_info, left_on='num', right_on='num')
    element_abund['abund'] = 0.0
    element_abund['abund_std_err'] = 0.0


    #calculate abundances for each element
    for i, num in enumerate(element_abund['num']): # beforer, instead of element_abund it was "element_list"
        print (num)

        line_abundances = abund_star[abund_star['num'] == num]

        if element_abund['nlines'][i] == 1:
            el_weighted_mean = line_abundances['abund']
            el_weighted_std = "NaN"
            element_abund['nlines'][i] = 1.0

        elif element_abund['nlines'][i] == 2:
            el_weighted_mean = np.mean(line_abundances['abund'])
            el_weighted_std = (np.max(line_abundances['abund']) - np.min(line_abundances['abund']))/2
            element_abund['nlines'][i] = 2.0

        elif element_abund['nlines'][i] > 2:
            #removing very outliers i.e. at 10sigma level
            removing_outliers = np.abs(line_abundances['abund'] - np.median(line_abundances['abund']))/(np.std(line_abundances['abund'])+1E-13) < 3
            line_best_abundances = line_abundances[removing_outliers]

            #calculating the weighted mean and its error
            el_weighted_mean, el_weighted_std = weighted_avg_and_std(np.array(line_best_abundances['abund']))
            element_abund['nlines'][i] = len(line_best_abundances['abund'])

        el_weighted_mean

        element_abund['abund'][i] =  round_float(el_weighted_mean) #np.round(el_weighted_mean, 3)
        element_abund['abund_std_err'][i] = round_float(el_weighted_std)#np.round(el_weighted_std,3)

    element_abund.to_csv('./moog_abundances/moog_abstar/%s_%s_abund.dat' % (starname, param), sep = '\t', index = False)

    os.system('mv outab.dat moog_abundances/moog_outab/%s_%s_outab.dat' % (starname, param))
    os.system('mv standard_out moog_abundances/moog_std_out/%s_%s_std_out.dat' % (starname, param))
    os.system('rm out.atm lines.dat fort* abfind.par blends.par outab_formated.tmp')

def abfind_elem_lt4lines(starname, param):
    """
    Derives abundances for a given star ('starname') and given 'parameter' (param).
    param = orig (the main abundances); teff (varying the teff) etc.
    """

    #create the abfind.par file
    abfind_par(1)

    #run MOOGSILENT2019

    os.system('cp atmospheres/%s_%s.atm out.atm' % (starname,param))
    os.system('cp ./ares_ew_corr/err_ew/%s_aresout_ew_err.dat lines.dat'  % starname)
    os.system('echo abfind.par | MOOGSILENT2019')

    #formatting the moog output
    ALPHA = string.ascii_letters
    with open('outab.dat', 'r') as fin, open('outab_formated.tmp', 'w') as fout:
        for line in fin:
            if line.lstrip():
                line = line.lstrip()
                if not line.startswith(tuple(ALPHA)):
                    new_line = re.sub(' +','\t', line)
                    fout.write(new_line)

    #open the formated moog output abudnances
    abund_star = pd.read_table("outab_formated.tmp", skiprows=1, usecols=(0,1,4,6,7), names=['wave', 'num', 'EWin',  'abund', "delavg"],
                          delimiter=r'\s+')
    abund_star.to_csv('./moog_abundances/moog_outab_formatted/%s_%s_outab_formated_ew_err.dat' % (starname, param), sep = '\t', index = False)

    #List of elements that have lines in the moog output
    element_list_tmp = abund_star.groupby('num')['num'].count()
    # element_list has two columns: num - element's atomic number and nlines - number of lines
    element_list = pd.DataFrame([element_list_tmp.index, element_list_tmp.values], index = ['num', 'nlines']).transpose()

    element_info = pd.read_table("elements.rdb")
    element_abund = pd.merge(left=element_list,right=element_info, left_on='num', right_on='num')
    element_abund['abund'] = 0.0
    element_abund['abund_std_err_ew'] = 0.0

    #calculate abundances for each element
    for i, num in enumerate(element_abund['num']):# beforer, instead of element_abund it was "element_list"

        line_abundances = abund_star[abund_star['num'] == num]

        if element_abund['nlines'][i] == 1:
            el_weighted_mean = line_abundances['abund']
            el_weighted_std = "NaN"
            element_abund['nlines'][i] = 1.0

        elif element_abund['nlines'][i] == 2:
            el_weighted_mean = np.mean(line_abundances['abund'])
            el_weighted_std = (np.max(line_abundances['abund']) - np.min(line_abundances['abund']))/2
            element_abund['nlines'][i] = 2.0

        elif element_abund['nlines'][i] > 2:
            #removing very outliers i.e. at 10sigma level
            removing_outliers = np.abs(line_abundances['abund'] - np.median(line_abundances['abund']))/(np.std(line_abundances['abund'])+1E-13) < 3
            line_best_abundances = line_abundances[removing_outliers]

            #calculating the weighted mean and its error
            el_weighted_mean, el_weighted_std = weighted_avg_and_std(np.array(line_best_abundances['abund']))
            element_abund['nlines'][i] = len(line_best_abundances['abund'])

        element_abund['abund'][i] = round_float(el_weighted_mean) #np.round(el_weighted_mean, 3)
        element_abund['abund_std_err_ew'][i] = round_float(el_weighted_std) #np.round(el_weighted_std,3)

    element_abund_selected_columns = element_abund[['num', 'abund_std_err_ew']]
    element_abund_selected_columns.to_csv('./moog_abundances/moog_abstar/%s_%s_abund_ew_err.dat' % (starname, param), sep = '\t', index = False)

    os.system('mv outab.dat moog_abundances/moog_outab/%s_%s_outab_ew_err.dat' % (starname, param))
    os.system('mv standard_out moog_abundances/moog_std_out/%s_%s_std_out_ew_err.dat' % (starname, param))
    os.system('rm out.atm lines.dat fort* abfind.par blends.par outab_formated.tmp')


def abfind_err_param(starname):

    """
    Deriving errors of the abundances
    """

    #Opening the table with abundances
    element_abund_error = pd.read_csv("./moog_abundances/moog_abstar/%s_orig_abund.dat" % starname,
    usecols=(2,8), skiprows=1, names=['element', 'abund_orig'], delimiter='\t')
    #Adding new columns for the abundance differences due to variation of parameters


    element_abund_error['abund_diff_erlogg'] = 0.0
    element_abund_error['abund_diff_erteff'] = 0.0
    element_abund_error['abund_diff_erfeh'] = 0.0
    element_abund_error['abund_diff_ervtur'] = 0.0

    parameters = ['erlogg', 'erteff', 'erfeh', 'ervtur']
    for i, parameter in enumerate(parameters):

        #Running abfind function to derive abundances for each parameter variation
        abfind(starname, parameter)

        #Opening the tbale with abundances and slecting only abundanes and elemnt name columns
        slected_columns = pd.read_csv("./moog_abundances/moog_abstar/%s_%s_abund.dat" % (starname, parameter),
        usecols=(2,8), skiprows=1, names=['element', 'abund_%s' % parameter], delimiter='\t')



        element_abund_error_merged = element_abund_error.merge(slected_columns, on='element')
        #Adding a new column for abundances difference due to parameter variation
        element_abund_error_merged['abund_diff_%s' % parameter] = 0.0
        element_abund_error_merged['abund_diff_%s' % parameter] = np.abs(element_abund_error_merged['abund_orig']-element_abund_error_merged['abund_%s' %parameter ])
        element_abund_error['abund_diff_%s' % parameter] = element_abund_error_merged['abund_diff_%s' % parameter]


        slected_columns.to_csv('./moog_abundances/moog_abstar/%s_%s_abund.dat' % (starname, parameter), sep = '\t', index = False)

    #Calculating total error due to parameter variation
    element_abund_error['abund_param_err'] = 0.0
    element_abund_error['abund_param_err'] = (element_abund_error['abund_diff_erlogg']**2 + element_abund_error['abund_diff_erteff']**2 +
    element_abund_error['abund_diff_erfeh']**2 + element_abund_error['abund_diff_ervtur']**2)**0.5
    element_abund_param_error = element_abund_error[['element', 'abund_param_err']]

    element_abund = pd.read_table("./moog_abundances/moog_abstar/%s_orig_abund.dat" % starname, delimiter='\t')
    element_abund_all = element_abund.merge(element_abund_param_error, on='element')

    element_abund_all.to_csv('./moog_abundances/moog_abstar/%s_abund_param_err.dat' % starname, sep = '\t', index = False)
    #Deleting unnesseary files
    os.system('rm ./moog_abundances/moog_abstar/*erlogg* ./moog_abundances/moog_abstar/*erteff* ./moog_abundances/moog_abstar/*erfeh* ./moog_abundances/moog_abstar/*ervtur*')
    os.system('rm ./moog_abundances/moog_outab_formatted/*erlogg* ./moog_abundances/moog_outab_formatted/*erteff* ./moog_abundances/moog_outab_formatted/*erfeh* ./moog_abundances/moog_outab_formatted/*ervtur*')
    os.system('rm ./moog_abundances/moog_outab/*erlogg* ./moog_abundances/moog_outab/*erteff* ./moog_abundances/moog_outab/*erfeh* ./moog_abundances/moog_outab/*ervtur*')
    os.system('rm ./moog_abundances/moog_std_out/*erlogg* ./moog_abundances/moog_std_out/*erteff* ./moog_abundances/moog_std_out/*erfeh* ./moog_abundances/moog_std_out/*ervtur*')


def abund_err(starname):
    abfind_elem_lt4lines('%s' %  starname, 'orig')

    element_abund = pd.read_table("./moog_abundances/moog_abstar/%s_abund_param_err.dat" % starname, delimiter='\t')
    element_abund_ew_err = pd.read_table("./moog_abundances/moog_abstar/%s_orig_abund_ew_err.dat" % starname, delimiter='\t')
    element_abund_all = element_abund.merge(element_abund_ew_err, on='num', how='outer')
    element_abund_all['total_err'] = 0


    element_abund_all.loc[element_abund_all.abund_std_err_ew > 0, 'total_err'] =  (element_abund_all['abund_std_err_ew']**2 + element_abund_all['abund_param_err']**2)**0.5
    element_abund_all.loc[element_abund_all.abund_std_err_ew != element_abund_all.abund_std_err_ew, 'total_err'] =  (element_abund_all['abund_std_err']**2 + element_abund_all['abund_param_err']**2)**0.5

    element_abund_all['abund_param_err'] =  np.round(element_abund_all['abund_param_err'], 3)
    element_abund_all['total_err'] = np.round(element_abund_all['total_err'], 3)
    element_abund_all.to_csv('./moog_abundances/moog_abstar/%s_abund_err.dat' % starname, sep = '\t', index = False)



def abundances(input_starlist):
    """
    Main program to derive abundances and errors for a list of stars
    """
    starlist = pd.read_csv(input_starlist, sep='\s+', skiprows=[1])

    for i in range(len(starlist)):
        print (starlist.star[i])

        abfind('%s' %  starlist.star[i],  'orig')
        abfind_err_param('%s' %  starlist.star[i])
        abund_err('%s' %  starlist.star[i])
        #abfind_hfs('%s' %  starlist.star[i],  'orig')
        #abfind_hfs_err_param('%s' %  starlist.star[i])
        #abund_hfs_err('%s' %  starlist.star[i])
        os.system('rm ./moog_abundances/moog_abstar/*abund_param_err.dat ./moog_abundances/moog_abstar/*orig_abund.dat ./moog_abundances/moog_abstar/*orig_abund_ew_err.dat')
        os.system('rm ./moog_abundances/abundances_hfs/moog_abstar/*abund_param_err.dat ./moog_abundances/abundances_hfs/moog_abstar/*orig_abund.dat ./moog_abundances/abundances_hfs/moog_abstar/*orig_abund_ew_err.dat')

abundances('input_param_error.rdb')


"""
if __name__ == '__main__':
    abundances(str(sys.argv[1]))
"""
