#!/usr/bin/python
#Imports
import numpy as np
import os
import pandas as pd
import decimal
import sys

class EWMeasurements:
    def __init__(self, input_params):
        self.input_df = input_params.copy()

    def round_up0(self, i):
        """Round up to 2nd decimal"""
        rounded = decimal.Decimal(str(i)).quantize(decimal.Decimal('1.11'), rounding=decimal.ROUND_HALF_UP)
        return float(rounded)


    def get_snr(self, fname='logARES.txt'):
        """Get the SNR from ARES

        Input
        -----
        fname : str
        Name of the logfile of ARES

        Output
        ------
        snr : int
        The SNR either provided by ARES or measured
        """

        with open(fname, 'r') as lines:
            for line in lines:
                if line.startswith('S/N'):
                    break
        line = line.strip('\n').split(':')
        snr = int(float(line[1]))
        return snr


    def mine_opt(self, spectra_name):
        '''Creates the mine.opt file'''
        fout = "specfits='./%s.fits\'\n" % spectra_name
        fout += "readlinedat= 'cdo.dat'\n"
        fout += "fileout='aresout.dat'\n"
        fout += "lambdai= 3000\n"
        fout += "lambdaf=8000\n"
        fout += "smoothder=14\n"
        fout += "space=3.0\n"
        fout += "rejt= 3;5764,5766,6047,6052,6068,6076\n"
        #fout += "rejt= 0.998258\n"
        fout += "lineresol=0.1\n"
        fout += "miniline=2\n"
        fout += "plots_flag=0\n"
        #fout += "rvmask='4,6767.78,6771.04,6772.32,6777.4,6771.04,20'\n"
        #fout += "rvmask='6,5352.5,5353.39,5361.62,5362.84,5364.88,5365.4,5362.84,20'\n"
        fout += "rvmask='3,6021.8,6024.06,6027.06,6024.06,20'\n"
        #fout += "rvmask='0,0.4'\n"
        with open('mine.opt', 'w') as f:
            f.writelines(fout)




    def make_linelist(self, starname):
        """Merging linelist with ares file
        """
        os.system("cut -f1 < linelist_damp.rdb > cdo.dat")
        os.system('ARES')

        #opening the linelist file
        linelist = pd.read_csv("linelist_damp.rdb", skiprows=2,
                            names=['WL', 'EP', 'loggf', 'element',  'num', "damp", 'EWsun'],
                            delimiter=r'\s+', converters={'WL': self.round_up0})
        linelist = linelist.sort_values(['WL'])

        #opening the ares output file
        aresout = pd.read_csv("aresout.dat", usecols=(0,2,3,4,5), names=['wave', 'depth', 'FWHM',  'EW', "EWerr"],
                            delimiter=r'\s+')
        aresout = aresout.sort_values(['wave'])

        #mearging the linelist and ares output files
        data_output = pd.merge(left=linelist, right=aresout, left_on='WL', right_on='wave', how='inner')
        #Remove lines that have depths smaller than 0 and larger than 1
        data_output = data_output[(data_output.depth >= 0) & (data_output.depth < 1.0)]
        data_output.to_csv('data_output.rdb', index=False, sep = '\t')


        #creating the table with EWs for MOOG input
        snr = self.get_snr()
        hdr = '%s - S/N: %s' % (starname, snr)
        ew_moog_input = data_output.loc[:, ['WL', 'num', 'EP', 'loggf', "damp", 'EW']]
        ew_moog_input = ew_moog_input.sort_values(['num'])
        fileout = "aresout_" + starname + ".dat"
        np.savetxt(fileout, ew_moog_input.values, fmt=('%9.3f', '%10.1f', '%9.2f', '%9.3f', '%9.2f', '%19.1f'), header=' %s' % hdr)


        #Selecting elements that have less than 4 lines
        element_list_tmp = data_output.groupby('num')['num'].count()
        element_list = pd.DataFrame([element_list_tmp.index, element_list_tmp.values], index = ['num', 'nlines']).transpose()
        elements_with_lt_4_lines = element_list[element_list.nlines < 4]

        lines_of_elements_with_lt_4_line = pd.merge(left=data_output,right=elements_with_lt_4_lines, left_on='num', right_on='num')

        #Calculating the error in EW with the Cayrel formulae
        lines_of_elements_with_lt_4_line['ew_plus_error'] = 0.0
        lines_of_elements_with_lt_4_line['ew_plus_error'] = lines_of_elements_with_lt_4_line['EW'] + 1000*2.45*((lines_of_elements_with_lt_4_line['FWHM']*0.035/2.3548)**0.5)/snr \
        + 6*1000*lines_of_elements_with_lt_4_line['FWHM']/snr/2.3548
        lines_of_elements_with_lt_4_line['ew_plus_error'] = np.round(lines_of_elements_with_lt_4_line['ew_plus_error'], 1)

        lines_of_elements_with_lt_4_line['ew_minus_error'] = 0.0
        lines_of_elements_with_lt_4_line['ew_minus_error'] = lines_of_elements_with_lt_4_line['EW'] - (1000*2.45*((lines_of_elements_with_lt_4_line['FWHM']*0.035/2.3548)**0.5)/snr \
        + 6*1000*lines_of_elements_with_lt_4_line['FWHM']/snr/2.3548)
        lines_of_elements_with_lt_4_line['ew_minus_error'] = np.round(lines_of_elements_with_lt_4_line['ew_minus_error'], 1)

        #creating the table with EWs (EW +/- error in EW) for MOOG input
        snr = self.get_snr()
        hdr = '%s - S/N: %s' % (starname, snr)
        ew_moog_input_plus_err_tmp = lines_of_elements_with_lt_4_line.loc[:, ['WL', 'num', 'EP', 'loggf', "damp", 'ew_plus_error']]
        ew_moog_input_plus_err_tmp = ew_moog_input_plus_err_tmp.rename(index=str, columns={"ew_plus_error": "ew_error"})
        ew_moog_input_minus_err_tmp = lines_of_elements_with_lt_4_line.loc[:, ['WL', 'num', 'EP', 'loggf', "damp", 'ew_minus_error']]
        ew_moog_input_minus_err_tmp = ew_moog_input_minus_err_tmp.rename(index=str, columns={"ew_minus_error": "ew_error"})
        ew_moog_input_plus_minus_err_tmp = pd.concat([ew_moog_input_plus_err_tmp, ew_moog_input_minus_err_tmp], ignore_index=True)
        #removing lines with EW - err_EW < = 0
        ew_moog_input_err = ew_moog_input_plus_minus_err_tmp[ew_moog_input_plus_minus_err_tmp['ew_error']>0]
        ew_moog_input = ew_moog_input_err.sort_values(['num'])
        fileout = "aresout_" + starname + "_ew_err.dat"
        np.savetxt(fileout, ew_moog_input.values, fmt=('%9.3f', '%10.1f', '%9.2f', '%9.3f', '%9.2f', '%19.1f'), header=' %s' % hdr)


        os.system("rm cdo.dat mine.opt aresout.dat")



    def ew_measurements(self):
        #The main program to measure EWs for a list of stars
        os.system("mkdir ares_log")
        os.system("mkdir ares_out")
        os.system('mkdir ares_ew')
        os.system('mkdir ares_ew/err_ew')

        starlist = self.input_df


        for i in range(len(starlist)):
            print (starlist.star[i])

            #Creating the Non-HFS linelist
            self.mine_opt("./%s" % starlist.star[i])
            self.make_linelist("%s" % starlist.star[i])

            os.system('mv  logARES.txt ./ares_log/%s.log' % starlist.star[i])
            os.system('mv  data_output.rdb ./ares_out/%s_aresout.dat' % starlist.star[i])
            os.system('mv  aresout_%s.dat ./ares_ew/%s_aresout.dat' % (starlist.star[i], starlist.star[i]))
            os.system('mv  aresout_%s_ew_err.dat ./ares_ew/err_ew/%s_aresout_ew_err.dat' % (starlist.star[i], starlist.star[i]))


        os.system('cp -r ares_ew ares_ew_corr')
        print ("Done")


if __name__ == '__main__':
    input_params = pd.read_table('input_param_error.rdb')
    ew_ares = EWMeasurements(input_params)
    ew_measurements = ew_ares.ew_measurements()

    
