#!/usr/bin/python
#Imports
import time
import numpy as np
import string
import os
import re
import math
import pandas as pd
import sys


class AtmosphereModeling():
    def __init__(self, input_params):
        self.input_df = input_params.copy()


    def atmo_kurucz(self, starname, teff, logg, feh, vtur, erteff, erlogg,erfeh,ervtur):
        """Creating models of atmospheres """

        os.system('echo %s %s %s | intermod_kurucz' % (teff, logg, feh))
        os.system('echo %s | transform_kurucz' % vtur)
        os.system('mv out.atm %s_orig.atm' %starname)

        os.system('echo %s %s %s | intermod_kurucz' % (teff, logg, feh + erfeh))
        os.system('echo %s | transform_kurucz' % vtur)
        os.system('mv out.atm %s_erfeh.atm' %starname)

        os.system('echo %s %s %s | intermod_kurucz' % (teff+erteff, logg, feh))
        os.system('echo %s | transform_kurucz' % vtur)
        os.system('mv out.atm %s_erteff.atm' %starname)

        os.system('echo %s %s %s | intermod_kurucz' % (teff, logg+erlogg, feh))
        os.system('echo %s | transform_kurucz' % vtur)
        os.system('mv out.atm %s_erlogg.atm' %starname)

        os.system('echo %s %s %s | intermod_kurucz' % (teff, logg, feh))
        os.system('echo %s | transform_kurucz' % (vtur + ervtur))
        os.system('mv out.atm %s_ervtur.atm' %starname)

        os.system('rm mod*')


    def model_atmospheres(self):
        #The main program to measure EWs for a list of stars
        os.system("mkdir atmospheres")

        starlist = self.input_df

        for i in range(len(starlist)):

            self.atmo_kurucz(starlist.star[i], starlist.teff[i], starlist.logg[i], starlist.feh[i], starlist.vtur[i],
            starlist.erteff[i], starlist.erlogg[i], starlist.erfeh[i], starlist.ervtur[i])

            os.system('mv %s_orig.atm atmospheres/%s_orig.atm' % (starlist.star[i], starlist.star[i]))
            os.system('mv %s_erfeh.atm atmospheres/%s_erfeh.atm' % (starlist.star[i], starlist.star[i]))
            os.system('mv %s_erteff.atm atmospheres/%s_erteff.atm' % (starlist.star[i], starlist.star[i]))
            os.system('mv %s_erlogg.atm atmospheres/%s_erlogg.atm' % (starlist.star[i], starlist.star[i]))
            os.system('mv %s_ervtur.atm atmospheres/%s_ervtur.atm' % (starlist.star[i], starlist.star[i]))


if __name__ == '__main__':
    input_params = pd.read_table('input_param_error.rdb')
    atmospheres = AtmosphereModeling(input_params)
    atmospheres.model_atmospheres()

