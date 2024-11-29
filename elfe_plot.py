#!/usr/bin/python

import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import math
import mpl_toolkits.axisartist as axisartist
import seaborn as sns
import pandas as pd
import statsmodels.formula.api as sm
from matplotlib.ticker import MultipleLocator


harps_sample = pd.read_table('abundances_harps1111_o_c.rdb', sep = '\t')
harps_solar_analogs = harps_sample[(harps_sample['teff'] > 5277) & (harps_sample['teff'] < 6277) & (harps_sample['feh'] < 1.4)]

target_name = ''






def N_columns_rows_for_plot(N_elements):
    n_rows = 0
    n_columns = 0

    if N_elements == 4:
        n_rows = 2
        n_columns = 2
    elif N_elements >=5 and N_elements <= 6:
        n_rows = 3
        n_columns = 2
    elif N_elements >= 7 and N_elements <= 9:
        n_rows = 3
        n_columns = 3
    elif N_elements >= 10 and N_elements <= 12:
        n_rows = 4
        n_columns = 3
    elif N_elements >= 13 and N_elements <= 15:
        n_rows = 5
        n_columns = 3
    elif N_elements == 16:
        n_rows = 4
        n_columns = 4
    elif N_elements >= 17 and N_elements <= 20:
        n_rows = 5
        n_columns = 4
    elif N_elements >= 21 and N_elements <= 24:
        n_rows = 6
        n_columns = 4
    elif N_elements >= 22 and N_elements <= 28:
        n_rows = 7
        n_columns = 4

    return n_rows, n_columns

def tick_parameters(n_rows, n_columns):
    top_left = ()
    for i in range(n_rows-1):
        top_left += +i*n_columns,
    top_left += 100,

    top_right = ()
    for j in range(n_columns-1):
        for i in range(n_rows-1):
            top_right += (1+i*n_columns)+(j),
    top_right += 100,

    bottom_left = ()
    bottom_left += (n_rows-1)*n_columns,
    bottom_left += 100,

    bottom_right = ()
    for j in range(n_columns-1):
        bottom_right += (n_rows-1)*n_columns + j+1,
    bottom_right += 100,

    #print bottom_left, top_left, top_right, bottom_right
    return bottom_left, top_left, top_right, bottom_right

def elfe_feh():
    abundances_init = pd.read_table("all_abundances_ordered_star.rdb")
    elements_df = pd.read_table("elements_to_plot.rdb")

    n_rows, n_columns = N_columns_rows_for_plot(len(elements_df.element))

    print (abundances_init)

    extreame_all_elfe = []
    for i, element in enumerate(elements_df.element):
        harps_solar_analogs_best = harps_solar_analogs[(harps_solar_analogs['%s' % element] < 1) & (harps_solar_analogs['%s' % element] > -2)]
        abundances_df = abundances_init[(abundances_init['%s_rel_sun' % element] > -2)]
        elfe = abundances_df['%s_rel_sun' % element] - abundances_df['feh'] + 0.02
        max_elfe = np.max(elfe)
        min_elfe = np.min(elfe)
        extreame_all_elfe.append(min_elfe)
        extreame_all_elfe.append(max_elfe)
    min_all_elfe = math.ceil(np.min(extreame_all_elfe)*10)/10-0.1
    max_all_elfe = math.ceil(np.max(extreame_all_elfe)*10)/10


    for i, element in enumerate(elements_df.element):
        print (element)
        harps_solar_analogs_best = harps_solar_analogs[(harps_solar_analogs['%s' % element] < 1) & (harps_solar_analogs['%s' % element] > -2)]
        abundances_df = abundances_init[(abundances_init['%s_rel_sun' % element] > -2)]
        elfe = abundances_df['%s_rel_sun' % element] - abundances_df['feh'] + 0.02
        elfe_err = (abundances_df['%s_err' % element]**2 + abundances_df['erfeh']**2)**0.5

        giants = abundances_df[abundances_df.teff < 4300]
        elfe_giants = giants['%s_rel_sun' % element] - giants['feh'] + 0.02
        elfe_err_giants = (giants['%s_err' % element]**2 + giants['erfeh']**2)**0.5

        dwarfs = abundances_df[abundances_df.logg >= 0.5]
        elfe_dwarfs = dwarfs['%s_rel_sun' % element] - dwarfs['feh'] + 0.02
        elfe_err_dwarfs = (dwarfs['%s_err' % element]**2 + dwarfs['erfeh']**2)**0.5

        target_abund = abundances_df[abundances_df.star == target_name]
        elfe_target_abund = target_abund['%s_rel_sun' % element] - target_abund['feh'] + 0.02
        elfe_err_target_abund = (target_abund['%s_err' % element]**2 + target_abund['erfeh']**2)**0.5



        plt.figure(1)
        sns.set_style("white")
        sns.set_style("ticks")
        sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
        sns.set_context("paper", font_scale=1.5)
        fig = plt.gcf()
        #Fig size depends on the number of columns and rows
        fig.set_size_inches(3*n_columns, 3*n_rows)
        fig.subplots_adjust(left=0.13, right=0.98, bottom = 0.09, top = 0.98, wspace=0.0, hspace=0.0)
        ax = plt.subplot(n_rows, n_columns, i+1)
        ax.text(.9, .85,'%s' % element, horizontalalignment='right', transform=ax.transAxes, fontsize = 14)
        plt.scatter(harps_solar_analogs_best.feh, (harps_solar_analogs_best['%s' % element]-harps_solar_analogs_best.feh), c = 'grey', s = 20, alpha = 0.5, edgecolors='none')
        plt.scatter(dwarfs['feh'], elfe_dwarfs, linewidth=1, color='black', alpha = 0.6, s = 70)
        plt.errorbar(dwarfs['feh'], elfe_dwarfs, yerr= elfe_err_dwarfs, xerr= dwarfs['erfeh'], color='black', fmt='o', markersize=1, linewidth=0.5, alpha = 0.6)
        plt.scatter(giants['feh'], elfe_giants, linewidth=1, color='red', alpha = 0.6, s = 90)
        #plt.errorbar(giants['feh'], elfe_giants, yerr= elfe_err_giants, xerr= giants['erfeh'], color='red', fmt='o', markersize=1, linewidth=0.5, alpha = 0.6)
        plt.scatter(target_abund['feh'], elfe_target_abund, linewidth=1, color='red', alpha = 0.6, s = 150)
        plt.errorbar(target_abund['feh'], elfe_target_abund, yerr= elfe_err_target_abund, xerr= target_abund['erfeh'], linewidth=2, color='red', alpha = 0.99, ms = 15)

        plt.axhline(y=0, color='blue', linestyle='--', linewidth=1.)
        plt.axvline(x=0, color='blue', linestyle='--', linewidth=1.)
        #plt.axhline(y=0.12, color='blue', linestyle='--', linewidth=1.)
        x_min = math.ceil(np.min(abundances_df['feh'])*10)/10-0.2
        x_max = math.ceil(np.max(abundances_df['feh'])*10)/10 + 0.1

        #print dwarfs[(elfe_dwarfs > 0.12) & (dwarfs['feh'] > -0.2)]
        print (dwarfs[(elfe_dwarfs > 0.6)])


        plt.ylim(-0.3, 0.6)
        plt.xlim(-0.9, 0.5)
        #plt.xticks(np.arange(-0.1,0.14, 0.05))

        #plt.yticks(np.arange(-0.2,0.4, 0.1))
        #plt.minorticks_on()
        plt.tick_params(labelsize=13)

        bottom_left, top_left, top_right, bottom_right = tick_parameters(n_rows, n_columns)

        if i in top_right:
            plt.tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            right='on',      # ticks along the bottom edge are off
            left='on',         # ticks along the top edge are off
            labelleft='off')
            plt.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            top='on',      # ticks along the bottom edge are off
            bottom='on',         # ticks along the top edge are off
            labelbottom='off', labelsize=13)


        if i in top_left:
            plt.tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            right='on',      # ticks along the bottom edge are off
            left='on',         # ticks along the top edge are off
            labelleft='on',labelsize=13)
            plt.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            top='on',      # ticks along the bottom edge are off
            bottom='on',         # ticks along the top edge are off
            labelbottom='off', labelsize=14)

        if i in bottom_right:
            plt.tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            right='on',      # ticks along the bottom edge are off
            left='on',         # ticks along the top edge are off
            labelleft='off',labelsize=14)
            plt.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            top='on',      # ticks along the bottom edge are off
            bottom='on',         # ticks along the top edge are off
            labelbottom='on',labelsize=12.5)

        if i in bottom_left:
            plt.tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            right='on',      # ticks along the bottom edge are off
            left='on',         # ticks along the top edge are off
            labelleft='on',labelsize=13)
            plt.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            top='on',      # ticks along the bottom edge are off
            bottom='on',         # ticks along the top edge are off
            labelbottom='on',labelsize=12.5)

        fig.text(0.55, 0.01, '[Fe/H]', ha='center', fontsize = 17)
        fig.text(0.01, 0.5, '[X/Fe]', va='center', rotation='vertical', fontsize = 17)
        plt.savefig('elfe_feh.pdf', format='pdf')
    plt.clf()

def elfe_teff():
     abundances_init = pd.read_table("all_abundances_ordered_star.rdb")
     elements_df = pd.read_table("elements_to_plot.rdb")

     n_rows, n_columns = N_columns_rows_for_plot(len(elements_df.element))

     extreame_all_elfe = []
     for i, element in enumerate(elements_df.element):
         harps_solar_analogs_best = harps_solar_analogs[(harps_solar_analogs['%s' % element] < 1) & (harps_solar_analogs['%s' % element] > -2)]
         abundances_df = abundances_init[(abundances_init['%s_rel_sun' % element] > -2)]
         elfe = abundances_df['%s_rel_sun' % element] - abundances_df['feh'] + 0.02
         max_elfe = np.max(elfe)
         min_elfe = np.min(elfe)
         extreame_all_elfe.append(min_elfe)
         extreame_all_elfe.append(max_elfe)
     min_all_elfe = math.ceil(np.min(extreame_all_elfe)*10)/10-0.1
     max_all_elfe = math.ceil(np.max(extreame_all_elfe)*10)/10


     for i, element in enumerate(elements_df.element):
         print (element)
         harps_best = harps_sample[(harps_sample['%s' % element] < 1) & (harps_sample['%s' % element] > -2)]
         abundances_df = abundances_init[(abundances_init['%s_rel_sun' % element] > -2)]
         elfe = abundances_df['%s_rel_sun' % element] - abundances_df['feh'] + 0.02
         elfe_err = (abundances_df['%s_err' % element]**2 + abundances_df['erfeh']**2)**0.5


         giants = abundances_df[abundances_df.logg < 3.5]
         elfe_giants = giants['%s_rel_sun' % element] - giants['feh'] + 0.02
         elfe_err_giants = (giants['%s_err' % element]**2 + giants['erfeh']**2)**0.5

         dwarfs = abundances_df[abundances_df.logg >= 3.5]
         elfe_dwarfs = dwarfs['%s_rel_sun' % element] - dwarfs['feh'] + 0.02
         elfe_err_dwarfs = (dwarfs['%s_err' % element]**2 + dwarfs['erfeh']**2)**0.5



         plt.figure(1)
         sns.set_style("white")
         sns.set_style("ticks")
         sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
         fig = plt.gcf()
         #Fig size depends on the number of columns and rows
         fig.set_size_inches(3*n_columns, 3*n_rows)
         fig.subplots_adjust(left=0.1, right=0.98, bottom = 0.07, top = 0.98, wspace=0.0, hspace=0.0)
         ax = plt.subplot(n_rows, n_columns, i+1)
         ax.text(.9, .85,'%s' % element, horizontalalignment='right', transform=ax.transAxes, fontsize = 14)
         plt.scatter(harps_best.teff, (harps_best['%s' % element]-harps_best.feh), c = 'grey', s = 20, alpha = 0.5, edgecolors='none')
         plt.scatter(dwarfs['teff'], elfe_dwarfs, linewidth=1, color='black', alpha = 0.6, s = 70)
         plt.errorbar(dwarfs['teff'], elfe_dwarfs, yerr= elfe_err_dwarfs, xerr= dwarfs['erfeh'], color='black', fmt='o', markersize=1, linewidth=0.5, alpha = 0.6)
         plt.scatter(giants['teff'], elfe_giants, linewidth=1, color='red', alpha = 0.6, s = 90)
         plt.errorbar(giants['teff'], elfe_giants, yerr= elfe_err_giants, xerr= giants['erfeh'], color='red', fmt='o', markersize=1, linewidth=0.5, alpha = 0.6)
         plt.axhline(y=0, color='blue', linestyle='--', linewidth=1.)
         plt.axvline(x=5777, color='blue', linestyle='--', linewidth=1.)
         x_min = math.ceil(np.min(abundances_df['teff'])/100)*100-100
         x_max = math.ceil(np.max(abundances_df['teff'])/100)*100+100

         plt.ylim(-0.5, 0.8)
         plt.xlim(4400, 6800)
         plt.xticks(np.arange(4500,6600, 500))
         #plt.yticks(np.arange(-0.2,0.4, 0.1))
         #plt.minorticks_on()
         plt.tick_params(labelsize=13)


         bottom_left, top_left, top_right, bottom_right = tick_parameters(n_rows, n_columns)

         if i in top_right:
             plt.tick_params(
             axis='y',          # changes apply to the x-axis
             which='both',      # both major and minor ticks are affected
             right='on',      # ticks along the bottom edge are off
             left='on',         # ticks along the top edge are off
             labelleft='off')
             plt.tick_params(
             axis='x',          # changes apply to the x-axis
             which='both',      # both major and minor ticks are affected
             top='on',      # ticks along the bottom edge are off
             bottom='on',         # ticks along the top edge are off
             labelbottom='off', labelsize=16)


         if i in top_left:
             plt.tick_params(
             axis='y',          # changes apply to the x-axis
             which='both',      # both major and minor ticks are affected
             right='on',      # ticks along the bottom edge are off
             left='on',         # ticks along the top edge are off
             labelleft='on',labelsize=13)
             plt.tick_params(
             axis='x',          # changes apply to the x-axis
             which='both',      # both major and minor ticks are affected
             top='on',      # ticks along the bottom edge are off
             bottom='on',         # ticks along the top edge are off
             labelbottom='off', labelsize=16)

         if i in bottom_right:
             plt.tick_params(
             axis='y',          # changes apply to the x-axis
             which='both',      # both major and minor ticks are affected
             right='on',      # ticks along the bottom edge are off
             left='on',         # ticks along the top edge are off
             labelleft='off',labelsize=16)
             plt.tick_params(
             axis='x',          # changes apply to the x-axis
             which='both',      # both major and minor ticks are affected
             top='on',      # ticks along the bottom edge are off
             bottom='on',         # ticks along the top edge are off
             labelbottom='on',labelsize=12.5)

         if i in bottom_left:
             plt.tick_params(
             axis='y',          # changes apply to the x-axis
             which='both',      # both major and minor ticks are affected
             right='on',      # ticks along the bottom edge are off
             left='on',         # ticks along the top edge are off
             labelleft='on',labelsize=13)
             plt.tick_params(
             axis='x',          # changes apply to the x-axis
             which='both',      # both major and minor ticks are affected
             top='on',      # ticks along the bottom edge are off
             bottom='on',         # ticks along the top edge are off
             labelbottom='on',labelsize=12.5)

         fig.text(0.55, 0.01, '$T_{eff}$', ha='center', fontsize = 17)
         fig.text(0.02, 0.5, '[X/Fe]', va='center', rotation='vertical', fontsize = 17)
         plt.savefig('elfe_teff.pdf', format='pdf')
     plt.clf()

def elfe_logg():
     abundances_init = pd.read_table("all_abundances_ordered_star.rdb")
     elements_df = pd.read_table("elements_to_plot.rdb")

     n_rows, n_columns = N_columns_rows_for_plot(len(elements_df.element))

     extreame_all_elfe = []
     for i, element in enumerate(elements_df.element):
         harps_solar_analogs_best = harps_solar_analogs[(harps_solar_analogs['%s' % element] < 1) & (harps_solar_analogs['%s' % element] > -2)]
         abundances_df = abundances_init[(abundances_init['%s_rel_sun' % element] > -2)]
         elfe = abundances_df['%s_rel_sun' % element] - abundances_df['feh'] + 0.02
         max_elfe = np.max(elfe)
         min_elfe = np.min(elfe)
         extreame_all_elfe.append(min_elfe)
         extreame_all_elfe.append(max_elfe)
     min_all_elfe = math.ceil(np.min(extreame_all_elfe)*10)/10-0.1
     max_all_elfe = math.ceil(np.max(extreame_all_elfe)*10)/10


     for i, element in enumerate(elements_df.element):
         print (element)
         harps_best = harps_sample[(harps_sample['%s' % element] < 1) & (harps_sample['%s' % element] > -2)]
         abundances_df = abundances_init[(abundances_init['%s_rel_sun' % element] > -2)]
         elfe = abundances_df['%s_rel_sun' % element] - abundances_df['feh'] + 0.02
         elfe_err = (abundances_df['%s_err' % element]**2 + abundances_df['erfeh']**2)**0.5

         giants = abundances_df[abundances_df.logg < 3.5]
         elfe_giants = giants['%s_rel_sun' % element] - giants['feh'] + 0.02
         elfe_err_giants = (giants['%s_err' % element]**2 + giants['erfeh']**2)**0.5

         dwarfs = abundances_df[abundances_df.logg >= 3.5]
         elfe_dwarfs = dwarfs['%s_rel_sun' % element] - dwarfs['feh'] + 0.02
         elfe_err_dwarfs = (dwarfs['%s_err' % element]**2 + dwarfs['erfeh']**2)**0.5

         plt.figure(1)
         sns.set_style("white")
         sns.set_style("ticks")
         sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
         fig = plt.gcf()
         #Fig size depends on the number of columns and rows
         fig.set_size_inches(3*n_columns, 3*n_rows)
         fig.subplots_adjust(left=0.1, right=0.98, bottom = 0.07, top = 0.98, wspace=0.0, hspace=0.0)
         ax = plt.subplot(n_rows, n_columns, i+1)
         ax.text(.8, .85,'%s' % element, horizontalalignment='right', transform=ax.transAxes, fontsize = 14)
         plt.scatter(harps_best.logg, (harps_best['%s' % element]-harps_best.feh), c = 'grey', s = 20, alpha = 0.5, edgecolors='none')
         plt.scatter(dwarfs['logg'], elfe_dwarfs, linewidth=1, color='black', alpha = 0.6, s = 70)
         plt.errorbar(dwarfs['logg'], elfe_dwarfs, yerr= elfe_err_dwarfs, xerr= dwarfs['erfeh'], color='black', fmt='o', markersize=1, linewidth=0.5, alpha = 0.6)
         plt.scatter(giants['logg'], elfe_giants, linewidth=1, color='red', alpha = 0.6, s = 90)
         plt.errorbar(giants['logg'], elfe_giants, yerr= elfe_err_giants, xerr= giants['erfeh'], color='red', fmt='o', markersize=1, linewidth=0.5, alpha = 0.6)
         plt.axhline(y=0, color='blue', linestyle='--', linewidth=1.)
         plt.axvline(x=4.44, color='blue', linestyle='--', linewidth=1.)
         x_min = math.ceil(np.min(abundances_df['logg'])*10)/10-0.2
         x_max = math.ceil(np.max(abundances_df['logg'])*10)/10+0.1

         plt.ylim(-0.5, 0.8)
         plt.xlim(3.25, 4.75)
         #plt.xticks(np.arange(2.5,4.6, 0.5))
         #plt.yticks(np.arange(-0.2,0.4, 0.1))
         #plt.minorticks_on()
         plt.tick_params(labelsize=13)


         bottom_left, top_left, top_right, bottom_right = tick_parameters(n_rows, n_columns)

         if i in top_right:
             plt.tick_params(
             axis='y',          # changes apply to the x-axis
             which='both',      # both major and minor ticks are affected
             right='on',      # ticks along the bottom edge are off
             left='on',         # ticks along the top edge are off
             labelleft='off')
             plt.tick_params(
             axis='x',          # changes apply to the x-axis
             which='both',      # both major and minor ticks are affected
             top='on',      # ticks along the bottom edge are off
             bottom='on',         # ticks along the top edge are off
             labelbottom='off')


         if i in top_left:
             plt.tick_params(
             axis='y',          # changes apply to the x-axis
             which='both',      # both major and minor ticks are affected
             right='on',      # ticks along the bottom edge are off
             left='on',         # ticks along the top edge are off
             labelleft='on', labelsize=13)
             plt.tick_params(
             axis='x',          # changes apply to the x-axis
             which='both',      # both major and minor ticks are affected
             top='on',      # ticks along the bottom edge are off
             bottom='on',         # ticks along the top edge are off
             labelbottom='off')

         if i in bottom_right:
             plt.tick_params(
             axis='y',          # changes apply to the x-axis
             which='both',      # both major and minor ticks are affected
             right='on',      # ticks along the bottom edge are off
             left='on',         # ticks along the top edge are off
             labelleft='off')
             plt.tick_params(
             axis='x',          # changes apply to the x-axis
             which='both',      # both major and minor ticks are affected
             top='on',      # ticks along the bottom edge are off
             bottom='on',         # ticks along the top edge are off
             labelbottom='on',labelsize=12)

         if i in bottom_left:
             plt.tick_params(
             axis='y',          # changes apply to the x-axis
             which='both',      # both major and minor ticks are affected
             right='on',      # ticks along the bottom edge are off
             left='on',         # ticks along the top edge are off
             labelleft='on',labelsize=13)
             plt.tick_params(
             axis='x',          # changes apply to the x-axis
             which='both',      # both major and minor ticks are affected
             top='on',      # ticks along the bottom edge are off
             bottom='on',         # ticks along the top edge are off
             labelbottom='on',labelsize=12.)

         plt.savefig('elfe_logg.pdf', format='pdf')
         fig.text(0.55, 0.01, '$\log g$', ha='center', fontsize = 17)
         fig.text(0.02, 0.5, '[X/Fe]', va='center', rotation='vertical', fontsize = 17)
     plt.clf()


def HR():
     abundances_init = pd.read_table("all_abundances_ordered_star.rdb")

     plt.figure(1)
     sns.set_style("white")
     sns.set_style("ticks")
     sns.set_style({"xtick.direction": "in","ytick.direction": "in"})
     fig = plt.gcf()

     plt.scatter(abundances_init['teff'], abundances_init['logg'], linewidth=1, color='black')
     plt.gca().invert_yaxis()
     plt.gca().invert_xaxis()
     plt.ylabel('log g')
     plt.xlabel('Teff')
     plt.savefig('HR.pdf', format='pdf')

#HR()
#elfe_teff()
elfe_feh()
#elfe_logg()
