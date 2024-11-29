#!/usr/bin/python
#Imports

import pandas as pd
import numpy as np

def find_loc(tablename, elem):
    """
    Find the last occurrence of `elem` in the 'element' column of `tablename`.
    """
    match = tablename[tablename['element'] == elem]
    if match.empty:
        raise ValueError(f"Element '{elem}' not found in the table.")
    return match.iloc[-1].name

def average_abund(tablename, elem):
    """
    Compute and update the average abundance and total error for an element.
    """
    try:
        elI = find_loc(tablename, f"{elem}I")
        el = find_loc(tablename, elem)
        elII = find_loc(tablename, f"{elem}II")

        # Compute average abundance and error
        avg_abund = np.round((tablename.at[elI, 'abund'] + tablename.at[elII, 'abund']) / 2, 2)
        avg_err = np.round(((tablename.at[elI, 'total_err']**2 + tablename.at[elII, 'total_err']**2)**0.5) / 2, 2)

        # Update values in tablename
        tablename.at[el, 'abund'] = avg_abund
        tablename.at[el, 'total_err'] = avg_err
    except KeyError as e:
        raise KeyError(f"Missing column in the table: {e}")
    except ValueError as e:
        raise ValueError(e)

def abundances_results_formatting_column(input_starlist, element_file, output_file):
    """
    Format abundances into columns with stars as separate entries.
    """
    input_data = pd.read_table(input_starlist)
    element_info = pd.read_table(element_file)


    for star in input_data['star']:
        print(f"Processing star: {star}")
        star_abund_file = f'./moog_abundances/moog_abstar/{star}_abund_err.dat'
        star_abund_data = pd.read_table(star_abund_file)[['element', 'abund', 'total_err']]
        merged_data = element_info.merge(star_abund_data, on='element', how='outer')
        merged_data = merged_data.set_index('element').reindex(element_info['element']).reset_index()   

        # Update columns in element_info
        element_info[f'{star}_abund'] = merged_data['abund']
        element_info[f'{star}_abund_err'] = merged_data['total_err']
        element_info[f'{star}_abund_rel_sun'] = merged_data['abund'] - merged_data['sun_vesta']

    element_info.to_csv(output_file, sep='\t', index=False, float_format='%.3f')

def abundances_results_formatting_line(input_starlist, column_file, output_file):
    """
    Format abundances into rows with elements as separate entries.
    """
    input_data = pd.read_table(input_starlist)
    abundances_formatted_column = pd.read_table(column_file)

    for j, element in enumerate(abundances_formatted_column.element):
        print (element)
        input_data['%s' % element] = 0.0
        input_data['%s_err' % element] = 0.0
        input_data['%s_rel_sun' % element] = 0.0
        for i, star in enumerate(input_data.star):
            input_data['%s' % element][i] = abundances_formatted_column['%s_abund' % star][j]
            input_data['%s_err' % element][i] = abundances_formatted_column['%s_abund_err' % star][j]
            input_data['%s_rel_sun' % element][i] = abundances_formatted_column['%s_abund_rel_sun' % star][j]

    
        print (input_data)

    input_data.to_csv(output_file, sep='\t', index=False, float_format='%.3f')


# Input files
input_starlist = 'input_param_error.rdb'
element_file = 'elements_formatting.rdb'
output_column_file = 'all_abundances_ordered_element.rdb'
output_line_file = 'all_abundances_ordered_star.rdb'

# Step 1: Format abundances into columns
abundances_results_formatting_column(input_starlist, element_file, output_column_file)

# Step 2: Format abundances into rows
abundances_results_formatting_line(input_starlist, output_column_file, output_line_file)

