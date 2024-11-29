#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd

from ares_ew import EWMeasurements
from atmospheres_kurucz import AtmosphereModeling
#from abundances_3sigma_outlier import AbundanceDetermination

class SpectralAnalysisPipeline:
    def __init__(self, stellar_params):
        self.stellar_params = stellar_params
        self.ew_measurements = EWMeasurements(stellar_params)
        self.atmosphere_modeling = None
        self.abundance_determination = None

    def measure_EWs(self):
        self.ew_measurements.ew_measurements()
        print("Equivalent Widths (EWs) measured.")

    
    def model_atmospheres(self):
        self.atmosphere_modeling = AtmosphereModeling(self.stellar_params)
        self.atmosphere_modeling.model_atmospheres()
        print("Atmosphere models created.")

    """
    def determine_abundances(self):
        self.abundance_determination = AbundanceDetermination(self.stellar_params)
        self.abundance_determination.calculate_abundances()
        print("Abundances determined.")
    """

    def run_pipeline(self):
        self.measure_EWs()
        self.model_atmospheres()
        #self.determine_abundances()
        print("Spectral analysis pipeline completed.")

""" 
Example usage
stellar_params = ...  # Load your spectra data
pipeline = SpectralAnalysisPipeline(stellar_params)
pipeline.run_pipeline()
"""
stellar_params = pd.read_table('input_param_error.rdb')
pipeline = SpectralAnalysisPipeline(stellar_params)
pipeline.run_pipeline()