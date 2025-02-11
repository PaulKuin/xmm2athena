import numpy as np
from astroML.time_series import lomb_scargle, lomb_scargle_bootstrap
from astroML.plotting import setup_text_plots
setup_text_plots(fontsize=8, usetex=True)
import matplotlib.pyplot as plt

def classify_light_curve(time, flux, err):

    # Compute Lomb-Scargle periodogram
    omega = 2 * np.pi * np.linspace(0.01, 10, 1000) # frequencies
    LS = lomb_scargle(time, flux, dy=err,omega=omega)

    # Bootstrap to estimate significance
    LS_bootstrap = lomb_scargle_bootstrap(time, flux, dy=err,
        #dy=np.ones_like(flux),
        N_bootstraps=100, random_state=0)

    # Compute the power and significance level at a few frequencies
    power = lomb_scargle(time, flux, dy=err, omega=omega)
    significance = lomb_scargle_bootstrap(time, flux, dy=err,
                                           N_bootstraps=100, random_state=0,
                                           omega=omega)

    # Classify light curve based on the maximum power and significance
    max_power_freq = omega[np.argmax(power)]
    max_significance_freq = omega[np.argmax(significance)]
    
    # Define thresholds for classification
    periodic_threshold = 0.05
    stochastic_threshold = 0.9
    
    # Classify based on thresholds
    if max_power_freq > periodic_threshold:
        classification = 'Periodic'
    elif max_significance_freq < stochastic_threshold:
        classification = 'Stochastic'
    else:
        classification = 'Eruptive'
    
    return classification

# Example usage:
# Replace 'time' and 'flux' with your actual time and flux data arrays
# time = ...
# flux = ...
# classification = classify_light_curve(time, flux)
# print("Classification:", classification)

'''
chatgpt -- 
In this code:

The lomb_scargle function computes the Lomb-Scargle periodogram of the light curve.

The lomb_scargle_bootstrap function performs bootstrap resampling to estimate the significance levels of the periodogram peaks.

The power and significance of the periodogram are computed at various frequencies.

The light curve is classified as periodic if the maximum power frequency exceeds 
a threshold, as stochastic if the maximum significance frequency is below a 
threshold, and as eruptive otherwise.

You'll need to replace 'time' and 'flux' with your actual time and flux data arrays. 

Additionally, you may need to adjust the threshold values for your specific 
dataset and classification criteria.

'''