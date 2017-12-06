# README #

Codes to compare ECCC model GEM-MACH-GHG output to TCCON columns.

Not "user friendly"

## smoothing.py

Reads TCCON and GEM-MACH-GHG data.

Homogenize TCCON a priori and GEM-MACH-GHG profiles (time coincidence and same atmospheric levels)

Computes TCCON a priori columns from the apriori profiles.

Computes GEM-MACH-GHG columns.

Use TCCON averaging kernel to smooth the data.

Save the data to be read by analysis.py for statistics and plots.

This version of the code uses the GGG 2014 data output/format for the TCCON measurements.

There are 3 "MODIFY HERE" parts in the code that should be the only needed modifications if the model data format does not change

## analysis.py

Reads the output of smoothing.py and creates plots of the data

The plotting part was based on bokeh 0.12.4 ; the latest version of bokeh is now 0.12.11 and the code was not tested with it

### Who do I talk to? ###

sebastien.roche@mail.utoronto.ca