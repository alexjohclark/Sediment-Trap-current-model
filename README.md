READ_ME
Model and datasets for submitted publication "Clark et al. (2024), Clumped isotope temperatures of coccolithophores from global sediment traps, Paleoceanography and Paleoclimatology, submitted".

Model files:
  d18O_surface.m
    d18O of the surface ocean txt file from NASA into matrix around set location
  findcarb.m
    Function to derive carbonate parameters around given location
    
  Run_model_file.m
    Converts nc datasets from Multi Observation Global Ocean to cell array
  curpred.m
    Function of a sinking particle until depth for given time period and year
  temppred.m
    Function deriving the sinking path of a particle and temperature from origin
  temppred2.m
    Function deriving the sinking path of a particle and temperature at place of sediment trap
  temppred_avg.m
    Function deriving the sinking path of a particle and temperature from origin with averaged paths
  temppred_month.m
    Function to find what months fit D47 temperatures

Data files:
  d18O_dataset.nc
    Dataset of d18O of the surface ocean
  Model_output_timing.xlsx
    Dataset containing the longitude and latitude of the approximate provenance location, timing of when the 
    sediment trap was open and would capture coccoliths for both without and with the current model
  Model_output_curpred_publ.xlsx
    Dataset documenting the movement in the x and y direction, sinking depths, 
    and total distance travelled for each sediment trap location for 75, 100, and 125 mday-1 sinking velocity 
  Model_output_tempred_75_publ.xlsx
    Dataset showing heterogeneity of temperatures for the sinking velocity of 75 mday-1 for the start, average, 
    and end of the opening of the sediment trap for each depth level for the capture cone around the location 
    of the sediment trap and the provenance location
  Model_output_tempred_100_publ.xlsx
    Dataset showing heterogeneity of temperatures for the sinking velocity of 100 mday-1 for the start, average, 
    and end of the opening of the sediment trap for each depth level for the capture cone around the location 
    of the sediment trap and the provenance location 
  Model_output_tempred_125_publ.xlsx
    Dataset showing heterogeneity of temperatures for the sinking velocity of 125 mday-1 for the start, average, 
    and end of the opening of the sediment trap for each depth level for the capture cone around the location 
    of the sediment trap and the provenance location
