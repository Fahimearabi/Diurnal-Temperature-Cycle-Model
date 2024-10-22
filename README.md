The code named "Soil_Heat_Flux_at_DAY" and "Soil_Heat_Flux_at_night" are for estimating soil heat flux during day and night, respectively. To estimate soil heat flux using MODIS sensor images, the product of day and night land surface temperature, emissivity, albedo and NDVI vegetation index is needed. These images should be resized and their spatial resolution should be the same. Using these two codes, you can calculate the heat flux of the soil 4 times every day.
	
 By using the code "Modeling_the_daily_cycle_of_soil_heat_flux" and four images of soil heat flux every day, the daily cycle of soil heat flux can be estimated and the hourly time series of soil heat flux can be obtained. In this method, sunrise, sunset, and day length are calculated for each date, and changes in soil heat flux are modeled in each day.
	
 In the code named "thermal_inertia_final", the thermal inertia estimation method is presented using Modis sensor images, which requires four images of the land surface temperature every day and night and albedo.
	
 In the "cod_G_INERTIA" code, the soil heat flux estimation method is presented using the harmonic analysis method.
