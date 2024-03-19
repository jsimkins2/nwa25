import xarray as xr
import matplotlib.pyplot as plt
import os
import pandas as pd
import cmocean


# Create the output directory if it doesn't exist
output_dir = './transects'
os.makedirs(output_dir, exist_ok=True)

# Load the netCDF file
#file_path = '/Users/james/Documents/Rutgers/output/nwa25/sim4.5/transects/test2.nc'
file_path = '19970101.ocean_5day.nc'
ds = xr.open_dataset(file_path)

ds = ds.drop({'average_T1','average_T2','average_DT','salt','volcello','dzRegrid'})
ds = ds.isel(time=-1)
# Function to plot and save vertical transects
def plot_transects(variable, var_name, lat_var, depth_var, output_dir):
    latitudes = ds[lat_var].values
    depths = ds[depth_var].values
    
    time = pd.to_datetime(str(ds.time.values))
    time_str = time.strftime('%Y%m%d')
    
    # Loop over the latitudes at 10-degree intervals
    for lat in range(6, 50, 10):
        # Find the nearest latitude index
        lat_idx = abs(latitudes - lat).argmin()
        
        # Select the data slice
        transect = ds[variable].isel(**{lat_var: lat_idx})
        
        # Plotting
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.set_facecolor('dimgray')  # Set the background color to grey
        cmap = cmocean.cm.balance if variable == 'v' else cmocean.cm.thermal
        transect.plot(y=depth_var, yincrease=False, ax=ax, cmap=cmap, add_colorbar=True)
        plt.title(f'{var_name} Transect at {lat}° Latitude')
        plt.xlabel('Longitude')
        plt.ylabel('Depth (m)')
        
        # Save the figure
        plt.savefig(f'{output_dir}/{variable}_transect_{lat}deg_{time_str}.png')
        plt.close()

# Plot and save transects for 'temp' (potential temperature) and 'v' (meridional velocity)
plot_transects('temp', 'Potential Temperature', 'yh', 'z_l', output_dir)
plot_transects('v', 'Meridional Velocity', 'yq', 'z_l', output_dir)

# Close the xarray dataset
ds.close()