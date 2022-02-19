# River Discharge Data Download

Currently, we are working with GloFAS-ERA5. This can be downloaded from: https://cds.climate.copernicus.eu/cdsapp#!/dataset/cems-glofas-historical?tab=form. An account with Copernicus ECMWF must be created before one can download. ECMWF has developed a Python API package to handle the downloading of the data. This package is called `cdsapi` and can be installed via `pip install cdsapi`. One can select the parameters at the link above, which will generate a python query for `cdsapi` that looks like;

```
import cdsapi

c = cdsapi.Client()

c.retrieve(
    'cems-glofas-historical',
    {
        'variable': 'river_discharge_in_the_last_24_hours',
        'format': 'grib',
        'system_version': 'version_3_1',
        'hyear': '1995',
        'hmonth': [
            'april', 'august', 'december',
            'february', 'january', 'july',
            'june', 'march', 'may',
            'november', 'october', 'september',
        ],
        'hydrological_model': 'lisflood',
        'hday': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
            '13', '14', '15',
            '16', '17', '18',
            '19', '20', '21',
            '22', '23', '24',
            '25', '26', '27',
            '28', '29', '30',
            '31',
        ],
        'area': [
            89.95, -100, -35,
            40,
        ],
        'product_type': 'consolidated',
    },
    'glofas-era5_1995.grib')
```

Once the file is download, Climate Data Operators (cdo) can be used to conver the file to netcdf.

```
cdo -f nc copy glofas-era5_1995.grib glofas-era5_1995.nc
```
