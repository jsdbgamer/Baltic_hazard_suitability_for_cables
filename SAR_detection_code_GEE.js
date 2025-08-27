// SAR detection layer for the month of November 2024
// Author: Shaaz
// Date: 10/8/2025

// Bands:
// R = VV/VH ratio (95th percentile)
// G = VV (95th percentile)
// B = VH (95th percentile)

// Notes:
// Sentinel-1 GRD values in GEE are linear sigma0.
// Exported in 8-bit RGB for computational efficiency.

// 1. Polygon shape
var aoi = ee.Geometry.Polygon(
  [[[16.6984, 56.5107],
    [16.7314, 55.2287],
    [20.8073, 55.1911],
    [20.9729, 55.9702],
    [21.2037, 56.2337],
    [21.4673, 57.5478],
    [17.2156, 57.5478]]]);

// 2. Load Sentinel-1 GRD with IW mode and VV+VH covering the 
// whole month of November 2024
var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(aoi)
  .filter(ee.Filter.eq('instrumentMode','IW'))
  .filterDate('2024-11-01','2024-11-30')
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation','VH'));

// 3. Computing the 95th-percentile (VV & VH)
var vv95 = s1.select('VV').reduce(ee.Reducer.percentile([95])).rename('VV_95');
var vh95 = s1.select('VH').reduce(ee.Reducer.percentile([95])).rename('VH_95');

// 4. Build a ratio band (VV/VH)
var ratio95 = vv95.divide(vh95).rename('ratio_95');

// 5. Stack into a 3-band image and clip
var rgb95 = ee.Image.cat([ratio95, vv95, vh95]).clip(aoi);

// 6. Visualisation parameters
var visParams = {
  bands: ['ratio_95','VV_95','VH_95'],
  min:   [0.5,      -25,   -25],
  max:   [5,         5,      0],
  gamma: [1,         1,      1]
};

// 7. Preview of the 8-bit RGB image
var rgb8 = rgb95.visualize(visParams);

// 8. Show on map
Map.centerObject(aoi, 7);
Map.addLayer(rgb8, {}, '95th-pct RGB (8-bit)', false);

// 9. Export and save (customised to be in UTM 33N)
Export.image.toDrive({
  image:         rgb8,                          
  description:   'S1_95pct_RGB_Baltic_Nov2024',
  folder:        'EarthEngineExports',
  fileNamePrefix:'S1_95pct_RGB_8bit',
  region:        aoi,
  scale:         10,
  crs:           'EPSG:32633',
  maxPixels:     1e10,
  fileFormat:    'GeoTIFF',
  formatOptions: {
  }
});

// End of code