## This code blends two workflows to generate a multiband raster of landsat-derived 
## predictor variables for spring, summer, and fall. User can define year or years of interest.

### Wesley Rancher, Hana Matsumoto
### 10 July 2024 // cleaned:15 August 2024

### Hurni, K., Würsch, L. & Heinimann, A. (2017): Google Earth Engine Image Pre-processing Tool.
### https://doi.org/10.1016/j.rse.2019.111225
### Massey et al., 2023 //  follow their code instructions for proper citation

# %% start GEE session
import ee
import geemap
import pandas as pd

#ee.Reset()
ee.Authenticate()
ee.Initialize()

# %% Read in some ROI files and sampling data
#AK_landscape = ee.FeatureCollection('projects/ee-vegshiftsalaska/assets/LandisModelRegion') 
AK_landscape = ee.FeatureCollection('projects/ee-vegshiftsalaska/assets/Dalton_Landis') 

Map = geemap.Map(center=[64, -152], zoom=5, basemap='Esri.WorldGrayCanvas')

# read in external predictors
soil = ee.Image('projects/ee-vegshiftsalaska/assets/soil_classes_clipped_250m')
topo = ee.Image('projects/ee-vegshiftsalaska/assets/dem_clipped_30m')

# rename the bands
soil_renamed = soil.rename(["soil_classes"])
topo_renamed = topo.rename(['elev', 'slope', 'aspect'])

# read in CAFI data
cafiPlots = ee.FeatureCollection('projects/ee-vegshiftsalaska/assets/CAFI_SITE_INTERIOR')

# read in climate and disturbance vars
Clim_vars = ee.Image("WORLDCLIM/V1/BIO").select('bio01', 'bio04', 'bio05', 'bio06', 'bio12', 
                                                'bio13', 'bio14', 'bio15')
Clim_vars_renamed = Clim_vars.rename('annual_mean_temp', 'temp_seasonality', 'max_temp_warmest_month', 
                                     'min_temp_coldest_month', 'annual_precip', 'precip_wettest_month', 
                                     'precip_driest_month', 'precip_seasonality')

# %% Define years and functions for calculating SVIs
years = ee.List.sequence(2000, 2023)
year_info = years.size().getInfo()

# functions for vegetation indices calculations
bands = ee.List(['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'ndvi', 'evi', 'mndwi', 'nbr', 'vari', 'savi', 'tcb', 'tcg', 'tcw'])

# normalized difference vegetation index
def ndvi_calc(img):
    ndvi = img.normalizedDifference(['nir', 'red']) \
              .rename('ndvi') \
              .multiply(10000.0) \
              .toInt16()
    return ndvi
# enhanced vegetation index
def evi_calc(img):
    # coefficients
    c1 = 2.5
    c2 = 6 
    c3 = 7.5
    c4 = 1.0

    # define bands
    red = img.select('red')
    nir = img.select('nir')
    blue = img.select('blue')

    # equation
    evi = nir.subtract(red) \
        .divide(nir.add(red.multiply(c2)).subtract(blue.multiply(c3)).add(c4)).multiply(c1) \
        .rename('evi') \
        .multiply(10000.0) \
        .toInt16()
    return evi
# modified normalized difference water index
def mndwi_calc(img):
    mndwi = img.normalizedDifference(['green', 'swir2']) \
              .rename('mndwi') \
              .multiply(10000.0) \
              .toInt16()
    return mndwi
# nbr
def nbr_calc(img):
    nbr = img.normalizedDifference(['nir', 'swir2']) \
            .rename('nbr') \
            .multiply(10000.0) \
            .toInt16()
    return nbr
# visible atmospherically resistant index
def vari_calc(img):
    vari = img.select('red') \
              .subtract(img.select('green')) \
              .divide(img.select('red').add(img.select('green')).subtract(img.select('blue'))) \
              .rename('vari') \
              .multiply(10000.0) \
              .toInt16()
    return vari
# soil adjusted vegetation index
def savi_calc(img):
    L = 0.5  # coefficient

    # bands
    nir = img.select('nir')
    red = img.select('red')

    # equation
    savi = nir.subtract(red).multiply(1+L) \
        .divide((nir).add(red).add(L)) \
        .rename('savi') \
        .multiply(10000.0) \
        .toInt16()
    return savi
# tasseled cap transformations // are these coefficients correct for our imagery?
def tasseled_cap(img):
    # coefficients
    brightness_coeff = ee.Image([0.2043, 0.4158, 0.5524, 0.5741, 0.3124, 0.2303])
    greenness_coeff = ee.Image([-0.1603, -0.2819, -0.4934, 0.7940, 0.0002, -0.1446])
    wetness_coeff = ee.Image([0.0315, 0.2021, 0.3102, 0.1594, -0.6806, -0.6109])

    # bands to apply coefficients to
    # img = img.select(['B2', 'B3', 'B4', 'B5', 'B6', 'B7'])
    img = img.select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2'])

    # math
    brightness = (img.multiply(brightness_coeff) \
        .reduce(ee.Reducer.sum()) \
        .rename('tcb') \
        .multiply(10000.0) \
        .toInt16())
    greenness = (img.multiply(greenness_coeff) \
        .reduce(ee.Reducer.sum()) \
        .rename('tcg') \
        .multiply(10000.0) \
        .toInt16())
    wetness = (img.multiply(wetness_coeff) \
        .reduce(ee.Reducer.sum()) \
        .rename('tcw') \
        .multiply(10000.0) \
        .toInt16())
    return ee.Image([brightness, greenness, wetness])

# %% Add indices and band names functions
def add_indices(image):
    scaled_image = image.toFloat().divide(10000.0)
    return image.addBands(ndvi_calc(scaled_image))\
                .addBands(evi_calc(scaled_image))\
                .addBands(mndwi_calc(scaled_image))\
                .addBands(nbr_calc(scaled_image))\
                .addBands(vari_calc(scaled_image))\
                .addBands(savi_calc(scaled_image))\
                .addBands(tasseled_cap(scaled_image))

# add suffix to all band names
def add_suffix(in_image, suffix_str):
    # convert band names to lowercase and add suffix
    def append_suffix(band_name):
        return ee.String(band_name).toLowerCase().cat('_').cat(suffix_str)
    # apply the suffix to all band names
    bandnames = in_image.bandNames().map(append_suffix)
    # the number of bands
    nb = bandnames.size()
    # select bands with the new names
    return in_image.select(ee.List.sequence(0, ee.Number(nb).subtract(1)), bandnames)

# %% Landsat sensor corrections and cloud masks

# Hurni harmonization technique (option in lieu of Massey approach)
def harmonization(img):
    slopes = ee.Image([0.8474, 0.8483, 0.9047, 0.8462, 0.8937, 0.9071])
    intercepts = ee.Image([0.0003, 0.0088, 0.0061, 0.0412, 0.0254, 0.0172])
    img_harm = img.select(['blue', 'green', 'red', 'nir', 'swir1', 'swir2']) \
                  .multiply(slopes) \
                  .add(intercepts.multiply(10000)) \
                  .int16()
    return img.select().addBands(img_harm).addBands(img.select('pixel_qa'))


# function for cloud mask images in collection // includes water mask
def cloud_mask_landsat8(img):
    img = ee.Image(img)
    quality_band = img.select('pixel_qa') #fmask for hls and QA_PIXEL for landsat
    water = quality_band.bitwiseAnd(1).neq(0) 
    shadow = quality_band.bitwiseAnd(8).neq(0)  
    cloud = quality_band.bitwiseAnd(32).neq(0) 
    cloud_confidence = quality_band.bitwiseAnd(64).add(quality_band.bitwiseAnd(128)).interpolate([0, 64, 128, 192], [0, 1, 2, 3], 'clamp').int()
    cloud_confidence_medium_high = cloud_confidence.gte(2)
    cloudM = water.Or(shadow).Or(cloud).Or(cloud_confidence_medium_high).select([0], ['cloudM'])
    
    # add cirrus confidence to cloud mask (cloudM) for Landsat 8
    cirrus_confidence = quality_band.bitwiseAnd(256).add(quality_band.bitwiseAnd(512)).interpolate([0, 256, 512, 768], [0, 1, 2, 3], 'clamp').int()
    cirrus_confidence_medium_high = cirrus_confidence.gte(2)
    cloudM = cloudM.Or(cirrus_confidence_medium_high)
    cloudM = cloudM.Not()
 
    # mask image with cloud mask and add as band
    image_cloud_masked = img.updateMask(cloudM).addBands(cloudM)
    return image_cloud_masked

# renaming bands
def get_landsat_collection_sr(sensor):
    if sensor in ['LC08', 'LC09']:
        bands = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'QA_PIXEL']
    else:  # for Landsat 5 and 7
        bands = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'QA_PIXEL']

    band_names_landsat = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'pixel_qa']
    cloud_threshold = 80 
    collection_filtered_without_date = ee.ImageCollection('LANDSAT/' + sensor + '/C02/T1_L2') \
        .filterBounds(AK_landscape)\
        .filterMetadata('CLOUD_COVER', 'less_than', cloud_threshold) \
        .select(bands, band_names_landsat)  #rename bands
    collection_filtered_without_date = collection_filtered_without_date.map(harmonization)
    return collection_filtered_without_date

# %% Iterate
pv_dataframe_landis = []
for i in range(year_info):
    year = years.get(i).getInfo()
    startDate = ee.Date.fromYMD(year, 1, 1)
    endDate =   ee.Date.fromYMD(year, 12, 31)

    #seasonal windows
    startSeason1 = ee.Date.fromYMD(year, 3, 1)
    endSeason1 = ee.Date.fromYMD(year, 5, 31)
    startSeason2 = ee.Date.fromYMD(year, 6, 1)
    endSeason2 = ee.Date.fromYMD(year, 8, 31)
    startSeason3 = ee.Date.fromYMD(year, 9, 1)
    endSeason3 = ee.Date.fromYMD(year, 11, 30)
        
    ################################################
    #get this to work for startSeason endSeason
    def get_landsat_Images(sensor, AK_landscape, startDate, endDate, startSeason, endSeason):  
        collection = get_landsat_collection_sr(sensor)
        cleaned_images = (collection
                            .filterBounds(AK_landscape)
                            .filterDate(startDate, endDate)
                            .filterDate(startSeason, endSeason)
                            .map(cloud_mask_landsat8)
                            .map(add_indices))
        return cleaned_images

    # function to bring everything together            
    def make_ls(AK_landscape):
        sensors = ['LC08', 'LC09', 'LE07', 'LT05']
        spring_images = ee.ImageCollection([])
        summer_images = ee.ImageCollection([])
        fall_images = ee.ImageCollection([])

        for sensor in sensors:
            spring_images = spring_images.merge(get_landsat_Images(sensor, AK_landscape, startDate, endDate, startSeason1, endSeason1))
            summer_images = summer_images.merge(get_landsat_Images(sensor, AK_landscape, startDate, endDate, startSeason2, endSeason2))
            fall_images = fall_images.merge(get_landsat_Images(sensor, AK_landscape, startDate, endDate, startSeason3, endSeason3))
        # Median composites by season
        springtime = add_suffix(spring_images.median().select(bands), '1').unmask(-9999)
        summertime = add_suffix(summer_images.median().select(bands), '2').unmask(-9999)
        falltime = add_suffix(fall_images.median().select(bands), '3').unmask(-9999)

        # add each composite as bands
        return springtime.addBands(summertime).addBands(falltime)


    # apply the function appending the image at each iteration to our list of images
    landsat_composite = make_ls(AK_landscape).clip(AK_landscape).reproject(crs='EPSG:3338', scale=30)
    full_composite = landsat_composite.addBands(soil_renamed)\
                                      .addBands(topo_renamed)\
                                      .addBands(Clim_vars_renamed)\
                                      .clip(AK_landscape)\
                                      .reproject(crs='EPSG:3338', scale=30)\
                                      .int16()
                                     

    #sampling and exporting
    img = full_composite
    def get_pixel_values(f, img):
        mean = full_composite.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=f.geometry(),
            scale=30,
            crs='EPSG:3338')
        return f.setMulti(mean)

    
    pv_sampling = cafiPlots.map(lambda f: get_pixel_values(f, img))
    

    ################ export tasks 
    # Image composite to G drive
    # Sampling data to local or g drive

    #export image to google drive
    geom = AK_landscape.geometry()
    task = ee.batch.Export.image.toDrive(
        image=full_composite,
        description=f'FullComp_Dalton_{year}',
        folder='Alaska_Proj',
        region=geom,
        scale=30,
        crs='EPSG:3338',
        maxPixels=1e13)
    task.start()

    # apply the sample function and export locally for interior
    #pv_values_interior = pv_sampling.getInfo()



    # apply the sample function and export locally for model region
    # for feature in pv_values_interior['features']:
    #     properties = feature['properties']
    #     pv_dataframe_landis.append(properties)
    # extracted_pv_vals_landis = pd.DataFrame(pv_dataframe_landis)
    # extracted_pv_vals_landis.to_csv(f'D:/MS-Research/alaska/data/output/csv/cleaned/GEE_Sampling_Data_BiomassLandisRegion/PixelVals_Interior_{year}.csv')
    
    # # or to g drive
    # landis_region_task = ee.batch.Export.table.toDrive(
    #     collection=pv_sampling_landis,
    #     description=f'PixelVals_Landis_{year}_V2',
    #     folder='Alaska_Proj',
    #     fileFormat='CSV')
    # landis_region_task.start()
    
    
# %%