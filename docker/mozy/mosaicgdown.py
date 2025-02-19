#%% import gdown
import dask.array as da
import rasterio
from dask.diagnostics import ProgressBar
from rasterio.merge import merge
import os
import gdown
#%% download images from google drive
def download_images(file_ids, dd):
    for file_id in file_ids:
        file_path = os.path.join(dd, f'{file_id}.tif')
        gdown.download(f'https://drive.google.com/uc?id={file_id}', file_path, quiet=False)
#%% file ids
file_ids = ['1Wd4aHf2zuVbmaNR-k0wDYdoVohDuLihl', '1cjEeRfrzAK5v5Hh1LHgGNDtNbBBNVrlH', 
            '1Fa6_vfksNCGLTB9f2dzqxCge0BH148fx', '1gE1OC99b0svG4v0qZYwFjawCNBiD5G2V',
            '1WW8q1wGf6_bORXgfiDd4OfzYWugbGahf', '1atsJN2J4OvWBPFHwe6R8vt0hkT_ApVpR']

dd = '/Users/wancher/Downloads/'
download_images(file_ids, dd)
src_files = [rasterio.open(os.path.join(dd, f"{file_id}.tif")) for file_id in file_ids]

#%% read rasters in chunks
dask_arrays = []
for src in src_files:
    block_shape = (src.height // 10, src.width // 10)
    chunks = (block_shape[0], block_shape[1])
    
    #read into dask array
    arr = da.from_array(src.read(1), chunks=chunks)
    dask_arrays.append(arr)

#%% combine using dask
mosaic = da.concatenate(dask_arrays, axis=0)
with ProgressBar():
    mosaic_result = mosaic.compute() 

#%% save
output_path = '/path/to/output/mosaic.tif'
with rasterio.open(output_path, 'w', driver='GTiff', count=1, dtype='float32',
                   crs=src.crs, transform=src.transform, width=mosaic_result.shape[1], height=mosaic_result.shape[0]) as dst:
    dst.write(mosaic_result, 1)

print(f"Mosaic saved as '{output_path}'")
