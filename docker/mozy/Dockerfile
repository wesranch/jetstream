# Start from a minimal Python 3 base image
FROM python:3.10-slim

#gdal
RUN apt-get update && apt-get install -y \
build-essential \
gdal-bin \
libgdal-dev \
&& rm -rf /var/lib/apt/lists/*
    
ENV GDAL_CONFIG=/usr/bin/gdal-config

#python libs
RUN pip install --no-cache-dir \
google-auth google-auth-oauthlib google-auth-httplib2 \
gdown \
dask[complete] \
#rasterio \
geopandas \
geemap \
earthengine-api

#can be done later
RUN mkdir /scripts

COPY ./mosaicgdown.py /scripts
#RUN earthengine authenticate --quiet
EXPOSE 8888
#if your script is ready
#CMD ["python", "mosaic_script.py"]
