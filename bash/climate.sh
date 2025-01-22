#Configure AWS for file transfer
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install --update
aws --version

echo "setting up model directory"
mkdir -p climate-processing
mkdir -p climate-processing/biomass-data climate-processing/scripts climate-processing/climate-data climate-processing/perma climate-processing/output


#think this is causing the disconnect
endpoint="https://usgs2.osn.mghpcc.org"
aws configure

# Downloading files
echo "Downloading biomass data from s3 to this instance of jetstream into "biomass-data/""
cd climate-processing
aws s3 cp s3://landis-ii-input/rs-data/biomass-data/ --endpoint-url "$endpoint" biomass-data --recursive

#Download climate vars from SNAP
echo 'fetching climate data and downloading into climate-data/...'
cd climate-data/
curl "http://data.snap.uaf.edu/data/Base/AK_771m/projected/AR5_CMIP5_models/Projected_Monthly_and_Derived_Precipitation_Products_771m_CMIP5_AR5/monthly/pr_total_mm_AK_AR5_GFDL-CM3_rcp85_01_2006-12_2100.zip" --output pr_GFDL-CM3_rcp85.zip
curl "http://data.snap.uaf.edu/data/Base/AK_771m/projected/AR5_CMIP5_models/Projected_Monthly_and_Derived_Precipitation_Products_771m_CMIP5_AR5/monthly/pr_total_mm_AK_AR5_CCSM4_rcp85_01_2006-12_2100.zip" --output pr_CCSM4_rcp85.zip
curl "http://data.snap.uaf.edu/data/Base/AK_771m/projected/AR5_CMIP5_models/Projected_Monthly_and_Derived_Temperature_Products_771m_CMIP5_AR5/monthly/tas_mean_C_AK_AR5_GFDL-CM3_rcp85_01_2006-12_2100.zip" --output tas_GFDL-CM3_rcp85.zip
curl "http://data.snap.uaf.edu/data/Base/AK_771m/projected/AR5_CMIP5_models/Projected_Monthly_and_Derived_Temperature_Products_771m_CMIP5_AR5/monthly/tas_mean_C_AK_AR5_CCSM4_rcp85_01_2006-12_2100.zip" --output tas_CCSM4_rcp85.zip

echo "unzipping files..."
unzip "*.zip"
cd ../perma
curl "http://data.snap.uaf.edu/data/Base/Other/LandCarbon/Permafrost/mu_permafrost_0_100_2.zip" --output permafrost.zip
unzip permafrost.zip
cd ../

#dependency for terra
echo "installing gdal so we can load terra in R"
sudo apt update
sudo apt install -y gdal-bin libgdal-dev
gdal-config --version