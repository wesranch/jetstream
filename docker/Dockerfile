FROM ubuntu:22.04

#need to get to build in js2

################################################################
# PREPARATIONS
################################################################

# Set the DEBIAN_FRONTEND environment variable to noninteractive
# Avoid prompts when installing packages
ENV DEBIAN_FRONTEND=noninteractive

# PREPARING PACKAGES AND UPGRADING, AND INSTALLING DEPENDENCIES OF LANDIS-II
# apt-get clean and rm -rf /var/lib/apt/lists/* are used to clear the packages folder caches
# to avoid putting them in the docker layerin system, freeing space
RUN apt-get update && apt-get upgrade -y \
&& apt-get install -y sudo apt-transport-https \
&& apt-get install -y ca-certificates gnupg curl software-properties-common \
&& apt-get install -y wget vim pip nano git python3 python-is-python3 \
&& apt-get install -y libjpeg62 \
&& apt-get install -y libpng16-16 \
&& apt-get install -y gdal-bin \
&& apt-get install -y libgdal-dev \
#&& apt-get install -y unzip \
# Installing libssl1 sometimes needed for Dotnet
&& wget http://archive.ubuntu.com/ubuntu/pool/main/o/openssl/libssl1.1_1.1.0g-2ubuntu4_amd64.deb \
&& dpkg -i libssl1.1_1.1.0g-2ubuntu4_amd64.deb \
&& apt-get clean && rm -rf /var/lib/apt/lists/*

ENV C_INCLUDE_PATH=/usr/include/gdal
ENV CPLUS_INCLUDE_PATH=/usr/include/gdal

# INSTALLING DOTNET SDK AND RUNTIME 8.0
RUN mkdir /bin/.dotnet/ && cd /bin/.dotnet/ \
&& wget https://dot.net/v1/dotnet-install.sh -O dotnet-install.sh \
&& chmod +x ./dotnet-install.sh \
&& ./dotnet-install.sh --channel 8.0 -InstallDir /bin/.dotnet/  \
&& ./dotnet-install.sh --channel 8.0 --runtime aspnetcore -InstallDir /bin/.dotnet/ \
&& apt-get clean && rm -rf /var/lib/apt/lists/*
ENV DOTNET_ROOT=/bin/.dotnet
ENV PATH=$PATH:$DOTNET_ROOT:$DOTNET_ROOT/tools

# PREPARING VARIABLES FOR DOWNLOADS
# WARNING : With a bad internet connection, git clone can sometimes fail too.
# These commands should help git on bad internet connections
RUN git config --global http.version HTTP/1.1 \
&& git config --global http.postBuffer 524288000 \
&& git config --global http.lowSpeedLimit 0 \
&& git config --global http.lowSpeedTime 999999

# PREPARING VARIABLE FOR MORE EASY COMMANDS 
ENV LANDIS_EXTENSIONS_TOOL="/bin/LANDIS_Linux/Core-Model-v8-LINUX/build/Release/Landis.Extensions.dll"
ENV LANDIS_CONSOLE="/bin/LANDIS_Linux/build/Release/Landis.Console.dll"
ENV LANDIS_FOLDER="/bin/LANDIS_Linux"

################################################################
# COMPILING LANDIS-II
################################################################

###### COMPILING THE CORE 
RUN mkdir /bin/LANDIS_Linux \
&& cd /bin/LANDIS_Linux && git clone https://github.com/LANDIS-II-Foundation/Core-Model-v8-LINUX.git \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Tool-Console/src && dotnet build -c Release 
# && mkdir /bin/LANDIS_Linux/build && mv /bin/LANDIS_Linux/Core-Model-v8-LINUX/build/* /bin/LANDIS_Linux/build 
# && rm -r /bin/LANDIS_Linux/Core-Model-v8-LINUX

# RUN mkdir /bin/LANDIS_Linux \
# && cd /bin/LANDIS_Linux && git clone https://github.com/tgestabrook/Core-Model-v8-LINUX.git \
# && cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Tool-Console/src && dotnet build -c Release \
# && mkdir /bin/LANDIS_Linux/build && mv /bin/LANDIS_Linux/Core-Model-v8-LINUX/build/* /bin/LANDIS_Linux/build \
# && rm -r /bin/LANDIS_Linux/Core-Model-v8-LINUX

# Downloading support libraries and installing them
RUN cd /bin/LANDIS_Linux && git clone https://github.com/LANDIS-II-Foundation/Support-Library-Dlls-v8.git \
&& mv /bin/LANDIS_Linux/Support-Library-Dlls-v8/* /bin/LANDIS_Linux/Core-Model-v8-LINUX/build/extensions \
&& rm -r /bin/LANDIS_Linux/Support-Library-Dlls-v8

# RUN mkdir /bin/LANDIS_Linux \
# && cd /bin/LANDIS_Linux && git clone https://github.com/tgestabrook/Core-Model-v8-LINUX.git \
# && git clone https://github.com/LANDIS-II-Foundation/Support-Library-Dlls-v8.git \
# && mv /bin/LANDIS_Linux/Support-Library-Dlls-v8/* /bin/LANDIS_Linux/build/extensions \
# && rm -r /bin/LANDIS_Linux/Support-Library-Dlls-v8

# Transfering the python script that helps the compilation of the extensions
COPY ./files_to_help_compilation/editing_csproj_LANDIS-II_files.py /bin/LANDIS_Linux/Core-Model-v8-LINUX

# Transfering the sh script that allows us to only download one folder from one commit from github
COPY ./files_to_help_compilation/downloadSpecificGitCommitAndFolder.sh /bin/LANDIS_Linux/Core-Model-v8-LINUX
RUN chmod +x /bin/LANDIS_Linux/Core-Model-v8-LINUX/downloadSpecificGitCommitAndFolder.sh

###### COMPILING AND REGISTERING EXTENSIONS 
# INFO: The files necessary for the compilation are downloaded through the script downloadSpecificGitCommitAndFolder.sh
# It is used to avoid downloading the entire repository as it's not needed.
#
# It takes three arguments :
# The URL of the repo
# The commit hash for the commit from which we want to download files
# The folder we want to download files from (excludes all of the others)
#
# The folder we want is almost always /src which contains the source code.
#
# We also always end up deleting the files downloaded for the build to avoid saving them in docker layering system (saves space).
#
# For the .txt file necessary to register the extension in the extension.xml file, I've decided to download them through
# their unique Github URL (corresponding to the commit). This is again to avoid downloading a lot of files for nothing.

#### SUCCESSION EXTENSIONS

# Biomass Succession
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit 58ad3673e02abe82f437a6b68c44220c51351091 is latest at time of writing.
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-Biomass-Succession.git 58ad3673e02abe82f437a6b68c44220c51351091 /src \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-Biomass-Succession/src \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./biomass-succession.csproj \
&& dotnet build -c Release \
&& wget -O Biomass_Succession_7.txt https://raw.githubusercontent.com/LANDIS-II-Foundation/Extension-Biomass-Succession/58ad3673e02abe82f437a6b68c44220c51351091/deploy/installer/Biomass%20Succession%207.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "Biomass_Succession_7.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-Biomass-Succession 

# NECN Succession
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit 37ce246c37bab3448e3db134373deb56063e14ac is latest at time of writing.
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-NECN-Succession.git 37ce246c37bab3448e3db134373deb56063e14ac /src \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-NECN-Succession/src \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./NECN-succession.csproj \
# The devellopers of NECN forgot to add a line referencing the location of one of the support libraries.
# I add it here manually at the right place in the file.
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-NECN-Succession/src && sed -i '39i    <HintPath>..\\..\\build\\extensions\\Landis.Library.Climate-v5.dll</HintPath>' NECN-succession.csproj \
&& dotnet build -c Release \
&& wget https://raw.githubusercontent.com/LANDIS-II-Foundation/Extension-NECN-Succession/37ce246c37bab3448e3db134373deb56063e14ac/deploy/installer/NECN_Succession8.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "NECN_Succession8.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-NECN-Succession

# # PnET Succession
# Currently left until https://github.com/LANDIS-II-Foundation/Extension-PnET-Succession/issues/17 is fixed.
# RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# # Commit 13ffece64f0b0180490aaa96cf91bd212e9e2e50 is the most recent at date of writing..
# && ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-PnET-Succession.git 13ffece64f0b0180490aaa96cf91bd212e9e2e50 /src \ 
# && cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-PnET-Succession/src \
# && python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./PnET-Succession.csproj \
# # A definition seems to be missing in Plugin.cs, even on newer versions; I add it here. Strange ?
# # && cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-PnET-Succession/src && sed -i '57i        public int ThreadCount;' PlugIn.cs \
# # We need to remove the .sln file; it's not useful to build, and is full of bad/missing references to other projects.
# && rm /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-PnET-Succession/src/PnET-Succession.sln
# RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-PnET-Succession/src \
# && dotnet build -c Release \
# && cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-PnET-Succession \
# # To finish, we have to transfer some "Default files" of PnET in the extension folder
# && ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-PnET-Succession.git 13ffece64f0b0180490aaa96cf91bd212e9e2e50 /deploy \ 
# && cp -r /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-PnET-Succession/deploy/Defaults /bin/LANDIS_Linux/Core-Model-v8-LINUX/build/extensions \
# && wget https://raw.githubusercontent.com/LANDIS-II-Foundation/Extension-PnET-Succession/13ffece64f0b0180490aaa96cf91bd212e9e2e50/deploy/installer/PnET-Succession.txt \
# && dotnet $LANDIS_EXTENSIONS_TOOL add "PnET-Succession.txt" \
# && cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-PnET-Succession

# # # Forest Carbon Succession (ForCS)
# Not yet ready for corev8 as of writing. Will have to wait.

#GIPL
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Library-GIPL.git master /src \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Library-GIPL/src \
#&& sed -i -e 's/lib\\Landis.Library/..\\..\\Support-Library-Dlls-v8-master\\Landis.Library/g' GiplDamm.csproj \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./GiplDamm.csproj \
&& dotnet build -c Release GiplDamm.csproj \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Library-GIPL

#SHAW
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Library-Shaw.git master /src \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Library-Shaw/src \
#&& sed -i -e 's/lib\\Landis.Library/..\\..\\Support-Library-Dlls-v8-master\\Landis.Library/g' ShawDammNewInputs.csproj \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./ShawDammNewInputs.csproj \
&& dotnet build -c Release ShawDammNewInputs.csproj \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Library-Shaw

# # # DGS Succession extension : Not yet ready for v8.
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# just grab the latest for v8
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-DGS-Succession.git v8 /src \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-DGS-Succession/src \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./DGS-succession.csproj \
&& dotnet build -c Release \
&& wget -O DGS_Succession_2.4.txt https://raw.githubusercontent.com/LANDIS-II-Foundation/Extension-DGS-Succession/v8/deploy/installer/DGS_Succession%202.4.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "DGS_Succession_2.4.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-DGS-Succession 

#Update hint paths for GIPL and SHAW

#### DISTURBANCE EXTENSIONS

## FIRE

# Base Fire
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit 40169a04b62be126e109c7990fa6a6d64a8ff895 is latest as of writing this.
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-Base-Fire.git 40169a04b62be126e109c7990fa6a6d64a8ff895 /src \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-Base-Fire/src \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./original-fire.csproj \
# We remove the .sln as it's not useful for the build
&& rm original-fire.sln \
# We got to add two lines to indicate where the support libraries are
&& sed -i '31i    <HintPath>..\\..\\build\\extensions\\Landis.Library.UniversalCohorts-v1.dll</HintPath>' original-fire.csproj \
&& sed -i '37i    <HintPath>..\\..\\build\\extensions\\Landis.Library.Parameters-v2.dll</HintPath>' original-fire.csproj \
&& dotnet build -c Release \
&& wget https://raw.githubusercontent.com/LANDIS-II-Foundation/Extension-Base-Fire/40169a04b62be126e109c7990fa6a6d64a8ff895/deploy/installer/Original%20Fire%205.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "Original Fire 5.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-Base-Fire

# # BFOLDS Fire - Not ready yet for v8, to check with Marc Ouelette.

# Climate-social Fire
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit b463ea378f1bcde4369907a408dfe64b9cc52c7a is latest as of writing this.
#try master?
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-Social-Climate-Fire.git master /src \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-Social-Climate-Fire/src \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./SocialClimateFire.csproj \
# We got to add two lines to indicate where the support libraries are
&& sed -i '39i    <HintPath>..\\..\\build\\extensions\\Landis.Library.Climate-v5.dll</HintPath>' SocialClimateFire.csproj \
&& dotnet build -c Release \
&& wget -O Scrapple_4.txt https://raw.githubusercontent.com/LANDIS-II-Foundation/Extension-Social-Climate-Fire/master/deploy/installer/Scrapple%204.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "Scrapple_4.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-Social-Climate-Fire

# # Dynamic Fuel system - Not referenced anymore in the LANDIS-II extension, but is only refered as "Dynamic Fuels & Fire System" with a link to 	http://landis-ii-foundation.github.io/Extension-Dynamic-Fire-System/. Will let this one be. Seems to have been replaced by dynamic biomass fuels ?

# # Dynamic Biomass Fuels
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit 06dd67482b20a74a0e075782e66b09b4fe42c248 is latest as of writing this.
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-Dynamic-Biomass-Fuels.git 06dd67482b20a74a0e075782e66b09b4fe42c248 /src \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-Dynamic-Biomass-Fuels/src \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./dynamic-fuels.csproj \
&& dotnet build -c Release \
&& wget https://raw.githubusercontent.com/LANDIS-II-Foundation/Extension-Dynamic-Biomass-Fuels/06dd67482b20a74a0e075782e66b09b4fe42c248/deploy/installer/Dynamic%20Fuels%204.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "Dynamic Fuels 4.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-Dynamic-Biomass-Fuels

# Dynamic Fire System
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit 4970f846bd2b22f014f201ba2af2278436aefd7c is latest as of writing this.
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-Dynamic-Fire-System.git 4970f846bd2b22f014f201ba2af2278436aefd7c /src \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-Dynamic-Fire-System/src \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./dynamic-fire.csproj \
# We got to add two lines to indicate where the support libraries are
&& sed -i '37i    <HintPath>..\\..\\build\\extensions\\Landis.Library.UniversalCohorts-v1.dll</HintPath>' dynamic-fire.csproj \
&& sed -i '40i    <HintPath>..\\..\\build\\extensions\\Landis.Library.Parameters-v2.dll</HintPath>' dynamic-fire.csproj \
&& dotnet build -c Release \
&& wget https://raw.githubusercontent.com/LANDIS-II-Foundation/Extension-Dynamic-Fire-System/4970f846bd2b22f014f201ba2af2278436aefd7c/deploy/installer/Dynamic%20Fire%20Component%204.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "Dynamic Fire Component 4.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-Dynamic-Fire-System

# ## HARVEST

# Biomass Harvest
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit 24b01fea5a90b05b2732c3e52e09a02fdb47db59 is latest as of writing this.
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-Biomass-Harvest.git 24b01fea5a90b05b2732c3e52e09a02fdb47db59 /src \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-Biomass-Harvest/src \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./biomass-harvest-ext.csproj \
&& dotnet build -c Release \
&& wget https://raw.githubusercontent.com/LANDIS-II-Foundation/Extension-Biomass-Harvest/24b01fea5a90b05b2732c3e52e09a02fdb47db59/deploy/installer/Biomass%20Harvest%206.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "Biomass Harvest 6.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-Biomass-Harvest

# # SOSIEL Harvest - Not ready for v8 as time of writing

# FRS Module
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit b2f1345f33bf0d089dabba75f8a1afe5de2d9d26 is latest as of writing this.
# Not using the download script here as the repo is organized differently (source code is in the root)
&& git clone https://github.com/Klemet/LANDIS-II-Forest-Roads-Simulation-extension.git \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/LANDIS-II-Forest-Roads-Simulation-extension && git checkout b2f1345f33bf0d089dabba75f8a1afe5de2d9d26 \
&& rm /bin/LANDIS_Linux/Core-Model-v8-LINUX/LANDIS-II-Forest-Roads-Simulation-extension/Forest-Roads-Extension.sln
# Sadly, the csproj file here is completly wrong for compiling on Linux. 
# I got to replace it entirely by something that is more in line with the other extensions.
COPY files_to_help_compilation/Forest-Roads-Extension.csproj /bin/LANDIS_Linux/Core-Model-v8-LINUX/LANDIS-II-Forest-Roads-Simulation-extension
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/LANDIS-II-Forest-Roads-Simulation-extension && dotnet build -c Release \
&& wget https://raw.githubusercontent.com/Klemet/LANDIS-II-Forest-Roads-Simulation-extension/b2f1345f33bf0d089dabba75f8a1afe5de2d9d26/Deploy/Installation%20Files/plug-ins-installer-files/Forest%20Roads%20Simulation%202.0.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "Forest Roads Simulation 2.0.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r LANDIS-II-Forest-Roads-Simulation-extension

# # Magic harvest
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit 4def8fc98e3d816367f8ba2ec6d2b8ea76d279f6 is latest as of writing this.
# Not using the download script here as the repo is organized differently (source code is in the root)
&& git clone https://github.com/Klemet/LANDIS-II-Magic-Harvest.git \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/LANDIS-II-Magic-Harvest/ && git checkout 4def8fc98e3d816367f8ba2ec6d2b8ea76d279f6 \
&& rm /bin/LANDIS_Linux/Core-Model-v8-LINUX/LANDIS-II-Magic-Harvest/"Magic Harvest.sln"
# Sadly, the csproj file here is completly wrong for compiling on Linux. 
# I got to replace it entirely by something that is more in line with the other extensions.
COPY ["./files_to_help_compilation/Magic Harvest.csproj", "/bin/LANDIS_Linux/Core-Model-v8-LINUX/LANDIS-II-Magic-Harvest"]  
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/LANDIS-II-Magic-Harvest && dotnet build -c Release \
&& wget https://raw.githubusercontent.com/Klemet/LANDIS-II-Magic-Harvest/4def8fc98e3d816367f8ba2ec6d2b8ea76d279f6/Deploy/Installation%20Files/plug-ins-installer-files/Magic%20harvest%20v2.1.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "Magic harvest v2.1.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r LANDIS-II-Magic-Harvest


## WIND

# Base Wind - Now called Original wind
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit e0796d577a6221c78b10502da4fc7dd0ddfe08e9 is latest as of writing this.
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-Base-Wind.git e0796d577a6221c78b10502da4fc7dd0ddfe08e9 /src \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-Base-Wind/src \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./original-wind.csproj \
# We need to add a line to the csproj.
&& sed -i '35i    <HintPath>..\\..\\build\\extensions\\Landis.Library.UniversalCohorts-v1.dll</HintPath>' original-wind.csproj \
&& dotnet build -c Release \
&& wget https://raw.githubusercontent.com/LANDIS-II-Foundation/Extension-Base-Wind/e0796d577a6221c78b10502da4fc7dd0ddfe08e9/deploy/installer/Original%20Wind%204.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "Original Wind 4.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-Base-Wind

# Linear Wind
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit b8efe5ca20ca386fc978db97670e50424941153e is latest as of writing this.
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-LinearWind.git b8efe5ca20ca386fc978db97670e50424941153e /src \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-LinearWind/src \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./linear-wind.csproj \
# We need to add a line to the csproj.
&& sed -i '37i    <HintPath>..\\..\\build\\extensions\\Landis.Library.UniversalCohorts-v1.dll</HintPath>' linear-wind.csproj \
&& dotnet build -c Release \
&& wget https://raw.githubusercontent.com/LANDIS-II-Foundation/Extension-LinearWind/b8efe5ca20ca386fc978db97670e50424941153e/deploy/installer/Linear%20Wind%203.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "Linear Wind 3.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-LinearWind

# # Hurricane
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit a12806d77d4b251d8800766f124c39adf90541be is latest as of writing this.
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-Biomass-Hurricane.git a12806d77d4b251d8800766f124c39adf90541be /src \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-Biomass-Hurricane/src \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./hurricane-extension.csproj \
# We need to add a line to the csproj.
&& sed -i '31i    <HintPath>..\\..\\build\\extensions\\Landis.Library.UniversalCohorts-v1.dll</HintPath>' hurricane-extension.csproj \
&& sed -i '34i    <HintPath>..\\..\\build\\extensions\\Landis.Library.Metadata-v2.dll</HintPath>' hurricane-extension.csproj \
&& dotnet build -c Release \
&& wget https://raw.githubusercontent.com/LANDIS-II-Foundation/Extension-Biomass-Hurricane/a12806d77d4b251d8800766f124c39adf90541be/deploy/current/Hurricane%203.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "Hurricane 3.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-Biomass-Hurricane


## DISEASES, INSECTS AND OTHERS

# Base BDA (now named Climate biological Disturbance agent)
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit eb1d998a14b7555ddd7c527dda797669b0c99546 is latest as of writing this.
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-Base-BDA.git eb1d998a14b7555ddd7c527dda797669b0c99546 /src \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-Base-BDA/src \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./BDA-Climate.csproj \
# We need to add a line to the csproj.
&& sed -i '38i    <HintPath>..\\..\\build\\extensions\\Landis.Library.UniversalCohorts-v1.dll</HintPath>' BDA-Climate.csproj \
&& sed -i '41i    <HintPath>..\\..\\build\\extensions\\Landis.Library.Climate-v5.dll</HintPath>' BDA-Climate.csproj \
&& sed -i '44i    <HintPath>..\\..\\build\\extensions\\Landis.Library.Metadata-v2.dll</HintPath>' BDA-Climate.csproj \
&& dotnet build -c Release \
&& wget https://raw.githubusercontent.com/LANDIS-II-Foundation/Extension-Base-BDA/eb1d998a14b7555ddd7c527dda797669b0c99546/deploy/installer/Climate%20BDA%205.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "Climate BDA 5.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-Base-BDA

# Root rot - Not yet available for v8

# # Biomass Browse - Not yet available for v8

# # Land-use plus
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit 574940aa6382ed9e5840b78b0544300bf5a40cd2 is latest as of writing this.
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-Land-Use-Plus.git 574940aa6382ed9e5840b78b0544300bf5a40cd2 /src \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-Land-Use-Plus/src \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./land-use.csproj \
&& dotnet build -c Release \
&& wget https://raw.githubusercontent.com/LANDIS-II-Foundation/Extension-Land-Use-Plus/574940aa6382ed9e5840b78b0544300bf5a40cd2/deploy/installer/Land%20Use%204.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "Land Use 4.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-Land-Use-Plus


### OUTPUT EXTENSIONS

# Max species Age
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit bba5b5a4879d0d6cbfbcce4867702bd4df3ac350 is latest as of writing this.
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-Output-Max-Species-Age.git bba5b5a4879d0d6cbfbcce4867702bd4df3ac350 /src \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-Output-Max-Species-Age/src \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./max-species-age.csproj \
# We need to add a line to the csproj.
&& sed -i '37i    <HintPath>..\\..\\build\\extensions\\Landis.Library.UniversalCohorts-v1.dll</HintPath>' max-species-age.csproj \
&& dotnet build -c Release \
&& wget https://raw.githubusercontent.com/LANDIS-II-Foundation/Extension-Output-Max-Species-Age/bba5b5a4879d0d6cbfbcce4867702bd4df3ac350/deploy/installer/Output%20MaxSpeciesAge%204.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "Output MaxSpeciesAge 4.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-Output-Max-Species-Age

# Biomass-by-age
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit 0419cc64634f57ad3660590408ded3aef88ecf9d is latest as of writing this.
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-Output-Biomass-By-Age.git 0419cc64634f57ad3660590408ded3aef88ecf9d /src \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-Output-Biomass-By-Age/src \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./output-biomass-by-age.csproj \
# We need to add a line to the csproj.
&& sed -i '37i    <HintPath>..\\..\\build\\extensions\\Landis.Library.UniversalCohorts-v1.dll</HintPath>' output-biomass-by-age.csproj \
&& dotnet build -c Release \
&& wget https://raw.githubusercontent.com/LANDIS-II-Foundation/Extension-Output-Biomass-By-Age/0419cc64634f57ad3660590408ded3aef88ecf9d/deploy/installer/Output%20Biomass-by-Age%204.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "Output Biomass-by-Age 4.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-Output-Biomass-By-Age

# Biomass Community
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit 58252f441cc393cc1e63ea6c36175e15bba93916 is latest as of writing this.
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-Output-Biomass-Community.git 58252f441cc393cc1e63ea6c36175e15bba93916 /src \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-Output-Biomass-Community/src \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./output-biomass-community.csproj \
# We need to add a line to the csproj.
&& sed -i '36i    <HintPath>..\\..\\build\\extensions\\Landis.Library.UniversalCohorts-v1.dll</HintPath>' output-biomass-community.csproj \
&& dotnet build -c Release \
&& wget https://raw.githubusercontent.com/LANDIS-II-Foundation/Extension-Output-Biomass-Community/58252f441cc393cc1e63ea6c36175e15bba93916/deploy/installer/Output%20Biomass%20Community%203.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "Output Biomass Community 3.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-Output-Biomass-Community

# Biomass Output
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit d5cb256f7669df36a76d9337c779cdc7f1cdbd0b is latest as of writing this.
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-Output-Biomass.git d5cb256f7669df36a76d9337c779cdc7f1cdbd0b /src \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-Output-Biomass/src \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./output-biomass.csproj \
&& dotnet build -c Release \
&& wget https://raw.githubusercontent.com/LANDIS-II-Foundation/Extension-Output-Biomass/d5cb256f7669df36a76d9337c779cdc7f1cdbd0b/deploy/installer/Output%20Biomass%204.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "Output Biomass 4.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-Output-Biomass

# Biomass reclassification
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit fad7e9f7e39b9cf72e1e55210cb0e8cd09082671 is latest as of writing this.
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-Output-Biomass-Reclass.git fad7e9f7e39b9cf72e1e55210cb0e8cd09082671 /src \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-Output-Biomass-Reclass/src \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./output-biomass.csproj \
&& dotnet build -c Release \
&& wget https://raw.githubusercontent.com/LANDIS-II-Foundation/Extension-Output-Biomass-Reclass/fad7e9f7e39b9cf72e1e55210cb0e8cd09082671/deploy/installer/Output%20Biomass%20Reclass%204.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "Output Biomass Reclass 4.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-Output-Biomass-Reclass

# # Cohort statistics outputs
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit 045272850c77b8b5e8c36ba1fe8c5041b7a523c2 is latest as of writing this.
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-Output-Cohort-Statistics.git 045272850c77b8b5e8c36ba1fe8c5041b7a523c2 /src \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-Output-Cohort-Statistics/src \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./output-cohort-stats.csproj \
# We need to add a line to the csproj.
&& sed -i '37i    <HintPath>..\\..\\build\\extensions\\Landis.Library.UniversalCohorts-v1.dll</HintPath>' output-cohort-stats.csproj \
&& dotnet build -c Release \
&& wget https://raw.githubusercontent.com/LANDIS-II-Foundation/Extension-Output-Cohort-Statistics/045272850c77b8b5e8c36ba1fe8c5041b7a523c2/deploy/installer/Output%20Cohort%20Statistics%204.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "Output Cohort Statistics 4.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-Output-Cohort-Statistics

# # Landscape habitat Output - Not yet available for v8

# Local habitat suitability
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit 1366a092625e0a26fff870e16529c1fe3e071c14 is latest as of writing this.
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-Local-Habitat-Suitability-Output.git 1366a092625e0a26fff870e16529c1fe3e071c14 /src \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-Local-Habitat-Suitability-Output/src \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./local-habitat-output.csproj \
&& dotnet build -c Release \
&& wget https://raw.githubusercontent.com/LANDIS-II-Foundation/Extension-Local-Habitat-Suitability-Output/1366a092625e0a26fff870e16529c1fe3e071c14/deploy/installer/Local%20Habitat%20Output.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "Local Habitat Output.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-Local-Habitat-Suitability-Output

# # PnET Output - Not yet available for v8

# Wildlife habitat Output
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit 695a03ba11a21d8a12eb714a6c8759a8284290f2 is latest as of writing this.
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-Output-Wildlife-Habitat.git 695a03ba11a21d8a12eb714a6c8759a8284290f2 /src \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-Output-Wildlife-Habitat/src \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./wildlife-habitat.csproj \
&& dotnet build -c Release \
&& wget https://raw.githubusercontent.com/LANDIS-II-Foundation/Extension-Output-Wildlife-Habitat/695a03ba11a21d8a12eb714a6c8759a8284290f2/deploy/installer/Wildlife%20Habitat%20Output%203.txt \
&& dotnet $LANDIS_EXTENSIONS_TOOL add "Wildlife Habitat Output 3.txt" \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Extension-Output-Wildlife-Habitat

# ################################################################
# # RECOMPILING SOME THINGS TO SOLVE METADATALIBRARY ERROR
# ################################################################

# # Metadata library
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# Commit 72b8caa14cdd6af81c6e1e4541c4c3e18bc63eca is latest at time of writing.
&& ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Library-Metadata.git 72b8caa14cdd6af81c6e1e4541c4c3e18bc63eca /src \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Library-Metadata/src \
&& python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./Metadata.csproj \
&& dotnet build -c Release \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Library-Metadata

# Universal cohorts library // needs v2??
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
&& rm /bin/LANDIS_Linux/Core-Model-v8-LINUX/build/extensions/Landis.Library.UniversalCohorts-v1.dll \
# Commit a1e84adc8073fb7c89ee078a38349b5d578d4179 is latest at time of writing.
&& mkdir /bin/LANDIS_Linux/Core-Model-v8-LINUX/Library-Universal-Cohort/ && cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Library-Universal-Cohort/ \
&& git clone https://github.com/LANDIS-II-Foundation/Library-Universal-Cohort.git \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Library-Universal-Cohort/Library-Universal-Cohort/ && git checkout a1e84adc8073fb7c89ee078a38349b5d578d4179 \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Library-Universal-Cohort/Library-Universal-Cohort/ && python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./universal-cohort.csproj \
&& dotnet build -c Release \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Library-Universal-Cohort

# Universal Initial community library
RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
&& rm /bin/LANDIS_Linux/Core-Model-v8-LINUX/build/extensions/Landis.Library.InitialCommunity.Universal.dll \
# Commit 5dc6dd299eef88ded1c88871470d58c26c1a4093 is latest at time of writing.
&& mkdir /bin/LANDIS_Linux/Core-Model-v8-LINUX/Library-Initial-Community/ && cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Library-Initial-Community/ \
&& git clone https://github.com/LANDIS-II-Foundation/Library-Initial-Community.git \ 
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Library-Initial-Community/Library-Initial-Community/ && git checkout 5dc6dd299eef88ded1c88871470d58c26c1a4093 \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Library-Initial-Community/Library-Initial-Community/ && python /bin/LANDIS_Linux/Core-Model-v8-LINUX/editing_csproj_LANDIS-II_files.py ./initial-community.csproj \
&& dotnet build -c Release \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/ && rm -r Library-Initial-Community

# Recompiling the core as per https://github.com/CU-ESIIL/ExtremeWildfire/blob/main/docker/landis2/Dockerfile
RUN sed -i '/<\/Project>/i \
<ItemGroup>\n\
    <Reference Include="Landis.Extension.Succession.Biomass-v7.dll">\n\
        <HintPath>../../build/extensions/Landis.Extension.Succession.Biomass-v7.dll</HintPath>\n\
    </Reference>\n\
    <Reference Include="Landis.Library.Metadata-v2">\n\
        <HintPath>../../build/extensions/Landis.Library.Metadata-v2.dll</HintPath>\n\
    </Reference>\n\
    <Reference Include="Landis.Library.Succession-v9.dll">\n\
        <HintPath>../../build/extensions/Landis.Library.Succession-v9.dll</HintPath>\n\
    </Reference>\n\
    <Reference Include="Landis.Extension.Succession.NECN-v8.dll">\n\
        <HintPath>../../build/extensions/Landis.Extension.Succession.NECN-v8.dll</HintPath>\n\
    </Reference>\n\
    <Reference Include="Landis.Extension.GiplDamm.dll">\n\
        <HintPath>../../build/extensions/Landis.Extension.GiplDamm.dll</HintPath>\n\
    </Reference>\n\
        <Reference Include="Landis.Extension.ShawDamm.dll">\n\
        <HintPath>../../build/extensions/Landis.Extension.ShawDamm.dll</HintPath>\n\
    </Reference>\n\
    <Reference Include="Landis.Extension.Succession.DGS.dll">\n\
        <HintPath>../../build/extensions/Landis.Extension.Succession.DGS.dll</HintPath>\n\
    </Reference>\n\
    <Reference Include="Landis.Extension.SocialClimateFire-v4.dll">\n\
        <HintPath>../../build/extensions/Landis.Extension.SocialClimateFire-v4.dll</HintPath>\n\
    </Reference>\n\
    <Reference Include="Landis.Extension.Output.Biomass-v4.dll">\n\
        <HintPath>../../build/extensions/Landis.Extension.Output.Biomass-v4.dll</HintPath>\n\
    </Reference>\n\
    <Reference Include="Landis.Extension.Output.BiomassCommunity-v3.dll">\n\
        <HintPath>../../build/extensions/Landis.Extension.Output.BiomassCommunity-v3.dll</HintPath>\n\
    </Reference>\n\
    <Reference Include="Landis.Extension.Output.BiomassReclass-v4.dll">\n\
        <HintPath>../../build/extensions/Landis.Extension.Output.BiomassReclass-v4.dll</HintPath>\n\
    </Reference>\n\
    <Reference Include="Landis.Library.BiomassHarvest-v3.dll">\n\
        <HintPath>../../build/extensions/Landis.Library.BiomassHarvest-v3.dll</HintPath>\n\
    </Reference>\n\
    <Reference Include="Landis.Library.Climate-v5.dll">\n\
        <HintPath>../../build/extensions/Landis.Library.Climate-v5.dll</HintPath>\n\
    </Reference>\n\
    <Reference Include="Landis.Library.HarvestManagement-v4.dll">\n\
        <HintPath>../../build/extensions/Landis.Library.HarvestManagement-v4.dll</HintPath>\n\
    </Reference>\n\
    <Reference Include="Landis.Library.HarvestManagement-v4.dll">\n\
        <HintPath>../../build/extensions/Landis.Library.HarvestManagement-v4.dll</HintPath>\n\
    </Reference>\n\
    <Reference Include="log4net.dll">\n\
        <HintPath>../../build/extensions/log4net.dll</HintPath>\n\
    </Reference>\n\
    <Reference Include="Landis.Library.UniversalCohorts-v1.dll">\n\
        <HintPath>../../build/extensions/Landis.Library.UniversalCohorts-v1.dll</HintPath>\n\
    </Reference>\n\
</ItemGroup>' /bin/LANDIS_Linux/Core-Model-v8-LINUX/Tool-Console/src/Console.csproj \
# && cp /bin/LANDIS_Linux/build/extensions/extensions.xml /bin/LANDIS_Linux \
&& cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Tool-Console/src && dotnet build -c Release \
&& mkdir /bin/LANDIS_Linux/build && cp -r /bin/LANDIS_Linux/Core-Model-v8-LINUX/build/* /bin/LANDIS_Linux/build 
# && cp /bin/LANDIS_Linux/extensions.xml /bin/LANDIS_Linux/build/extensions/ \
# && rm -r /bin/LANDIS_Linux/Core-Model-v8-LINUX


#something for sudo
RUN adduser --disabled-password --gecos "VICE_User" --uid 1000 user  && \
    usermod -aG sudo user && \
    echo "$LOCAL_USER ALL=NOPASSWD: $PRIV_CMDS" >> /etc/sudoers

# ################################################################
# # TESTING THE EXTENSIONS
# # (Optional; you can comment this whole section).
# ################################################################


# # We create a temp file that is used as a time reference to delete the output files created by the tests,
# # in case we don't want to keep them.
# # RUN touch /tmp/tempFileForTestFileRemoval.txt
# # ENV temp_file="/tmp/tempFileForTestFileRemoval.txt"

# # # The version here is 6.1; and there are test files for v6.1. But these return an error where launched.
# # # Seems like it's the precision of the biomass of a cohort in the initial community files that is the issue.
# # # The test files for v6 works perfectly.

# # Biomass Succession
# RUN cd /bin/LANDIS_Linux/Core-Model-v8-LINUX \
# # Commit 58ad3673e02abe82f437a6b68c44220c51351091 is latest at time of writing.
# && ./downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-Biomass-Succession.git 58ad3673e02abe82f437a6b68c44220c51351091 /testings \
# && cd /bin/LANDIS_Linux/Core-Model-v8-LINUX/Extension-Biomass-Succession/testings/CoreV8.0-BiomassSuccession7.0 \
# # We have to rename the initial communities file of the testing folder because LANDIS-II of Linux is case sensitive, and there is a mismatch between how its written in the parameter file and its name.
# && mv "initial-communities.tif" "Initial-Communities.tif" \
# && dotnet $LANDIS_CONSOLE scenario.txt

# # # RUN cd $LANDIS_FOLDER/Extension-Biomass-Succession/testings/CoreV7.0-BiomassSuccession6.0 \
# # # && dotnet $LANDIS_CONSOLE scenario.txt
# RUN cd /bin/LANDIS_Linux \
# && /bin/LANDIS_Linux/Core-Model-v8-LINUX/downloadSpecificGitCommitAndFolder.sh https://github.com/LANDIS-II-Foundation/Extension-NECN-Succession.git 37ce246c37bab3448e3db134373deb56063e14ac /testing \
# && cd $LANDIS_FOLDER/Extension-NECN-Succession/testing/Core8-NECN8-Landscape \
# && dotnet $LANDIS_CONSOLE scenario_LTB_landscape.txt

# # We delete all of the files generating during testing (more recent than the temporary file we created before.
# # Comment this line to keep the files generated during testing.
# # RUN find $LANDIS_FOLDER//Extension-Age-Only-Succession-Archived/testing/Corev7 -type f -newer "/tmp/tempFileForTestFileRemoval.txt" -exec rm {} \;

# ################################################################
# # FINISHING
# ################################################################

# # Re-configure git for latest version of HTTP protocol
RUN git config --global --unset http.version

# install ttyd
RUN wget "https://github.com/tsl0922/ttyd/releases/download/1.6.3/ttyd.x86_64" > ttyd \
&& chmod a+x ttyd \
&& mv ttyd /usr/local/bin/ttyd

RUN apt-get update \
&& apt-get install -y curl grep sed dpkg \
&& curl -L "https://github.com/krallin/tini/releases/download/v0.19.0/tini_0.19.0-amd64.deb" > tini.deb \
&& dpkg -i tini.deb \
&& rm tini.deb \
&& apt-get clean \
&& rm -rf /var/lib/apt/lists/*

# Install a few dependencies for iCommands, text editing, and monitoring instances
RUN apt-get update \
&& apt-get install -y lsb-release apt-transport-https curl gnupg2 libfuse2 gcc less nodejs software-properties-common apt-utils glances htop nano \
&& apt-get clean \
&& rm -rf /var/lib/apt/lists/*

RUN apt-get update \
&& apt-get install -y vim-nox \
&& apt-get clean \
&& rm -rf /var/lib/apt/lists/*

RUN apt-get update \
&& apt-get install -y emacs-nox \
&& apt-get clean \
&& rm -rf /var/lib/apt/lists/*

RUN wget -qO - https://packages.irods.org/irods-signing-key.asc | apt-key add - \
&& echo "deb [arch=amd64] https://packages.irods.org/apt/ $(lsb_release -sc) main" >> /etc/apt/sources.list.d/renci-irods.list \
&& apt-get update \
&& apt install -y irods-icommands 

USER root

# set shell as bash and terminal as linux
ENV SHELL=bash
ENV TERM=xterm
EXPOSE 7681

# changes tmux layout while running
COPY entry.sh /bin

RUN chmod +x /bin/entry.sh

RUN echo 'set-option -g status off' >> ~/.tmux.conf

# add iRODS iCommands to user profile as JSON
RUN mkdir .irods 

ENTRYPOINT ["/bin/bash"]
#CMD ["ttyd", "bash"]