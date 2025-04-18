FROM ubuntu:20.04

USER root

# Set non-interactive status for builds
ENV DEBIAN_FRONTEND=noninteractive

# Install basic dependencies including sudo, ca-certificates, and gnupg for key management
RUN apt-get update && \
    apt-get install -y sudo apt-transport-https \
    ca-certificates gnupg curl software-properties-common && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Add GitHub CLI GPG key and repository
RUN curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg | dd of=/etc/apt/trusted.gpg.d/githubcli-archive-keyring.gpg && \
    echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/trusted.gpg.d/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" > /etc/apt/sources.list.d/github-cli.list && \
    apt-get update && \
    apt-get install -y gh && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install additional dependencies
RUN apt-get update && \
    apt-get install -y git \
    libjpeg62 \
    libpng16-16 \
    wget \
    unzip && \
    add-apt-repository universe && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /opt

RUN wget -q https://packages.microsoft.com/config/ubuntu/20.04/packages-microsoft-prod.deb && \
    dpkg -i packages-microsoft-prod.deb && \
    apt-get update && \
    apt-get install -y dotnet-sdk-8.0 && \
    apt-get clean && \
    rm packages-microsoft-prod.deb && \
    rm -rf /var/lib/apt/lists/*

# Install Extensions
# Install LANDIS
RUN wget -O v8.zip https://github.com/LANDIS-II-Foundation/Core-Model-v8-LINUX/archive/refs/heads/main.zip && \
    unzip v8.zip && \
    wget -O main.zip https://github.com/LANDIS-II-Foundation/Support-Library-Dlls-v8/archive/refs/heads/main.zip && \
    unzip -n main.zip -d /opt/Core-Model-v8-LINUX-main/build/extensions/ && \
    rm main.zip && \
    cp /opt/Core-Model-v8-LINUX-main/build/extensions/Support-Library-Dlls-v8-main/*.dll /opt/Core-Model-v8-LINUX-main/Tool-Console/src/lib/ && \
    cd /opt/Core-Model-v8-LINUX-main/Tool-Extension-Admin/src && dotnet build -c Release && \
    cd /opt/Core-Model-v8-LINUX-main/Tool-Console/src && dotnet build -c Release

#SCRPPLE test DGS compliant (this is the version within ML's personal repo- WORKS)
RUN wget https://github.com/LANDIS-II-Foundation/Extension-Social-Climate-Fire/archive/refs/heads/master.zip && \
    unzip master.zip -d /opt/Core-Model-v8-LINUX-main/build/extensions/ && \
    rm master.zip && \
    cp /opt/Core-Model-v8-LINUX-main/build/extensions/Support-Library-Dlls-v8-main/*.dll /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Social-Climate-Fire-master/src/lib/ && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Social-Climate-Fire-master/src; sed -i -e 's/lib\\Landis.Library/..\\..\\Support-Library-Dlls-v8-master\\Landis.Library/g' SocialClimateFire.csproj && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Social-Climate-Fire-master/src; sed -i -e 's/lib\\Ether./..\\..\\Support-Library-Dlls-v8-master\\Ether./g' SocialClimateFire.csproj && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Social-Climate-Fire-master/src; dotnet build -c Release SocialClimateFire.csproj && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Social-Climate-Fire-master/deploy/installer; dotnet /opt/Core-Model-v8-LINUX-main/build/Release/Landis.Extensions.dll add "/opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Social-Climate-Fire-master/deploy/installer/Scrapple 4.txt" && \
    ls /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Social-Climate-Fire-master/src/obj/Release/ && \
    cp /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Social-Climate-Fire-master/src/obj/Release/netstandard2.0/Landis.Extension.SocialClimateFire-v4.dll /opt/Core-Model-v8-LINUX-main/build/extensions/.

# Output Biomass (WORKS)
RUN wget -O master.zip https://github.com/LANDIS-II-Foundation/Extension-Output-Biomass/archive/refs/heads/master.zip && \
    unzip master.zip -d /opt/Core-Model-v8-LINUX-main/build/extensions/ && \
    rm master.zip && \
    cp /opt/Core-Model-v8-LINUX-main/build/extensions/Support-Library-Dlls-v8-main/*.dll /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Output-Biomass-master/src/lib/ && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Output-Biomass-master/src; sed -i -e 's/lib\\Landis.Library/..\\..\\Support-Library-Dlls-v8-master\\Landis.Library/g' output-biomass.csproj && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Output-Biomass-master/src; dotnet build -c Release output-biomass.csproj && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Output-Biomass-master/deploy/installer; dotnet /opt/Core-Model-v8-LINUX-main/build/Release/Landis.Extensions.dll add "Output Biomass 4.txt" && \
    cp /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Output-Biomass-master/src/obj/Release/netstandard2.0/Landis.Extension.Output.Biomass-v4.dll /opt/Core-Model-v8-LINUX-main/build/extensions/.

# Output Biomass Community 
RUN wget -O master.zip https://github.com/LANDIS-II-Foundation/Extension-Output-Biomass-Community/archive/refs/heads/master.zip && \
    unzip master.zip -d /opt/Core-Model-v8-LINUX-main/build/extensions/ && \
    rm master.zip && \
    cp /opt/Core-Model-v8-LINUX-main/build/extensions/Support-Library-Dlls-v8-main/*.dll /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Output-Biomass-Community-master/src/lib/ && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Output-Biomass-Community-master/src; sed -i -e 's/lib\\Landis.Library/..\\..\\Support-Library-Dlls-v8-master\\Landis.Library/g' output-biomass-community.csproj && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Output-Biomass-Community-master/src; dotnet build -c Release output-biomass-community.csproj && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Output-Biomass-Community-master/deploy/installer; dotnet /opt/Core-Model-v8-LINUX-main/build/Release/Landis.Extensions.dll add "Output Biomass Community 3.txt" && \
    cp /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Output-Biomass-Community-master/src/obj/Release/netstandard2.0/Landis.Extension.Output.BiomassCommunity-v3.dll /opt/Core-Model-v8-LINUX-main/build/extensions/.

# Output Biomass Succession
RUN wget -O master.zip https://github.com/LANDIS-II-Foundation/Extension-Biomass-Succession/archive/refs/heads/master.zip && \
    unzip master.zip -d /opt/Core-Model-v8-LINUX-main/build/extensions/ && \
    rm master.zip && \
    cp /opt/Core-Model-v8-LINUX-main/build/extensions/Support-Library-Dlls-v8-main/*.dll /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Biomass-Succession-master/src/lib/ && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Biomass-Succession-master/src; sed -i -e 's/lib\\Landis.Library/..\\..\\Support-Library-Dlls-v8-master\\Landis.Library/g' biomass-succession.csproj && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Biomass-Succession-master/src; dotnet build -c Release biomass-succession.csproj && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Biomass-Succession-master/deploy/installer; dotnet /opt/Core-Model-v8-LINUX-main/build/Release/Landis.Extensions.dll add "/opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Biomass-Succession-master/deploy/installer/Biomass Succession 7.txt" && \
    cp /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Biomass-Succession-master/src/obj/Release/netstandard2.0/Landis.Extension.Succession.Biomass-v7.dll /opt/Core-Model-v8-LINUX-main/build/extensions/.

# Output Biomass Reclass
RUN wget -O master.zip https://github.com/LANDIS-II-Foundation/Extension-Output-Biomass-Reclass/archive/refs/heads/master.zip && \
    unzip master.zip -d /opt/Core-Model-v8-LINUX-main/build/extensions/ && \
    rm master.zip && \
    cp /opt/Core-Model-v8-LINUX-main/build/extensions/Support-Library-Dlls-v8-main/*.dll /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Output-Biomass-Reclass-master/src/lib/ && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Output-Biomass-Reclass-master/src; sed -i -e 's/lib\\Landis.Library/..\\..\\Support-Library-Dlls-v8-master\\Landis.Library/g' output-biomass.csproj && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Output-Biomass-Reclass-master/src; dotnet build -c Release output-biomass.csproj && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Output-Biomass-Reclass-master/deploy/installer; dotnet /opt/Core-Model-v8-LINUX-main/build/Release/Landis.Extensions.dll add "Output Biomass Reclass 4.txt" && \
    cp /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-Output-Biomass-Reclass-master/src/obj/Release/netstandard2.0/Landis.Extension.Output.BiomassReclass-v4.dll /opt/Core-Model-v8-LINUX-main/build/extensions/.

# GIPL (WORKS)
RUN wget https://github.com/LANDIS-II-Foundation/Library-GIPL/archive/master.zip && \
    unzip master.zip -d /opt/Core-Model-v8-LINUX-main/build/extensions/ && \
    rm master.zip && \
    mv /opt/Core-Model-v8-LINUX-main/build/extensions/Library-GIPL-master/ /opt/Core-Model-v8-LINUX-main/build/extensions/Library-GIPL/ && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Library-GIPL/src; sed -i -e 's/lib\\Landis.Library/..\\..\\Support-Library-Dlls-v8-master\\Landis.Library/g' GiplDamm.csproj && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Library-GIPL/src; dotnet build -c Release GiplDamm.csproj && \
    cp /opt/Core-Model-v8-LINUX-main/build/extensions/Library-GIPL/src/bin/Debug/netstandard2.0/Landis.Extension.GiplDamm.dll /opt/Core-Model-v8-LINUX-main/build/extensions/.  
 
# SHAW (WORKS)
RUN wget https://github.com/LANDIS-II-Foundation/Library-Shaw/archive/refs/heads/master.zip && \
    unzip master.zip -d /opt/Core-Model-v8-LINUX-main/build/extensions/ && \
    rm master.zip && \
    mv /opt/Core-Model-v8-LINUX-main/build/extensions/Library-Shaw-master/ /opt/Core-Model-v8-LINUX-main/build/extensions/Library-Shaw/ && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Library-Shaw/src; sed -i -e 's/lib\\Landis.Library/..\\..\\Support-Library-Dlls-v8-master\\Landis.Library/g' ShawDammNewInputs.csproj && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Library-Shaw/src; dotnet build -c Release ShawDammNewInputs.csproj && \
    cp /opt/Core-Model-v8-LINUX-main/build/extensions/Library-Shaw/src/bin/Debug/netstandard2.0/Landis.Extension.ShawDamm.dll /opt/Core-Model-v8-LINUX-main/build/extensions/. 

# DGS (WORKS)
# RUN wget https://github.com/LANDIS-II-Foundation/Extension-DGS-Succession/archive/refs/heads/master.zip && \
#     unzip master.zip -d /opt/Core-Model-v8-LINUX-main/build/extensions/ && \
#     rm master.zip && \
#     mv /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-DGS-Succession-master/ /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-DGS-Succession/ && \
#     cd /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-DGS-Succession/src; sed -i -e 's/lib\\Landis.Library/..\\..\\Support-Library-Dlls-v8-master\\Landis.Library/g' DGS-succession.csproj && \
#     cd /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-DGS-Succession/src; dotnet build -c Release DGS-succession.csproj && \
#     cd /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-DGS-Succession/deploy/installer; dotnet /opt/Core-Model-v8-LINUX-main/build/Release/Landis.Extensions.dll add "DGS_Succession 1.102.txt" && \
#     cp /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-DGS-Succession/src/bin/Debug/netstandard2.0/Landis.Extension.Succession.DGS.dll /opt/Core-Model-v8-LINUX-main/build/extensions/. 

# NECN (WORKS)
RUN wget https://github.com/LANDIS-II-Foundation/Extension-NECN-Succession/archive/refs/heads/master.zip && \
    unzip master.zip -d /opt/Core-Model-v8-LINUX-main/build/extensions/ && \
    rm master.zip && \
    cp /opt/Core-Model-v8-LINUX-main/build/extensions/Support-Library-Dlls-v8-main/*.dll /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-NECN-Succession-master/src/lib/ && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-NECN-Succession-master/src; sed -i -e 's/lib\\Landis.Library/..\\..\\Support-Library-Dlls-v8-master\\Landis.Library/g' NECN-succession.csproj && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-NECN-Succession-master/src; dotnet build -c Release NECN-succession.csproj && \
    cd /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-NECN-Succession-master/deploy/installer; dotnet /opt/Core-Model-v8-LINUX-main/build/Release/Landis.Extensions.dll add "/opt/Core-Model-v8-LINUX-main/build/extensions/Extension-NECN-Succession-master/deploy/installer/NECN_Succession8.txt" && \
    cp /opt/Core-Model-v8-LINUX-main/build/extensions/Extension-NECN-Succession-master/src/obj/Release/netstandard2.0/Landis.Extension.Succession.NECN-v8.dll /opt/Core-Model-v8-LINUX-main/build/extensions/.

# Add a hint path for all of the extensions to the Tool-Console csproj file and then rebuild it 
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
    <Reference Include="Landis.Extension.SocialClimateFire-v4.dll">\n\
        <HintPath>../../build/extensions/Landis.Extension.SocialClimateFire-v4.dll</HintPath>\n\
    </Reference>\n\
    <Reference Include="Landis.Extension.GiplDamm.dll">\n\
        <HintPath>../../build/extensions/Landis.Extension.GiplDamm.dll</HintPath>\n\
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
    <Reference Include="Landis.Extension.ShawDamm.dll">\n\
        <HintPath>../../build/extensions/Landis.Extension.ShawDamm.dll</HintPath>\n\
    </Reference>\n\
    <Reference Include="Landis.Extension.Succession.DGS.dll">\n\
        <HintPath>../../build/extensions/Landis.Extension.Succession.DGS.dll</HintPath>\n\
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
</ItemGroup>' /opt/Core-Model-v8-LINUX-main/Tool-Console/src/Console.csproj

# Move Support Libraries to the correct place
RUN mv /opt/Core-Model-v8-LINUX-main/build/extensions/Support-Library-Dlls-v8-main/*.dll /opt/Core-Model-v8-LINUX-main/build/extensions/

RUN chown -R 1000:1000 /opt/Core-Model-v8-LINUX-main
    
ARG LOCAL_USER=user
ARG PRIV_CMDS='/bin/ch*,/bin/cat,/bin/gunzip,/bin/tar,/bin/mkdir,/bin/ps,/bin/mv,/bin/cp,/usr/bin/apt*,/usr/bin/pip*,/bin/yum'

RUN adduser --disabled-password --gecos "VICE_User" --uid 1000 user  && \
    usermod -aG sudo user && \
    echo "$LOCAL_USER ALL=NOPASSWD: $PRIV_CMDS" >> /etc/sudoers

RUN apt-get update && \
    apt-get install -y wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 \
    gettext-base git mercurial subversion \
    tmux && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> /home/user/.bashrc && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> /home/user/.zshrc && \
    chown -R user:user /opt/conda

# install ttyd
RUN curl -L "https://github.com/tsl0922/ttyd/releases/download/1.6.3/ttyd.x86_64" > ttyd && \
    chmod a+x ttyd && \
    mv ttyd /usr/local/bin/ttyd

RUN apt-get update && \
    apt-get install -y curl grep sed dpkg && \
    curl -L "https://github.com/krallin/tini/releases/download/v0.19.0/tini_0.19.0-amd64.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install a few dependencies for iCommands, text editing, and monitoring instances
RUN apt-get update && \
    apt-get install -y lsb-release apt-transport-https curl gnupg2 libfuse2 gcc less nodejs software-properties-common apt-utils glances htop nano  && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN apt-get update && \
    apt-get install -y vim-nox && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN apt-get update && \
    apt-get install -y emacs-nox && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget -qO - https://packages.irods.org/irods-signing-key.asc | apt-key add - && \
    echo "deb [arch=amd64] https://packages.irods.org/apt/ $(lsb_release -sc) main" >> /etc/apt/sources.list.d/renci-irods.list && \
    apt-get update && \
    apt install -y irods-icommands 

USER user

# set shell as bash and terminal as linux
ENV SHELL=bash
ENV TERM=xterm

# open port 7681 for ttyd
EXPOSE 7681

# changes tmux layout while running
COPY entry.sh /bin
RUN echo 'set-option -g status off' >> ~/.tmux.conf

# add iRODS iCommands to user profile as JSON
RUN mkdir /home/user/.irods 

RUN cd /opt/Core-Model-v8-LINUX-main/Tool-Extension-Admin/src && dotnet build -c Release
RUN echo "alias console='dotnet /opt/Core-Model-v8-LINUX-main/build/Release/Landis.Console.dll'" >> ~/.bashrc
RUN echo "cd /opt/Core-Model-v8-LINUX-main/Tool-Console/src && dotnet build -c Release && cd" >> ~/.bashrc
ENTRYPOINT ["bash", "/bin/entry.sh"]

CMD ["ttyd", "bash"]

WORKDIR /home/user