function CreateFolders()
global epname
if(~isfolder('ErrorPlots'))
    mkdir('ErrorPlots');
    cd ErrorPlots
    mkdir(epname)
    cd ../
else
    cd ErrorPlots
    if(~isfolder(epname))
        mkdir(epname)
    end
    cd ../
end
if(~isfolder('Phasespace'))
    mkdir('Phasespace');
    cd Phasespace
    mkdir(epname)
    cd ../
else
    cd Phasespace
    if(~isfolder(epname))
        mkdir(epname)
    end
    cd ../
end
if(~isfolder('Basisplots'))
    mkdir('Basisplots');
    cd Basisplots;
    mkdir(epname)
    cd ../
else
    cd Basisplots
    if(~isfolder(epname))
        mkdir(epname)
    end
    cd ../
end
end