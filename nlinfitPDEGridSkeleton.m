%Skeleton MATLAB scripts to illustrate the image processing and PDE model}
%fitting algorithms used in 'Quantitative analysis reveals that Actin and
%Src-family kinases regulate nuclear YAP1 and its export' by Ege et al.
%Please reference 'Quantitative analysis reveals that Actin and
%Src-family kinases regulate nuclear YAP1 and its export' by Ege et al. if
%using or copying any part of this code

function output = nlinfitPDEGridSkeleton(pars,Tnlinfit)

    %Reads in the gridpoint neighbours such that the correct boundary
    %conditions are imposed at each gridpoint.
    directory='C:\DirectoryOfCellNeighbours\';
    ExcelRead=[directory 'gridneighbours.xlsx'];
    [nucnuc,txt,raw]=xlsread(ExcelRead,'NucNucNeighbours');
    [nuccyto,txt,raw]=xlsread(ExcelRead,'NucCytoNeighbours');
    [cytonuc,txt,raw]=xlsread(ExcelRead,'CytoNucNeighbours');
    [cytocyto,txt,raw]=xlsread(ExcelRead,'CytoCytoNeighbours');

    %Change this to 'BleachpointIndex' output from 'PDEFitSkeleton.m' so
    %model knows where bleach is happening
    bleachid=10;

    [nonuc,~]=size(nucnuc);
    [nocyto,~]=size(cytocyto);
    tspanodedata=Tnlinfit(1:length(Tnlinfit)/(nonuc+nocyto));      

    %Example free parameters to be fitted (import/export between
    %compartments and immobile/mobile reactions in the cytoplasm are not
    %included here.)
    k1=pars(1);%Assoc rate in nucleus
    eta=pars(2);%Decay rate due to bleaching
    C0=pars(3);%Initial condition in cytoplasm
    %Example fixed parameters that must be input by the user. The rate of
    %dissociation from immobile to mobile state in the nucleus. Must take a
    %numerical value.
    km1=KM1;
    %The diffusion parameter used in the numerics where D1 is the true rate
    %of diffusion and h accounts for the gridsize in the numerics. Must
    %take a numerical value.
    D=D1/h^2;

    %Example generation of initial conditions: Generate initial conditions
    %for each gridpoint in the nucleus and cytoplasm based on assumed
    %transfer rates between each state and compartment.
    
    %e.g. initial mobile fraction in nucleus is a function, f, of
    %cytoplasmic initial concentration.
    M0=f(C0);
    %e.g. initial immobile fraction in nucleus is a function of mobile
    %fraction in nucleus and association and dissociation rates.
    I0=abs(k1/km1*M0); 
    %Set each gridpoint in nucleus and cytoplasm to these initial
    %conditions.
    for I=1:nonuc;
        z0(2*I-1)=I0;
        z0(2*I)=M0;
    end;
    for I=1:nocyto;
        z0(I+2*nonuc)=C0;
    end;

    %PDE reduced to ODE model to be solved using ode solver
    function dz = FlipFullModelODE(t,z)
        %Assume we are in a system where there is a mobile and immobile
        %fraction in the nucleus but the protein is all mobile in the
        %cytoplasm
        dz = zeros(2*nonuc+nocyto,1); 
        %Generating odes for nuclear gridpoints
        for I=1:nonuc
            I;
            %Determine the nuclear and cytoplasmic neighbours for each
            %nuclear gridpoint
            nucneighbours=nucnuc(I,:);
            nucneighbours=nucneighbours(nucneighbours>0);
            cytoneighbours=nuccyto(I,:);
            cytoneighbours=cytoneighbours(cytoneighbours>0);

            %If the gridpoint is the bleachpoint we need to include decay
            %due to bleaching (eta).
            if(I==bleachid)
                %ODE for immobile fraction. Incorporates reactions between
                %mobile and immobile states and decay due to bleaching.
                dz(2*I-1)=-abs(km1)*z(2*I-1)+abs(k1)*z(2*I)-eta*z(2*I-1);
                %ODE for mobile fraction. Incorporates transfer between
                %mobile and immobile fractions, decay due to bleaching and
                %diffusion with the rest of the nucleus.
                dz(2*I) = abs(km1)*z(2*I-1)-abs(k1)*z(2*I)...
                    +D*sum(z(2*nucneighbours))...
                    -D*length(nucneighbours)*z(2*I)-eta*z(2*I); 

            else
                %If the gridpoint is elsewhere in the nucleus include
                %conditions for boundaries. I.e. include cases for nuclear
                %only neighbours, nuclear and cytoplasmic neighbours or
                %cytoplasmic only neighbours (caused by coarse
                %discretization of cell). These will determine the
                %inclusion of diffusion and import/export functions etc. in
                %the numerics
            end
        end;

        for I=1:nocyto
            nucneighbours=cytonuc(I,:);
            nucneighbours=nucneighbours(nucneighbours>0);
            cytoneighbours=cytocyto(I,:);
            cytoneighbours=cytoneighbours(cytoneighbours>0);
            %Generate relevant and equivalent odes for each cytoplasmic
            %gridpoint here.

        end;
    end;

    %Output of ode solver that nlinfit attempts to fit to the experimental
    %data.
    [T,Z] = ode15s(@FlipFullModelODE,tspanodedata,z0);
    output=[];
    for I=1:nonuc
        output=[output;Z(:,2*I-1)+Z(:,2*I)];
    end;
    for I=1:nocyto
        output=[output;Z(:,2*nonuc+I)];
    end;
end
  
 
  
  
  