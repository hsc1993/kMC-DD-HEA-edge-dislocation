classdef simulation_parameters
    %Class to hold all the simulation parameters to be easily passed
    %accessed throughout
    
    properties

%Physical Constants%%%%
Kb = 1.3806*10E-23; %Boltzmann Constant
Na = 6.022*10E23; %Avogadro's number in atoms/mole

%%%%Material Constants%%%%
density = 7.86;  %Bulk density of material [g/cm^3]
molarMass = 55.845; %Molar mass [g/mol]
vacFormationEnergy = 2.7237*10E-19; %Energy per vacancy formation [Joules]

b = 0.25*10E-9; %Burger's vector [m]
G = 80*10E9; %Shear modulus in [Pa]
poissonRatio = 0.29; %Poisson's Ratio
Eo = 9.6*10E-20; %Vacancy Migration Energy [Joules]
v = 1*10E12; %Vibrational frequency [Hz]

atomicVolume = 0.77*((0.25*10E-9)^3); %Atomic volume [m^3] = 0.77(b^3)
volumetricStrain = 0.1; %Volumetric strain
nearestNeighbors = 8; %8 Nearest neighbors for bcc
correlationFactor = 0.78; %Correlation factor 

%%%%Simulation Constants%%%%
% int numNodes = 14; %Number of nodes to discretize across the straight dislocation that runs across the simulation box
% double fractZOverhang = .1; %Fraction of the total z-length that the z-length straight dislocaiton extends above and below the simulation box
% double fractXOverhang = .2; %Fraction of the total x-length that the dislocation loop extends beyond the simulation box
% int DisStructUpdateFreq = 5; %Number of runID trials between updating the DislocationStructure and absorbing vacancies
% int PrintStatsFreq = 100; %Frequency of outputting the test statistics
% int PrintEVLFreq = 100; %Frequency of outputting the evl file
% int PrintCPUTime = 0; %Whether or not to print the running CPU time for the simulation in results/stats.txt

kmcAcceptableError = 30; %Acceptable percentage error for self-guess of kMC timestep
outputGlobalTimeStep = 0; %1=true or 0=false for printing globaltimestep file -- printed in "test" folder
outputV = 0; %1 or 0 for printing V_#.txt files -- printed in "evl1" folder
outputEVL=1; %1 or 0 for printing evl_#.txt files -- printed in "evl1" folder

%%Do not change%%%
totalGlobalTime = 0; %Running counter for the total gloabl time
lastTotalGlobalTime = 0; %Running counter for the previous gloabl time to be used in velocity calculation
lastDistanceMoved = 0;

%%%%Mesh-Related%%%%
meshType = 0; % 0 = rectangular prism mesh, 1 = cylindrical mesh

L1=1632; %the side length of the cube, in units of Burgers vector
L2=1510; %the side length of the cube, in units of Burgers vector
L3=1520; %the side length of the cube, in units of Burgers vector

%If meshType==0
%L1 = x-length
%L2 = y-length
%L3 = z-length

%If meshType==0
%L1 = radius
%L2 = z-length

%%%%Physically Relevant Simulation Constants%%%%
Temp = 1250; %Temperature (Kelvin)
minVacJump = (0.25*10E-9)/100; %Minimum jump distance of vacancy for MC purposes = b/100
maxVacJump = 0.25*10E-9; %Maximum jump distance of vacancy for MC purposes = b

%%%%Vacancy Related%%%%
distToAbsorbption = 2; %Distance in Burger's units for a vacancy to be absorbed by a dislocation
vacancyConcentration =0; %2e21; %Volumetric Concentration of vacancies in vac/m^3 - overrides the calculation of the number of vacancies if not vacancyConcentration=0
useEmission = 1; %0=do not use emission of vacancaies -- randomly replace them each time they are absorbed, 1=do not automatically replace vacancies -- have segments emit vacancies and negatively climb
useDiscreteEmission = 1; %1 = use discrete emission from local concentration of vacancies around a dislocation segment, 0=use a global average of vacancy concentration for vacancy emission

%%DO NOT CHANGE%%
RunningVacAbsorbed=0; %Running counter of how many vacancies have been absorbed
RunningLastVacNumber = 0; %Running counter of how many vacancies were previously absorbed before last calculation of dislocation velocity
RunningVacIDnum = 0; %Running counter for the ID of vacancies
RunningVacEmitted=0; %Running counter of how many vacancies have been emitted
RunningLastVacEmitted=0; %Running counter of how many vacancies were previously emitted before last calculation of dislocation velocity
RunningBalancedVacs = 0; %Running counter of how many vancacies were added (+) or subtracted (-) from the system to keep the concentration constant
RunningCPUTime = 0; %Running counter of the total CPU time elapsed
vacNum; %Number of vacancies - to be calculated later

%%%%Random Generators%%%%
% random_device generator; %Generator to select random numbers throughout the simulation
% 
% %%%%For Vacancy Movement%%%%
% uniform_real_distribution<double> ZeroOnedistribution(0, 1); %Uniform distribution to be used throughout program
% uniform_real_distribution<double> Stepdistribution(minVacJump, maxVacJump); %Uniform distribution for step distances to be used throughout the program
int
%%%%For Random Vacancy Placement In Mesh%%%% -- had to explicitly redefine each generator within the vacancies/defectiveCrystal file due to dynamicBoxResizing 
%std::uniform_real_distribution<double> L1Dist(-L1/2,L1/2);
%std::uniform_real_distribution<double> L2Dist(-L3/2,L2/2);
%std::uniform_real_distribution<double> L3Dist(-L3/2,L3/2);

%%%%Parametric Testing of Climb Parameters%%%%
% int useParametricStudy = 1; %0=run non-parametric study, 1 = run parametric study by varying the variables listed below
DDSimulationText=0; %0=minimal output from the DD module, 1=full DD and vacancy simulation text output

dynamicBoxResizing = 1; %0=do not dynamically resize the box to maintain a constant vacancy number, 1 = resize box as to keep a constant vacancy number throughout the parametric study
constVacNumber = 100; %Number of vacancies to consider for each trial of the parametric stdy IF dynamicBoxResizing==1

requireConstantVacs = 1; %0=do not add/remove vacancies to keep the pre-set vacancy concentration, 1=ensure that the pre-set vacancy concentration is kept by adding/removing vacancies
% 
% int tempTrials = 4; %Number of intermediary temperature trials to do within min and max bounds
% int concTrials = 4; %Number of intermediary concentration multiplier trials to do within min and max bounds
% int stressTrials = 4; %Number of intermediary stress trials to do within min and max bounds
% int trials = 4; %Number of trials to do within the min and max bounds. For example, if minTemp=100, maxTemp=200, and trials = 3, simulations would be run at 100, 150, and 300.
% int repeatTrials = 1; %How many times to repeat each simulation
% 
% double minTemp = 1000; %Minimum temperature to use in parametric study [Kelvins]
% double maxTemp = 1500; %Maximum temperature to use in parametric study [Kelvins]
% 
% double minVacConcentration = 1; %Minimum vacancy concentration to use in parametric study  - must be a multiplier of the thermal vacancy concentration[vacs/m^3] 
% double maxVacConcentration = 100; %Maximum vacancy concentration to use in parametric study - must be a multiplier of the thermal vacancy concentration [vacs/m^3]
% %If the thermal vacancy concentration is 1E10 and maxVacConcentration=5, the maximum vacancy concentration will be 5E10.
% 
% int useStress = 1;  %0 = do not use applied stress in the parametric study, 1 = use applied stress in the parametric study
% double minPressure = -100*E6); %Minimum stress [pressure] to use in parametric study [Pa]
% double maxPressure = 100*E6);  %Maximum stress [pressure] to use in parametric study [Pa]
% double appliedPressure = 0; %DO NOT CHANGE -- this is the default (thermal) pressure used to calculate the enthalpy of vacancy formation
%       
          
    end
    
    

end

