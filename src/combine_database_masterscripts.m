function combine_database_masterscripts()

% unfinished

%% it works with haining's 2D database
% global vr_set theta_set X_set Y_set color_set targetpar_set all_Y all_X all_phases

utility = Utility.getInstance();
filesep = utility.getFilesep();

compType = computer;
if strcmp(compType, 'PCWIN') == 1
    outputPath =  'D:\Dropbox\VIV\src\ILCFprediction\databaseComb\';
    dataPath_JMD =  'D:\Dropbox\VIV\src\ILCFprediction\databaseJMD\';
    dataPath_HZ =  'D:\Dropbox\VIV\src\ILCFprediction\databaseHZ\'; 
else
    outputPath =  '/Users/haining/VIV/src/ILCFprediction/databaseComb/';
        dataPath_JMD =  '/Users/haining/VIV/src/ILCFprediction/databaseJMD/';
    dataPath_HZ =  '/Users/haining/VIV/src/ILCFprediction/databaseHZ/'; 
end

utility.conditionalMkdir(outputPath);

color_set = {'r','y','g','c','b','m','k',[255/255 140/255 0],'r'};
targetpar_set = {'CL1','CL3','CD2','Cmy', 'Cmx','pow','CDv'}; 

%% Dataset JMD
% 6*6*8*8 = 2304  6*6*9*8 = 2592
X_set_JMD = {'0','0.15','0.3','0.45','0.6','0.75'};
Y_set_JMD = {'0.25','0.5','0.75','1', '1.25', '1.5'};
theta_set_JMD = {'180n' '135n' '90n' '45n' '0' '45' '90' '135'};
vr_set_JMD = {'4p5' '5' '5p5' '6' '6p5' '7' '7p5' '8'}; 

all_Y_JMD = [0.25:0.25:1.5];

X_vector_JMD = [0:0.15:0.75 [0:0.15:0.75] [0:0.15:0.75] [0:0.15:0.75] [0:0.15:0.75] [0:0.15:0.75]];
Y_vector_JMD = [0.25*ones(1,6) 0.5*ones(1,6) 0.75*ones(1,6) 1*ones(1,6) 1.25*ones(1,6) 1.5*ones(1,6)];
EXP_matrix_JMD = [X_vector_JMD; Y_vector_JMD];

vrIndex_JMD = [4];
AmpXIndex_JMD = [1:6];
AmpYIndex_JMD = [1:6];
targetpar_JMD = 'Cmy';

%% Dataset 2
% 5*6*18*9 = 4860 5*6*19*9 = 5130
% 21*18*9 = 3402 21*19*9 =3591
vr_set_HZ = {'4' '4.5' '5' '5.5' '6' '6.5' '7' '7.5' '8'}; 
X_set_HZ = {'0.05','0.1','0.15','0.2','0.25','0.3'};
Y_set_HZ = {'0.15','0.25','0.5','0.75', '1'};
theta_set_HZ = {'180n' '157.5n' '135n' '112.5n' '90n' '67.5n' '45n' '30n' '15n' '0'...
     '15' '30' '45' '67.5' '90' '112.5' '135' '157.5'};
 
all_Y_HZ = [0.15,0.25,0.5,0.75, 1 ];

%% 21
X_vector_HZ = [0.05 0.1 0.05 0.1 0.15 0.05 0.1 0.15 0.2 0.05 0.1 0.15 0.2 0.25 0.3 0.05 0.1 0.15 0.2 0.25 0.3];
Y_vector_HZ = [0.15 0.15 0.25 0.25 0.25 0.5 0.5 0.5 0.5 0.75 0.75 0.75 0.75 0.75 0.75 1 1 1 1 1 1];
EXP_matrix_HZ = [X_vector_HZ; Y_vector_HZ];

vrIndex_HZ = [5];
AmpXIndex_HZ = [1:6];
AmpYIndex_HZ = [1:5];
targetpar_HZ = 'Cmy';

%% combined Database


% 7*10*9*18 

[rawDatabaseStrucHZ, rawDatabaseDSHZ] = loadDatabaseHZ();

[rawDatabaseStrucJMD, rawDatabaseDSJMD] = loadDatabaseJMD();

rawDatabaseStruc_comb = combineDatabase(rawDatabaseStrucHZ, rawDatabaseStrucJMD);

JoinDatabase = cat(1, rawDatabaseDSJMD,rawDatabaseDSHZ);

close all
