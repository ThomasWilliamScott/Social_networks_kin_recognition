% INFINITE ISLAND MODEL AGENT BASED SIMULATION (TESTING FOR BALANCING 
% SELECTION IN A FINITE POPULATION). 

% This script generates data from a version of the agent based simulation
% in which tag mutation is absent. In this version of the model, tag
% diversity is inevitably lost, and we record the amount of time this takes
% to happen under different parameterisations (tag fixation time), as this 
% information can be used to examine whether balancing selection is in
% operation. 

% The data generated using this script is saved in the 'Saved_data' folder 
% as 'Balancing_sel_finite_pop_data.mat' and was used to generate Fig. 4f 
% & S8.

% This script runs the agent based simulation model, for a given set of 
% parameter values (specified above), and records the time taken for tag 
% diversity to be lost (tag fixation time). It then records 'tag fixation
% time' in a results matrix called 'resfix'. Each run of the simulation
% generates a 'tag fixation time', which is stored in a unique row of
% resfix. Specifically, resfix is a 6 column matrix, with each row
% corresponding to a run of the simulation. Column 1 records the 'c' value
% for that run; column 2 records the population size for that run (D*N);
% column 3 records the 'alpha' value for that run; column 4 records the
% 'csearch' value for that run; column 6 records the trait mutation rate
% for that run; column 5 records the output - i.e. the tag fixation time -
% for that run.

% Running this script will generate a results matrix (resfix) with some
% results entries corresponding to certain parameter regimes. In order to 
% generate all the data featured in our 'Balancing_sel_finite_pop_data.mat'
% file, this script must be run many times, for different parameter 
% regimes, and for multiple iterations of each parameter regime. This can 
% be acheived by adding in 'for loops' to this script, so that the results 
% matrix (resfix) acquires data for various different parameter 
% combinations.

clearvars
clc

TagMu = 0; % Tag mutation rate. This needs to be set to zero in this script (should not be changed from zero).
count=1; % This count variable will be updated according to count=count+1 every time a new result is obtained. This will allow us to sequentially add results to a results matrix. 

% PARAMETERS TO BE SET AS INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=0.01; % Migration rate
r=0.01; % Recombination rate
TraitMu = 0.005; % Trait mutation rate
alphaR = 1; % Partner-reassociation probability
L=2; % Maximum segregating tags (N.B. This is labelled as 'L_{max}' in the main text, and is often labelled as 'tag' in other scripts).
T=1000000; % This gives an upper bound on the length of time (generations) that each simulation is run for. This should be set to be sufficiently high that the upper bound is never reached (reaching it would generate artificial tag fixation times).
N=7; % Deme size 
DR=[30 35 40 45 50]; % Number of demes in the population.
alpha=1; % Encounter parameter.
c=0.5; % Cost of cooperation
b=4.5; % Benefit of cooperation
csearch=0; % Cost of abandononing a partner and re-associating for a new social encounter.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for cur_D=1:numel(DR) % This for loop is used to generate results for different population sizes. Similar for loops can be added for other variables to generate data for other parameter combinations.
    
    D = DR(cur_D);
    nInd = N*D;% Overall population size.
    
% The population is tracked using a 2 dimensional matrix called 'pop'. 
% This matrix is updated each generation. Each row in 'pop' represents an 
% individual in the population (there are therefore nInd rows in total). 
% Column 1 denotes an individual's ID, which is a unique number ranging from
% 1 to nInd. IDs are attributed to each individual because it helps with
% indexing. Column 2 denotes an individual's deme ID. This ranges from 1 to
% D. Column 3 denotes an individual's trait ID. This is 0 if the individual
% is a defector, and 1 if the individual is a cooperator. Column  4 denotes
% an individual's tag ID. This ranges from 1 to L. Column 5 denotes the
% number of same-tag defectors within an individual's deme. Column 6 
% denotes the number of same-tag cooperators within an individual's deme. 
% Column 7 denotes the fitness of the individual (fitness depends on 
% success in social interactions).

% SETTING INITIAL POPULATION STATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pop = zeros(nInd,7);
pop(:,1) = 1:nInd ; % Column 1: individual ID
pop(:,2) = repelem(1:D,1,N) ; % Column 2: deme ID
pop(:,3) = randsample(2,nInd,true,[0.5 0.5]) -1; % Column 3: Trait ID. Conditional helper (1) or defector (0). Cooperator frequency is assumed to be ~0.5 initially.
maj = 1/L;
abc = rand(1,L-1);
abc = [maj , (abc ./ sum(abc)) .* (1-maj)];
abc(1) = maj;
abc = abc ./ sum(abc);
pop(:,4) = randsample(L,nInd,true,abc) ; % Column 4: tag ID. Assume tag freqs initially roughly equal in each deme.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following 'help' and 'def' matrices record the number of cooperators 
% and defectors bearing each tag. This will be updated as generations pass.
for l=1:L
help(l,1) = sum(pop(:,3)==1 & pop(:,4)==l) ./nInd; 
def(l,1) = sum(pop(:,3)==0 & pop(:,4)==l) ./nInd; 
end


for t=1:T % This loops over generations, and updates 'pop' each generation.
   clear W parentID parent2ID mutpop

% The following three for loops are used to update columns 5 & 6 of 'pop'.    
 for l = 1:L
    for d = 1:D  
        for j=0:1
            pop( pop(:,2)==d & pop(:,3)==j & pop(:,4)==l ,5) = sum( pop(:,2)==d & pop(:,3)==0 & pop(:,4)==l ) ; % This enters the number of same-tag defectors in each individual's deme (including the individual itself).
            pop( pop(:,2)==d & pop(:,3)==j & pop(:,4)==l ,6) = sum( pop(:,2)==d & pop(:,3)==1 & pop(:,4)==l ) ; % This enters the number of same-tag cooperators in each individual's deme (including the individual itself).       
        end
    end
 end

 % The following two lines are used to update columns 5 & 6 again. Before,
 % column 5 (6) gave the number of defectors (cooperators) in the
 % individual's deme, including the individual itself. But we need to
 % update this so that column 5 (6) gives the number of defectors 
 % (cooperators) in the individual's deme, *discounting* the individual 
 % itself ('others-only').
 pop(:,5) = pop(:,5) - (1-pop(:,3)) ; 
 pop(:,6) = pop(:,6) - pop(:,3) ; 
 
 % The following line enters individual fitness into column 7.
 pop(:,7) = 1 + ( b .* (pop(:,6)./(N-1)) - c .* pop(:,3) .* ((pop(:,5)+pop(:,6))./(N-1)) -csearch.*alpha.*((N-1-pop(:,5)-pop(:,6))./(N-1))) ./ (1 - alpha .* ((N-1-pop(:,5)-pop(:,6))./(N-1)) )  ;
 
 % The following line converts NaN fitness values into 1. NaN values are
 % obtained when an individual has no tag-mates in its deme (so can't
 % socially interact). We therefore need to convert these NaN values into 1
 % (baseline fitness).
 pop(isnan(pop))=1; 
 
for d = 1:D   
    % For a given deme (d), we define an array 'W', which gives the
    % contribution of gametes from each individual to that deme. Each
    % individual contributes a proportion 1-m gametes to its local deme,
    % and a proportion m gametes to other demes. 'W' is defined (for a
    % given deme, d) using the following two lines of code.
    W(pop(:,2)==d) = pop(pop(:,2)==d,7).*(1-m); 
    W(pop(:,2)~=d) = pop(pop(:,2)~=d,7).*(m/(D-1));     
    
    % Having obtained the fecundity array 'W' for a given deme, d, we use
    % this to choose 'parents', which seed the next generation.
    % Specifically, for each deme, we choose N (deme size) 'parent 1s', and
    % N 'parent 2s'. The probability of being chosen as a parent is
    % proportional to the individual's deme-specific fecundity ('W'). By
    % looping over all demes (d), we obtain two arrays, 'parent2ID' & 
    % 'parentID', which give the IDs of all parent 1s and parent 2s across 
    % the population.
    parent2ID((d-1)*N+1:d*N) = randsample(pop(:,1),N,true,W) ;
    parentID((d-1)*N+1:d*N) = randsample(pop(:,1),N,true,W) ;  
end

% We assume that parent 1 and parent 2 mate together to produce a zygote.
% The probability of a recombination event occuring in the zygote is given
% by the following array:
 rec = rand(nInd,1)>(1-r) ; 

% The trait & tag ID of the zygote (after recombination) is obtained, and 
% entered into columns 3 & 4 of 'pop' (matrix is updated to account for the
% next generation), using the following three lines of code.
pop(:,3:4) = pop(parentID,3:4);  
recomb = pop(parent2ID,4) ;
pop(:,4) = pop(:,4).*(1-rec) + recomb.*rec ;
 
 
% The probability, for each individual, of a mutation event occuring at the
% trait locus is given by the following array:
mutpop = rand(nInd,1)>1-TraitMu ;
% We update column 3 accordingly to account for trait mutation.
pop(:,3) = pop(:,3)+mutpop ;
pop( pop(:,3) == 2 , 3) = 0;

% The probability, for each individual, of a mutation event occuring at the
% tag locus is given by the following array:
mutpop = rand(nInd,1)>1-TagMu ;
% We update column 4 accordingly to account for tag mutation.
pop(:,4) = pop(:,4) .* (1-mutpop) + mutpop .* randi(L,[nInd,1]) ;

% We can then update our 'help' and 'def' matrices, which keep track of the
% number of cooperators and defectors bearing each tag, across generations.
for l=1:L
help(l,t+1) = sum(pop(:,3)==1 & pop(:,4)==l) ./nInd ; 
def(l,t+1) = sum(pop(:,3)==0 & pop(:,4)==l) ./nInd ; 
end

% The following calculation and 'if statement' serves to calculate tag
% diversity, and to stop the simulation either if tag diversity = 1 
% (i.e 1 tag has gone all the way to fixation) or if the simulation has
% been running for some upper length of time given by T (note that T should
% be set to be sufficiently high than runs never reach this upper bound - 
% reaching this upper bound generates artificial tag fixation times).
AAA(:,t) = (help(:,t) + def(:,t)).^2 ;
if 1/sum(mean(AAA(:,t),2)) == 1 | t==T  
    resfix(count,1) = c; % cost of cooperation is entered into column 1
    resfix(count,2) = D * N; % population size is entered into column 2
    resfix(count,3) = alpha; % encounter parameter is entered into column 3
    resfix(count,4) = csearch; % cost of partner search is entered into column 4
    resfix(count,5) = t; % tag fixation time is entered into column 5
    resfix(count,6) = TraitMu; % trait mutation is entered into column 6
    count=count+1; % This updating of the count parameter ensures that the result of each successive run of the simulation is entered into the next row of the results matrix (resfix).
    break
end
end
end


