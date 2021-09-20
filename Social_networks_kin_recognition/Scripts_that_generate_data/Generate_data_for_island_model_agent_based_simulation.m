% INFINITE ISLAND MODEL AGENT BASED SIMULATION: STRONGER SELECTION & 
% SMALLER DEMES (SUPP. INFO. 4D).

% This script generates data from the agent based simulation of the island
% model described in Supp. Info. 4d. The data generated from this script 
% is saved in the 'Saved_data' folder as 'Island_strong_sel_b=3c.mat' and 
% 'Island_strong_sel_b=9c.mat' and was used to generate Supp. Fig. 13.

clearvars
clc
close all

% PARAMETERS TO BE SET AS INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mR=[0:0.1:1]; % Range of migration rates
rR=[0:0.05:0.5]; % Range of recombination rates
TraitMu = 0.005; % Trait mutation rate
TagMu = 0.005; % Tag mutation rate (per-generation probability that a given tag  mutates into a randomly chosen tag (this includes mutation back to itself).
alpha = 1; % Partner-reassociation probability
L=10; % Maximum segregating tags (N.B. This is labelled as 'L_{max}' in the main text, and is often labelled as 'tag' in other scripts).
T=5000; % Generations
c=0.3; % Cost of cooperation
b=2.7; % Benefit of cooperation
N=3; % Deme size 
D=167; % Number of demes in the population.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nInd = N*D;% Overall population size.

for cur_r=1:numel(rR)
    
    r = rR(cur_r);
 
for cur_m=1:numel(mR)
    
    m = mR(cur_m);
    
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
 pop(:,7) = 1 + ( b .* (pop(:,6)./(N-1)) - c .* pop(:,3) .* ((pop(:,5)+pop(:,6))./(N-1)) ) ./ (1 - alpha .* ((N-1-pop(:,5)-pop(:,6))./(N-1)) )  ;
 
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

end

XXX = sum(help,1);  % This array gives the population cooperator freqeuncy across generations.
resalt(cur_r,cur_m) = mean(XXX(round(T/4):T)); % This array gives the mean population cooperator freqeuncy across generations, discounting the first T/4 generations (to allow time for the equilibrium to be approached).

% The following two lines calculate the tag diversity (number of tags
% maintained) at evolutionary equilibrium. (The first T/4 generations are
% discounted in this calculation to allow time for the equilibrium to be 
% approached).
AAA = (help + def).^2 ;
resmark(cur_r,cur_m) = 1/sum(mean(AAA(:,round(T/4):T),2));

    %end
end
clear  pop 
end

% After running this script (for a given input of parameters), we obtain
% two summary statistics (matrices), 'resalt' and 'resmark'. These matrices
% respectively give the equilibrium cooperator frequency and number of tags
% maintained, for different r & m combinations.