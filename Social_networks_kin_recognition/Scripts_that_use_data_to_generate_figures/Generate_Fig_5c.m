% This figure generates data and uses this data to plot Fig. 5c (Main text). 
% This plot records the probability of social interaction, for varying 
% alpha, and for varying tag frequencies (assuming all tags in the 
% population are at equal frequency). The probability of social interaction
% is taken as an average over all migration and recombination probabilities
% within the region of parameter space where kin discrimination is favoured
% by kin selection (Equation 1 in main text is satisfied).

clearvars
close all
clc

N=30; % Deme size
b=0.3; % Benefit of cooperation
c=0.1; % Cost of cooperation
alphR = [0:0.01:1]; % This array plots a range of alpha values (partner recycling probabilities).
tagR = [ 5 10 1000000000]; % This array plots three different values of L_{max} (maximum number of tags). The inverse of this (1/L_{max}) gives tag frequency.
rR=[0:0.0025:0.5]; % This array plots a range of recombination probabilities.
mR=[0:0.002:1]; % This array plots a range of migration probabilities.
int = zeros(length(alphR),length(tagR)); % This empty matrix will be populated with results.
listres=NaN; % This will be populated with results.

for cur_tag = 1:length(tagR)

    tag = tagR(cur_tag);

for cur_alpha = 1:length(alphR)

    alpha = alphR(cur_alpha);
    
    count=1;
    
for cur_r = 1:length(rR)
    
    r = rR(cur_r);
    
for cur_m = 1:length(mR)
    
    m = mR(cur_m);
    
% gamma (below) is the probability that 2 pairs of genes across 3 
% individuals coalesce (see Supp. Info. 4b).
gamma = ((-1+m)^3*(-12*(-1+r)^3-16*m^7*(-3+N)*(-2+N)*(-1+N)^3*(-1+r)^3+2*m^8*(-3+N)*(-2+N)*(-1+N)^3*(-1+r)^3+2*N^5*(-2+r)*r^2+4*N*(-1+r)*(-2+3*r)*(-5+4*r)+N^4*r*(-26+(44-17*r)*r)+N^3*(-1+r)*(18+13*r*(-7+4*r))-N^2*(-1+r)*(45+r*(-134+73*r))+m^6*(-2+N)*(-1+N)^2*(-1+r)*(168*(-1+r)^2+N*(-221+(448-225*r)*r)+N^2*(55+56*(-2+r)*r))-m^5*(-1+N)^2*(-1+r)*(-672*(-1+r)^2+N^3*(106+r*(-225+113*r))+4*N*(299+r*(-616+311*r))-2*N^2*(321+r*(-671+338*r)))+m^4*(-1+N)*(-1+r)*(840*(-1+r)^2+N^3*(-883+(1951-989*r)*r)+N^4*(122+r*(-281+143*r))-10*N*(229+r*(-476+241*r))+N^2*(2209+r*(-4750+2417*r)))+m^2*(-1+N)*(336*(-1+r)^3-N^3*(-1+r)*(260+r*(-756+397*r))-2*N*(-1+r)*(431+r*(-952+491*r))+N^2*(-1+r)*(751+r*(-1884+989*r))+N^4*(-24+r*(122+r*(-158+57*r))))-m^3*(-1+N)*(672*(-1+r)^3-4*N^3*(-1+r)*(159+2*r*(-193+99*r))-8*N*(-1+r)*(223+r*(-476+243*r))+4*N^2*(-1+r)*(415+r*(-947+488*r))+N^4*(-80+r*(298+r*(-334+115*r))))+m*(96*(-1+r)^3-4*N^5*r*(5+2*r*(-5+2*r))+22*N^2*(-1+r)*(19+r*(-49+26*r))-4*N*(-1+r)*(83+r*(-184+95*r))-N^3*(-1+r)*(230+r*(-743+403*r))+N^4*(-48+r*(278+r*(-364+131*r))))))/((-1+(-2+m)*m*(-1+N))*(-12*(-1+r)^3-9*m^8*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+m^9*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+2*N^5*(-2+r)*r^2+4*N*(-1+r)*(-2+3*r)*(-5+4*r)+N^4*r*(-26+(44-17*r)*r)+N^3*(-1+r)*(18+13*r*(-7+4*r))-N^2*(-1+r)*(45+r*(-134+73*r))+m^7*(-2+N)*(-1+N)^2*(-1+r)*(-216*(-1+r)^2+N^3*(-7+6*r)*(-5+6*r)+N^2*(-211+430*r-216*r^2)+N*(390-786*r+394*r^2))-m^6*(-2+N)*(-1+N)^2*(-1+r)*(-504*(-1+r)^2+N^2*(-466+(987-500*r)*r)+14*N*(63+r*(-129+65*r))+N^3*(76+r*(-167+84*r)))+m^5*(-2+N)*(-1+N)^2*(-1+r)*(-756*(-1+r)^2+42*N*(30+r*(-63+32*r))-3*N^2*(211+4*r*(-119+61*r))+N^3*(98+r*(-244+125*r)))+m*(-1+N)*(-108*(-1+r)^3+4*N^5*(-2+r)*r^2-2*N^4*r*(-5+3*r)*(-6+7*r)+3*N^3*(-1+r)*(22+r*(-108+61*r))+4*N*(-1+r)*(69+r*(-159+83*r))-N^2*(-1+r)*(227+r*(-680+369*r)))-m^2*(-1+N)*(-432*(-1+r)^3+2*N^5*r*(16+r*(-32+13*r))-6*N^2*(-1+r)*(193+5*r*(-100+53*r))+12*N*(-1+r)*(99+r*(-219+113*r))+N^3*(-1+r)*(484+r*(-1607+874*r))+6*N^4*(12+r*(-78+(106-39*r)*r)))+m^3*(-1+N)*(-1008*(-1+r)^3+140*N*(-1+r)*(21+r*(-45+23*r))+N^3*(-1+r)*(1583+6*r*(-716+379*r))-2*N^2*(-1+r)*(1592+r*(-3759+1957*r))+2*N^4*(174+r*(-764+(915-323*r)*r))+2*N^5*(-12+r*(74+r*(-101+37*r))))-m^4*(-1+N)*(-1512*(-1+r)^3+140*N*(-1+r)*(33+r*(-69+35*r))-4*N^2*(-1+r)*(1339+4*r*(-742+381*r))+N^3*(-1+r)*(2932+r*(-7029+3646*r))+N^4*(754+r*(-2771+(3077-1059*r)*r))+N^5*(-72+r*(298+r*(-348+121*r))))));
% phi (below) is the probability that 2 pairs of genes across 2 individuals
% coalesce at 2 loci (see Supp. Info. 4b).
phi = -(((-1+m)^2*(-12*(-1+r)^3-9*m^8*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+m^9*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+2*N^5*(-2+r)*r^2+4*N*(-1+r)*(-2+3*r)*(-5+4*r)+N^4*r*(-26+(44-17*r)*r)+N^3*(-1+r)*(18+13*r*(-7+4*r))-N^2*(-1+r)*(45+r*(-134+73*r))+m^7*(-2+N)*(-1+N)^2*(-1+r)*(-216*(-1+r)^2+35*N^3*(-1+r)^2+N^2*(-211-212*(-2+r)*r)+N*(390-786*r+394*r^2))-m^6*(-2+N)*(-1+N)^2*(-1+r)*(-504*(-1+r)^2+N^3*(-1+r)*(-76+77*r)+N^2*(-466+(945-472*r)*r)+14*N*(63+r*(-129+65*r)))+m^5*(-2+N)*(-1+N)*(-1+r)*(756*(-1+r)^2+2*N^4*(-1+r)*(-49+52*r)+N^3*(-731-753*(-2+r)*r)-42*N*(48+r*(-99+50*r))+N^2*(1893+2*r*(-1977+998*r)))-m^4*(-1+N)*(-1512*(-1+r)^3+2*N^5*(-1+r)^2*(-36+43*r)+140*N*(-1+r)*(33+r*(-69+35*r))-4*N^2*(-1+r)*(1339+2*r*(-1439+732*r))+N^3*(-1+r)*(2932+r*(-6321+3194*r))+N^4*(754+r*(-2371+(2441-823*r)*r)))+m^3*(-1+N)*(-1008*(-1+r)^3+8*N^5*(-1+r)^2*(-3+5*r)+140*N*(-1+r)*(21+r*(-45+23*r))-2*N^2*(-1+r)*(1592+3*r*(-1203+619*r))+N^3*(-1+r)*(1583+2*r*(-1837+938*r))+N^4*(348+r*(-1168+(1253-429*r)*r)))+m^2*(-432*(-1+r)^3-8*N^6*(-1+r)^2*r+N^3*(-1+r)*(1642+15*r*(-279+146*r))+12*N*(-1+r)*(135+r*(-291+149*r))-2*N^2*(-1+r)*(1173+r*(-2748+1429*r))+N^4*(556+r*(-2097+(2362-815*r)*r))+N^5*(-72+r*(304+r*(-374+135*r))))+m*(108*(-1+r)^3-N^3*(-1+r)*(293-924*r+500*r^2)-2*N^5*r*(12+5*r*(-5+2*r))-4*N*(-1+r)*(96+r*(-213+110*r))+N^2*(-1+r)*(503+r*(-1292+685*r))+N^4*(-66+r*(360+r*(-465+167*r))))))/((-1+(-2+m)*m*(-1+N))*(-12*(-1+r)^3-9*m^8*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+m^9*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+2*N^5*(-2+r)*r^2+4*N*(-1+r)*(-2+3*r)*(-5+4*r)+N^4*r*(-26+(44-17*r)*r)+N^3*(-1+r)*(18+13*r*(-7+4*r))-N^2*(-1+r)*(45+r*(-134+73*r))+m^7*(-2+N)*(-1+N)^2*(-1+r)*(-216*(-1+r)^2+N^3*(-7+6*r)*(-5+6*r)+N^2*(-211+430*r-216*r^2)+N*(390-786*r+394*r^2))-m^6*(-2+N)*(-1+N)^2*(-1+r)*(-504*(-1+r)^2+N^2*(-466+(987-500*r)*r)+14*N*(63+r*(-129+65*r))+N^3*(76+r*(-167+84*r)))+m^5*(-2+N)*(-1+N)^2*(-1+r)*(-756*(-1+r)^2+42*N*(30+r*(-63+32*r))-3*N^2*(211+4*r*(-119+61*r))+N^3*(98+r*(-244+125*r)))+m*(-1+N)*(-108*(-1+r)^3+4*N^5*(-2+r)*r^2-2*N^4*r*(-5+3*r)*(-6+7*r)+3*N^3*(-1+r)*(22+r*(-108+61*r))+4*N*(-1+r)*(69+r*(-159+83*r))-N^2*(-1+r)*(227+r*(-680+369*r)))-m^2*(-1+N)*(-432*(-1+r)^3+2*N^5*r*(16+r*(-32+13*r))-6*N^2*(-1+r)*(193+5*r*(-100+53*r))+12*N*(-1+r)*(99+r*(-219+113*r))+N^3*(-1+r)*(484+r*(-1607+874*r))+6*N^4*(12+r*(-78+(106-39*r)*r)))+m^3*(-1+N)*(-1008*(-1+r)^3+140*N*(-1+r)*(21+r*(-45+23*r))+N^3*(-1+r)*(1583+6*r*(-716+379*r))-2*N^2*(-1+r)*(1592+r*(-3759+1957*r))+2*N^4*(174+r*(-764+(915-323*r)*r))+2*N^5*(-12+r*(74+r*(-101+37*r))))-m^4*(-1+N)*(-1512*(-1+r)^3+140*N*(-1+r)*(33+r*(-69+35*r))-4*N^2*(-1+r)*(1339+4*r*(-742+381*r))+N^3*(-1+r)*(2932+r*(-7029+3646*r))+N^4*(754+r*(-2771+(3077-1059*r)*r))+N^5*(-72+r*(298+r*(-348+121*r)))))));
% F (below) is the probability that 1 pair of genes across 2 individuals 
% coalesces (see Supp. Info. 4b).
F = (1-m)^2 / ( (1-m)^2 + N*(1-(1-m)^2)) ;

denom = (1 - alpha .* (1-1/tag) .* (1-F)) ; % This variable (short for denominator) is defined as it is used in the following two calculations.
M1 = ((F+(1-F).* 1/tag) ./ denom); % This array gives the per-generation probability of obtaining a social interaction (assuming all tags are equal frequency).
M2 = ( phi + (F-phi).* 1/tag)   ./ denom; % This array gives the per-generation probability of receiving help (assuming all tags are equal frequency).

% Satisfaction of the below condition indicates that kin discrimination is
% favoured by kin selection (Hamilton's Rule).
if b .*  (M2/M1)  - c    - (1-m) .* (b-c)  .* ( 1/N  + (M2/M1)/N + ((N-2)/N) .* (gamma  + (F-gamma).* ( (1/tag + (1-1/tag).*alpha.*M1)) )./M1) > 0            
listres(count) = M1; % This array records per-generation interaction probability. Values are only recorded if kin discrimination is favoured by kin selection (above condition satisfied).
count=count+1; % This counting is used to ensure that all M1 values are recorded across the full range of m & r within which kin discrimination is favoured by kin selection.
end
      
end
end

int(cur_alpha,cur_tag) = nanmean(listres); % This records, for each alpha & tag combination, the mean social interaction probability, where the mean is taken across all r & m values where kin recognition is favoured by kin selection (Equation 1 in main text is satisfied).

end
end

% The for loop below plots the average social interaction probability for
% different alpha values (X axis) and different tag frequencies
% (1/L_{max}; solid lines). 
for i = 1:numel(tagR)
plot(alphR,int(:,i),'Color','k','LineWidth',2)
hold on
end
hold off
box off
set(gcf,'color','white')
ylim([0 1])