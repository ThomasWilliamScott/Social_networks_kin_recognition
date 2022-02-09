% This script generates the 'relatedness as regression' illustrative plots 
% (Fig. S3 & S5).

close all
clearvars
clc

% Parameter values are chosen for illustration.

% PARAMETER VALUES %%%%%%%%%%%%%%
m=0.1; % migration rate
r=0.1; % recombination rate
N=10; % deme size
tag=20; % maximum tag number (denoted by L_{max} in the text)
alpha=0.4; % encounter parameter
xi=1/20; % tag popualtion frequency (assumed same across tags, for simplicity)
pii=0.3; % cooperator proportion (assumed same across all tags, for simplicity). It is labelled as pii rather than pi, as we have labelled it elsewhere, because we will need to use the value of the mathematical constant pi (i.e. 3.14159...) in the following calculations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% coalescence probabilities %%%%%%%%%
gamma = ((-1+m)^3*(-12*(-1+r)^3-16*m^7*(-3+N)*(-2+N)*(-1+N)^3*(-1+r)^3+2*m^8*(-3+N)*(-2+N)*(-1+N)^3*(-1+r)^3+2*N^5*(-2+r)*r^2+4*N*(-1+r)*(-2+3*r)*(-5+4*r)+N^4*r*(-26+(44-17*r)*r)+N^3*(-1+r)*(18+13*r*(-7+4*r))-N^2*(-1+r)*(45+r*(-134+73*r))+m^6*(-2+N)*(-1+N)^2*(-1+r)*(168*(-1+r)^2+N*(-221+(448-225*r)*r)+N^2*(55+56*(-2+r)*r))-m^5*(-1+N)^2*(-1+r)*(-672*(-1+r)^2+N^3*(106+r*(-225+113*r))+4*N*(299+r*(-616+311*r))-2*N^2*(321+r*(-671+338*r)))+m^4*(-1+N)*(-1+r)*(840*(-1+r)^2+N^3*(-883+(1951-989*r)*r)+N^4*(122+r*(-281+143*r))-10*N*(229+r*(-476+241*r))+N^2*(2209+r*(-4750+2417*r)))+m^2*(-1+N)*(336*(-1+r)^3-N^3*(-1+r)*(260+r*(-756+397*r))-2*N*(-1+r)*(431+r*(-952+491*r))+N^2*(-1+r)*(751+r*(-1884+989*r))+N^4*(-24+r*(122+r*(-158+57*r))))-m^3*(-1+N)*(672*(-1+r)^3-4*N^3*(-1+r)*(159+2*r*(-193+99*r))-8*N*(-1+r)*(223+r*(-476+243*r))+4*N^2*(-1+r)*(415+r*(-947+488*r))+N^4*(-80+r*(298+r*(-334+115*r))))+m*(96*(-1+r)^3-4*N^5*r*(5+2*r*(-5+2*r))+22*N^2*(-1+r)*(19+r*(-49+26*r))-4*N*(-1+r)*(83+r*(-184+95*r))-N^3*(-1+r)*(230+r*(-743+403*r))+N^4*(-48+r*(278+r*(-364+131*r))))))/((-1+(-2+m)*m*(-1+N))*(-12*(-1+r)^3-9*m^8*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+m^9*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+2*N^5*(-2+r)*r^2+4*N*(-1+r)*(-2+3*r)*(-5+4*r)+N^4*r*(-26+(44-17*r)*r)+N^3*(-1+r)*(18+13*r*(-7+4*r))-N^2*(-1+r)*(45+r*(-134+73*r))+m^7*(-2+N)*(-1+N)^2*(-1+r)*(-216*(-1+r)^2+N^3*(-7+6*r)*(-5+6*r)+N^2*(-211+430*r-216*r^2)+N*(390-786*r+394*r^2))-m^6*(-2+N)*(-1+N)^2*(-1+r)*(-504*(-1+r)^2+N^2*(-466+(987-500*r)*r)+14*N*(63+r*(-129+65*r))+N^3*(76+r*(-167+84*r)))+m^5*(-2+N)*(-1+N)^2*(-1+r)*(-756*(-1+r)^2+42*N*(30+r*(-63+32*r))-3*N^2*(211+4*r*(-119+61*r))+N^3*(98+r*(-244+125*r)))+m*(-1+N)*(-108*(-1+r)^3+4*N^5*(-2+r)*r^2-2*N^4*r*(-5+3*r)*(-6+7*r)+3*N^3*(-1+r)*(22+r*(-108+61*r))+4*N*(-1+r)*(69+r*(-159+83*r))-N^2*(-1+r)*(227+r*(-680+369*r)))-m^2*(-1+N)*(-432*(-1+r)^3+2*N^5*r*(16+r*(-32+13*r))-6*N^2*(-1+r)*(193+5*r*(-100+53*r))+12*N*(-1+r)*(99+r*(-219+113*r))+N^3*(-1+r)*(484+r*(-1607+874*r))+6*N^4*(12+r*(-78+(106-39*r)*r)))+m^3*(-1+N)*(-1008*(-1+r)^3+140*N*(-1+r)*(21+r*(-45+23*r))+N^3*(-1+r)*(1583+6*r*(-716+379*r))-2*N^2*(-1+r)*(1592+r*(-3759+1957*r))+2*N^4*(174+r*(-764+(915-323*r)*r))+2*N^5*(-12+r*(74+r*(-101+37*r))))-m^4*(-1+N)*(-1512*(-1+r)^3+140*N*(-1+r)*(33+r*(-69+35*r))-4*N^2*(-1+r)*(1339+4*r*(-742+381*r))+N^3*(-1+r)*(2932+r*(-7029+3646*r))+N^4*(754+r*(-2771+(3077-1059*r)*r))+N^5*(-72+r*(298+r*(-348+121*r))))));
phi = -(((-1+m)^2*(-12*(-1+r)^3-9*m^8*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+m^9*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+2*N^5*(-2+r)*r^2+4*N*(-1+r)*(-2+3*r)*(-5+4*r)+N^4*r*(-26+(44-17*r)*r)+N^3*(-1+r)*(18+13*r*(-7+4*r))-N^2*(-1+r)*(45+r*(-134+73*r))+m^7*(-2+N)*(-1+N)^2*(-1+r)*(-216*(-1+r)^2+35*N^3*(-1+r)^2+N^2*(-211-212*(-2+r)*r)+N*(390-786*r+394*r^2))-m^6*(-2+N)*(-1+N)^2*(-1+r)*(-504*(-1+r)^2+N^3*(-1+r)*(-76+77*r)+N^2*(-466+(945-472*r)*r)+14*N*(63+r*(-129+65*r)))+m^5*(-2+N)*(-1+N)*(-1+r)*(756*(-1+r)^2+2*N^4*(-1+r)*(-49+52*r)+N^3*(-731-753*(-2+r)*r)-42*N*(48+r*(-99+50*r))+N^2*(1893+2*r*(-1977+998*r)))-m^4*(-1+N)*(-1512*(-1+r)^3+2*N^5*(-1+r)^2*(-36+43*r)+140*N*(-1+r)*(33+r*(-69+35*r))-4*N^2*(-1+r)*(1339+2*r*(-1439+732*r))+N^3*(-1+r)*(2932+r*(-6321+3194*r))+N^4*(754+r*(-2371+(2441-823*r)*r)))+m^3*(-1+N)*(-1008*(-1+r)^3+8*N^5*(-1+r)^2*(-3+5*r)+140*N*(-1+r)*(21+r*(-45+23*r))-2*N^2*(-1+r)*(1592+3*r*(-1203+619*r))+N^3*(-1+r)*(1583+2*r*(-1837+938*r))+N^4*(348+r*(-1168+(1253-429*r)*r)))+m^2*(-432*(-1+r)^3-8*N^6*(-1+r)^2*r+N^3*(-1+r)*(1642+15*r*(-279+146*r))+12*N*(-1+r)*(135+r*(-291+149*r))-2*N^2*(-1+r)*(1173+r*(-2748+1429*r))+N^4*(556+r*(-2097+(2362-815*r)*r))+N^5*(-72+r*(304+r*(-374+135*r))))+m*(108*(-1+r)^3-N^3*(-1+r)*(293-924*r+500*r^2)-2*N^5*r*(12+5*r*(-5+2*r))-4*N*(-1+r)*(96+r*(-213+110*r))+N^2*(-1+r)*(503+r*(-1292+685*r))+N^4*(-66+r*(360+r*(-465+167*r))))))/((-1+(-2+m)*m*(-1+N))*(-12*(-1+r)^3-9*m^8*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+m^9*(-3+N)*(-2+N)^2*(-1+N)^3*(-1+r)^3+2*N^5*(-2+r)*r^2+4*N*(-1+r)*(-2+3*r)*(-5+4*r)+N^4*r*(-26+(44-17*r)*r)+N^3*(-1+r)*(18+13*r*(-7+4*r))-N^2*(-1+r)*(45+r*(-134+73*r))+m^7*(-2+N)*(-1+N)^2*(-1+r)*(-216*(-1+r)^2+N^3*(-7+6*r)*(-5+6*r)+N^2*(-211+430*r-216*r^2)+N*(390-786*r+394*r^2))-m^6*(-2+N)*(-1+N)^2*(-1+r)*(-504*(-1+r)^2+N^2*(-466+(987-500*r)*r)+14*N*(63+r*(-129+65*r))+N^3*(76+r*(-167+84*r)))+m^5*(-2+N)*(-1+N)^2*(-1+r)*(-756*(-1+r)^2+42*N*(30+r*(-63+32*r))-3*N^2*(211+4*r*(-119+61*r))+N^3*(98+r*(-244+125*r)))+m*(-1+N)*(-108*(-1+r)^3+4*N^5*(-2+r)*r^2-2*N^4*r*(-5+3*r)*(-6+7*r)+3*N^3*(-1+r)*(22+r*(-108+61*r))+4*N*(-1+r)*(69+r*(-159+83*r))-N^2*(-1+r)*(227+r*(-680+369*r)))-m^2*(-1+N)*(-432*(-1+r)^3+2*N^5*r*(16+r*(-32+13*r))-6*N^2*(-1+r)*(193+5*r*(-100+53*r))+12*N*(-1+r)*(99+r*(-219+113*r))+N^3*(-1+r)*(484+r*(-1607+874*r))+6*N^4*(12+r*(-78+(106-39*r)*r)))+m^3*(-1+N)*(-1008*(-1+r)^3+140*N*(-1+r)*(21+r*(-45+23*r))+N^3*(-1+r)*(1583+6*r*(-716+379*r))-2*N^2*(-1+r)*(1592+r*(-3759+1957*r))+2*N^4*(174+r*(-764+(915-323*r)*r))+2*N^5*(-12+r*(74+r*(-101+37*r))))-m^4*(-1+N)*(-1512*(-1+r)^3+140*N*(-1+r)*(33+r*(-69+35*r))-4*N^2*(-1+r)*(1339+4*r*(-742+381*r))+N^3*(-1+r)*(2932+r*(-7029+3646*r))+N^4*(754+r*(-2771+(3077-1059*r)*r))+N^5*(-72+r*(298+r*(-348+121*r)))))));
F = (1-m)^2 / ( (1-m)^2 + N*(1-(1-m)^2)) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% relatedness between social interactants: 
Rtag = (phi + (F-phi)*(1/tag)) / (phi + (F-phi)*(1+1/tag)+(1-2*F+phi)*(1/tag)) 

% relatedness between competitors is calculated in the following 2 lines.
partner = (1/tag + ((tag-1)/tag)*alpha*((F+(1-F)*(1/tag))/(1-alpha*((tag-1)/tag)*(1-F))));
Rdisplaced = (1-m)*(1/N + Rtag/N + ((N-2)/N) * ( ( gamma + (F-gamma) * partner) / ( gamma + (F-gamma)*(1+partner) + (1-2*F+gamma)*partner ) ))

denom = (1 - alpha .* (1-xi) .* (1-F));
Mint = (F+(1-F).* xi) ./ denom; % probability per-generation of interacting
Mhelped1 = (phi+(F-phi)*pii+(F-phi)*xi+(1-2*F+phi)*pii*xi)/denom; % probability per-generation of being helped, for a cooperator.
Mhelped0 = ((F-phi)*pii+(1-2*F+phi)*pii *xi)/denom; % probability per-generation of being helped, for a defector.

pointAx = Mhelped0 / Mint; % This gives the social partner's expected genic value when the focal individual's genic value is zero (defector).
pointBx = Mhelped1 / Mint; % This gives the social partner's expected genic value when the focal individual's genic value is one (helper).

% The following 4 lines of code will be used to generate circles with areas
% corresponding to the frequencies of different types of interaction.
pointA = 2*sqrt(pointAx / pi); 
pointB = 2*sqrt(pointBx / pi);
pointC = 2*sqrt((1-pointAx) / pi);
pointD = 2*sqrt((1-pointBx) / pi);

% The following 4 lines of code store different colours.
col1 = [0, 0.4470, 0.7410];
col2 = 	[0.8500, 0.3250, 0.0980];
col3 = [0.9290, 0.6940, 0.1250, 0.3]	;
col4 = [0.9290, 0.6940, 0.1250, 1]	;

% The following lines of code plot and format Fig. S3

hold on

plot(0,0,'r*','MarkerSize',20,'Marker','x','Color',col3,'LineWidth',5);
plot(0,1,'r*','MarkerSize',20,'Marker','x','Color',col3,'LineWidth',5);
plot(1,0,'r*','MarkerSize',20,'Marker','x','Color',col3,'LineWidth',5);
plot(1,1,'r*','MarkerSize',20,'Marker','x','Color',col3,'LineWidth',5);

pos = [0 - 0.5*(pointC) 0 - 0.5*(pointC) (pointC) (pointC)];
rectangle('Position',pos,'Curvature',[1 1],'FaceColor',col3,'EdgeColor',col4)
axis equal
   
pos = [0-0.5*pointA 1-0.5*pointA (pointA) (pointA)];
rectangle('Position',pos,'Curvature',[1 1],'FaceColor',col3,'EdgeColor',col4)
axis equal
   
pos = [1-0.5*(pointD) 0-0.5*(pointD) (pointD) (pointD)];
rectangle('Position',pos,'Curvature',[1 1],'FaceColor',col3,'EdgeColor',col4)
axis equal
    
pos = [1-0.5*pointB 1-0.5*pointB (pointB) (pointB)];
rectangle('Position',pos,'Curvature',[1 1],'FaceColor',col3,'EdgeColor',col4)
axis equal
    
plot([0,1],[pointAx,pointBx],'Color',col1,'LineWidth',2,'LineStyle','--')

H=gca;
H.LineWidth=1.5;
set(gca,'FontWeight','bold')
xlim([0,1])
ylim([0,1])

set(gcf,'color','white')
 
 
 
 
 
figure % new panel to generate Fig. S5
 
Mdisplaced1 = m*pii + (1-m)*(( (gamma / Mint) + ((F-gamma)*pii) / Mint + (F-gamma) * ( (xi+(1-xi)*alpha*Mint) / Mint) + (1-2*F+gamma)*( (xi+(1-xi)*alpha*Mint)/Mint)*pii)*((N-2)/N) + (Mhelped1/Mint)*(1/N) + (1/N)); 
Mdisplaced0 = m*pii + (1-m)*(( ((F-gamma)*pii) / Mint + (1-2*F+gamma)*( (xi+(1-xi)*alpha*Mint)/Mint)*pii)*((N-2)/N) + (Mhelped0/Mint)*(1/N)); 
pointAx = Mdisplaced0; % This gives the competitor's expected genic value when the focal individual's genic value is zero (defector).
pointBx = Mdisplaced1; % This gives the competitor's expected genic value when the focal individual's genic value is one (helper).

% The following 4 lines of code will be used to generate circles with areas
% corresponding to the frequencies of different types of interaction.
pointA = 2*sqrt(pointAx / pi);
pointB = 2*sqrt(pointBx / pi);
pointC = 2*sqrt((1-pointAx) / pi);
pointD = 2*sqrt((1-pointBx) / pi);

% The following 4 lines of code store different colours.
col1 = [0, 0.4470, 0.7410];
col2 = 	[0.8500, 0.3250, 0.0980];
col3 = [0.9290, 0.6940, 0.1250, 0.3]	;
col4 = [0.9290, 0.6940, 0.1250, 1]	;

% The following lines of code plot and format Fig. S5
 
hold on

plot(0,0,'r*','MarkerSize',20,'Marker','x','Color',col3,'LineWidth',5);
plot(0,1,'r*','MarkerSize',20,'Marker','x','Color',col3,'LineWidth',5);
plot(1,0,'r*','MarkerSize',20,'Marker','x','Color',col3,'LineWidth',5);
plot(1,1,'r*','MarkerSize',20,'Marker','x','Color',col3,'LineWidth',5);

pos = [0 - 0.5*(pointC) 0 - 0.5*(pointC) (pointC) (pointC)];
rectangle('Position',pos,'Curvature',[1 1],'FaceColor',col3,'EdgeColor',col4)
axis equal
   
pos = [0-0.5*pointA 1-0.5*pointA (pointA) (pointA)];
rectangle('Position',pos,'Curvature',[1 1],'FaceColor',col3,'EdgeColor',col4)
axis equal
   
pos = [1-0.5*(pointD) 0-0.5*(pointD) (pointD) (pointD)];
rectangle('Position',pos,'Curvature',[1 1],'FaceColor',col3,'EdgeColor',col4)
axis equal
    
pos = [1-0.5*pointB 1-0.5*pointB (pointB) (pointB)];
rectangle('Position',pos,'Curvature',[1 1],'FaceColor',col3,'EdgeColor',col4)
axis equal
    
plot([0,1],[pointAx,pointBx],'Color',col1,'LineWidth',2,'LineStyle','--')

H=gca;
H.LineWidth=1.5;
set(gca,'FontWeight','bold')
xlim([0,1])
ylim([0,1])

set(gcf,'color','white')
