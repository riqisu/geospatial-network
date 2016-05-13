#!/bin/bash

####################################################################
#                     Basic Code Description                       #
# ---------------------------------------------------------------- #
# OS: Linux Ubuntu 14.04 with gcc version 4.9.2, and 			   #
#	Matlab 8.4.0.150421 (R2014b)                                   #
# Author: Riqi Su                                                  #
# Last updated: Thu Nov 19 2015	                                   #

#define matlab path; replace with your own path;
matlab='/usr/local/MATLAB/R2014b/bin/matlab'

#define parameter sets
Stamp=28
Degree=5
Nodes=30
DataDef=500
NmDef='[0.60]'
EpsDef='[1e-3]'
RepTimeDef=1

# compile source codes and generate networks and time series;
if [ ! -e BANet.cpp ]
then
	echo 'Missing BANet.cpp !'
	exit -1
fi

if [ ! -e DDEv3.cpp ]
then
	echo 'Missing DDEv3.cpp !'
	exit -1
fi

g++ BANet.cpp -o BANet.x 
g++ DDEv3.cpp -o DDEv3.x -O3

./BANet.x ${Nodes} ${Degree} $Stamp
./DDEv3.x BAN${Nodes}k${Degree}Stam$Stamp.txt $DataDef

# Generate matlab code to reconstruct geospatial networks.
tempFile='matlab_BANet.m'

# Generate matlab code
cat > ${tempFile} << EOF
clear all; close all; clc;
Stamp=$Stamp; Nodes=$Nodes;  Deg=$Degree;
RepTime = $RepTimeDef;
%specify locations for sensing package;
path(path, './Optimization');

NmRatio= ${NmDef};
EpsList=${EpsDef};
EpsNum= size( EpsList,2);

%load time series and network files (used to evaluate reconstruction)
outputName= ['Node' num2str(Nodes) 'ID' num2str(Stamp) 'NeedRmEps1e-3'];
NetworkName=['BAN' num2str(Nodes) 'k' num2str(Deg) 'Stam' ...
    num2str(Stamp)];
load( ['TimeSer' NetworkName, '.mat'], 'Xt', 'Yt', 'Zt', 'dXdt',  'dYdt', 'dZdt', ...
    'RK_H', 'iniWmatrix', 'DelayMat', 'Speed');

Posit= load(['PosN' num2str(Nodes) 'k' num2str(Deg) 'Stam' ...
    num2str(Stamp) '.txt']);

%set the threshold for nzWei terms;
minCoup = min( iniWmatrix( abs(iniWmatrix) > 1e-8)) / 2;
minDelay= min( DelayMat(abs(DelayMat) > 1e-8)) * minCoup / 2;

Start =100;
End = size( Xt, 1);

Xt = Xt(Start:End, :);
Yt = Yt(Start:End , :);
Zt = Zt(Start:End , :);

dXdt = dXdt(Start:End , :);
dYdt = dYdt(Start:End , :);
dZdt = dZdt(Start:End , :);

Wmatrix = iniWmatrix;

maxPow = 3; %max expand Taylor order;
%%%%%%%%%%%%%%%%%%%%%%%%%%
term = 0;
for powSum =0:maxPow
    for k =0: powSum
        for j= 0: powSum-k
            term = term +1;
            xPower(term)= powSum - k - j;
            yPower(term)= j;
            zPower(term)= k;
        end
    end
end

Nodes = size(dXdt,2);
xOffset = term;
yOffset = term +  1 * (Nodes -1);
zOffset = term +  2 * (Nodes -1);

dxOffset = term + 3 * (Nodes-1);
dyOffset = term + 4 * (Nodes-1);
dzOffset = term + 5 * (Nodes-1);

xIndex = ( xPower== 1) .* ( yPower == 0) .* ( zPower == 0);
xTerm = find( xIndex == 1);
yIndex = ( xPower== 0) .* ( yPower == 1) .* ( zPower == 0);
yTerm = find( yIndex == 1);
zIndex = ( xPower== 0) .* ( yPower == 0) .* ( zPower == 1);
zTerm = find( zIndex == 1);
xyIndex = ( xPower== 1) .* ( yPower == 1) .* ( zPower == 0);
xyTerm = find( xyIndex == 1);
xzIndex = ( xPower== 1) .* ( yPower == 0) .* ( zPower == 1);
xzTerm = find( xzIndex == 1);

% required sampling;
NeedPropList = round( NmRatio * ( term + 6*(Nodes-1)));
PropNum = size(NeedPropList,2);


% allT is max available time series;
[allT,~] = size(dXdt);

nzCoupMat = zeros( EpsNum, Nodes* RepTime);
nzDelayMat = zeros( EpsNum, Nodes* RepTime);

zCoupMat = zeros( EpsNum, Nodes* RepTime);
zDelayMat = zeros( EpsNum, Nodes* RepTime);

StoreCoup = zeros( RepTime, EpsNum);
StoreDelay = zeros( RepTime, EpsNum);

zStoreCoup = zeros( RepTime, EpsNum);
zStoreDelay = zeros( RepTime, EpsNum);

usedT = 0;

% save Matrix
CoupMat_Pred = zeros(Nodes,Nodes-1);
DelayMat_Pred = zeros(Nodes,Nodes-1);

for r = 1: RepTime
    %         r =1;
    for p = 1:PropNum
        %for p = 15:15
        NeedData = NeedPropList(p);
        %         NeedData = 120;
        
        tempL = randperm(allT);
        Select = tempL(1:NeedData);
        
        partXt = Xt( Select,:);
        partYt = Yt( Select,:);
        partZt = Zt( Select,:);
        
        partDX = dXdt( Select,:);
        partDY = dYdt( Select,:);
        partDZ = dZdt( Select,:);
        for epsID=1:EpsNum
            for NodeID = 1:Nodes
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Amatrix = DataProcessDL( NodeID, partXt, partYt, partZt,partDX, ...
                		partDY, partDZ, xPower, yPower, zPower);
                zerosList = zeros(NeedData, 1);
                [Xpred, NeedTime, redo] = SensingP2( Amatrix, partDX(:,NodeID), ...
                		EpsList(epsID));
                
                usedT = usedT + NeedTime;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                temp = Wmatrix(NodeID,:);
                temp(NodeID)= [];
                Xreal_Coup = temp;
                
                tempDey = DelayMat(NodeID,:);
                tempDey(NodeID) =[];
                Xreal_Delay = - tempDey .* temp;
                
                nzCoup = find( abs(Xreal_Coup) > minCoup);
                nzDelay= find( abs(Xreal_Delay) > minDelay);
                
                zCoup= 1:(Nodes-1); zCoup(nzCoup)=[];
                zDelay=1:(Nodes-1); zDelay(nzDelay)=[];
                
                Xpred_Coup = Xpred(zOffset+1: zOffset+ Nodes-1)';
                Xpred_Delay= Xpred(dzOffset+1: dzOffset+ Nodes-1)';
                
                CoupMat_Pred(NodeID, :) = Xpred_Coup;
                DelayMat_Pred(NodeID,:) = Xpred_Delay;
                
                
                nzCoupMat( p, (r-1)*Nodes + NodeID ) = ...
                    mean( abs(  (Xreal_Coup( nzCoup) - Xpred_Coup( nzCoup) ) ./ Xreal_Coup( nzCoup) ) ) ;
                nzDelayMat(p, (r-1)*Nodes + NodeID )= ...
                    mean( abs(  (Xreal_Delay( nzDelay) - Xpred_Delay( nzDelay) )./ Xreal_Delay( nzDelay)));
                zCoupMat(p, (r-1)*Nodes + NodeID ) = ...
                    mean( abs(  (Xreal_Coup( zCoup) - Xpred_Coup( zCoup) )));
                zDelayMat(p, (r-1)*Nodes + NodeID )= ...
                    mean( abs(  (Xreal_Delay( zDelay) - Xpred_Delay( zDelay))));
                
            end
            StoreCoup(r, p)= mean( nzCoupMat( p, (r-1)*Nodes+1: r*Nodes));
            StoreDelay(r, p)= mean( nzDelayMat( p, (r-1)*Nodes+1: r*Nodes));
            zStoreCoup(r, p)= mean( zCoupMat( p, (r-1)*Nodes+1: r*Nodes));
            zStoreDelay(r, p)= mean( zDelayMat( p, (r-1)*Nodes+1: r*Nodes));
        end
    end
end


% reconstruct adjaciancy matrix from sensing results;
a = diag(ones(1,Nodes));
Index = find(a==0);

MatCoup = zeros(Nodes, Nodes);
MatCoup(Index) = CoupMat_Pred';
MatCoup = MatCoup';

ListDeley= find( MatCoup > minCoup);

% calculate the estimated delays: tau_{ij}= bij/wij;
MatDelay = zeros(Nodes, Nodes);
MatDelay(Index) = DelayMat_Pred';
MatDelay = MatDelay';

MatDelay(ListDeley)= -MatDelay(ListDeley)./MatCoup(ListDeley);


minD=min( DelayMat(abs(DelayMat) > 1e-8))/2;
maxD=max( DelayMat(abs(DelayMat) > 1e-8))*2;

% set the distance matrix 
distMatAA=zeros(Nodes, Nodes);
indexMat= find(( MatDelay >minD).*(MatDelay <maxD));
distMatAA(indexMat)= MatDelay(indexMat);

distMat= 0.5.*(distMatAA+ distMatAA')*100/Speed;

% locating node position from delay matrix;
beaList=1:4;
estPosit= triLocation( distMat, Posit(beaList,:), beaList);
% estPosit= triLocation( DelayMat*100/Speed, Posit(1:12,:), 1:12);

preList=[1:Nodes];
preList(beaList)=[];

plot( Posit(preList, 1), Posit(preList, 2), 'ko', 'MarkerSize', 10); hold on;
plot( Posit(beaList, 1), Posit(beaList, 2), 'rs','MarkerSize', 10); hold on;

for i=1:Nodes
    plot( [Posit(i, 1) estPosit(i,1)]', ...
        [Posit(i, 2) estPosit(i,2)]', ...
        '-', 'Color', [0.6 0.6 0.6], 'lineWidth', 2); hold on;
end
hold off;
xlim([0 1]); ylim([ 0 1]);
set(gca, 'xTick', [], 'yTick', [])

set(1, 'PaperPositionMode', 'auto');
print( 1, '-painters', 'fig1Triangle.pdf', '-depsc');

% save results;

save( [ outputName '.mat'],'Nodes', 'usedT', 'PropNum', 'NmRatio','NeedPropList', ...
			'RepTime', 'Wmatrix', 'DelayMat','Nodes', 'maxPow');
save( [ outputName '.mat'],'-append', 'StoreCoup', 'StoreDelay','nzCoupMat', 'nzDelayMat');
save( [ outputName '.mat'],'-append', 'zStoreCoup', 'zStoreDelay','zCoupMat', 'zDelayMat');
save( [ outputName '.mat'],'-append', 'MatCoup', 'MatDelay', 'Posit', 'estPosit');

EOF

$matlab  -r "run $tempFile;"
