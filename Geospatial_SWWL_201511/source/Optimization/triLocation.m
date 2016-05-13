function pose=triLocation(Dij, Posit, BeaconList)

% GivenMat= distMat(BeaconList, BeaconList);
% nBeacon= length(BeaconList);
% xMat= repmat(Posit(:,1), 1, nBeacon);
% yMat= repmat(Posit(:,2), 1, nBeacon);
% 
% estDist= sqrt((xMat-xMat').^2 + (yMat-yMat').^2);
% 
% nzList= (GivenMat>0);
% aveTau= mean( estDist(nzList)./GivenMat(nzList));
% 
% disp(aveTau)

% Dij= max(distMat, distMat');

AdjMat= double( Dij> 0);

Nodes= size(AdjMat, 1);

LogiBeacon= zeros(1, Nodes);
LogiBeacon(BeaconList)= 1;
LogiBeacon=logical(LogiBeacon);

NeibList= (sum( AdjMat(:, LogiBeacon),2) >=3)';
AvailList= (NeibList- LogiBeacon)>0;

pose= zeros(Nodes, 6);% store final results;
pose(LogiBeacon, 1:2)= Posit;

DoneBeacon=LogiBeacon;
% plot(Posit(LogiBeacon, 1), Posit(LogiBeacon, 2), 'x'); hold on;
while nnz(AvailList)
    
    numAvail= nnz(AvailList);
    idAvail= find(AvailList);
    for ID=1:numAvail
        Target=idAvail(ID);
        BeaconTar= (LogiBeacon.*AdjMat(Target, :) > 0);
        [estPosit, valid] = localize2d( Dij(BeaconTar, Target),...
            pose(BeaconTar, :));
        
        pose(Target, :)= estPosit;
        PosiMat= repmat(estPosit, nnz(BeaconTar), 1);
        estDist= sqrt(sum((PosiMat-pose(BeaconTar, :)).^2,2));
        
        errList= abs(estDist- Dij(BeaconTar, Target));
        if mean(errList(1:3) )<0.08
            LogiBeacon(Target)= true;
        end
        
    end
    
    DoneBeacon(AvailList)=true;
    
    NeibList= (sum( AdjMat(:, LogiBeacon),2) >=3)';
    AvailList= (NeibList- DoneBeacon)>0;
    
    if nnz(AvailList)==0
        LogiBeacon=DoneBeacon;
        NeibList= (sum( AdjMat(:, LogiBeacon),2) >=3)';
        AvailList= (NeibList- LogiBeacon)>0;
    end
end

end
%