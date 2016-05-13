%Process input time series data, get base Amatrix and the nomalize row revNormA;

function  Amatrix = DataProcessSP( NodeID, NodeID2, Xt, Yt, Zt, xPower, yPower, zPower)

% Prepare A matrix;
count = size( xPower, 2);

useT= size( Xt, 1);

TserMatX  = repmat( Xt(:, NodeID), 1 , count);
TserMatY  = repmat( Yt(:, NodeID), 1 , count);
TserMatZ  = repmat( Zt(:, NodeID), 1, count);

xPowerMat  = repmat( xPower,  useT, 1);
yPowerMat =  repmat( yPower, useT, 1);
zPowerMat  = repmat( zPower,  useT, 1);
% dxPowerMat = repmat( dxPower, useT, 1);
% dzPowerMat = repmat( dzPower, useT, 1);

Amatrix_temp = (TserMatX.^xPowerMat) .* ( TserMatY .^ yPowerMat) .* ( TserMatZ .^ zPowerMat);

TserMatX  = repmat( Xt(:, NodeID2), 1 , count);
TserMatY  = repmat( Yt(:, NodeID2), 1 , count);
TserMatZ  = repmat( Zt(:, NodeID2), 1, count);

Amatrix_temp2 = (TserMatX.^xPowerMat) .* ( TserMatY .^ yPowerMat) .* ( TserMatZ .^ zPowerMat);

Xt_temp= Xt;
Xt_temp(:, NodeID)=[];

Yt_temp= Yt;
Yt_temp(:, NodeID)=[];

Zt_temp= Zt;
Zt_temp(:, NodeID)=[];

Amatrix = [Amatrix_temp Amatrix_temp2 Xt_temp Yt_temp Zt_temp];

end
