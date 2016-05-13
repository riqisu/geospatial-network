%Process input time series data, get base Amatrix and the nomalize row revNormA;

function  Amatrix = DataProcessHD( NodeID, Xt, Zt, DX, DZ, xPower, zPower, dxPower, dzPower)

% Prepare A matrix;
count = size( xPower, 2);

useT= size( Xt, 1);

TserMatX  = repmat( Xt(:, NodeID), 1 , count);
TserMatZ  = repmat( Zt(:, NodeID), 1, count);
TserMatDX = repmat( DX(:, NodeID), 1, count);
TserMatDZ = repmat( DZ(:, NodeID), 1, count);

xPowerMat  = repmat( xPower,  useT, 1);
% yPowerMat = repmat( yPower, useT, 1);
zPowerMat  = repmat( zPower,  useT, 1);
dxPowerMat = repmat( dxPower, useT, 1);
dzPowerMat = repmat( dzPower, useT, 1);

Amatrix_temp = (TserMatX.^xPowerMat) .* ( TserMatZ .^ zPowerMat) .* ( TserMatDX .^ dxPowerMat) .* ( TserMatDZ .^ dzPowerMat);

Xt_temp= Xt;
Xt_temp(:, NodeID)=[];

Zt_temp= Zt;
Zt_temp(:, NodeID)=[];

DX_temp= DX;
DX_temp(:, NodeID)=[];

DZ_temp= DZ;
DZ_temp(:, NodeID)=[];

Amatrix = [Amatrix_temp Xt_temp Zt_temp DX_temp DZ_temp ];

end
