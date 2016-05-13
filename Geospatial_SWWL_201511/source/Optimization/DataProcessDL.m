%Process input time series data, get base Amatrix and the nomalize row revNormA;

function  Amatrix = DataProcessDL( NodeID, Xt, Yt, Zt, dXdt, dYdt, dZdt, xPower, yPower, zPower)

% Prepare A matrix;
count = size( xPower, 2);

[ useT Nodes] = size( Xt);

TserMatX = repmat( Xt(:, NodeID), 1 , count);
TserMatY = repmat( Yt(:, NodeID), 1 , count);
TserMatZ = repmat( Zt(:, NodeID), 1 , count);

xPowerMat = repmat( xPower, useT, 1);
yPowerMat = repmat( yPower, useT, 1);
zPowerMat = repmat( zPower, useT, 1);

Amatrix_temp = (TserMatX.^xPowerMat) .* ( TserMatY .^ yPowerMat) .* ( TserMatZ .^ zPowerMat);

Xt(:, NodeID) = [];
Yt(:, NodeID) = [];
Zt(:, NodeID) = [];

dXdt(:, NodeID) = [];
dYdt(:, NodeID) = [];
dZdt(:, NodeID) = [];

% d2Xdt(:, NodeID) = [];
% d2Ydt(:, NodeID) = [];
% d2Zdt(:, NodeID) = [];

% Amatrix = [ones(useT, 1) Amatrix_temp Xt Yt Zt dXdt dYdt dZdt d2Xdt d2Ydt d2Zdt];
Amatrix = [Amatrix_temp Xt Yt Zt dXdt dYdt dZdt];

end
