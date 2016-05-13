%Process input time series data, get base Amatrix and the nomalize row revNormA;

function  Amatrix = DataProcessHD(  Xt, Yt, Zt, xPower, yPower, zPower)

% Prepare A matrix;
count = size( xPower, 2);
[useT, Nodes] = size(Xt);

xPowerMat  = repmat( xPower,  useT, 1);
yPowerMat =  repmat( yPower,  useT, 1);
zPowerMat  = repmat( zPower,  useT, 1);

Amatrix_Pow = zeros( useT, Nodes * count);

for i= 1:Nodes

    TserMatX  = repmat( Xt(:, i), 1 , count);
    TserMatY  = repmat( Yt(:, i), 1 , count);
    TserMatZ  = repmat( Zt(:, i), 1 , count);

    Amatrix_temp = (TserMatX.^xPowerMat) .* ( TserMatY .^ yPowerMat) .* ( TserMatZ .^ zPowerMat);

    Amatrix_Pow(:, (i-1)*count+1:i*count ) = Amatrix_temp;
end

Amatrix = [ones(useT,1) Amatrix_Pow];

end
