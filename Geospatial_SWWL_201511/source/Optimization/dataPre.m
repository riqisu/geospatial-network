%Noise;
function [Xt, Yt, Zt, dXdt, dYdt, dZdt] = dataPre( TimeSer, NoiStr, RK_H)

    [allT, allNodes] = size( TimeSer);
    
    TimeSerNoi = TimeSer(:, 2:allNodes) + ( rand(allT, allNodes-1) - 0.5)* NoiStr;
    
    Xt = TimeSerNoi( 2:3:allT, 2:3:allNodes);
    Yt = TimeSerNoi( 2:3:allT, 3:3:allNodes);
    Zt = TimeSerNoi( 2:3:allT, 4:3:allNodes);
    
    dXdt = ( TimeSerNoi( 3:3:allT, 2:3:allNodes) - TimeSerNoi( 1:3:allT, 2:3:allNodes) ) ./ ( 2* RK_H);
    dYdt = ( TimeSerNoi( 3:3:allT, 3:3:allNodes) - TimeSerNoi( 1:3:allT, 3:3:allNodes) ) ./ ( 2* RK_H);
    dZdt = ( TimeSerNoi( 3:3:allT, 4:3:allNodes) - TimeSerNoi( 1:3:allT, 4:3:allNodes) ) ./ ( 2* RK_H);
    
    