function [ nzErr zErr] = StatisticCoup( Nodes, Xpred, Wei)
    
    termNum = size( Xpred, 1);
    termNode = round( termNum / Nodes );

    Xpredwei = Xpred;
    Xpredwei( (NodeID-1)*termNode : NodeID * termNode, 1) = 0; %see the dynamic term as 0;

    Xreal = zeros( termNum, 1 );
    Xreal( IndexWei) = Wei;

    nzList = find( Xreal > 1e-4);
    zList = find( Xreal < 1e-4);

    nzErr = (Xpredwei( nzList) - Xreal( nzList)) ./ Xreal( nzList);
    zErr = Xpredwei( zList) - Xreal( zList);
    
end