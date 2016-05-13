function [Xpred, NeedTime, redo] = SensingP2( Amatrix, dXdt, epsion)	

    

	partAmat = Amatrix;
    partDX = dXdt;
    NeedData = size(dXdt, 1);
    
%     A= rand(NeedData, NeedData);
%     
%     partAmat= A* partAmat;
%     partDX= A*partDX;
    
     % Normalize;
    NormA = sqrt ( sum( partAmat .* partAmat) );
    revNormA = 1 ./ NormA;
    partAmat = partAmat .*  repmat( revNormA, NeedData, 1);
    
	tic;
	X0 = partAmat' * partDX;
    
%     [Xpred redo]= l1eq_pd( X0, partAmat , [],  partDX, 1e-10);
     Xpred = l1qc_logbarrier(X0,  partAmat, [],  partDX, epsion, 1e-10);

%     x0 = inv( partAmat'* partAmat)* partAmat'* partDX;
%    Xpred = l1decode_pd(X0, partAmat, [],  partDX, 1e-3);
	NeedTime = toc;

	Xpred = Xpred .* revNormA';

    redo=0;
    
    if ( sum( isnan(Xpred) ) > 0 )
        redo = 1;
    end

end
