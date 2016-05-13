function [target count] = CoupCof( maxPow, Comp)

	count =0;
	
	for i=0:maxPow
	    for j =0:maxPow-i
		for k= 0:maxPow -i -j
		        if ( i+j + k)>0
		            count = count +1;
		            xPower(count)= k;
		            yPower(count)= j;
		            zPower(count)= i;
		        end
		end
	    end
	end
	
	if Comp=='X'
		TarRow = (xPower == 1) .* ( yPower == 0) .* ( zPower == 0);
	elseif Comp == 'Y'
		TarRow = (xPower == 0) .* ( yPower == 1) .* ( zPower == 0);
	elseif Comp == 'Z'
		TarRow = (xPower == 0) .* ( yPower == 0) .* ( zPower == 1);
	else
		fprintf( 'No such component');
	end
	
	Posit = find( TarRow == 1 );
	target = Posit(1);
end
