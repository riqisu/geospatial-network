function OutputName = readData( NetworkName, RK_H)

TimeSerialName = ['TimeSer' NetworkName];
WeiTau = [ NetworkName];

%delay list and as the name of output;
OutputName = ['TD_' NetworkName];

combine = load( [ TimeSerialName '.txt']);
% combine = load( 'combine.txt');

[allT Nodes3] = size( combine);

Nodes = Nodes3 / 3;
DimX = Nodes;
DimY = Nodes;

allT = floor( allT / 3) * 3;

DataX = combine(:, 1:3: Nodes3);
DataY = combine(:, 2:3: Nodes3);
DataZ = combine(:, 3:3: Nodes3);

Xt = DataX(2:3:allT, :);
Yt = DataY(2:3:allT, :);
Zt = DataZ(2:3:allT, :);

dXdt = ( DataX(3:3:allT, :) - DataX(1:3:allT, :) ) ./ ( 2* RK_H);
dYdt = ( DataY(3:3:allT, :) - DataY(1:3:allT, :) ) ./ ( 2* RK_H);
dZdt = ( DataZ(3:3:allT, :) - DataZ(1:3:allT, :) ) ./ ( 2* RK_H);

d2Xdt = ( DataX(3:3:allT, :) + DataX(1:3:allT, :) -2* DataX(2:3:allT, :)) ./ ( RK_H * RK_H);
d2Ydt = ( DataY(3:3:allT, :) + DataY(1:3:allT, :) -2* DataY(2:3:allT, :)) ./ ( RK_H * RK_H);
d2Zdt = ( DataZ(3:3:allT, :) + DataZ(1:3:allT, :) -2* DataZ(2:3:allT, :)) ./ ( RK_H * RK_H);
% 
% 
tempSparse = load([ WeiTau '.txt']);

[ links btemp] = size( tempSparse);
matTemp1 = sparse( tempSparse(:,1)+1, tempSparse(:,2)+1, tempSparse(:,3), DimX, DimY);
matTemp2 = sparse( tempSparse(:,2)+1, tempSparse(:,1)+1, tempSparse(:,3), DimX, DimY);
iniWmatrix = full(matTemp1 + matTemp2);

DelayMat1 = sparse( tempSparse(:,1)+1, tempSparse(:,2)+1, tempSparse(:,4), DimX, DimY);
% DelayMat2 = sparse( tempSparse(:,2)+1, tempSparse(:,1)+1, tempSparse(:,4), DimX, DimY);
DelayMat = full( DelayMat1);

save( [OutputName '.mat'], 'Xt', 'Yt', 'Zt', 'dXdt', 'dYdt', 'dZdt', 'd2Xdt', 'd2Ydt', 'd2Zdt', 'Nodes', 'RK_H');
save( [OutputName '.mat'],'-append', 'iniWmatrix', 'DelayMat');

% 
% figure( 1 )
% plot3(DataX2(:,8), DataY2(:,8), DataZ2(:,8), 'b.');
% hold on;
% plot3(DataX1(:,8), DataY1(:,8), DataZ1(:,8), 'r.');
% hold on;
% plot3(DataX3(:,8), DataY3(:,8), DataZ3(:,8), 'k.');
% hold on;
% print( 1,'-depsc2', ['Node' num2str(1) '.eps']);
% hold off;
