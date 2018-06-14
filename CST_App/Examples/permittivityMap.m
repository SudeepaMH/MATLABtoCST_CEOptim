%Read in data for a permittivity map and create the irregular grid profile
%in CST. This example uses the permittivity map from a pre-defined
%quasi-optical coordinate transformation used to disguise a bump in a PEC
%ground as described in "Hao et. al. 'Design of a Carpet Cloak to Conceal an Antenna
%Located Underneath', IEEE Trans. Antennas + Propagation, 2012, 
%DOI: 10.1109/TAP.2012.2207058"

%Edit path as required
M = dlmread('permittivityMap.txt');

X = M(:,1);
Y = M(:,2);
Er = M(:,3);
map = jet(128);

X = reshape(X,21,41);
Y = reshape(Y,21,41);
Er = reshape(Er,21,41);

figure;
surf(X,Y,Er)
colormap(map)
view([0 90])

CST = CST_MicrowaveStudio(cd,'test.cst');

ni = sqrt(Er);
nAv = mapAverage(ni);
minn = min(nAv(:));
maxn = max(nAv(:));
            
map = jet(16);

freq = 10;
lamda = 3e8/10e9*1e3;
height = lamda/2;

ii = addGradedIndexMaterialCST_iregular(CST,X,Y,nAv,map,1,'scaleFactor',[1,1,height],'CRange',[minn maxn],'MatchedImpedance',false);
CST.setBackgroundLimits([0 0],[50 0],[0 0])
CST.addFieldMonitor('e-field',f)
CST.addFieldMonitor('farfield',f)
CST.setFreq(10*0.95,10*1.15);
CST.addSymmetryPlane('Z','Electric')

Xmin = min(min(X))-50;
Xmax = max(max(X))+50;
Ymin = min(min(Y));
CST.addBrick([Xmin,Xmax],[Ymin Ymin],[0 height],'Ground1','Component2','PEC');
CST.addCylinder(56.6,55.6,'z',0,-40,[0 height],'cylinder 1','Component2','PEC');

CST.addWaveguidePort('ymax',[-20 20],200, [0 7.1]);
CST.selectSolver('frequency')

% The user should now rotate the waveguide port around the Z axis at X=0, Y=0 by an
% anle of 45 degrees (make sure the 'Rotate axis aligned box' is
% unchecked), before running the simulations and viewing the results. The
% simulation will be much quicker to run if all materials of the same
% permittivity are added together using the booleon operations in CST
% 
%

