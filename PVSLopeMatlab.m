statefile='ocean_his_p5.nc';
diagfile='ocean_dia_p5.nc';
slice = {0, 0, 0, 0};
nx = 402; ny=5; nz=60; nt=240;
% Q = GetVarROMSRutgers(statefile, statefile, {'PV', '(1)'}, slice);
pm = ncread(statefile, 'pm');
dx = 1./pm;
pn = ncread(statefile, 'pn');
dy = 1./pn;
zw = depths(statefile, statefile, 5, 0, 1);
dz = diff(zw,1,3);

Q = NaN(nx, ny, nz, nt);
W = Q; T = Q;
JDDIV = Q;
JFDIV = Q;

z = depths(statefile, statefile, 1, 0, 1); % XX- SHOULD BE IN FOR LOOP...BUT SLOW..

for i=1:nt
    disp(num2str(i));
    sliceT = {slice{1}, slice{2}, slice{3}, [i i]};

U = GetVarROMSRutgers(statefile, statefile, {'u', '(1)'}, sliceT);
V = GetVarROMSRutgers(statefile, statefile, {'v', '(1)'}, sliceT);
% B = GetVarROMSRutgers(statefile, statefile, {'b', '(1)'}, sliceT);
T(:,:,:,i) = GetVarROMSRutgers(statefile, statefile, {'temp', '(1)'}, sliceT);
B = 9.8*2e-4.*squeeze(T(:,:,:,i));
W(:,:,:,i) = GetVarROMSRutgers(statefile, statefile, {'w', '(1)'}, sliceT);

Wx = DrvS(pm, z, squeeze(W(:,:,:,i)), 'x');
Wy = DrvS(pn, z, squeeze(W(:,:,:,i)), 'y');
Uy = DrvS(pn, z, U, 'y');
Vx = DrvS(pm, z, V, 'x');
Uz = DrvROMS(z,U, 'z');
Vz = DrvROMS(z, V, 'z');
Bx = DrvS(pm, z, B, 'x');
By = DrvS(pn, z, B, 'y');
Bz = DrvROMS(z, B, 'z');

f = ncread(statefile, 'f');

OMEGAZ = (repmat(f,[1 1 nz])+Vx-Uy);
OMEGAX = -Vz + Wy;
OMEGAY = Uz - Wx;
Q(:,:,:,i) = OMEGAZ.*Bz + OMEGAY.*By +OMEGAX.*Bx;

D = 9.8*2e-4.*GetVarROMSRutgers(statefile, diagfile, {'temp_vdiff', 'temp_hdiff', '(1)+(2)'}, sliceT);
JDz = -OMEGAZ.*D;
JDx = -OMEGAX.*D;
JDy = -OMEGAY.*D;
JDDIV(:,:,:,i) = DrvS(pm, z,JDx, 'x') + DrvS(pn, z, JDy, 'y') + DrvROMS(z,JDz,'z');
% JDDIV(:,:,:,i) = DrvROMS(z,JDz,'z')+DrvS(pm, z,JDx, 'x');


Fx = GetVarROMSRutgers(statefile, diagfile, {'u_vvisc', 'u_hvisc', '(1)+(2)'}, sliceT);
Fy = GetVarROMSRutgers(statefile, diagfile, {'v_vvisc', 'v_hvisc', '(1)+(2)'}, sliceT);
JFz = Bx.*Fy - By.*Fx;
JFx = -Bz.*Fy;
JFy = Bz.*Fx;
JFDIV(:,:,:,i) = DrvS(pm, z,JFx, 'x') + DrvS(pn, z, JFy, 'y') + DrvROMS(z,JFz,'z');
% JFDIV(:,:,:,i) =  DrvROMS(z,JFz,'z')+DrvS(pm, z,JFx, 'x');


end
%%
[nx ny nz nt] = size(Q);
dxf = repmat(dx, [1 1 nz nt]);
dyf = repmat(dy, [1 1 nz nt]);
dzf = repmat(dz, [1 1 1 nt]);
dys = squeeze(dy(1,1));
dxs = squeeze(dx(1,1));

Qa = squeeze(nansum(nansum(nansum(Q.*dxf.*dyf.*dzf))));
Qa = Qa-Qa(1);
% [~, ~, ~, Qt] = gradient(Q, 3*3600
% Qta = squeeze(nansum(nansum(nansum(Qt.*dxf.*dyf.*dzf))));
Qta = gradient(Qa, 3600*3);
JTOTFdt = squeeze(nansum(nansum(nansum((JFDIV).*dxf.*dyf.*dzf))));
JTOTDdt = squeeze(nansum(nansum(nansum((JDDIV).*dxf.*dyf.*dzf))));

JTOTF = cumtrapz(JTOTFdt).*3600*3;
JTOTD = cumtrapz(JTOTDdt).*3600*3;

%%
figure
% plotyy(1:nt, Qta, 1:nt, -JTOTdt)
subplot(2,1,1)
plot(Qta, 'LineWidth', 2); hold on; 
plot(-JTOTFdt, 'LineWidth', 2); 
plot(-JTOTDdt, 'LineWidth', 2);
plot(-(JTOTFdt + JTOTDdt),'--', 'LineWidth', 2);
hold off
l = legend('$dq/dt$', '$J_F$', '$J_D$', 'SUM');
set(l,  'Interpreter', 'Latex')

grid on
subplot(2,1,2)
plot(Qa, 'LineWidth', 2); hold on; 
plot(-JTOTF, 'LineWidth', 2); 
plot(-JTOTD, 'LineWidth', 2);
plot(-(JTOTF + JTOTD),'--', 'LineWidth', 2);
hold off
l = legend('$\Delta q$', '$\int_t J_F$', '$\int_t J_D$', 'SUM');
set(l,  'Interpreter', 'Latex')
grid on

set(gcf, 'Color','w')
% scatter(Qta, -JTOTdt)
%  plot(Qa-Qa(10));
%  hold on
%  plot((JTOT-JTOT(10)));
%  hold off
% plotyy(1:nt, Qa, 1:nt, JTOT/1021);
%%
psi = NaN(nx, nz, nt);
for i=1:nz
    for t=1:nt
psi(:,i,t) = cumtrapz(squeeze(W(:,3,i,t))).*squeeze(dxf(:,3,i,t));
    end
end
%%
cl = [-1 1].*1e-4;
% cl = [-1 1].*5e-3;
Z = squeeze(z(:,3,:));
X = repmat((1:nx).', [1 60]);
figure
for i=1:9
    ti = 1+floor(nt./9).*i;
    subplot(3,3,i);
            pcolor(X.', Z.', squeeze(W(:,3,:,ti)).'); shading interp

    title(num2str(ti./8, 2))
    hold on
    contour(X.', Z.', squeeze(psi(:,:,ti)).', 20,'Color', [.5 .5 .5]); shading interp

    contour(X.', Z.', squeeze(T(:,3,:,ti)).', -10:.25:40, 'k');
    contour(X.', Z.', squeeze(Q(:,3,:,ti)).', [0 0], 'k', 'LineWidth', 2);
    hold off
        set(gca, 'clim', cl);

end
colormap(cptcmap('balance'));
set(gcf, 'Color','w', 'Position',[675   216   895   758]);

%%


cl = [-1 1].*1e-3;
% cl = [-1 1].*5e-3;
Z = squeeze(z(:,3,:));
X = repmat((1:nx).', [1 60]);
for i=1:9
    ti = 1+floor(nt./18).*i;
    subplot(3,3,i);
    pcolor(X.', Z.', squeeze(W(:,3,:,ti)).'); shading interp
    title(num2str(ti./8, 2))
    hold on
    contour(X.', Z.', squeeze(T(:,3,:,ti)).', -10:.25:40, 'k');
    contour(X.', Z.', squeeze(Q(:,3,:,ti)).', [0 0], 'k', 'LineWidth', 2);
    hold off
        set(gca, 'clim', cl);

end

