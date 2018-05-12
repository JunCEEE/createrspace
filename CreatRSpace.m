latticeConstants = {
    {1,'Si','Diamond',1,[5.4305,5.4305,5.4305]}
    {2,'Cu','FCC',2,[3.615,3.615,3.615]}
    {3,'Pt','FCC',2,[3.9242,3.9242,3.9242]}
    {4,'KTaO3','Simple Cubic',4,[3.988,3.988,3.988]}
    }; % No.,Material,Crystal System, Unitcell Type Number, Lattice Constants 

unitcells = {
    {1,'Diamond',[90,90,90],[[0.0,0.0,0.0];[0.5,0.5,0.0];[0.5,0.0,0.5];[0.0,0.5,0.5];[0.25,0.25,0.25];[0.75,0.75,0.25];[0.75,0.25,0.75];[0.25,0.75,0.75]]}
    {2,'FCC',[90,90,90],[[0.0,0.0,0.0];[0.5,0.5,0.0];[0.5,0.0,0.5];[0.0,0.5,0.5]]}
    {3,'BCC',[90,90,90],[[0.0,0.0,0.0];[0.5,0.5,0.5]]}
    {4,'Simple Cubic',[90,90,90],[0.0,0.0,0.0]}
    {5,'HCP',[90,90,120],[[0.0,0.0,0.0];[1.0/3, 2.0/3, 0.5]]}
    }; % Unitcell Type Number,Crystal System, [alpha,beta,gamma],Atom Fractional Coordinates 

%%%% Fiber Calculation %%%%
lcn = 3;
lc = latticeConstants{lcn}{5}; %Lattice constants: a, b, c; unit: angstrom
unitcell = latticeConstants{lcn}{4}; %1 for Diamond, 2 for FCC ...
hkl = getHKL(4,4,4);
sample_axis_x = [1,-1,0];
sample_axis_z = [1,1,1];
fiber_axis = [1,1,1];
[fp,real_hkl1]  = calculationFiber(hkl,unitcells,unitcell,lc,sample_axis_x,sample_axis_z,fiber_axis);
fp = 2*pi*fp;


%%%% Single crystal Calculation %%%%
lcn = 1;
lc = latticeConstants{lcn}{5}; %Lattice constants: a, b, c; unit: angstrom
unitcell = latticeConstants{lcn}{4}; %1 for Diamond, 2 for FCC ...
hkl = getHKL(10,10,10);
sample_axis_x = [1,0,0];
sample_axis_z = [0,0,1];
[q,real_hkl] = calculationSingle(hkl,unitcells,unitcell,lc,sample_axis_x,sample_axis_z);
q = 2*pi*q;

%%%% Polycrystalline Calculation %%%%
lcn = 4;
lc = latticeConstants{lcn}{5}; %Lattice constants: a, b, c; unit: angstrom
unitcell = latticeConstants{lcn}{4}; %1 for Diamond, 2 for FCC ...
hkl = getHKL(3,2,1);
[q_r,real_hkl2] = calculationPoly(hkl,unitcells,unitcell,lc);
q_r = 2*pi*q_r;
%%%%%%%% 3D Polycrystalline Calculation %%%%%%%%
% [qp,real_hkl3,d_q] = calculation3DPoly(hkl,unitcells,unitcell,lc);
% qp = 2*pi*qp;
% d_q = 2*pi*d_q;

%%%% BB Calculation
v = [0,1,0]; %Transverse vector
wavelength = 1.54;
tth_scanrange = [0,0.1,90];% two_theta range (degrees)
[k_BB,tth_BB] = BB(v,wavelength,tth_scanrange);
k_BB = 2*pi*k_BB;

%%%% PB Calculation
v = [0,1,0]; %Transverse vector
wavelength = 1.54;
alpha = 1.0; %Incident angle
tth_scanrange = [0,0.1,90];% two_theta range (degrees)
[k_PB,tth_PB] = PB(v,wavelength,alpha,tth_scanrange);
k_PB = 2*pi*k_PB;    

%%%% Single crystal/Fiber Plot %%%%
reciprocal_axis_x = [0,1,0]; 
reciprocal_axis_y = [0,0,1];
fig1 = figure;
axis equal;

pq = proj(q,reciprocal_axis_x, reciprocal_axis_y);
scatter(pq(:,1),pq(:,2),'r');
hold on
pfp = proj(fp,reciprocal_axis_x, reciprocal_axis_y);
scatter(pfp(:,1),pfp(:,2),'g');
hold on
pk_BB = proj(k_BB,reciprocal_axis_x, reciprocal_axis_y);
scatter(pk_BB(:,1),pk_BB(:,2),'b');
hold on
pk_PB = proj(k_PB,reciprocal_axis_x, reciprocal_axis_y);
scatter(pk_PB(:,1),pk_PB(:,2),'m');
hold on

% pqp = proj(qp,reciprocal_axis_x, reciprocal_axis_y);
% scatter(pqp(:,1),pqp(:,2),'b');
% hold on

%%%% Poly crystal Plot %%%%
for i = 1:length(q_r)
    circle(0,0,q_r(i));
end

%%%% Save %%%%
save_proj(real_hkl,pq,'Projection.dat');
save_rspace(real_hkl1,fp,'FiberReciprocalSpace.dat');
save_rspace(real_hkl,q,'SiReciprocalSpace.dat');
% save_3DPoly(real_hkl3,qp,d_q,'PolyReciprocalSpace.dat');
save_ewald(k_BB,tth_BB,'BB.dat');
save_ewald(k_PB,tth_PB,'PB.dat');


%%%%%%%%%%%%%%%%%%%
%%%% Functions %%%%
%%%%%%%%%%%%%%%%%%%
function [k_BB,tth] = BB(v,wavelength,tth_range)
r = 1/wavelength;
k_BB = [];
tth = [];
for theta = tth_range(1):tth_range(2)/2:tth_range(3)/2
    k_0 = rodrigues_rot(v,[-1,0,0],deg2rad(theta));
    c_0 = -k_0*r;
    k_1 = c_0 + rodrigues_rot(k_0,[1,0,0],deg2rad(theta*2))*r;
    k_BB = [k_BB;k_1];
    tth = [tth;theta*2];
end
end

function [k_PB,tth] = PB(v,wavelength,alpha,tth_range)
alpha = deg2rad(alpha);
r = 1/wavelength;
k_0 = rodrigues_rot(v,[-1,0,0],alpha);
c_0 = -k_0*r;
k_PB = [];
tth = [];
for twotheta = tth_range(1):tth_range(2):tth_range(3)
    k_1 = c_0 + rodrigues_rot(k_0,[1,0,0],deg2rad(twotheta))*r;
    k_PB = [k_PB;k_1];
    tth = [tth;twotheta];
end
end

function save_ewald(k,tth,filename)
fdata = [k,tth];
fileID = fopen(filename,'w');
fprintf(fileID,'%i\n',length(fdata));
fprintf(fileID,'kx ky kz tth\n');
for i = 1:length(fdata)
fprintf(fileID,'%f %f %f %f\n',fdata(i,:));
end
fclose(fileID);
end



function save_proj(hkl,q,filename)
fdata = [hkl,q];
fileID = fopen(filename,'w');
fprintf(fileID,'%i\n',length(fdata));
fprintf(fileID,'h k l x y\n');
for i = 1:length(fdata)
fprintf(fileID,'%f %f %f %f %f\n',fdata(i,:));
end
fclose(fileID);
end

function save_rspace(hkl,q,filename)
fdata = [hkl,q];
fileID = fopen(filename,'w');
fprintf(fileID,'%i\n',length(fdata));
fprintf(fileID,'h k l qx qy qz\n');
for i = 1:length(fdata)
fprintf(fileID,'%f %f %f %f %f %f\n',fdata(i,:));
end
fclose(fileID);
end


function save_3DPoly(hkl,q,d_r,filename)
fdata = [hkl,q,d_r];
fileID = fopen(filename,'w');
fprintf(fileID,'%i\n',length(fdata));
fprintf(fileID,'h k l qx qy qz d_r\n');
for i = 1:length(fdata)
fprintf(fileID,'%f %f %f %f %f %f %f\n',fdata(i,:));
end
fclose(fileID);
end



function [fp,real_hkl1] = calculationFiber(hkl,unitcells,unitcell,lc,sample_axis_x,sample_axis_z,fiber_axis)
RUnitCellVec = getRUnitCellVec(unitcells{unitcell}{3},lc,sample_axis_x,sample_axis_z);
q_fiber = fiber_axis*RUnitCellVec;
fp = [];
real_hkl1 = [];
[p,real_hkl00] = calculationSingle(hkl,unitcells,unitcell,lc,sample_axis_x,sample_axis_z);
for theta = 0:2/3*pi/100:2/3*pi
    for i = 1:length(p)
        p(i,:) = rodrigues_rot(p(i,:),q_fiber,theta);
    end
fp = [fp;p];
real_hkl1 = [real_hkl1;real_hkl00];
end
end


function fp = calculationFiber1(hkl,unitcells,unitcell,lc,sample_axis_x,sample_axis_z,fiber_axis,pOnFiber)
atoms = atomPosition(unitcells{unitcell}{3},unitcells{unitcell}{4});
real_hkl = getReal_hkl(hkl,atoms);
RUnitCellVec = getRUnitCellVec(unitcells{unitcell}{3},lc,sample_axis_x,sample_axis_z);
q = real_hkl*RUnitCellVec;
q_fiber = fiber_axis*RUnitCellVec;
q_fiber = q_fiber/norm(q_fiber);
q_pOnFiber = pOnFiber*RUnitCellVec;
AB = q - q_pOnFiber;
n_BC = AB * q_fiber';
BC = n_BC.*q_fiber;
CA = AB - BC;
r = zeros(length(q),1);
for i = 1:length(q)
r(i) = norm(CA(i,:));
end
% A = [int64(n_BC*1000000),int64(r*1000000)];
% [AA,ia] = unique(A,'rows');
% r = double(AA(:,2))/1000000;
% r = zeros(length(q),1);
% for i = 1:length(q)
%     r(i) = norm(cross(q_fiber,q(i,:)-q_fiber))/norm(q_fiber);
% end
% A = [real_hkl,int64(r*1000000)];
% r = double(r)/1000000;
% real_hkl = A(ia,1:3);
% q = q(ia,:);
fp = getFiberCircle(q,q_pOnFiber,q_fiber,r);
end

function fp = getFiberCircle(q,q_pOnFiber,q_fiber,r)
q_fiber = q_fiber/norm(q_fiber);
AB = q - q_pOnFiber;
BC = (AB * q_fiber').*q_fiber;
CA = AB - BC;
a = zeros(length(q),3);
for i = 1:length(q)
a0 = cross(q_fiber,CA(i,:));
a(i,:) = a0/norm(a0)*r(i);
end
gamma = 0:pi/50:2*pi;
s_gamma = size(gamma);
fp = zeros(length(q)*s_gamma(2),4);
nn = 0;
for i = 1:length(q)
    for j = 1:s_gamma(2)
        nn = nn+1;        
        fp(nn,1:3) = CA(i,:).*cos(gamma(j)) + a(i,:).*sin(gamma(j))+BC(i,:);
        fp(nn,4) = r(i);
    end
end
end


function [q,real_hkl] = calculationSingle(hkl,unitcells,unitcell,lc,sample_axis_x,sample_axis_z)
atoms = atomPosition(unitcells{unitcell}{3},unitcells{unitcell}{4});
real_hkl = getReal_hkl(hkl,atoms);
RUnitCellVec = getRUnitCellVec(unitcells{unitcell}{3},lc,sample_axis_x,sample_axis_z);
q = real_hkl*RUnitCellVec;
end

function [d_q,real_hkl] = calculationPoly(hkl,unitcells,unitcell,lc)
atoms = atomPosition(unitcells{unitcell}{3},unitcells{unitcell}{4});
real_hkl = getReal_hkl(hkl,atoms);
RUnitCellVec = getRUnitCellVec(unitcells{unitcell}{3},lc,[1,0,0],[0,0,1]);
q = real_hkl*RUnitCellVec;
d_q = sqrt(sum(q.^2, 2));
A = [real_hkl,int64(d_q*1000000)];
[d_q,ia] = unique(A(:,4),'rows','last');
d_q = double(d_q)/1000000.0;
real_hkl = double(A(ia,1:3));
end

function [q,real_hkl2,d_q2] = calculation3DPoly(hkl,unitcells,unitcell,lc)
atoms = atomPosition(unitcells{unitcell}{3},unitcells{unitcell}{4});
real_hkl = getReal_hkl(hkl,atoms);
RUnitCellVec = getRUnitCellVec(unitcells{unitcell}{3},lc,[1,0,0],[0,0,1]);
q = real_hkl*RUnitCellVec;
d_q = sqrt(sum(q.^2, 2));
A = [real_hkl,int64(d_q*1000000)];
[d_q,ia] = unique(A(:,4),'rows','last');
d_q = double(d_q)/1000000.0;
real_hkl = A(ia,1:3);
q = [];
real_hkl2 = [];
d_q2 = [];
for theta = 0:pi/100:pi
    for phi = 0:2*pi/100:2*pi
        x = d_q*sin(theta)*cos(phi);
        y = d_q*sin(theta)*sin(phi);
        z = d_q*cos(theta);
        q = [q;[x,y,z]];
        real_hkl2 = [real_hkl2;real_hkl];
        d_q2 = [d_q2;d_q];
    end
end
real_hkl2 = double(real_hkl2);
end



function RUnitCellVec = getRUnitCellVec(abg,lc,sample_axis_x,sample_axis_z)
a = atomPosition(abg,[lc(1),0,0]);
a = rot_all(a,sample_axis_x,sample_axis_z);
b = atomPosition(abg,[0,lc(2),0]);
b = rot_all(b,sample_axis_x,sample_axis_z);
c = atomPosition(abg,[0,0,lc(3)]);
c = rot_all(c,sample_axis_x,sample_axis_z);
V = dot(a,cross(b,c));
a_r = cross(b,c)/V; 
b_r = cross(c,a)/V;
c_r = cross(a,b)/V;
RUnitCellVec = [a_r;b_r;c_r];
end

function atoms = atomPosition(abg,unitcell)
v1 = [1,0,0];
v2 = [cosd(abg(3)),sind(abg(3)),0];
cz = sqrt(1-cosd(abg(1))*cosd(abg(1))-cosd(abg(2))*cosd(abg(2))-cosd(abg(3))*cosd(abg(3))+2*cosd(abg(1))*cosd(abg(2))*cosd(abg(3)))/sind(abg(3));
v3 = [cosd(abg(2)),(cosd(abg(1))-cosd(abg(2))*cosd(abg(3)))/sind(abg(3)),cz];
v3 = v3/norm(v3);
% vx = unitcell(:,1).*v1;
% vy = unitcell(:,2).*v2;
% vz = unitcell(:,3).*v3;
vv = [v1;v2;v3];
atoms = unitcell*vv;
end



function hkl = getHKL(a,b,c)
H = -a:a;
K = -b:b;
L = -c:c;
nq = length(H)*length(K)*length(L);
hkl = zeros(nq,3);
nl = 0;
for i = H
    for j = K
        for k = L
            nl = nl+1;
            hkl(nl,:) = [i,j,k];
        end
    end
end
end

function pp = proj(vecs,axis_x, axis_y)
axis_x = axis_x/norm(axis_x);
axis_y = axis_y/norm(axis_y);
proj_x = vecs*axis_x';
proj_y = vecs*axis_y';
pp = [proj_x,proj_y];
end

function real_hkl = getReal_hkl(hkl,atoms)
real_hkl = [];
for i = 1:length(hkl)
    a = hkl(i,:);
    kr = atoms*a';
    sf0 = exp(2i*pi.*kr);
    sf = sum(sf0);
    if (abs(sf) > 1e-7)
        real_hkl = [real_hkl ;hkl(i,:)];
    end
end
end

function after_rot = rot_all(vecs,sample_axis_x,sample_axis_z)
r1 = vrrotvec(sample_axis_x,[1,0,0]);
% a1 = rodrigues_rot(sample_axis_x,r1(1:3),r1(4));
b1 = rodrigues_rot(sample_axis_z,r1(1:3),r1(4));
nl = size(vecs);
for i = 1:nl(1)
    vecs(i,:) = rodrigues_rot(vecs(i,:),r1(1:3),r1(4));
end
r2 = vrrotvec(b1,[0,0,1]);
% a = rodrigues_rot(a1,r2(1:3),r2(4));
% b = rodrigues_rot(b1,r2(1:3),r2(4));
for i = 1:nl(1)
    vecs(i,:) = rodrigues_rot(vecs(i,:),r2(1:3),r2(4));
end
after_rot = vecs;
end

function h = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
hold off
end

function v_rot = rodrigues_rot(v,k,theta)
    [m,n] = size(v);
    if (m ~= 3 && n ~= 3)
        error('input vector is/are not three dimensional'), end
    if (size(v) ~= size(k)) 
        error('rotation vector v and axis k have different dimensions'),end
    
    k = k/sqrt(k(1)^2 + k(2)^2 + k(3)^2); % normalize rotation axis
    No = numel(v)/3; % number of vectors in array
    v_rot = v; % initialize rotated vector array
    if ( n == 3 )
        crosskv = v(1,:); % initialize cross product k and v with right dim.
        for i = 1:No
            crosskv(1) = k(2)*v(i,3) - k(3)*v(i,2);
            crosskv(2) = k(3)*v(i,1) - k(1)*v(i,3); 
            crosskv(3) = k(1)*v(i,2) - k(2)*v(i,1);
            v_rot(i,:) = cos(theta)*v(i,:) + (crosskv)*sin(theta)...
                            + k*(dot(k,v(i,:)))*(1 - cos(theta));
        end
    else % if m == 3 && n ~= 3
        crosskv = v(:,1); % initialize cross product k and v with right dim.
        for i = 1:No
            crosskv(1) = k(2)*v(3,i) - k(3)*v(2,i);
            crosskv(2) = k(3)*v(1,i) - k(1)*v(3,i); 
            crosskv(3) = k(1)*v(2,i) - k(2)*v(1,i);
            v_rot(:,i) = cos(theta)*v(:,i) + (crosskv)*sin(theta)...
                            + k*(dot(k,v(:,i)))*(1 - cos(theta));
        end
    end
end
