% Test script for the Class: DiracGaussian
% Test script fot the function: constructGaussianCart
% Test script fot the function: constructGaussianPol

addpath([pwd, '\..'])


%%% Setup szeanrio %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = 0;
c = 1;
pot = 0;

k0x = 3;
k0y = 0;
b = 0.05; 

% set up spacial domain of area 4
x0 = -1;
x1 = 1;
x = linspace(x0, x1, 1000);
y0 = -1;
y1 = 1;
y = linspace(y0, y1, 1000);
[xx, yy] = meshgrid(x,y);


%% diracGaussianCart
disp('testint: DiracGaussianCart')


% construct wave packet
gwp = diracEq2D.constructGaussianCart(...
    k0x, ...  %kx0
    k0y, ...  %ky0
    b, ... %b
    'nk', 17, ...
    'mk', 8*pi, ...
    't0',  0, ...
    'x0',  -0.5, ...
    'y0',  -0.5, ...
    'potential', pot, ...
    'mass', m, ...
    'c', c, ...
    'solution', 1, ...
    'volumen', (x1-x0)*(y1-y0));


%%% test if solutions are correct %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_result = zeros(1, length(gwp.energies));
prec = 1e-13;
for i_k = 1:length(gwp.energies)
    
    k = gwp.wavenumbers(:,i_k);
    M = 1i*c*[(pot+m)/(1i*c), k(1)-1i*k(2);  -k(1)-1i*k(2), (pot-m)/(1i*c)];
    E = gwp.energies(i_k) ;
    A = gwp.eigenspinors(:,i_k);
    test_result(i_k) = norm(M*A - E*A) <= prec;

end

if all(test_result)
    disp('Test 1: passed')
else
    disp('Test 1: failed')
    disp(find(~test_result))
end

%%% test2: test if mean is correct %%%%%%%%%%%%
[k_avg, ~] = gwp.estimateEnergeticProperties();
prec = 1e-13;
if all(abs(k_avg - [k0x;k0y])<= prec) 
    disp('Test 2: passed')
else
    disp('Test 2: failed')
    disp(abs(k_avg - [k0x;k0y]))
end

%%% test if normalisation is correct %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[u, v] = gwp.getComponent(xx,yy,0);
w = abs(u).^2 + abs(v).^2;

p = trapz(y,  trapz(x , w, 2)); % numerical integrate over domain

prec = 1e-3;
if abs(p-1)<prec
        disp('Test 2: passed')
else
    disp('Test 2: failed')
    disp(abs(p))

end

%%% visualize wave packet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
surf(xx,yy, w);
shading interp
view(0,90)
axis square


%% diracGaussianPol
disp('testint: DiracGaussianPol')

% construct wave packet
gwp = diracEq2D.constructGaussianPol(...
    k0x, ...  %kx0
    k0y, ...  %ky0
    b , ... %b
    8*pi, ...
    1*pi/(x1-x0), ...
    1*pi/(y1-y0), ...
    't0',  0, ...
    'x0',  -0.5, ...
    'y0',  -0.5, ...
    'potential', pot, ...
    'mass', m, ...
    'c', c, ...
    'solution', 1, ...
    'volumen', (x1-x0)*(y1-y0));

%%% test if solutions are correct %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_result = zeros(1, length(gwp.energies));
prec = 1e-13;
for i_k = 1:length(gwp.energies)
    
    k = gwp.wavenumbers(:,i_k);
    M = 1i*c*[(pot+m)/(1i*c), k(1)-1i*k(2);  -k(1)-1i*k(2), (pot-m)/(1i*c)];
    E = gwp.energies(i_k) ;
    A = gwp.eigenspinors(:,i_k);
    test_result(i_k) = norm(M*A - E*A) <= prec;

end

if all(test_result)
    disp('Test 1: passed')
else
    disp('Test 1: failed')
    disp(find(~test_result))
end

%%% test2: test if mean is correct %%%%%%%%%%%%
[k_avg, ~] = gwp.estimateEnergeticProperties();
prec = 1e-13;
if all(abs(k_avg - [k0x;k0y])<= prec) 
    disp('Test 2: passed')
else
    disp('Test 2: failed')
    disp(abs(k_avg - [k0x;k0y]))
end

%%% test if normalisation is correct %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[u, v] = gwp.getComponent(xx, yy, 0);
w = abs(u).^2 + abs(v).^2;

p = trapz(y,  trapz(x , w, 2)); % numerical integrate over domain

prec = 1e-3;
if abs(p-1)<prec
        disp('Test 2: passed')
else
    disp('Test 2: failed')
    disp(abs(p))

end

%%% visualize wave packet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
surf(xx,yy, w);
shading interp
view(0,90)
axis square

%% dirtyGaussian
disp('testint: DiracGaussianDirty')

% construct wave packet
normalize = true;

[u, v] = diracEq2D. constructGaussianDirty(...
    xx, ... %xx
    yy, ... %yy
    k0x, ...%kx0
    k0y, ...%ky0
    0.2 , ... %bx
    0.2 , ... %by
    'x0',  -0.5, ...
    'y0',  -0.5, ...
    'potential', pot, ...
    'mass', m, ...
    'c', c, ...
    'solution', 1, ...
    'normalize', normalize);
w = abs(u).^2 + abs(v).^2;


%%% test if normalisation is correct %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = trapz(y,  trapz(x , w, 2)); % numerical integrate over domain

prec = 1e-3;
if abs(p-1)<prec
    disp('Test 1: passed')
elseif normalize
    disp('Test 1: failed')
    disp(abs(p))
else 
    disp('Test 1: passed')
end
    
%%% visualize wave packet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
surf(xx,yy, w);
shading interp
view(0,90)
axis square
