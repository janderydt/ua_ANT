function TrainANN

geometry = 2000;
cycle = 1;
net = feedforwardnet([10],'trainlm');

% set early stopping parameters
net.divideFcn= 'dividerand';
net.divideParam.trainRatio = 0.8; % training set [%]
net.divideParam.valRatio = 0.1; % validation set [%]
net.divideParam.testRatio = 0.1; % test set [%]
 
load("inversiondata.mat");

startgeometry=[data(:).startgeometry];
Ind = find(startgeomtry==geometry);
gaA=[data(Ind).gaA];
gaC=[data(Ind).gaC];
gsA=[data(Ind).gsA];
gsC=[data(Ind).gsC];
m=[data(Ind).m];
n=[data(Ind).n];

for ii=1:numel(Ind)
    qGL(ii) = data(Ind).qGL(1,cycle);
end

[X(1,:),XC(1),XS(1)] = normalize(gaA);
[X(2,:),XC(2),XS(2)] = normalize(gaC); 
[X(3,:),XC(3),XS(3)] = normalize(gsA); 
[X(4,:),XC(4),XS(4)] = normalize(gsC); 
[X(5,:),XC(5),XS(5)] = normalize(m); 
[X(6,:),XC(6),XS(6)] = normalize(n); 

X = gpuArray(double(X));            

% training data
[Vtmp,VC,VS] = normalize(qGL(:)');
V = gpuArray(double(Vtmp));

net = configure(net,X,V);
%view(net);

trainedNet = train(net,X,V,'useParallel','yes','showResources','yes','useGPU','only');
Y = trainedNet(X);
perf = perform(trainedNet,Y,V)

figure; plot(Y,V,'*k'); hold on;
plot([-1.5 3],[-1.5 3]);