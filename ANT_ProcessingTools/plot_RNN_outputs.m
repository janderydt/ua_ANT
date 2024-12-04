function plot_RNN_outputs

perturbation='Calv_dh';
cycle=1;

pct = string(99);

load("Delta_u_AS_Weertman.mat");
T = Delta_u.(perturbation).map(:,:,cycle);
T_mean=mean(T,1);
T = T-repmat(T_mean(:)',size(T,1),1);

for ii=1:numel(pct)

    net=importNetworkFromTensorFlow("./RNN/TF_files/tuned_model_forAll_N"+pct(ii));
    netUpdated=initialize(net);
    load("./RNN/mat_files/data_N0k"+pct(ii)+".mat");
    load("./RNN/mat_files/SVD_N0k"+pct(ii)+".mat");
    X_full = [X_train; X_val; X_test];
    T_full = [T_train; T_val; T_test];
    
    Y_train=double(predict(netUpdated,X_train));
    Y_val=double(predict(netUpdated,X_val));
    Y_test=double(predict(netUpdated,X_test));

    for m=1:size(T_full,2)
        figure; hold on;
        plot(T_train(:,m),Y_train(:,m),'ok');
        plot(T_val(:,m),Y_val(:,m),'xk');
        plot(T_test(:,m),Y_test(:,m),'dk');
        plot(T_test(10,m),Y_test(10,m),'dr');
        %plotregression(T_full,Y(:,m));
        title("mode "+num2str(m));
    end

    % reconstruct spatial maps
    Y_test = double(predict(netUpdated,X_test(50,:)));
    Y_test = Y_test.*T_train_S+T_train_C; % undo normalization
    Y_reproj = Y_test*B_trunc;% reproject onto nodal basis
    %Y_reproj = Y_test(:)\T_reproj;

    figure; tiledlayout(1,3,"TileSpacing","tight");
    nexttile; hold on;
    PlotMeshScalarVariable([],MUA_2018,Y_reproj'+T_mean');
    caxis([-1000 1000]); axis equal; axis tight;

    Y_test_orig = T_test(50,:).*T_train_S+T_train_C; % undo normalization
    Y_reproj_orig = Y_test_orig*B_trunc;% reproject onto nodal basis
    %Y_reproj_orig = (T_reproj'\Y_reproj_orig')';
    nexttile; hold on;
    PlotMeshScalarVariable([],MUA_2018,Y_reproj_orig'+T_mean');
    caxis([-1000 1000]); axis equal; axis tight;

    Ind_orig = seq(626+78+50);
    T_orig=T(Ind_orig,:);
    nexttile; hold on;
    PlotMeshScalarVariable([],MUA_2018,T_orig'+T_mean');
    caxis([-1000 1000]); axis equal; axis tight;

end