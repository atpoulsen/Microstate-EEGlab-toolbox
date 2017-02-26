function topo_micro(A,chanlocs)
%% Initialisation
K = size(A,2);
% scrsz = get(groot,'ScreenSize');
% figure('Position',[scrsz(3)/10 scrsz(4)/3 scrsz(3)/2+scrsz(4)/3*K scrsz(4)/3])


%% Plot
for k = 1:K
    subplot(1,K,k)
    topoplot(A(:,k)',chanlocs);
end