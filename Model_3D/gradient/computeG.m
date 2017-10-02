function G = computeG(ind_cont,ind_slip,Fc_mc,P_mc,Zdes)
%%% It is ok
%     G = zeros(6,3*length(ind_cont)+3*length(ind_slip));
%     G11 = repmat(eye(3),1,length(ind_cont));
%     G21 = zeros(3,3*length(ind_cont));
%     for i=1:length(ind_cont)
%         G21(:,(3*i)-2:3*i) = [0 0 P_mc(2,i)-Zdes(2);0 0 -(P_mc(1,i)-Zdes(1));-P_mc(2,i)+Zdes(2) P_mc(1,i)-Zdes(1) 0];
%     end
%     G22 = zeros(3,3*length(ind_slip));
%     I = [];
%     for i=1:length(ind_slip)
%         I(i) = find(ind_cont==ind_slip(i));
%     end
%     Fc_mc_splip = Fc_mc(:,I);
%     for i=1:length(ind_slip)
%         G22(:,(3*i)-2:3*i) = [0 Fc_mc_splip(3,i) 0;Fc_mc_splip(3,i) 0 0;Fc_mc_splip(2,i) -Fc_mc_splip(1,i) 0];
%     end
%     G = [G11 zeros(3,3*length(ind_slip));G21 G22];
    G11 = repmat(eye(3),1,length(ind_cont));
    G21 = zeros(3,3*length(ind_cont));
    for i=1:length(ind_cont)
        G21(:,(3*i)-2:3*i) = [0 0 P_mc(2,i)-Zdes(2);0 0 -(P_mc(1,i)-Zdes(1));-P_mc(2,i)+Zdes(2) P_mc(1,i)-Zdes(1) 0];
    end
    G22 = zeros(3,2*length(ind_slip));
    if ~isempty(ind_slip)
        I = [];
        for i=1:length(ind_slip)
            I(i) = find(ind_cont==ind_slip(i));
        end
        Fc_mc_splip = Fc_mc(:,I);
        for i=1:length(ind_slip)
            G22(:,(2*i)-1:2*i) = [0 Fc_mc_splip(3,i);Fc_mc_splip(3,i) 0;Fc_mc_splip(2,i) -Fc_mc_splip(1,i)];
        end
    end
    G = [G11 zeros(3,2*length(ind_slip));G21 G22];
end