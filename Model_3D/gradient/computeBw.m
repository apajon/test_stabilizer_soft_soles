function Bw = computeBw(Fc_mc3,ind_cont,ind_slip,R,Wc,Pfree_c)
    %%% It is ok
%     delta_mc3 = Wc * Fc_mc3;
%     Fc_mc_hat = zeros(3*length(ind_cont),3);
%     delta_mc_hat = zeros(3*length(ind_cont),3);
%     
%     RPfree_c3 = R * Pfree_c; 
%     for i=1:length(ind_cont)
%         Fc_mc_hat((3*i)-2:3*i,:) = [0, -Fc_mc3(3*i), Fc_mc3(3*i-1); 
%                                  Fc_mc3(3*i), 0, -Fc_mc3(3*i-2);
%                                 -Fc_mc3(3*i-1), Fc_mc3(3*i-2), 0];
%         delta_mc_hat((3*i)-2:3*i,:) = [0, -delta_mc3(3*i), delta_mc3(3*i-1); 
%                                  delta_mc3(3*i), 0, -delta_mc3(3*i-2);
%                                 -delta_mc3(3*i-1), delta_mc3(3*i-2), 0];
%     end
%     Bw12 = Wc * Fc_mc_hat - delta_mc_hat;
%     for i=1:length(ind_cont)
%         Bw12((3*i)-2,2) = Bw12((3*i)-2,2) + RPfree_c3(3*i);
%         Bw12((3*i)-2,3) = Bw12((3*i)-2,3) - RPfree_c3((3*i)-1);
%         Bw12((3*i)-1,1) = Bw12((3*i)-1,1) - RPfree_c3(3*i);
%         Bw12((3*i)-1,3) = Bw12((3*i)-1,3) + RPfree_c3((3*i)-2);
%         Bw12((3*i),1) = Bw12((3*i),1) + RPfree_c3((3*i)-1);
%         Bw12((3*i),2) = Bw12((3*i),2) - RPfree_c3((3*i)-2);        
%     end
%     Bw11 = repmat(eye(3),length(ind_cont),1);
%     Bw = [Bw11 Bw12;zeros(3*length(ind_slip),3) zeros(3*length(ind_slip),3)];
    delta_mc3 = Wc * Fc_mc3;
    Fc_mc_hat = zeros(3*length(ind_cont),3);
    delta_mc_hat = zeros(3*length(ind_cont),3);
    
    RPfree_c3 = R * Pfree_c; 
    for i=1:length(ind_cont)
        Fc_mc_hat((3*i)-2:3*i,:) = [0, -Fc_mc3(3*i), Fc_mc3(3*i-1); 
                                 Fc_mc3(3*i), 0, -Fc_mc3(3*i-2);
                                -Fc_mc3(3*i-1), Fc_mc3(3*i-2), 0];
        delta_mc_hat((3*i)-2:3*i,:) = [0, -delta_mc3(3*i), delta_mc3(3*i-1); 
                                 delta_mc3(3*i), 0, -delta_mc3(3*i-2);
                                -delta_mc3(3*i-1), delta_mc3(3*i-2), 0];
    end
    Bw12 = Wc * Fc_mc_hat - delta_mc_hat;
    for i=1:length(ind_cont)
        Bw12((3*i)-2,2) = Bw12((3*i)-2,2) + RPfree_c3(3*i);
        Bw12((3*i)-2,3) = Bw12((3*i)-2,3) - RPfree_c3((3*i)-1);
        Bw12((3*i)-1,1) = Bw12((3*i)-1,1) - RPfree_c3(3*i);
        Bw12((3*i)-1,3) = Bw12((3*i)-1,3) + RPfree_c3((3*i)-2);
        Bw12((3*i),1) = Bw12((3*i),1) + RPfree_c3((3*i)-1);
        Bw12((3*i),2) = Bw12((3*i),2) - RPfree_c3((3*i)-2);        
    end
    Bw11 = repmat(eye(3),length(ind_cont),1);
    Bw = [Bw11 Bw12;zeros(2*length(ind_slip),3) zeros(2*length(ind_slip),3)];
end