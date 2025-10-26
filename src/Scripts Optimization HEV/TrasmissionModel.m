%% Modello di trasmissione
% Si occupa di trasferire la potenza generata dai motori (ICE ed EM) alle ruote attraverso il cambio/frizione e il differenziale.
% Di seguito é definito il rapporto di trasmissione meccanica in funzione della marcia innestata.

function [G,SelGear] = TrasmissionModel(i_fd,i_gb,n_gb,vel_cycle)

G = ones(size(vel_cycle));
SelGear = zeros(size(vel_cycle));
for i=1:numel(vel_cycle)-1
    if vel_cycle(i) > 0 && vel_cycle(i) <= 15/3.6
        SelGear(i) = n_gb(1);
        G(i) = i_fd*i_gb(1);
    elseif vel_cycle(i) > 15/3.6 && vel_cycle(i) <= 35/3.6
        SelGear(i) = n_gb(2);
        G(i) = i_fd*i_gb(2);
    elseif vel_cycle(i) > 35/3.6 && vel_cycle(i) <= 50/3.6
        SelGear(i) = n_gb(3);
        G(i) = i_fd*i_gb(3);
    elseif vel_cycle(i) > 50/3.6 && vel_cycle(i) <= 70/3.6
        SelGear(i) = n_gb(4);
        G(i) = i_fd*i_gb(4);
    elseif vel_cycle(i) > 70/3.6 && vel_cycle(i) <= 90/3.6
        SelGear(i) = n_gb(5);
        G(i) = i_fd*i_gb(5);
    elseif vel_cycle(i) > 90/3.6
        SelGear(i) = n_gb(6);
        G(i) = i_fd*i_gb(6);
    else
        SelGear(i) = 0;
        G(i) = 1;
    end
end

end


