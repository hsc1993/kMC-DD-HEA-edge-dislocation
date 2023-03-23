function [R_total,eventlist] = calTotalRate(rn,links,typelist,mu_v,Emv_list,kb,T,work_list,Eelastic_list)   
R_total = 0;
eventlist = [];
length_eventlist = 0;
q = [];
flag_plot = 0;
iter_jog = 0;

for idx_seg = 1:size(links,1)
%% jog segment
    if typelist(idx_seg)==2
        iter_jog = iter_jog+1;
        Emv = Emv_list(iter_jog);

        % forward jump
        E_forward = Emv+work_list(idx_seg,1)+Eelastic_list(idx_seg,1)+Eelastic_list(idx_seg,2);

        q_jog = mu_v*exp(-(E_forward)/(kb*T));               % eq (5)
        q = [q;q_jog];
        event_jog = event;
        event_jog.segtype = 'jog';
        event_jog.rate = q_jog;
        event_jog.idx_seg = idx_seg;
        event_jog.action = 1;
        eventlist = [eventlist;event_jog];

        E_backward = Emv+work_list(idx_seg,2)+Eelastic_list(idx_seg,3)+Eelastic_list(idx_seg,4);
        q_jog = mu_v*exp(-(E_backward)/(kb*T));               % eq (5)
        q = [q;q_jog];
        event_jog = event;
        event_jog.segtype = 'jog';
        event_jog.rate = q_jog;
        event_jog.idx_seg = idx_seg;
        event_jog.action = 2;
        eventlist = [eventlist;event_jog];

    end

end
length_eventlist = size(eventlist,1);
% calculate total rate using rate associated with each segment
for idx_event = 1:length_eventlist
    R_total = R_total+q(idx_event);
end
% Caution: when using R_total to sample events, it is fine, but when
% calculating total rates and dt_KMC, use half of R_total becauser of
% double counting 
