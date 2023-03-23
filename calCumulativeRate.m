function cumulative_ratelist = calCumulativeRate(eventlist)

length_eventlist = size(eventlist,1);
cumulative_ratelist = zeros(length_eventlist,1);

% calculate total rate using rate associated with each segment
for idx_event = 1:length_eventlist
    if idx_event==1
        cumulative_ratelist(idx_event) = eventlist(idx_event).rate;
    else
        cumulative_ratelist(idx_event) = cumulative_ratelist(idx_event-1)+eventlist(idx_event).rate;
    end
end

