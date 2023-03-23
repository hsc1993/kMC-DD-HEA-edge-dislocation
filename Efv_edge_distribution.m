function Efv_edge = Efv_edge_distribution(T,flag_plot)

% hypsecant probability
%     hypsecant.pdf(y) / scale with y = (x - loc) / scale
    loc = 0.9532944346350616;
    scale = 0.22462371533233122;
    syms x
    eqn = 1/pi*sech(x-loc)/scale == 0.01;
    S = solve(eqn,x);
    Efv_span_limit = [eval(S(1)),eval(S(2))];
    Efv_span = Efv_span_limit(2):0.01:Efv_span_limit(1);

    distribution = 1/pi*sech(Efv_span-loc)/scale;
    sum_distribution = sum(distribution);
    randnum = rand(1,1)*sum_distribution;
    sampling_probability = 0;
    i = 1;
    while randnum>sampling_probability
        sampling_probability = sampling_probability+distribution(i);
        i = i+1;
    end
    sampled_Efv_idx = i;
    Efv_edge = Efv_span(sampled_Efv_idx);

    Efv_mean = 0;
    for i = 1:size(distribution,2)
        Efv_mean = Efv_mean+Efv_span(i)*distribution(i);
    end
    Efv_mean = Efv_mean/sum_distribution;


    if flag_plot == 1
        subplot(1,3,3)

        hold on
        plot(Efv_span,distribution,LineWidth=2)
        scatter(Efv_edge,distribution(sampled_Efv_idx),LineWidth=2)
        xlabel('Efv [eV]', 'FontSize', 20)
        ylabel('counts', 'FontSize', 20)
        title('Efv edge distribution', 'FontSize', 20)
        l = legend(strcat('Efv edge mean = ',num2str(Efv_mean),'eV'),strcat('Efv edge sampled = ',num2str(Efv_edge),'eV'));
        l.FontSize = 14;
        hold off
    end



end