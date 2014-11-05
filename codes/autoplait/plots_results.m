function [  ] = plots_results( dir )
%PLOTS_RESULTS.M 
% 
% written by Yasuko Matsubara 
%


    %% load original data and segment information
    fid = fopen([dir, 'input']);
    fn = textscan(fid, '%s');
    fclose(fid);
    orgfn=fn{1}{1}
    datfn=[dir, 'segment.'];

    %% load org data
    org=load(orgfn);
    org=normalize(org)*10;
    st=1;
    ed=length(org)-1;
    miny=min(min(org));
    maxy=max(max(org));

    %% load label list
    fid = fopen([datfn, 'labels']);
    labels = textscan(fid, '%s %s %s %s %s');
    G=size(labels{1},1);
    fclose(fid);

    %% plot segments
    figure;
    subplot(2,1,1);
    %largefont;
    plot(org);
    
    plot_segments(datfn, G, st, ed, miny, maxy);

end

function plot_segments(datfn, G, st, ed, miny, maxy)
    w=0;
    %% for each segment set 's'
    for i=0: G-1

        cols=get(0,'DefaultAxesColorOrder');
        col=cols(mod(i,size(cols,1))+1,:);
        
        %% plot original data
        X=load([datfn, num2str(i)]);
        subplot(2,1,1)
        plot_lines(X(:,1));
        xlim([st,ed])
        ylim([miny, maxy*1.2])

        %% plot colorful segments
        subplot(G*2,1,G+i+1)
        plot_box(X, col, 1)
        xlim([st,ed])
        ylim([0 1])
        box on
        ylabel(num2str(i+1));       
        set(gca,'XTickLabel',{''})
        set(gca,'YTickLabel',{''})
        w=w+size(X,1);
    end

    disp(['w=', num2str(w)])

end


function [  ] = plot_box( X, color, alpha )
    n=size(X,1);
    for i=1:n
        hold on;
        plot_box_aux(X(i,1), X(i,2), color, alpha);
        hold off;
    end
end

function [  ] = plot_box_aux( st, ed, color, alpha )
    btm=-100000;
    top=100000;
    p=patch([st,ed, ed, st],[btm,btm, top, top], color);
    set(p,'EdgeColor','none')
    set(p,'FaceAlpha',alpha);
    ylim([btm, top])
end


function [  ] = plot_lines( X )
    btm=-1000;
    top=1000;
    hold on;    
    for i=1: length(X)
        plot([X(i), X(i)], [btm, top], 'black-')
    end
    ylim([btm, top])
end

function [ Xfull ] = normalize( X )
    Xfull=X;
    for i=1: size(Xfull,2)
        X=Xfull(:,i);
        Xn=(X-min(X))/(max(X)-min(X));
        Xfull(:,i) = Xn;
    end
end


