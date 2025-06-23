function gr = movgrowthrate2(t,a,badSegmentations,K,fittype,debugflag)
%  movgrowthrate   Moving cell growth rate value.
%     gr = movgrowthrate(c,K) for an MCell object c and odd positive integer 
%     scalar K computes a centered moving growth rate by sliding a window 
%     of length K along cell detection frames. Each element of gr is the 
%     local growth rate of the corresponding values of X inside
%     the window, with gr of the length c.lifeTime.  The sliding
%     window is truncated at the endpoints where there are fewer than K
%     elements to fill the window.
%
% Output:
%   grs - array of moving growth rates.
%   frames - array of frames.
%   a - array of cell areas
%% Parse parameters
if length(K)==1
    NB = ceil((K-1)/2);
    NF = floor((K-1)/2);
elseif length(K)==2
    NB = K(1);
    NF = K(2);
else
    error('Parameter K must have at most two values');
end
if nargin<5
    fittype='exp1';
end
if nargin<6
    debugflag=0;
end
a=a(:);
t=t(:);
%t=t-t(1);
gr = nan(size(a));
a(badSegmentations>0) = nan;
if debugflag==0
    for f=1:length(a)
        fw = max(f-NB,1) : min(f+NF,length(a));
        fa = a(fw);
        ft = t(fw);
        ind = ~isnan(fa);
        fa = fa(ind);
        ft = ft(ind);
        if numel(fa) >= max(3,(NF+NB)/2)
            try
                if strcmp(fittype,'exp1')
                    fitObj = fit(ft-ft(1),fa,'exp1');
                    gr(f) = fitObj.b;
%                     fitObj = fit(ft-ft(1),log(fa),'poly1');
%                     gr(f) = fitObj.p1;
                else
                    fitObj = fit(ft-ft(1),fa,'poly1');
                    gr(f) = fitObj.p1/fitObj.p2;
                end
            catch me
                disp(me);
            end
        end
    end
else
    figure,tiledlayout(1,2)
    ax1=nexttile(1);
    ax2=nexttile(2);
    for f=1:length(a)
        fw = max(f-NB,1) : min(f+NF,length(a));
        fa = a(fw);
        ft = t(fw);
        ind = ~isnan(fa);
        fa = fa(ind);
        ft = ft(ind);
        if numel(fa) >= 3%max(3,(NF+NB)/2)
            try
                if strcmp(fittype,'exp1')
                    fitObj = fit(ft,fa,'exp1');
                    gr(f) = fitObj.b;
                else
                    fitObj = fit(ft,fa,'poly1');
                    gr(f) = fitObj.p1/(fitObj.p2+fitObj.p1*ft(1));
                end
                cla(ax1)
                hold(ax1,'on')
                axes(ax1);
                plot(fitObj,t,a)
                plot(ft,fa,'bo');
                plot(ax2,t,gr,'-*')
                legend(ax1,'Location','best')
                pause()
            catch me
                disp(me);
            end
        end
    end
end