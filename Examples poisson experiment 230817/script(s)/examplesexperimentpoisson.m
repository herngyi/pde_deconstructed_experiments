%EXPERIMENT I - (DATE OF UPLOAD)
%ANALYSIS OF PROFICIENCY AND CONVERGENCE OF METHODS ON DIFFERENT MESHES
%
%   The following scripts runs all our previously coded methods for the
%resolution of the Poisson equation on deconstructed domains on various
%examples of meshes (see example.m function) and compares them using L2
%norm against the solution one would get from remeshing. It outputs seven
%different log/log plots showing the error and the rate at which this
%converges.
%
%   In order to test convergence, we will take maximum triangle areas which
%decrease as inverse powers of two, and then take the average edge length 
%as a measure of resolution, which corresponds to the horizontal axis of the
%plots.
%
%
%
%
%
%
%
exampless={'two disks';'snake on tube';'two rectangles';'disk on disk';'disk on annulus';'venus flytrap';'pair of pants'};
examples={'two_disks';'snake_on_tube';'two_rectangles';'disk_on_disk';'disk_on_annulus';'venus_flytrap';'pair_of_pants'};
for j=1:length(examples)
    l2naive=[];
    l2none=[];
    %l2smooth=[];
    l2delta=[];
    l2neumann=[];
    h=[];
    for i=1:8
        [TVA,TFA,TNA,na,TVB,TFB,TNB,nb,TH]=example('disk_on_disk',.1/(2^i));
        [l2naiv,l2delt,l2non,l2neuman]=comparison_poisson(TVA,TFA,TNA,na,TVB,TFB,TNB,nb,@(x,y) x.*0,@(x,y) [x,x],TH);
        l2naive=[l2naive,l2naiv];
        l2none=[l2none,l2non];
        %l2smooth=[l2smooth,l2smoot];
        l2neumann=[l2neumann,l2neuman];
        l2delta=[l2delta,l2delt];
        h=[h,avgedge(TVA,TFA)];
    end
    figure(j)
    disp(j)
    plot(log(h),log(l2naive),log(h),log(l2delta),log(h),log(l2neumann),log(h),log(l2none));
    axis equal
    legend('OSC','DSC','NSC','none')
    xlabel('log avg edge length')
    ylabel('log l2 error')
    str=['Error vs resolution in ',exampless(j,:)];
    title(str);
end
