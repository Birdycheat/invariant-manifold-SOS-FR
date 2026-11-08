function result = cm_readme_figures(root, outDir)
%CM_README_FIGURES Generate the README illustrations from MATLAB and reference data.
if ~isfolder(outDir), mkdir(outDir); end
files = strings(0,1);
fig = figure('Visible','off'); style(fig,24,5.5);
titles = {'Exact inputs','Invariant jet','Reduced SOS','Rational tuple','Exact checks'};
details = {'system coefficients','restricted basis','numerical centering','affine projection','identities + positivity'};
left = linspace(0.035,0.825,5); width = 0.14; bottom = 0.39; height = 0.30;
for k=1:5
    annotation(fig,'rectangle',[left(k),bottom,width,height], ...
        'Color',[0.15,0.25,0.35],'FaceColor',[0.94,0.97,0.99],'LineWidth',0.8);
    annotation(fig,'textbox',[left(k),bottom+0.09,width,0.12], ...
        'String',titles{k},'EdgeColor','none','HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
    annotation(fig,'textbox',[left(k)-0.014,bottom-0.17,width+0.028,0.13], ...
        'String',details{k},'EdgeColor','none','HorizontalAlignment','center','FontSize',8.5);
    if k<5
        annotation(fig,'arrow',[left(k)+width,left(k+1)],[bottom+height/2,bottom+height/2], ...
            'LineWidth',0.8,'HeadLength',5,'HeadWidth',5,'Color',[0.15,0.25,0.35]);
    end
end
annotation(fig,'textbox',[0.18,0.04,0.64,0.12],'String','Solver success is separate from exact certification.', ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',9);
files = [files;saveFigure(fig,outDir,'workflow')]; close(fig);

fig=figure('Visible','off');style(fig,22,9.2);
% Compare the same invariant slice in original and jet-adapted coordinates.
layout=tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
original=nexttile(layout);adapted=nexttile(layout);
[x,y]=meshgrid(linspace(-0.92,0.92,17),linspace(-0.5,1.9,17));
z=y-2*x.^2;fx=-x.^5+x.*z;fy=-4*x.^6+(4*x.^2-1).*z;
normalize=sqrt(fx.^2+fy.^2);normalize=max(normalize,eps);
quiver(original,x,y,fx./normalize,fy./normalize,0.45,'Color',[0.72,0.74,0.76],'LineWidth',0.5);
hold(original,'on');xx=linspace(-0.96,0.96,300);
curve=plot(original,xx,2*xx.^2,'k--','LineWidth',1.25);
[x,z]=meshgrid(linspace(-0.92,0.92,17),linspace(-0.85,0.85,17));
fx=-x.^5+x.*z;fz=-z;normalize=max(sqrt(fx.^2+fz.^2),eps);
quiver(adapted,x,z,fx./normalize,fz./normalize,0.45,'Color',[0.72,0.74,0.76],'LineWidth',0.5);
hold(adapted,'on');manifoldLine=plot(adapted,xx,zeros(size(xx)),'k--','LineWidth',1.25);
initial=[-0.75,0.6;0.75,-0.5;-0.4,-0.4;0.45,0.6];
% Trajectories illustrate the coordinate change; they are not certificate tests.
for k=1:size(initial,1)
    [~,trajectory]=ode45(@(~,s)[-s(1)^5+s(1)*s(2);-s(2)],[0,8],initial(k,:).', ...
        odeset('RelTol',1e-9,'AbsTol',1e-11));
    plot(original,trajectory(:,1),trajectory(:,2)+2*trajectory(:,1).^2,'Color',[0.08,0.35,0.62],'LineWidth',1.0);
    plot(adapted,trajectory(:,1),trajectory(:,2),'Color',[0.08,0.35,0.62],'LineWidth',1.0);
    plot(original,initial(k,1),initial(k,2)+2*initial(k,1)^2,'o', ...
        'Color',[0.08,0.35,0.62],'MarkerFaceColor','white','MarkerSize',4);
    plot(adapted,initial(k,1),initial(k,2),'o','Color',[0.08,0.35,0.62], ...
        'MarkerFaceColor','white','MarkerSize',4);
end
xlabel(original,'$x_1$','Interpreter','latex');ylabel(original,'$y$','Interpreter','latex');
xlabel(adapted,'$x_1$','Interpreter','latex');ylabel(adapted,'$z=y-2x_1^2$','Interpreter','latex');
title(original,'Original coordinates','FontWeight','normal');title(adapted,'Jet-adapted coordinates','FontWeight','normal');
xlim(original,[-0.98,0.98]);ylim(original,[-0.6,2.05]);
xlim(adapted,[-0.98,0.98]);ylim(adapted,[-0.95,0.95]);
legend(original,curve,{'$y=2x_1^2$'},'Interpreter','latex','Location','southoutside');
legend(adapted,manifoldLine,{'$z=0$'},'Interpreter','latex','Location','southoutside');
style(fig,22,9.2);files=[files;saveFigure(fig,outDir,'jet_geometry')];close(fig);

counts=zeros(4,2);
data=dir(fullfile(root,'data','benchmark','*.json'));
for k=1:numel(data)
    s=jsondecode(fileread(fullfile(data(k).folder,data(k).name)));
    if ismember(s.m,[5,7]),counts(1,1)=counts(1,1)+1;end
end
for pair={{'coupled',2},{'triple',3}}
    info=pair{1};data=dir(fullfile(root,'data',info{1},'instances','*.json'));
    for k=1:numel(data)
        s=jsondecode(fileread(fullfile(data(k).folder,data(k).name)));
        if isfield(s,'leading'),classification=s.leading.class;else,classification=s.class;end
        counts(info{2},1+strcmp(classification,'unstable'))=counts(info{2},1+strcmp(classification,'unstable'))+1;
    end
end
data=dir(fullfile(root,'data','beg','instances','*.json'));
for k=1:numel(data)
    s=jsondecode(fileread(fullfile(data(k).folder,data(k).name)));
    counts(4,1+strcmp(s.class,'unstable'))=counts(4,1+strcmp(s.class,'unstable'))+1;
end
assert(isequal(counts,[80,0;34,6;32,8;11,4]),'cm:Corpus','Unexpected benchmark class counts.');
fig=figure('Visible','off');style(fig,19,8.5);ax=axes(fig);
b=bar(ax,counts,'stacked','BarWidth',0.60);b(1).FaceColor=[0.18,0.38,0.56];b(2).FaceColor=[0.82,0.84,0.86];
for k=1:4
    text(ax,k,counts(k,1)/2,sprintf('%d',counts(k,1)),'Color','white','HorizontalAlignment','center','FontSize',10);
    if counts(k,2)>0
        text(ax,k,counts(k,1)+counts(k,2)/2,sprintf('%d',counts(k,2)), ...
            'Color','black','HorizontalAlignment','center','FontSize',10);
    end
end
xticks(ax,1:4);xticklabels(ax,{'Scalar','Two centers','Three centers','BEG'});
ylabel(ax,'System instances');ylim(ax,[0,90]);
legend(ax,{'Stable','Unstable'},'Location','northeast','Interpreter','none');
style(fig,19,8.5);files=[files;saveFigure(fig,outDir,'benchmark_cases')];close(fig);
result=struct('status','passed','files',{cellstr(files)},'main_benchmark_class_counts',counts, ...
    'geometry_scope','Exact invariant slice x2=0 of the two-center example. Arrows are normalized; paths are numerical illustrations, not a certified attraction region.');
fid=fopen(fullfile(outDir,'figures.json'),'w','n','UTF-8');
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(result,PrettyPrint=true));
end

function style(fig,width,height)
fig.Units='centimeters';fig.Position=[2,2,width,height];fig.Color='white';fig.PaperPositionMode='auto';
set(fig,'DefaultAxesFontSize',9,'DefaultAxesLineWidth',0.7,'DefaultAxesTickLabelInterpreter','latex');
for ax=findall(fig,'Type','axes').'
    ax.FontSize=9;ax.LineWidth=0.7;ax.Box='on';ax.Layer='top';
    ax.XGrid='on';ax.YGrid='on';ax.GridAlpha=0.12;ax.TickLabelInterpreter='latex';
end
end

function files=saveFigure(fig,folder,name)
svg=fullfile(folder,[name '.svg']);png=fullfile(folder,[name '.png']);
print(fig,svg,'-dsvg');
print(fig,png,'-dpng','-r300');
files=string({svg;png});
end
