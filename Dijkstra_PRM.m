clear all
close all
global N NumObs
N = 100;
NumObs = 15;
nObsX = 15;
nObsY = 10;
NumNodes = 100;
Thrs = 50;                      %thresold
Cobs  = ones(nObsX,nObsY);

[G]= GeneratePRM(Thrs,NumNodes,Cobs);
[PathOrder,Flag_Good] = Dijkstra(G);
if Flag_Good == 0
    disp('No Path Found')
    return
end

Goal = find(PathOrder(:,1) == G.Xg);
ii = 1;
Prev(1) = 0;
while(Prev~=G.Xs)
    Prev(ii) = PathOrder(Goal,3);
    Goal = Prev(ii);
    ii=ii+1;
end
NodesOrder = flip(Prev);
NodesOrder(end+1) = G.Xg;

for ii =1:length(NodesOrder)
    X(ii) = G.Node(1,1,NodesOrder(ii));
    Y(ii) = G.Node(1,2,NodesOrder(ii));
end
figure(1)
plot(X,Y,'--r','LineWidth',4)
title('Dijkstra - Shortest Path')
xlabel('X')
ylabel('Y')
saveas(gcf,'Dijkstra.fig')

%%
% clear all
% clc
% load('Geasy3')
% [PathOrder] = Dijkstra(G);
% % for ii =1:length(PathOrder)
% %     X(ii) = G.Node(1,1,PathOrder(ii));
% %     Y(ii) = G.Node(1,2,PathOrder(ii));
% % end
% figure(1)
% % plot(X,Y,'--r','LineWidth',2)


function[Closed,Flag_Good] = Dijkstra(G)
%%%INIT
%[Node Idx,Dist from start,PrevIdx]
Xs = [G.Xs 0 G.Xs];
Xg = G.Xg;
Flag_Good = 0;
OpenArray           = [];
Closed              = [];
OpenArray(1,:)      = Xs; %Idx of first node
while(~isempty(OpenArray))
    [MinDist2Start,I] = min((OpenArray(:,2)));
    BestXidx = OpenArray(I,1);
    Closed(BestXidx,:) = OpenArray(I,:);
    OpenArray(I,:) = [];
    Neighbers = G.DistMat(BestXidx,:);                         %finding nieghbores
    NeighbersIdx = find( (Neighbers~=0)&(Neighbers~=-1) );
    for zz = 1:length(NeighbersIdx)
        if any(Closed(:,1)== NeighbersIdx(zz))
            continue
        end
        %  if ~ismember(NeighbersIdx(zz),OpenArray) && Visited(NeighbersIdx(zz)) == false
        NewDist = G.DistMat(BestXidx,NeighbersIdx(zz)) + MinDist2Start;
        if ismember(NeighbersIdx(zz),OpenArray(:,1))
            II = find(OpenArray(:,1)==NeighbersIdx(zz));
            if NewDist < OpenArray(II,2)
                OpenArray(II,2) = NewDist;
                OpenArray(II,3) = BestXidx;
            end
        else
            OpenArray(end+1,:) = [NeighbersIdx(zz),NewDist,BestXidx];
        end
        
    end
    if BestXidx == Xg
        Flag_Good = 1;
        break
    end
end

if Flag_Good
    disp('Yay')
else
    disp('Fail')
end

end

function [G]= GeneratePRM(Thrs,NumNodes,Cobs)
global N NumObs
Map = zeros(N,N);
[nObsX,nObsY] = size(Cobs);
n = 0;
ObsX = zeros(NumObs,4);
ObsY = zeros(NumObs,4);
while(n<NumObs)
    ObsVertx   = [randi([1 N-nObsX],1,1) ; randi([1 N-nObsY],1,1)];
    if nnz(Map(ObsVertx(1):(ObsVertx(1)-1+nObsX),ObsVertx(2):(ObsVertx(2)-1+nObsY))) == 0
        Map(ObsVertx(1):(ObsVertx(1)-1+nObsX),ObsVertx(2):(ObsVertx(2)-1+nObsY)) = Cobs;
        n = n+1;
        ObsX(n,:) = [ObsVertx(1)-1+nObsX,ObsVertx(1),ObsVertx(1),ObsVertx(1)-1+nObsX];
        ObsY(n,:) = [ObsVertx(2),ObsVertx(2),ObsVertx(2)-1+nObsY,ObsVertx(2)-1+nObsY];
    end
end

n = 0;
while(n<NumNodes)
    NodesVertx = [randi([1 N],1,1) ; randi([1 N],1,1)];
    if Map(NodesVertx(1),NodesVertx(2)) == 0
        Map(NodesVertx(1),NodesVertx(2)) = 2;
        n = n + 1;
        Node(:,:,n) = NodesVertx';
    end
    
end

DistMat = -ones(NumNodes,NumNodes);
for jj = 1:NumNodes
    refNode = Node(:,:,jj);
    for zz = 1:NumNodes
        dist = norm(refNode - Node(:,:,zz));
        if dist<=Thrs
            for qq = 1:NumObs
                [xi,yi] = polyxpoly([refNode(1),Node(1,1,zz)],[refNode(2),Node(1,2,zz)],[ObsX(qq,:) ObsX(qq,1)],[ObsY(qq,:) ObsY(qq,1)]);
                %[a,b] = polyxpoly([refNode(1),Node(1,1,zz)],[refNode(2),Node(1,2,zz)],[ObsX(qq,:) ObsX(qq,1)],[ObsY(qq,:) ObsY(qq,1)]);
                if isempty(xi)~=1
                    break
                end
            end
            if isempty(xi) && isempty(yi)
                DistMat(jj,zz) = dist;
            end
        end
    end
end
%%PLOT%%
figure(1)
hold on
scatter(Node(1,1,:),Node(1,2,:))
% for ii=1:NumNodes
%     text(Node(1,1,ii),Node(1,2,ii),num2str(ii),'FontSize',15)
% end
for pp = 1:NumObs
    p = patch(ObsX(pp,:),ObsY(pp,:),'red');
end
for ii=1:NumNodes
    for jj=1:NumNodes
        if DistMat(ii,jj)~=-1
            plot([Node(1,1,ii),Node(1,1,jj)],[Node(1,2,ii),Node(1,2,jj)],'k');
        end
    end
end
for ii = 1:NumNodes
    DistCurr(ii) = norm([Node(1,1,ii) Node(1,2,ii)]);
end
[~,Xs] = min(DistCurr);
[~,Xg] = max(DistCurr);
figure(1)
plot(Node(1,1,Xs),Node(1,2,Xs),'x','LineWidth',5,'MarkerSize',10)
hold on
plot(Node(1,1,Xg),Node(1,2,Xg),'x','LineWidth',5,'MarkerSize',10)

G = struct;
G.Map     = Map;
G.DistMat = DistMat;
G.Node    = Node;
G.ObsX    = ObsX;
G.ObsY    = ObsY;
G.Xs      = Xs;
G.Xg      = Xg;
end

