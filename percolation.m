%% Advanced Engineering Analysis - ME 711 %%
%  Project #1  %
clc
clear all
close all

%% Model Parameters

L = 60;                     %Side length of the large square
t = 5;                      %Side length of the mini squares
r = 0.5*sqrt(2)*t;          %Radius of circle circumscribing mini square
Nt = 50;                    %Number of mini-squares 

%% Initialize Variables
l = 0;                      %Initialize length
finished = 0;               %Initialize finished
N = 0;                      %Initialize number of squares
path = 0;                   %Initialize path
ax = ones(1000,6);          %Preallocate ax
ay = ones(1000,6);          %Preallocate ay
dist = zeros(1000,1000);    %Preallocate dist
intersect = dist;           %Preallocate intersect
cluster = zeros(200,200);   %Preallocate cluster
rowcol2 = zeros(1,2);       %Preallocate rowcol2
rc = [];                    %Preallocate rc
Lrc = 0;                    %Preallocate Lrc

%% Generate Large Square Domain
rectangle('Position', [0,0,L,L])    %Plot large square
axis([-5,L+5,-5,L+5])               %Define axis limits
hold on                             %Keep displaying large square

%% Generate Mini Squares And Check For Path
while path == 0 %N <= Nt            %Loop until path exists (switch to commented while if a specific Nt is desired)
    N = N+1;                        %Iterate (increase number of mini squares)
    disp(N)                         %Display the number of squares
    row = [];                       %Initialize row
    col = [];                       %Initialize col
    cx = L*rand;                    %X coordinate of mini square center
    cy = L*rand;                    %Y coordinate of mini square center
    th = 90*rand;                   %Orientation of mini square center
    x = r*cosd(th);                 %X component of mini square's Q1 point
    y = r*sind(th);                 %Y component of mini square's Q1 point
    ax(N,:) = cx+[x,-y,-x,y,x,0];   %X coordinates of mini square vertices
    ay(N,:) = cy+[y,x,-y,-x,y,0];   %Y coordinates of mini square vertices
    %% Compute Distances Between Squares And If They Intersect
    for n = 1:N                     %For every square
        dist(N,n) = sqrt((ax(n,6)-ax(N,6))^2+(ay(n,6)-ay(N,6))^2);  %Compute their center's distance to the newest square's center
        dist(n,N) = dist(N,n);                                      %Distance from a to b = distance from b to a
        if dist(N,n) > sqrt(2)*t || dist(N,n) == 0  %If their distance is sufficiently far apart                             
            intersect(N,n) = 0;                     %They do not intersect    
            intersect(n,N) = 0;                     %a to b = b to a
        elseif dist(N,n) <= t && dist(N,n) ~= 0     %If their distance is sufficiently close
            intersect(N,n) = 1;                     %They do intersect
            intersect(n,N) = 1;                     %a to b = b to a
        else                                        %If their distance is "in the middle"
            [xint,yint]  = polyxpoly(ax(N,1:5),ay(N,1:5),ax(n,1:5),ay(n,1:5));  %Display intersection points (if they exist)
            if isempty(xint)                        %If there are no intersection points
                intersect(N,n) = 0;                 %They do not intersect
                intersect(n,N) = 0;                 %a to b = b to a
            else                                    %If there are intersection points
                intersect(N,n) = 1;                 %They do intersect
                intersect(n,N) = 1;                 %a to b = b to a
            end
        end
    end
    if ismember(1,intersect(:,N))                   %If the newest square intersects with any other squares
        [row,col] = find(intersect(:,N));           %Find the square(s) it intersects with (its pair(s))
        col(:) = N;                                 %Columns represent the newest square's number
    end
    if isempty(rc)                                  %If rc has not had any pairs added to it yet
        rc = [row,col];                             %Add the first pair(s)
    else                                            %If rc has existing pairs in it    
        rc = [rc;[row,col]];                        %Add newly found pairs to it
    end
    rowcol = rc';                                   %Transpose rc into a 2 x ... matrix
    %% Sort Squares Into Clusters
    if numel(rowcol) >= Lrc & Lrc > 0                           %If more pairs have been added
        for n = Lrc:2:numel(rowcol)                             %For new entries of pairs to rowcol
            if isempty(find(cluster == rowcol(n),1)) && ...
                    isempty(find(cluster == rowcol(n-1),1))     %If square and its pair are not in a cluster yet
                r0 = find(cluster(:,1) == 0,1);                 %Identify an empty row for a new cluster
                cluster(r0,1) = rowcol(n);                      %Place square in first column of new cluster
                cluster(r0,2) = rowcol(n-1);                    %Add its pair to the cluster
                continue
            elseif isempty(find(cluster == rowcol(n),1)) && ...
                    ismember(rowcol(n-1),cluster)               %If only the square is not in a cluster yet
                [r0,~] = find(cluster == rowcol(n-1));          %Row(s) of the pair
                r0 = unique(r0)';                               %Get rid of duplicate rows and transpose to row vector
                for n0 = r0(1:end)                              %For all rows the pair is in
                    c0 = find(cluster(n0,:) == 0,1);            %First empty column in that row
                    cluster(n0,c0) = rowcol(n);                 %Add the square to the cluster
                end
                continue
            elseif isempty(find(cluster == rowcol(n-1),1)) && ...
                    ismember(rowcol(n),cluster)                 %If only the pair is not in a cluster yet
                [r0,~] = find(cluster == rowcol(n));            %Row(s) of the square
                r0 = unique(r0)';                               %Get rid of duplicate rows and transpose to row vector
                for n0 = r0(1:end)                              %For all rows the square is in
                    c0 = find(cluster(n0,:) == 0,1);            %First empty column in that row
                    cluster(n0,c0) = rowcol(n-1);               %Add its pair to the cluster
                end
                continue
            else                                                %If both are already in a cluster
                [r0,~] = find(cluster == rowcol(n));            %Row(s) of the square
                [r1,~] = find(cluster == rowcol(n-1));          %Row(s) of the pair
                r0 = unique(r0)';                               %Get rid of duplicate rows and transpose to row vector
                r1 = unique(r1)';                               %Get rid of duplicate rows and transpose to row vector
                for n0 = r0(1:end)                              %For all rows the square is in
                    c0 = find(cluster(n0,:) == 0,1);            %First empty column in that row
                    cluster(n0,c0) = rowcol(n-1);               %Add its pair to the cluster
                end
                for n0 = r1(1:end)                              %For all rows the pair is in
                    c0 = find(cluster(n0,:) == 0,1);            %First empty column in that row
                    cluster(n0,c0) = rowcol(n);                 %Add the square to the cluster
                end
            end
        end
    end
    Lrc = numel(rowcol)+2;                                      %Initialize length for next loop
    
    %% Combine Clusters
    for n = 1:N  
        [rows,~] = find(cluster == n);                      %Rows of the duplicate entries
        if numel(rows) > 1                                  %If the square has duplicate entries in the cluster matrix
            rows = unique(rows)';                           %List every cluster square n is in, in order
            if numel(rows) > 1                              %If in multiple clusters
                radd = rows(1);                             %Lowest # cluster that it is in
                for n0 = rows(2:end)                        %For every other cluster it's in
                    Lrow = find(cluster(n0,:) == 0,1)-1;    %Length of the cluster
                    if numel(Lrow) == 0                     %If length is not defined this way
                        Lrow = size(cluster,2);             %Length is the columns in cluster
                    end
                    Lradd = find(cluster(radd,:) == 0,1)-1; %Length of the lowest # cluster
                    if numel(Lradd) == 0                    %If length is not defined this way
                        Lradd = size(cluster,2);            %Length is the columns in cluster
                    end
                    if Lradd+Lrow >= size(cluster,2)
                        add = Lradd+Lrow-size(cluster,2)+10;
                        zr = zeros(size(cluster,1),add);
                        cluster = [cluster,zr];
                    end
                    cluster(radd,Lradd+1:Lradd+Lrow) = ...          
                        cluster(n0,1:Lrow);                 %Combine the cluster to the lowest # cluster
                    cluster(n0,:) = zeros;                  %Overwrite the cluster with zeros
                end
            end
        end
    end 
    
    %% Remove No Longer Existing Clusters
    ze = find(cluster(:,1) == 0);                           %Rows in cluster matrix with only zeros
    ze = flipud(ze)';                                       %Reverse the order of ze
    for n = ze(1:end)                                       %For all rows of only zeros
        cluster(n,:) = [];                                  %Remove the row
        z = zeros(1,numel(cluster(1,:)));                   %Create a row vector of zeros
        cluster = [cluster;z];                              %Add the zero vector to preserve the size of the cluster matrix
    end
    
    %% Remove Duplicate Squares And Sort Clusters
    Nclusters = numel(unique(cluster(:,1)));                %Number of clusters
    if ismember(0,unique(cluster(:,1))) == 1                %If there is a zero in the unique matrix
        Nclusters = Nclusters-1;                            %Number of clusters is one less (zero not included)
    end
    for n = 1:Nclusters                                     %For every cluster
        ucluster = unique(cluster(n,:));                    %Matrix with all squares of the cluster sorted and not repeating
        if ismember(0,ucluster)                             %If there is a zero in the unique matrix
            ucluster(ucluster == 0) = [];                   %Remove the zero (it doesn't represent a square)
        end
        cluster(n,:) = zeros;                               %Remove values from the cluster
        cluster(n,1:numel(ucluster)) = ucluster;            %Replace with sorted unique values
    end
    
    %% Determine If Path Exists
    ymax = ones(1,N);                                       %Initialize ymax
    ymin = ones(1,N);                                       %Initialize ymin
    if N > L/(sqrt(2)*t)                                    %If there are enough squares to theoretically span L
        for n = 1:Nclusters                                 %For every cluster
            Lcluster = numel(unique(cluster(n,:)));         %Number of squares in the cluster
            if ismember(0,unique(cluster(n,:)))             %If there is a zero in the unique matrix
                Lcluster = Lcluster-1;                      %Length is one less (zero not included)
            end
            if Lcluster > 1                                 %If there are multiple squares in the cluster
                ymax(n) = max(max(ay(cluster(n,1:Lcluster),:)))';   %Fill in max value for each cluster
                if ymax(n) > L                              %If max value is outside the domain
                    ymax(n) = L;                            %Make the max value L
                end
                ymin(n) = min(min(ay(cluster(n,1:Lcluster),:)))';   %Fill in min value for each cluster
                if ymin(n) < 0                              %If min value is outside the domain
                    ymin(n) = 0;                            %Make the min value 0
                end
            end
        end
        LengthClusters = ymax-ymin;                         %Length of each cluster
        Length = max(LengthClusters);                       %Largest length of a cluster
        if Length >= L                                      %If the length spans the domain
            BlueCluster = find(LengthClusters == L,2);      %Identify the percolating cluster
            path = 1;                                       %A path exists
        end
    end
    %% Plot Squares While Code Runs
    for ny = 1:4
        plot([ax(N,ny), ax(N,ny+1)],[ay(N,ny), ay(N,ny+1)],'r')
        %text(ax(N,1),ay(N,1),num2str(N))
        pause(0.01)
        hold on
    end
    
end


%% Plot Squares, Highlight Path In Blue, Highlight Final Link In Green
for nx = 1:N
    for ny = 1:4
        plot([ax(nx,ny), ax(nx,ny+1)],[ay(nx,ny), ay(nx,ny+1)],'r')
    end
end
hold on

NelementsBlue = numel(unique(cluster(BlueCluster,:)));
if ismember(0,unique(cluster(BlueCluster,:)))
    NelementsBlue = NelementsBlue-1;
end
for n = 1:NelementsBlue
    nx = cluster(BlueCluster,n);
    for ny = 1:4
        plot([ax(nx,ny), ax(nx,ny+1)],[ay(nx,ny), ay(nx,ny+1)],'b')
    end
end
hold on

for ny = 1:4
    plot([ax(N,ny), ax(N,ny+1)],[ay(N,ny), ay(N,ny+1)],'g')
end
hold on

patch([L,L+t,L+t,L],[-t,-t,L+t,L+t],'w','EdgeColor','w')
patch([-t,L+t,L+t,-t],[L,L,L+t,L+t],'w','EdgeColor','w')
patch([-t,0,0,-t],[-t,-t,L+t,L+t],'w','EdgeColor','w')
patch([-t,L+t,L+t,-t],[-t,-t,0,0],'w','EdgeColor','w')
rectangle('Position', [0,0,L,L])