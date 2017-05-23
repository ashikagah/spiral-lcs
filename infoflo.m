function [uo,vo] = infoflo(mat)

% Calculates information flow based on binary time series
% 
% INPUT:    
%   mat         ... 2-D binary time series [N x M x time] (e.g. bmt)
%
% OUTPUT:
%   [uo,vo]     ... 2-D time series of information flow vector field (each pixel as a source)

javaaddpath('./JIDT/infodynamics-dist-1.3/infodynamics.jar');   % Change path as needed

uo = zeros(size(mat,1),size(mat,2));
vo = zeros(size(mat,1),size(mat,2));

% Define vectors from the center as information source to Moore neighbors 

y1 = [ 0; 1]; % from X to Y1            % Moore neighbors numbered counterclockwise
y2 = [ 1; 1]; % from X to Y2    
y3 = [ 1; 0]; % from X to Y3            %   Y4  -  Y3  -  Y2
y4 = [ 1;-1]; % from X to Y4            %   |      |      |
y5 = [ 0;-1]; % from X to Y5            %   Y5  -  X   -  Y1
y6 = [-1;-1]; % from X to Y6            %   |      |      |
y7 = [-1; 0]; % from X to Y7            %   Y6  -  Y7  -  Y8
y8 = [-1; 1]; % from X to Y8

for m = 2:size(mat,1)-1                 % m = row axis
	for n = 2:size(mat,2)-1             % n = column axis
% Define arrays
        X = squeeze(mat(m,n,:));        % Information source
        Y1 = squeeze(mat(m,n+1,:));     % Destination 1
        Y2 = squeeze(mat(m+1,n+1,:));   % Destination 2
        Y3 = squeeze(mat(m+1,n,:));     % Destination 3
        Y4 = squeeze(mat(m+1,n-1,:));   % Destination 4
        Y5 = squeeze(mat(m,n-1,:));     % Destination 5
        Y6 = squeeze(mat(m-1,n-1,:));   % Destination 6
        Y7 = squeeze(mat(m-1,n,:));     % Destination 7
        Y8 = squeeze(mat(m-1,n+1,:));   % Destination 8

% Calculate transfer entropy (TE) from information source to each destination
        teCalc=javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete',2,1);
        teCalc.initialise(); teCalc.addObservations(X,Y1); teXY1 = teCalc.computeAverageLocalOfObservations(); % TE X->Y1
        teCalc.initialise(); teCalc.addObservations(X,Y2); teXY2 = teCalc.computeAverageLocalOfObservations(); % TE X->Y2
        teCalc.initialise(); teCalc.addObservations(X,Y3); teXY3 = teCalc.computeAverageLocalOfObservations(); % TE X->Y3
        teCalc.initialise(); teCalc.addObservations(X,Y4); teXY4 = teCalc.computeAverageLocalOfObservations(); % TE X->Y4
        teCalc.initialise(); teCalc.addObservations(X,Y5); teXY5 = teCalc.computeAverageLocalOfObservations(); % TE X->Y5
        teCalc.initialise(); teCalc.addObservations(X,Y6); teXY6 = teCalc.computeAverageLocalOfObservations(); % TE X->Y6
        teCalc.initialise(); teCalc.addObservations(X,Y7); teXY7 = teCalc.computeAverageLocalOfObservations(); % TE X->Y7
        teCalc.initialise(); teCalc.addObservations(X,Y8); teXY8 = teCalc.computeAverageLocalOfObservations(); % TE X->Y8
        vec = teXY1*y1 + teXY2*y2 + teXY3*y3 + teXY4*y4 + teXY5*y5 + teXY6*y6 + teXY7*y7 + teXY8*y8;
        uo(m,n) = vec(2);
        vo(m,n) = vec(1);
    end                                 
end

m = 1; n = 1;                           % Left top corner of 2-D lattice
X = squeeze(mat(m,n,:));                % Information source
Y1 = squeeze(mat(m,n+1,:));             % Destination 1
Y2 = squeeze(mat(m+1,n+1,:));           % Destination 2
Y3 = squeeze(mat(m+1,n,:));             % Destination 3
teCalc=javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete',2,1);
teCalc.initialise(); teCalc.addObservations(X,Y1); teXY1 = teCalc.computeAverageLocalOfObservations(); % TE X->Y1
teCalc.initialise(); teCalc.addObservations(X,Y2); teXY2 = teCalc.computeAverageLocalOfObservations(); % TE X->Y2
teCalc.initialise(); teCalc.addObservations(X,Y3); teXY3 = teCalc.computeAverageLocalOfObservations(); % TE X->Y3
vec = teXY1*y1 + teXY2*y2 + teXY3*y3; 
uo(m,n) = vec(2);
vo(m,n) = vec(1);

m = 1; n = size(mat,2);                 % Right top corner of 2-D lattice
X = squeeze(mat(m,n,:));                % Information source
Y3 = squeeze(mat(m+1,n,:));             % Destination 3
Y4 = squeeze(mat(m+1,n-1,:));           % Destination 4
Y5 = squeeze(mat(m,n-1,:));             % Destination 5
teCalc=javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete',2,1);
teCalc.initialise(); teCalc.addObservations(X,Y3); teXY3 = teCalc.computeAverageLocalOfObservations(); % TE X->Y3
teCalc.initialise(); teCalc.addObservations(X,Y4); teXY4 = teCalc.computeAverageLocalOfObservations(); % TE X->Y4
teCalc.initialise(); teCalc.addObservations(X,Y5); teXY5 = teCalc.computeAverageLocalOfObservations(); % TE X->Y5
vec = teXY3*y3 + teXY4*y4 + teXY5*y5; 
uo(m,n) = vec(2);
vo(m,n) = vec(1);

m = size(mat,1); n = 1;                 % Left bottom corner of 2-D lattice
X = squeeze(mat(m,n,:));                % Information source
Y1 = squeeze(mat(m,n+1,:));             % Destination 1
Y7 = squeeze(mat(m-1,n,:));             % Destination 7
Y8 = squeeze(mat(m-1,n+1,:));           % Destination 8
teCalc=javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete',2,1);
teCalc.initialise(); teCalc.addObservations(X,Y1); teXY1 = teCalc.computeAverageLocalOfObservations(); % TE X->Y1
teCalc.initialise(); teCalc.addObservations(X,Y7); teXY7 = teCalc.computeAverageLocalOfObservations(); % TE X->Y7
teCalc.initialise(); teCalc.addObservations(X,Y8); teXY8 = teCalc.computeAverageLocalOfObservations(); % TE X->Y8
vec = teXY1*y1 + teXY7*y7 + teXY8*y8; 
uo(m,n) = vec(2);
vo(m,n) = vec(1);

m = size(mat,1); n = size(mat,2);       % Right bottom corner of 2-D lattice
X = squeeze(mat(m,n,:));                % Information source 
Y5 = squeeze(mat(m,n-1,:));             % Destination 5
Y6 = squeeze(mat(m-1,n-1,:));           % Destination 6
Y7 = squeeze(mat(m-1,n,:));             % Destination 7
teCalc=javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete',2,1);
teCalc.initialise(); teCalc.addObservations(X,Y5); teXY5 = teCalc.computeAverageLocalOfObservations(); % TE X->Y5
teCalc.initialise(); teCalc.addObservations(X,Y6); teXY6 = teCalc.computeAverageLocalOfObservations(); % TE X->Y6
teCalc.initialise(); teCalc.addObservations(X,Y7); teXY7 = teCalc.computeAverageLocalOfObservations(); % TE X->Y7
vec = teXY5*y5 + teXY6*y6 + teXY7*y7; 
uo(m,n) = vec(2);
vo(m,n) = vec(1);

m = 1;                                  % m = row axis; Top border of 2-D lattice
for n=2:size(mat,2)-1                   % n = column axis
    X = squeeze(mat(m,n,:));            % Information source
    Y1 = squeeze(mat(m,n+1,:));         % Destination 1
    Y2 = squeeze(mat(m+1,n+1,:));       % Destination 2
    Y3 = squeeze(mat(m+1,n,:));         % Destination 3
    Y4 = squeeze(mat(m+1,n-1,:));       % Destination 4
    Y5 = squeeze(mat(m,n-1,:));         % Destination 5
    teCalc=javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete',2,1);
    teCalc.initialise(); teCalc.addObservations(X,Y1); teXY1 = teCalc.computeAverageLocalOfObservations(); % TE X->Y1
    teCalc.initialise(); teCalc.addObservations(X,Y2); teXY2 = teCalc.computeAverageLocalOfObservations(); % TE X->Y2
    teCalc.initialise(); teCalc.addObservations(X,Y3); teXY3 = teCalc.computeAverageLocalOfObservations(); % TE X->Y3
    teCalc.initialise(); teCalc.addObservations(X,Y4); teXY4 = teCalc.computeAverageLocalOfObservations(); % TE X->Y4
    teCalc.initialise(); teCalc.addObservations(X,Y5); teXY5 = teCalc.computeAverageLocalOfObservations(); % TE X->Y5
    vec = teXY1*y1 + teXY2*y2 + teXY3*y3 + teXY4*y4 + teXY5*y5; 
    uo(m,n) = vec(2);
    vo(m,n) = vec(1);
end

n = 1;                                  % n = column axis; Left border of 2-D lattice
for m = 2:size(mat,1)-1                 % m = row axis
    X = squeeze(mat(m,n,:));            % Information source
    Y1 = squeeze(mat(m,n+1,:));         % Destination 1
    Y2 = squeeze(mat(m+1,n+1,:));       % Destination 2
    Y3 = squeeze(mat(m+1,n,:));         % Destination 3
    Y7 = squeeze(mat(m-1,n,:));         % Destination 7
    Y8 = squeeze(mat(m-1,n+1,:));       % Destination 8
    teCalc=javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete',2,1);
    teCalc.initialise(); teCalc.addObservations(X,Y1); teXY1 = teCalc.computeAverageLocalOfObservations(); % TE X->Y1
    teCalc.initialise(); teCalc.addObservations(X,Y2); teXY2 = teCalc.computeAverageLocalOfObservations(); % TE X->Y2
    teCalc.initialise(); teCalc.addObservations(X,Y3); teXY3 = teCalc.computeAverageLocalOfObservations(); % TE X->Y3
    teCalc.initialise(); teCalc.addObservations(X,Y7); teXY7 = teCalc.computeAverageLocalOfObservations(); % TE X->Y7
    teCalc.initialise(); teCalc.addObservations(X,Y8); teXY8 = teCalc.computeAverageLocalOfObservations(); % TE X->Y8
    vec = teXY1*y1 + teXY2*y2 + teXY3*y3 + teXY7*y7 + teXY8*y8; 
    uo(m,n) = vec(2);
    vo(m,n) = vec(1);
end

m = size(mat,1);                        % m = row axis; Bottom border of 2-D lattice
for n = 2:size(mat,2)-1                 % n = column axis
    X = squeeze(mat(m,n,:));            % Information source 
    Y1 = squeeze(mat(m,n+1,:));         % Destination 1
    Y5 = squeeze(mat(m,n-1,:));         % Destination 5
    Y6 = squeeze(mat(m-1,n-1,:));       % Destination 6
    Y7 = squeeze(mat(m-1,n,:));         % Destination 7
    Y8 = squeeze(mat(m-1,n+1,:));       % Destination 8
    teCalc=javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete',2,1);
    teCalc.initialise(); teCalc.addObservations(X,Y1); teXY1 = teCalc.computeAverageLocalOfObservations(); % TE X->Y1
    teCalc.initialise(); teCalc.addObservations(X,Y5); teXY5 = teCalc.computeAverageLocalOfObservations(); % TE X->Y5
    teCalc.initialise(); teCalc.addObservations(X,Y6); teXY6 = teCalc.computeAverageLocalOfObservations(); % TE X->Y6
    teCalc.initialise(); teCalc.addObservations(X,Y7); teXY7 = teCalc.computeAverageLocalOfObservations(); % TE X->Y7
    teCalc.initialise(); teCalc.addObservations(X,Y8); teXY8 = teCalc.computeAverageLocalOfObservations(); % TE X->Y8
    vec = teXY1*y1 + teXY5*y5 + teXY6*y6 + teXY7*y7 + teXY8*y8; 
    uo(m,n) = vec(2);
    vo(m,n) = vec(1);
end

n = size(mat,2);                        % n = column axis; Right border of 2-D lattice
for m = 2:size(mat,1)-1                 % m = row axis
    X = squeeze(mat(m,n,:));            % Information source
    Y3 = squeeze(mat(m+1,n,:));         % Destination 3
    Y4 = squeeze(mat(m+1,n-1,:));       % Destination 4
    Y5 = squeeze(mat(m,n-1,:));         % Destination 5
    Y6 = squeeze(mat(m-1,n-1,:));       % Destination 6
    Y7 = squeeze(mat(m-1,n,:));         % Destination 7
    teCalc=javaObject('infodynamics.measures.discrete.TransferEntropyCalculatorDiscrete',2,1);
    teCalc.initialise(); teCalc.addObservations(X,Y3); teXY3 = teCalc.computeAverageLocalOfObservations(); % TE X->Y3
    teCalc.initialise(); teCalc.addObservations(X,Y4); teXY4 = teCalc.computeAverageLocalOfObservations(); % TE X->Y4
    teCalc.initialise(); teCalc.addObservations(X,Y5); teXY5 = teCalc.computeAverageLocalOfObservations(); % TE X->Y5
    teCalc.initialise(); teCalc.addObservations(X,Y6); teXY6 = teCalc.computeAverageLocalOfObservations(); % TE X->Y6
    teCalc.initialise(); teCalc.addObservations(X,Y7); teXY7 = teCalc.computeAverageLocalOfObservations(); % TE X->Y7
    vec = teXY3*y3 + teXY4*y4 + teXY5*y5 + teXY6*y6 + teXY7*y7; 
    uo(m,n) = vec(2);
    vo(m,n) = vec(1);
end
           