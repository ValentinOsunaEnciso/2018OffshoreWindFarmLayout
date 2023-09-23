function E = DE_13(X,Nt,maxx,maxy)
% Output arguments:
% bestmem              : parameter vector with f_best solution
% EvalFO               : number of function evaluations
% Run DE minimization; strategy: DE/rand/1 ; ui = pm3 + F*(pm1 - pm2);
    rng('shuffle');
    Np=30; %Nt=46; 
    
    D=Nt*2; F=0.86; Cr=0.15; 
    maxIter=200; k = 0;
    CONVERGE=ones(1,maxIter);
    %maxx=2000; maxy=2000; 
%     [X,Y] = meshgrid(0:300:2000, 0:300:2000); % TURBINES POSITIONS.
%     rx=[]; ry=[]; tam=15;
%     for i1=1:tam
%         rx=[rx;X(1,:),X(2,:),X(3,:),X(4,:),X(5,:),X(6,:),X(7,:)];
%         ry=[ry;Y(1,:),Y(2,:),Y(3,:),Y(4,:),Y(5,:),Y(6,:),Y(7,:)]; 
%     end
%     rx=rx(:,1:46).*rand(tam,Nt); ry=ry(:,1:46).*rand(tam,Nt);
%     for i1=tam+1:Np
%         rx=[rx;maxx.*rand(1,Nt)];
%         ry=[ry;maxy.*rand(1,Nt)];
%     end
%     X=[rx,ry];
    Xold    = zeros(size(X));     % toggle population
    f       = zeros(1,Np);          % Create and reset the "cost array"
    bestmem   = zeros(1,D);           % f_best population member ever
    best = zeros(1,D);           % f_best population member in kation
    EvalFO    = 0;                    % number of function evaluations
    ibest   = 1;                      % start with first population member
    f(1)  = CaseD(Nt,maxx,maxy,X(ibest,1:Nt),X(ibest,Nt+1:end))*...
        (1+fitness8(X(ibest,1:Nt),X(ibest,Nt+1:end),40,Nt));%100*(X(ibest,2)-X(ibest,1)^2)^2 + (1-X(ibest,1))^2;
    f_best = f(1);                 % f_best objective function value so far
    EvalFO  = EvalFO + 1;
    for i1=2:Np                        % check the remaining members
        f(i1) = CaseD(Nt,maxx,maxy,X(i1,1:Nt),X(i1,Nt+1:end))*...
        (1+fitness8(X(i1,1:Nt),X(i1,Nt+1:end),40,Nt));
        EvalFO  = EvalFO + 1;
        if (f(i1) < f_best)           % if member is better
            ibest   = i1;                 % save its location
            f_best = f(i1);
        end   
    end
    best = X(ibest,:);         % f_best member of current kation
    bestmem = best;              % f_best member ever
    pm1 = zeros(Np,D);              % initialize population matrix 1
    pm2 = zeros(Np,D);              % initialize population matrix 2
    pm3 = zeros(Np,D);              % initialize population matrix 3
    pm4 = zeros(Np,D);              % initialize population matrix 4
    pm5 = zeros(Np,D);              % initialize population matrix 5
    a3  = zeros(Np);                % index array
    a4  = zeros(Np);                % index array
    a5  = zeros(Np);                % index array
    bm  = zeros(Np,D);              % initialize bestmember  matrix
    u  = zeros(Np,D);              % intermediate population of perturbed vectors
    mui = zeros(Np,D);              % mask for intermediate population
    mpo = zeros(Np,D);              % mask for old population
    rot = (0:1:Np-1);               % rotating index array
    rt  = zeros(Np);                % another rotating index array
    a1  = zeros(Np);                % index array
    a2  = zeros(Np);                % index array
    ind = zeros(4);    
    while (k < maxIter)
        Xold = X;                   % save the old population
        ind = randperm(4);              % index pointer array
        a1  = randperm(Np);             % shuffle locations of vectors
        rt = rem(rot+ind(1),Np);        % rotate indices by ind(1) positions
        a2  = a1(rt+1);                 % rotate vector locations   
        rt = rem(rot+ind(2),Np);
        a3  = a2(rt+1);                
        rt = rem(rot+ind(3),Np);
        a4  = a3(rt+1);               
        rt = rem(rot+ind(4),Np);
        a5  = a4(rt+1);   
        pm1 = Xold(a1,:);             % shuffled population 1
        pm2 = Xold(a2,:);             % shuffled population 2
        pm3 = Xold(a3,:);             % shuffled population 3
        pm4 = Xold(a4,:);             % shuffled population 4
        pm5 = Xold(a5,:);             % shuffled population 5
        for i1=1:Np                      % population filled with the f_best member
            bm(i1,:) = best;          % of the last iteration
        end
        mui = rand(Np,D) < Cr;          % all random numbers < Cr are 1, 0 otherwise
        mpo = mui < 0.5;                % inverse mask to mui
        u = pm3 + F*(pm1 - pm2);       % differential variation
        u = Xold.*mpo + u.*mui;     % binomial Crossover               
        for i1=1:Np % VECTORS ENTERING TO THE NEW POPULATION:     %%%%%%%%%
            for i2=1:Nt                        % X^k+1                                      
                if (u(i1,i2)<0)
                    u(i1,i2)=0;                    
                elseif u(i1,i2)>maxx
                    u(i1,i2)=maxx;
                end
             end
             for i2=Nt+1:2*Nt                   % Y^k+1                                         
                if (u(i1,i2)<0)
                    u(i1,i2)=0;                    
                elseif u(i1,i2)>maxy
                    u(i1,i2)=maxy;
                end
             end
            temp = CaseD(Nt,maxx,maxy,u(i1,1:Nt),u(i1,Nt+1:end))*...
                    (1+fitness8(u(i1,1:Nt),u(i1,Nt+1:end),40,Nt));%100*(u(i1,2)-u(i1,1)^2)^2 + (1-u(i1,1))^2;   % check cost of competitor
            EvalFO  = EvalFO + 1;
            if (temp <= f(i1))  % if competitor is better than value in "cost array"
                X(i1,:) = u(i1,:);  % replace old vector with new one (for new kation)
                f(i1)   = temp;  % save value in "cost array"
                if (temp < f_best)     % if competitor better than the f_best one ever
                    f_best = temp;     % new f_best value
                    bestmem = u(i1,:);     % new f_best parameter vector ever
                end
            end
        end 
        best = bestmem;       % freeze the f_best member of this iteration 
        fprintf(1,'Iteration: %d,  Best: %.6f,  F: %.2f,  Cr: %.2f,  Np: %d\n',k,f_best,F,Cr,Np);
        k = k + 1;
        CONVERGE(k)=f_best;
    end %---end while ((k < maxIter)
%     figure
%     plot(best(1,1:Nt),best(1,Nt+1:end),'gs','LineWidth',3)
%     grid on
%     figure
%     plot(CONVERGE)
    E=CaseC(Nt,maxx,maxy,best(1,1:Nt),best(1,Nt+1:end));
end
