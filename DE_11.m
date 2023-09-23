%% Differential Evolution, best/1/bin, Valentin Osuna-Enciso, CUCEI-UDG,
% CUTONALA, Junio, 2014
function f_best=DE_11(X,rx,ry)
    rng('shuffle');
    maxx=2000; maxy=2000;             % WINDFARM LIMITS (X,Y).
    Np=30;                            % POPULATION SIZE.
    Nt=46;                            % TURBINES NUMBER.
    F=0.6;                            % ESCALING FACTOR. .5
    Cr=0.6;                           % CROSS FACTOR..7
    maxIter=400;                       % MAX ITERATION NUMBER.
    k=0;                              % ITERATION NUMBER.
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
    %% INITIALIZATION: 
    f=zeros(Np,1);                    % POPULATION'S FITNESSES.
    for i1=1:Np
        f(i1,1)=CaseD(Nt,maxx,maxy,X(i1,1:Nt),X(i1,Nt+1:end));
    end
    EvalFO=Np; CONVERGE=ones(1,maxIter);
    [f_best,in1]= min(f); best=X(i1,:); 
    disp(sprintf('Iterac=%d, Fitnes=%.8f, Eval:%d\n',k,f_best,EvalFO));  
    while k<maxIter  
        r1=randperm(Np); r2=randperm(Np); r3=randperm(Np);
        for i1=1:Np
            v(1,:) = best(1,:) + F.*(X(r1(1,i1),:)-X(r2(1,i1),:));        % MUTANT VECTOR X
            u = X(i1,:);                               % TRIAL VECTOR X
            j_rand=randi(Nt*2);                % AT LEAST ONE GENOTYPE IS CHANGED
            for i2=1:Nt*2                         
                if(rand()<Cr || i2==j_rand) && v(1,i2)>=0 && v(1,i2)<=maxx %CROSS
                    u(1,i2)=v(1,i2);
                end
            end
            temp = CaseD(Nt,maxx,maxy,u(1,1:Nt),u(1,Nt+1:end));
            if(temp < f(i1,1))          
                X(i1,:) = u;f(i1,1) = temp; 
                if (temp < f_best)
                    f_best = temp; best = u;
                end
            end                             
        end   
        k=k+1; EvalFO=EvalFO+Np+1;
        CONVERGE(k)=f_best;  
        [f_worst,in2] = max(f); X(in2,:) = [rx(1,1:46).*rand(1,Nt),ry(1,1:46).*rand(1,Nt)];
        f(in2,1)=CaseD(Nt,maxx,maxy,X(in2,1:Nt),X(in2,Nt+1:end));
        disp(sprintf('Iterac=%d, Fitnes=%.8f, Eval:%d\n',k,f_best,EvalFO));   
    end
    %% Muestra resultados de la evolucion de DE:
    figure
    plot(best(1,1:Nt),best(1,Nt+1:end),'ks','LineWidth',3)
    grid on
    f_best=CaseC(Nt,maxx,maxy,best(1,1:Nt),best(1,Nt+1:end));
end
