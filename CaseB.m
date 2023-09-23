%% CASE B: UNIFORM WIND SPEED AND VARIABLE DIRECTION. %%%%%%%%%%%%%%%%%%%%%
%  VALENTIN OSUNA-ENCISO, CUTONALA, MARCH, 2016. %%%%%%%%%%%%%%%%%%%%%%%%%%
function PT=CaseB(Nt,maxx,maxy,rx,ry)
    % rx=maxx.*rand(1,Nt); ry=maxy.*rand(1,Nt); % TURBINE POSITIONS.
    % Nt=30;                                    % TURBINE NUMBER.
    % maxx=2000; maxy=2000;                     % WINDFARM SIZES (X,Y).
    fk=[0.02778,0.02778,0.02778,0.02778,0.02778,0.02778,0.02778,0.02778,...
    0.02778,0.02778,0.02778,0.02778,0.02778,0.02778,0.02778,0.02778,...
    0.02778,0.02778,0.02778,0.02778,0.02778,0.02778,0.02778,0.02778,...
    0.02778,0.02778,0.02778,0.02778,0.02778,0.02778,0.02778,0.02778,...
    0.02778,0.02778,0.02778,0.02778];                     
    v0=12;                                       % WIND FREE VELOCITY. 
    Vj=[]; Pii=zeros(1,Nt); i1=1;
    for teta=0:10:350
        Vj=[Vj;JENSEN_NG(Nt,v0,maxx,maxy,rx,ry,teta)];% VELOCITIES. 
        Pt=0;
        for i2=1:Nt  
            if Vj(i1,i2) > 2.3 && Vj(i1,i2) <= 12.8
                Pt=Pt+0.3*Vj(i1,i2)^3;
            elseif Vj(i1,i2) > 12.8 && Vj(i1,i2) <= 18
                Pt=Pt+630;
            end        
        end
        Pii(1,i1)=fk(1,i1)*Pt;
        i1=i1+1;
    end
    PT=sum(Pii);                                  % TOTAL GENERATED POWER.
end