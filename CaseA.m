%% CASE A: UNIFORM WIND SPEED AND UNIFORM DIRECTION. %%%%%%%%%%%%%%%%%%%%%%
%  VALENTIN OSUNA-ENCISO, CUTONALA, MARCH, 2016. %%%%%%%%%%%%%%%%%%%%%%%%%%
function PT=CaseA(Nt,maxx,maxy,rx,ry)
    % rx=maxx.*rand(1,Nt); ry=maxy.*rand(1,Nt); % TURBINE POSITIONS.
    % Nt=30;                                    % TURBINE NUMBER.
    % maxx=2000; maxy=2000;                     % WINDFARM SIZES (X,Y).
    fk=1;                                       % WIND PROBABILITY.
    teta=0;                                     % WIND DIRECTION.    
    v0=12;                                      % WIND FREE VELOCITY.    
    Vj=JENSEN_NG(Nt,v0,maxx,maxy,rx,ry,teta);   % VELOCITIES.
    PT=0;                                       % TOTAL POWER.
    for i1=1:Nt  
        if Vj(1,i1) > 2.3 && Vj(1,i1) <= 12.8
            PT=PT+0.3*Vj(1,i1)^3;
        elseif Vj(1,i1) > 12.8 && Vj(1,i1) <= 18
            PT=PT+630;
        end            
    end
    PT=PT*fk;                                   % TOTAL GENERATED POWER.
end