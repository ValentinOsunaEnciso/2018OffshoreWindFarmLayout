%% CASE C: NON-UNIFORM WIND SPEED AND VARIABLE DIRECTION. %%%%%%%%%%%%%%%%%
%  VALENTIN OSUNA-ENCISO, CUTONALA, MARCH, 2016. %%%%%%%%%%%%%%%%%%%%%%%%%%
function PT=CaseC(Nt,maxx,maxy,rx,ry)
    % rx=maxx.*rand(1,Nt); ry=maxy.*rand(1,Nt); % TURBINE POSITIONS.
    % Nt=30;                                    % TURBINE NUMBER.
    % maxx=2000; maxy=2000;                     % WINDFARM SIZES (X,Y).
    fk1=[0.0042,0.0042,0.0042,0.0042,0.0042,0.0042,0.0042,0.0042,0.0042,...
    0.0042,0.0042,0.0042,0.0042,0.0042,0.0042,0.0042,0.0042,0.0042,...
    0.0042,0.0042,0.0042,0.0042,0.0042,0.0042,0.0042,0.0042,0.0042,...
    0.0042,0.0042,0.0042,0.0042,0.0042,0.0042,0.0042,0.0042,0.0042]; %u0=8 
    fk2=[0.0084,0.0084,0.0084,0.0084,0.0084,0.0084,0.0084,0.0084,0.0084,...
    0.0084,0.0084,0.0084,0.0084,0.0084,0.0084,0.0084,0.0084,0.0084, ...
    0.0084,0.0084,0.0084,0.0084,0.0084,0.0084,0.0084,0.0084,0.0084, ...
    0.00107,0.0126,0.0149,0.0149,0.0195,0.0149,0.0149,0.0126,0.0102];%u0=12
    fk3=[0.0112,0.0112,0.0112,0.0112,0.0112,0.0112,0.0112,0.0112,0.0112,...
    0.0112,0.0112,0.0112,0.0112,0.0112,0.0112,0.0112,0.0112,0.0112, ...
    0.0112,0.0112,0.0112,0.0112,0.0112,0.0112,0.0112,0.0112,0.0112, ...
    0.0135,0.0163,0.0191,0.0302,0.0358,0.0307,0.0191,0.0163,0.0135]; %u0=17 
    v01=8; v02=12;  v03=17;                      % WIND FREE VELOCITIES.
    %% PRIMER CALCULO:
    Vj1=[];  Vj2=[]; Vj3=[]; Pt1=zeros(36,Nt); Pt2=zeros(36,Nt);
    Pt3=zeros(36,Nt); 
    for teta=0:10:350
        Vj1=[Vj1;JENSEN_NG(Nt,v01,maxx,maxy,rx,ry,teta)];% VELOCITIES.. 
        Vj2=[Vj2;JENSEN_NG(Nt,v02,maxx,maxy,rx,ry,teta)];% VELOCITIES.
        Vj3=[Vj3;JENSEN_NG(Nt,v03,maxx,maxy,rx,ry,teta)];% VELOCITIES.
    end
    for i1=1:36
        for i2=1:Nt 
            if Vj1(i1,i2) > 2.3 && Vj1(i1,i2) <= 12.8
                Pt1(i1,i2)=Pt1(i1,i2)+0.3*Vj1(i1,i2)^3;
            elseif Vj1(i1,i2) > 12.8 && Vj1(i1,i2) <= 18
                Pt1(i1,i2)=Pt1(i1,i2)+630;
            end 
            if Vj2(i1,i2) > 2.3 && Vj2(i1,i2) <= 12.8
                Pt2(i1,i2)=Pt2(i1,i2)+0.3*Vj2(i1,i2)^3;
            elseif Vj2(i1,i2) > 12.8 && Vj2(i1,i2) <= 18
                Pt2(i1,i2)=Pt2(i1,i2)+630;
            end 
            if Vj3(i1,i2) > 2.3 && Vj3(i1,i2) <= 12.8
                Pt3(i1,i2)=Pt3(i1,i2)+0.3*Vj3(i1,i2)^3;
            elseif Vj3(i1,i2) > 12.8 && Vj3(i1,i2) <= 18
                Pt3(i1,i2)=Pt3(i1,i2)+630;
            end         
        end  
        Pt1(i1,:)=Pt1(i1,:).*fk1(1,i1);
        Pt2(i1,:)=Pt2(i1,:).*(fk2(1,i1)+fk1(1,i1));
        Pt3(i1,:)=Pt3(i1,:).*(fk3(1,i1)+fk2(1,i1)+fk1(1,i1));
    end
    PT1=ones(1,Nt); PT2=ones(1,Nt); PT3=ones(1,Nt);
    for i1=1:Nt
       PT1(1,i1)=sum(Pt1(:,i1)); 
       PT2(1,i1)=sum(Pt2(:,i1));
       PT3(1,i1)=sum(Pt3(:,i1));
    end
    
    PT=sum(PT1)+sum(PT2)+sum(PT3);                % TOTAL GENERATED POWER.
end