% It calculates the far wake effect of Nt wind turbines (JENSEN), without 
% plotting. Valentín Osuna-Enciso, January,February,March, 2016.
function A=SHADOW(Nt,v0,maxx,maxy,rx,ry,teta)
    % INPUT DATA: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     rx=[158.739887574718,158.214499573145,108.781654711410,125.520923007868,446.461202642989,351.611612278146,277.868971359693,92.2168338788266,106.015421266160,38.6734040563384];
%     ry=[456.900205389784,353.357608848465,278.894483377438,156.714494968296,83.1017814510753,311.248629639948,493.967367476248,85.2160115284417,128.896125286007,198.399659316572];
%     Nt=10; v0=12;             % NUMBER OF TURBINES
%     maxx=500; maxy=500;       % TERRAIN LIMITS    
%     teta=180;                 % WIND ENTRANCE ANGLE
%     rx=maxx.*rand(1,Nt); ry=maxy.*rand(1,Nt);
    tetaAlterna=teta+90;    % EXPANSION OF THE WAKE
    z0=0.3;                % SURFACE ROUGHNESS
    z=60;                   % HUB HEIGHT
    r0=40;                  % ROTOR RADIUS
    alfa=(0.5)/(log(z/z0));    
    Ct=0.88;                % THRUST COEFFICIENT FOR OFFSHORE WINDFARM
    % DATA BEING CALCULATED: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Vj=v0.*ones(1,Nt);      % ENTRANCE VELOCITY TO EACH TURBINE
    M = zeros(2,Nt); B = zeros(2,Nt);    
    for i1=1:Nt 
        % ECUACION DE LA RECTA DEL CENTRO DE CADA TURBINA: %%%%%%%%%%%%%%%%
        x1=rx(1,i1);
        y1=ry(1,i1);
        m=0; b=0; %0o, 90o, 180o, 270o
        x2=x1; %Casos 90o y 270o
        if teta==270 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            y2=maxy; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif teta==90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
            y2=0;    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        else 
            m=tand(teta);
            b=y1-m*x1; %b, interseccion con el eje y 
            if teta<90
                x2=0;
                y2=b; 
                if y2<0
                    y2=0;
                    x2=(y2-b)/m;
                end
            elseif (teta>90 && teta <=180)
                x2=maxx;
                y2=m*x2 + b;
                if y2<0
                    y2=0;
                    x2=(y2-b)/m;
                end
            elseif (teta>180 && teta <270)
                x2=maxx;   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
                y2=m*x2 + b; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if y2>maxy
                    y2=maxy;
                    x2=(y2-b)/m; 
                end
            elseif (teta>270 && teta <360)
                x2 = 0;
                y2 = m*x2 + b;
                if y2>maxy
                    y2=maxy;
                    x2=(y2-b)/m; 
                end
            end                                    
        end
        M(1,i1)=m; B(1,i1)=b;
        % ECUACIÓN DE LA RECTA A 90o DE CADA TURBINA: %%%%%%%%%%%%%%%%%%%%%
        x1=rx(1,i1);
        y1=ry(1,i1);
        x2=x1; %Casos 90o y 270o
        if tetaAlterna>=360
           tetaAlterna=tetaAlterna-360; 
        end
        if tetaAlterna==270 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            y2=maxy; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif tetaAlterna==90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
            y2=0;    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        else 
            m=tand(tetaAlterna);
            b=y1-m*x1; %b, interseccion con el eje y 
            if tetaAlterna<90
                x2=0;
                y2=b; 
                if y2<0
                    y2=0;
                    x2=(y2-b)/m;
                end
            elseif (tetaAlterna>90 && tetaAlterna <=180)
                x2=maxx;
                y2=m*x2 + b;
                if y2<0
                    y2=0;
                    x2=(y2-b)/m;
                end
            elseif (tetaAlterna>180 && tetaAlterna <270)
                x2=maxx;   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
                y2=m*x2 + b; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if y2>maxy
                    y2=maxy;
                    x2=(y2-b)/m; 
                end
            elseif (tetaAlterna > 270 && tetaAlterna <360)
                x2 = 0;
                y2 = m*x2 + b;
                if y2>maxy
                    y2=maxy;
                    x2=(y2-b)/m; 
                end
            end                                    
        end
        M(2,i1)=m; B(2,i1)=b;        
        %dibujo aspas de la turbina:
        angAspa=180-tetaAlterna; %para 10o funciona
        y3=r0*sind(angAspa)+y1; y2=y1-r0*sind(angAspa);
        x2=r0*cosd(angAspa)+x1; x3=x1-r0*cosd(angAspa);
    end
    %CON LAS ECUACIONES DE LAS RECTAS, CALCULAR LOS PUNTOS DE INTERSECCIÓN:
    Xij=zeros(Nt,Nt); rij=Xij; Dij=Xij; Ashadowi=Xij; A0=pi*(r0^2);    
    J=zeros(Nt,Nt+1); %Matriz para calcular las velocidades
    I=zeros(Nt,Nt+1);
    for i1=1:Nt
        for i2=1:Nt
         if i1~=i2   
            x=(B(2,i2)-B(1,i1))/(M(1,i1)-M(2,i2));
            y=B(1,i1)+M(1,i1)*x;
            %%% 190o - 260o; 280o - 350o; 10o - 80o; 100o - 170o;
           if (((teta >180 && teta <=260)||(teta >270 && teta <=350)) && ((y-M(2,i1)*x)>=B(2,i1)) && y>=0 && x>=0 && y<=maxy && x<=maxx)||(((teta >0 && teta <=80)||(teta>90 && teta<=170)) && ((y-M(2,i1)*x)<=B(2,i1)) && y>=0 && x>=0 && y<=maxy && x<=maxx)
                Xij(i1,i2)=sqrt((x-rx(i1))^2+(y-ry(i1))^2);
                Dij(i1,i2)=sqrt((x-rx(i2))^2+(y-ry(i2))^2);
                rij(i1,i2)=r0+alfa*Xij(i1,i2);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if rij(i1,i2)>=Dij(i1,i2)+r0 %Afecta totalmente:
                    Ashadowi(i1,i2)=A0;
                    J(i2,end)=J(i2,end)+1;
                    I(i1,end)=I(i1,end)+1;
                    J(i2,J(i2,end))=i1;
                    I(i1,I(i1,end))=i2;
                    %dibujo afectación de la estela:
                    angAspa=180-tetaAlterna; %para 10o funciona
                    y3=rij(i1,i2)*sind(angAspa)+y; y2=y-rij(i1,i2)*sind(angAspa);
                    x2=rij(i1,i2)*cosd(angAspa)+x; x3=x-rij(i1,i2)*cosd(angAspa);
                    %afecta parcialmente:
                elseif (rij(i1,i2)<Dij(i1,i2)+r0) && (rij(i1,i2)>=Dij(i1,i2)-r0)
                    tetaw=acosd((Dij(i1,i2)^2+rij(i1,i2)^2-r0^2)/(2*Dij(i1,i2)*rij(i1,i2)));
                    Lij=rij(i1,i2)*cosd(tetaw);
                    zij=2*rij(i1,i2)*sind(tetaw);
                    Ashadowi(i1,i2)=rij(i1,i2)^2*cosd(Lij/rij(i1,i2))+...
                        r0^2*cosd((Dij(i1,i2)-Lij)/rij(i1,i2))-...
                        Dij(i1,i2)*zij;
                    J(i2,end)=J(i2,end)+1;
                    I(i1,end)=I(i1,end)+1;
                    J(i2,J(i2,end))=i1;
                    I(i1,I(i1,end))=i2;
                    %dibujo afectación de la estela:
                    angAspa=180-tetaAlterna; %para 10o funciona
                    y3=rij(i1,i2)*sind(angAspa)+y; y2=y-rij(i1,i2)*sind(angAspa);
                    x2=rij(i1,i2)*cosd(angAspa)+x; x3=x-rij(i1,i2)*cosd(angAspa);
                elseif (rij(i1,i2)<Dij(i1,i2)-r0) % No afecta
                    Ashadowi(i1,i2)=0;                    
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           elseif (teta==180 && rx(i2)>rx(i1) && rx(i2)>=0 && rx(i2)<=maxx) || (teta==0 && rx(i2)<rx(i1) && rx(i2)>=0 && rx(i2)<=maxx)
                y=ry(i1); %OKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKO
                x=rx(i2); %OKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKO
                Xij(i1,i2)=sqrt((x-rx(i1))^2+(y-ry(i1))^2);
                Dij(i1,i2)=sqrt((x-rx(i2))^2+(y-ry(i2))^2);
                rij(i1,i2)=r0+alfa*Xij(i1,i2); 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if rij(i1,i2)>=Dij(i1,i2)+r0 %Afecta totalmente:
                    Ashadowi(i1,i2)=A0;
                    J(i2,end)=J(i2,end)+1;
                    I(i1,end)=I(i1,end)+1;
                    J(i2,J(i2,end))=i1;
                    I(i1,I(i1,end))=i2;
                    %dibujo afectación de la estela:
                    angAspa=180-tetaAlterna; %para 10o funciona
                    y3=rij(i1,i2)*sind(angAspa)+y; y2=y-rij(i1,i2)*sind(angAspa);
                    x2=rij(i1,i2)*cosd(angAspa)+x; x3=x-rij(i1,i2)*cosd(angAspa);
                    %afecta parcialmente:
                elseif (rij(i1,i2)<Dij(i1,i2)+r0) && (rij(i1,i2)>=Dij(i1,i2)-r0)
                    tetaw=acosd((Dij(i1,i2)^2+rij(i1,i2)^2-r0^2)/(2*Dij(i1,i2)*rij(i1,i2)));
                    Lij=rij(i1,i2)*cosd(tetaw);
                    zij=2*rij(i1,i2)*sind(tetaw);
                    Ashadowi(i1,i2)=rij(i1,i2)^2*cosd(Lij/rij(i1,i2))+...
                        r0^2*cosd((Dij(i1,i2)-Lij)/rij(i1,i2))-...
                        Dij(i1,i2)*zij;
                    J(i2,end)=J(i2,end)+1;
                    I(i1,end)=I(i1,end)+1;
                    J(i2,J(i2,end))=i1;
                    I(i1,I(i1,end))=i2;
                    %dibujo afectación de la estela:
                    angAspa=180-tetaAlterna; %para 10o funciona
                    y3=rij(i1,i2)*sind(angAspa)+y; y2=y-rij(i1,i2)*sind(angAspa);
                    x2=rij(i1,i2)*cosd(angAspa)+x; x3=x-rij(i1,i2)*cosd(angAspa);
                elseif (rij(i1,i2)<Dij(i1,i2)-r0) % No afecta
                    Ashadowi(i1,i2)=0;                   
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           elseif ((teta==270 && ry(i2)>ry(i1) && ry(i2)>=0 && ry(i2)<=maxy) || (teta==90 && ry(i2)<ry(i1) && ry(i2)>=0 && ry(i2)<=maxy))
                y=ry(i2); %OKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKO
                x=rx(i1); %OKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKO
                Xij(i1,i2)=sqrt((x-rx(i1))^2+(y-ry(i1))^2);
                Dij(i1,i2)=sqrt((x-rx(i2))^2+(y-ry(i2))^2);
                rij(i1,i2)=r0+alfa*Xij(i1,i2);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if rij(i1,i2)>=Dij(i1,i2)+r0 %Afecta totalmente:
                    Ashadowi(i1,i2)=A0;
                    J(i2,end)=J(i2,end)+1;
                    I(i1,end)=I(i1,end)+1;
                    J(i2,J(i2,end))=i1;
                    I(i1,I(i1,end))=i2;
                    %dibujo afectación de la estela:
                    angAspa=180-tetaAlterna; %para 10o funciona
                    y3=rij(i1,i2)*sind(angAspa)+y; y2=y-rij(i1,i2)*sind(angAspa);
                    x2=rij(i1,i2)*cosd(angAspa)+x; x3=x-rij(i1,i2)*cosd(angAspa);
                    %afecta parcialmente:
                elseif (rij(i1,i2)<Dij(i1,i2)+r0) && (rij(i1,i2)>=Dij(i1,i2)-r0)
                    tetaw=acosd((Dij(i1,i2)^2+rij(i1,i2)^2-r0^2)/(2*Dij(i1,i2)*rij(i1,i2)));
                    Lij=rij(i1,i2)*cosd(tetaw);
                    zij=2*rij(i1,i2)*sind(tetaw);
                    Ashadowi(i1,i2)=rij(i1,i2)^2*cosd(Lij/rij(i1,i2))+...
                        r0^2*cosd((Dij(i1,i2)-Lij)/rij(i1,i2))-...
                        Dij(i1,i2)*zij;
                    J(i2,end)=J(i2,end)+1;
                    I(i1,end)=I(i1,end)+1;
                    J(i2,J(i2,end))=i1;
                    I(i1,I(i1,end))=i2;
                    %dibujo afectación de la estela:
                    angAspa=180-tetaAlterna; %para 10o funciona
                    y3=rij(i1,i2)*sind(angAspa)+y; y2=y-rij(i1,i2)*sind(angAspa);
                    x2=rij(i1,i2)*cosd(angAspa)+x; x3=x-rij(i1,i2)*cosd(angAspa);
                elseif (rij(i1,i2)<Dij(i1,i2)-r0) % No afecta
                    Ashadowi(i1,i2)=0;                  
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           end
         end
        end        
    end
    %%% CALCULUS OF THE VELOCITY ENTERING TO EACH TURBINE, CASE 0:
%     Cr=1-sqrt(1-Ct);
%     for i1=1:Nt         % AFFECTED TURBINE.
%        deficit=0;
%        for i2=1:Nt      % TURBINES AFFECTING.
%            if (rij(i2,i1) && Ashadowi(i2,i1))
%                 deficit=deficit+(Cr)*((r0/rij(i2,i1))^2)*(Ashadowi(i2,i1)/A0);               
%            end
%        end
%        if deficit>1 CASO 0
%           deficit=1;
%        end
%        Vj(1,i1)=Vj(1,i1)*(1-deficit);
%     end
    %%% CALCULUS OF THE VELOCITY ENTERING TO EACH TURBINE, CASE 1:
%     Cr=1-sqrt(1-Ct);
%     for i1=1:Nt         % AFFECTED TURBINE.
%        for i2=1:Nt      % TURBINES AFFECTING.
%            if (rij(i2,i1) && Ashadowi(i2,i1))
%                 deficit=(Cr)*((r0/rij(i2,i1))^2)*(Ashadowi(i2,i1)/A0);
%                 Vj(1,i1)=Vj(1,i1)*(1-deficit);
%            end
%        end
%     end
%     A=sum(sum(Ashadowi));
      %%% CALCULUS OF THE VELOCITY ENTERING TO EACH TURBINE, CASE 2:
        Cr=1-sqrt(1-Ct);
    for i1=1:Nt         % AFFECTED TURBINE.
       deficit=0;
       cont=0;
       for i2=1:Nt      % TURBINES AFFECTING.
           if (rij(i2,i1) && Ashadowi(i2,i1))
                deficit=deficit+(Cr)*((r0/rij(i2,i1))^2)*(Ashadowi(i2,i1)/A0)+cont; 
                cont=cont+0.01;
           end
       end
%        if deficit>1 
%           deficit=1;
%        end
       Vj(1,i1)=(deficit);
    end
    A=sum(Vj);
end   