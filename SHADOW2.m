% It calculates the far wake effect of Nt wind turbines (JENSEN), without 
% plotting. Valentín Osuna-Enciso, January,February,March, 2016.
function A=SHADOW2(Nt,v0,maxx,maxy,rx,ry,teta)
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
            y2=1.5*maxy; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif teta==90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
            y2=-0.5*maxy;    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        else 
            m=tand(teta);
            b=y1-m*x1; %b, interseccion con el eje y 
            if teta<90
                x2=-0.5*maxx;
                y2=m*x2+b; 
            elseif (teta>90 && teta <=180)
                x2=1.5*maxx;
                y2=m*x2 + b;
            elseif (teta>180 && teta <270)
                x2=1.5*maxx;   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
                y2=m*x2 + b; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif (teta>270 && teta <360)
                x2 = -0.5*maxx;
                y2 = m*x2 + b;
            end                                    
        end
        M(1,i1)=m; B(1,i1)=b;
%         plot([x1,x2],[y1,y2],'b-')%en vez de dibujarlo, saco la distancia  
        % ECUACIÓN DE LA RECTA A 90o DE CADA TURBINA: %%%%%%%%%%%%%%%%%%%%%
        x3=rx(1,i1);
        y3=ry(1,i1);
        x4=x3; %Casos 90o y 270o
        if tetaAlterna>=360
           tetaAlterna=tetaAlterna-360; 
        end
        if tetaAlterna==270 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            y4=maxy; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elseif tetaAlterna==90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
            y4=0;    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        else 
            m=tand(tetaAlterna);
            b=y3-m*x3; %b, interseccion con el eje y 
            y4=m*x4 + b;                                              
        end
        M(2,i1)=m; B(2,i1)=b;        
        %plot([x1,x2],[y1,y2],'r-')%en vez de dibujarlo, saco la distancia
        %dibujo aspas de la turbina y estela:
        angAspa=180-tetaAlterna; %para 10o funciona
        y5=r0*sind(angAspa)+y3; y4=y3-r0*sind(angAspa);
        x4=r0*cosd(angAspa)+x3; x5=x3-r0*cosd(angAspa); 
        %plot([x4,x5],[y4,y5],'r-','linewidth',6) %AEROGENERADOR
        xij=sqrt((x1-x2)^2+(y1-y2)^2);        
        rij=r0+alfa*xij;
        y7=rij*sind(angAspa)+y2; y6=y2-rij*sind(angAspa);
        x6=rij*cosd(angAspa)+x2; x7=x2-rij*cosd(angAspa); 
%         plot([x6,x7],[y6,y7],'k-','linewidth',2) %BASE ESTELA
%         plot([x4,x6],[y4,y6],'k-','linewidth',2) %LADO 1 ESTELA
%         plot([x5,x7],[y5,y7],'k-','linewidth',2) %LADO 2 ESTELA        
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
           if(((teta >180 && teta <=260)||(teta >270 && teta <=350))...
                   && ((y-M(2,i1)*x)>=B(2,i1)))||(((teta >0 && teta <=80)...
                   ||(teta>90 && teta<=170)) && ((y-M(2,i1)*x)<=B(2,i1)))    
                Xij(i1,i2)=sqrt((x-rx(i1))^2+(y-ry(i1))^2);
                Dij(i1,i2)=sqrt((x-rx(i2))^2+(y-ry(i2))^2);
                rij(i1,i2)=r0+alfa*Xij(i1,i2);
%                 plot(x,y,'ro','lineWidth',2)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if rij(i1,i2)>=Dij(i1,i2)+r0 %Afecta totalmente:
                    Ashadowi(i1,i2)=A0;
                    J(i2,end)=J(i2,end)+1;
                    I(i1,end)=I(i1,end)+1;
                    J(i2,J(i2,end))=i1;
                    I(i1,I(i1,end))=i2;
%                     plot(x,y,'g*','lineWidth',7)
                    %dibujo afectación de la estela:
                    angAspa=180-tetaAlterna; %para 10o funciona
                    y3=rij(i1,i2)*sind(angAspa)+y; y2=y-rij(i1,i2)*sind(angAspa);
                    x2=rij(i1,i2)*cosd(angAspa)+x; x3=x-rij(i1,i2)*cosd(angAspa);
%                     plot([x2,x3],[y2,y3],'g-','linewidth',2)
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
%                     plot(x,y,'m*','lineWidth',7)
                    %dibujo afectación de la estela:
                    angAspa=180-tetaAlterna; %para 10o funciona
                    y3=rij(i1,i2)*sind(angAspa)+y; y2=y-rij(i1,i2)*sind(angAspa);
                    x2=rij(i1,i2)*cosd(angAspa)+x; x3=x-rij(i1,i2)*cosd(angAspa);
%                     plot([x2,x3],[y2,y3],'c-','linewidth',2)
                elseif (rij(i1,i2)<Dij(i1,i2)-r0) % No afecta
                    Ashadowi(i1,i2)=0;
%                     plot(x,y,'w*','lineWidth',7)                    
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           elseif (teta==180 && rx(i2)>rx(i1) && rx(i2)>=0 && rx(i2)<=maxx) || (teta==0 && rx(i2)<rx(i1) && rx(i2)>=0 && rx(i2)<=maxx)
                y=ry(i1); %OKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKO
                x=rx(i2); %OKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKO
                Xij(i1,i2)=sqrt((x-rx(i1))^2+(y-ry(i1))^2);
                Dij(i1,i2)=sqrt((x-rx(i2))^2+(y-ry(i2))^2);
                rij(i1,i2)=r0+alfa*Xij(i1,i2);
%                 plot(x,y,'gp','lineWidth',2)   
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if rij(i1,i2)>=Dij(i1,i2)+r0 %Afecta totalmente:
                    Ashadowi(i1,i2)=A0;
                    J(i2,end)=J(i2,end)+1;
                    I(i1,end)=I(i1,end)+1;
                    J(i2,J(i2,end))=i1;
                    I(i1,I(i1,end))=i2;
%                     plot(x,y,'g*','lineWidth',7)
                    %dibujo afectación de la estela:
                    angAspa=180-tetaAlterna; %para 10o funciona
                    y3=rij(i1,i2)*sind(angAspa)+y; y2=y-rij(i1,i2)*sind(angAspa);
                    x2=rij(i1,i2)*cosd(angAspa)+x; x3=x-rij(i1,i2)*cosd(angAspa);
%                     plot([x2,x3],[y2,y3],'g-','linewidth',2)
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
%                     plot(x,y,'m*','lineWidth',7)
                    %dibujo afectación de la estela:
                    angAspa=180-tetaAlterna; %para 10o funciona
                    y3=rij(i1,i2)*sind(angAspa)+y; y2=y-rij(i1,i2)*sind(angAspa);
                    x2=rij(i1,i2)*cosd(angAspa)+x; x3=x-rij(i1,i2)*cosd(angAspa);
%                     plot([x2,x3],[y2,y3],'c-','linewidth',2)
                elseif (rij(i1,i2)<Dij(i1,i2)-r0) % No afecta
                    Ashadowi(i1,i2)=0;
%                     plot(x,y,'w*','lineWidth',7)                    
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           elseif ((teta==270 && ry(i2)>ry(i1) && ry(i2)>=0 && ry(i2)<=maxy) || (teta==90 && ry(i2)<ry(i1) && ry(i2)>=0 && ry(i2)<=maxy))
                y=ry(i2); %OKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKO
                x=rx(i1); %OKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKOKO
                Xij(i1,i2)=sqrt((x-rx(i1))^2+(y-ry(i1))^2);
                Dij(i1,i2)=sqrt((x-rx(i2))^2+(y-ry(i2))^2);
                rij(i1,i2)=r0+alfa*Xij(i1,i2);
%                 plot(x,y,'bh','lineWidth',2)  
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if rij(i1,i2)>=Dij(i1,i2)+r0 %Afecta totalmente:
                    Ashadowi(i1,i2)=A0;
                    J(i2,end)=J(i2,end)+1;
                    I(i1,end)=I(i1,end)+1;
                    J(i2,J(i2,end))=i1;
                    I(i1,I(i1,end))=i2;
%                     plot(x,y,'g*','lineWidth',7)
                    %dibujo afectación de la estela:
                    angAspa=180-tetaAlterna; %para 10o funciona
                    y3=rij(i1,i2)*sind(angAspa)+y; y2=y-rij(i1,i2)*sind(angAspa);
                    x2=rij(i1,i2)*cosd(angAspa)+x; x3=x-rij(i1,i2)*cosd(angAspa);
%                     plot([x2,x3],[y2,y3],'g-','linewidth',2)
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
%                     plot(x,y,'m*','lineWidth',7)
                    %dibujo afectación de la estela:
                    angAspa=180-tetaAlterna; %para 10o funciona
                    y3=rij(i1,i2)*sind(angAspa)+y; y2=y-rij(i1,i2)*sind(angAspa);
                    x2=rij(i1,i2)*cosd(angAspa)+x; x3=x-rij(i1,i2)*cosd(angAspa);
%                     plot([x2,x3],[y2,y3],'c-','linewidth',2)
                elseif (rij(i1,i2)<Dij(i1,i2)-r0) % No afecta
                    Ashadowi(i1,i2)=0;
%                     plot(x,y,'w*','lineWidth',7)                    
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