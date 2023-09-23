%Statistical analysis: DE_12, 13, 14, 15, 16
rng('shuffle');
Nt=46;
Np=30;
maxx=2000; maxy=2000; 
[X,Y] = meshgrid(0:300:2000, 0:300:2000); % TURBINES POSITIONS.
    rx=[]; ry=[]; tam=15;
    for i1=1:tam
        rx=[rx;X(1,:),X(2,:),X(3,:),X(4,:),X(5,:),X(6,:),X(7,:)];
        ry=[ry;Y(1,:),Y(2,:),Y(3,:),Y(4,:),Y(5,:),Y(6,:),Y(7,:)]; 
    end
    rx=rx(:,1:46).*rand(tam,Nt); ry=ry(:,1:46).*rand(tam,Nt);
    for i1=tam+1:Np
        rx=[rx;maxx.*rand(1,Nt)];
        ry=[ry;maxy.*rand(1,Nt)];
    end
    X=[rx,ry];
estadistica=zeros(30,6);
for i1=1:3
    estadistica(i1,:)=[i1,DE_12(X),DE_13(X),DE_14(X),DE_15(X),DE_16(X)];    
end
% DE_12(X);
% DE_13(X);
% DE_14(X);
% DE_15(X);
% DE_16(X);
i1=0;