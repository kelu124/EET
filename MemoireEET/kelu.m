function [time,ANTS,WORLD,PH]=kelu(N,X,Y,T,D,P,A)

%N = nombre de fourmis
%X,Y = Taille du monde
%D, décroissance des phéromones
%P, Probabilité d'aller tout droit.
%A attractiveness des pheromones
Thrs = 200; % max des pheromones
Diff = 0.01;
[WORLD,PH,ANTS]=init(N,X,Y,T);
tic;
po = 1;
figure;
for t=1:T
[WORLD,PH,ANTS]=MovingAnts(N,X,Y,ANTS,P,t,PH,WORLD,Thrs,A);
[PH] = Decay(PH,D,X,Y,Diff);
%imagesc(PH(:,:,1));
    if (mod(t,floor(T/5)) == 0)
        subplot(5,1,po);
        imagesc(PH(:,:,1)+WORLD);
        po = po + 1;
    end
end;
figure;
subplot(2,1,1);
imagesc(PH(:,:,1)+WORLD);
subplot(2,1,2);
imagesc(PH(:,:,2)+WORLD);

fourmis = zeros(X,Y);

for m=1:N
    if(ANTS(m,5)== 1)
        fourmis(ANTS(m,2),ANTS(m,3)) = fourmis(ANTS(m,2),ANTS(m,3)) + 1;
    end
end

figure;
imagesc(fourmis);

time = toc;
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           1
%       2   X   3
%           4
%
%

function [WORLD,PH,ANTS]=MovingAnts(N,X,Y,ANTS,P,t,PH,WORLD,Thrs,A)

for m=1:N
    x = ANTS(m,2);
    y = ANTS(m,3);
    
    FormerPosition = ANTS(m,4);
    AntState = ANTS(m,5);
    
    if(AntState == 1)
        PH(x,y,1)= PH(x,y,1) + 1;
        if (PH(x,y,1)>Thrs)
             PH(x,y,1)=Thrs;   
        end
    end
    if(AntState == -1)
            PH(x,y,2)= PH(x,y,2) + 1;
        if (PH(x,y,2)>Thrs)
             PH(x,y,2)=Thrs;   
        end
    end
    
    if (WORLD(x,y) == -5)
        if(AntState ~= -1)
        ANTS(m,5) = -1;
            ANTS(m,4) = mod(FormerPosition -1 + 2 ,4) + 1;
        end
    end
    
    if (WORLD(x,y) == -2)
        if(AntState ~= 1)
            ANTS(m,5) = 1;
            ANTS(m,4) = mod(FormerPosition -1 + 2 ,4) + 1;
        end
    end
    
    N_PH = zeros (4,2);
    Directivity = zeros (4,1);
    
    if (AntState == 1)
        N_PH(1,1) = PH(x,Check(Y,y+1),2);
        N_PH(2,1) = PH(Check(X,x-1),y,2);
        N_PH(3,1) = PH(x,Check(Y,y-1),2);
        N_PH(4,1) = PH(Check(X,x+1),y,2);
%     elseif(AntState == 1)
%         N_PH(1,1) = PH(x,Check(Y,y+1),2);
%         N_PH(2,1) = PH(Check(X,x-1),y,2);
%         N_PH(3,1) = PH(x,Check(Y,y-1),2);
%         N_PH(4,1) = PH(Check(X,x+1),y,2);
    end
    
    choice = rand(1,2);
    PP = 1/2 - P/2;
    

    
    Directivity(     FormerPosition               ) = P;
    Directivity( mod(FormerPosition -1 + 2, 4) + 1) = 0;
    Directivity( mod(FormerPosition -1 + 1, 4) + 1) = PP;
    Directivity( mod(FormerPosition -1 + 3, 4) + 1) = PP;
    
    if (AntState == 1)
        PhInfluence = Influence(FormerPosition,N_PH,PP,Thrs);
        Directivity = Directivity + A*PhInfluence;
    end;
  
  %% Prise de décision
  r = Decide ( choice(1) , Directivity(1), Directivity(2) ,Directivity(3) ,Directivity(4) );
  
  
  
  %% Déplacement effectif
    if    (r==1)
        ANTS(m,3) = Check(Y,y + 1); 
    elseif(r==2)
        ANTS(m,2) = Check(X,x - 1); 
    elseif(r==3)
        ANTS(m,3) = Check(Y,y - 1); 
    elseif(r==4)
        ANTS(m,2) = Check(X,x + 1); 
            elseif(r==0)
        ANTS(m,3) = Check(Y,ANTS(m,3)); 
    end 

    ANTS(m,4) = r;
    
end



return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [WORLD,PH,ANTS]=init(N,X,Y,T)
% INIT initialise les variables
WORLD = zeros(X,Y);
x = X - 50;
y = Y - 50;
WORLD(floor(2*x/3),floor(y/3)-1)=-2;WORLD(floor(2*x/3),floor(y/3))=-2;WORLD(floor(2*x/3),floor(y/3)+1)=-2;WORLD(floor(2*x/3),floor(y/3)+2)=-2;
WORLD(floor(x/3),floor(y/3)-1)=-5;WORLD(floor(x/3),floor(y/3))=-5;WORLD(floor(x/3),floor(y/3)+1)=-5;WORLD(floor(x/3),floor(y/3)+2)=-5;
% N fourmis, avec id, x, y, ancienne position, state, compteur, NpasX,
% NpasY

ANTS=zeros(N,8);
% deux phéromones
PH = zeros(X,Y,2);
ANTS(:,2)=mod(floor(rand(N,1)*X),X)+1;
ANTS(:,3)=mod(floor(rand(N,1)*Y),Y)+1;
ANTS(:,4)=mod(floor(rand(N,1)*4),4)+1;
ANTS(:,5)=0;
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [i]=Check(L,l)
i=mod(l,L);
if i==0
    i=L;
end
return;
%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [i]=Decide(d,p1,p2,p3,p4)

i=2;

if ( d <= p1 )
    i = 1;
elseif ( (d > p1) & (d <= p1+p2))
    i = 2;
elseif ( (d > p1+p2) & (d <= p1+p2+p3))
    i = 3;    
elseif ( d > p1 + p2 + p3)
    i = 4;
end;

return;
%%%%%%%%%%%%%%%%%%%%

function [PH]=Decay(PH,D,X,Y,Diff)

for i=1:X
    for j=1:Y
        if (PH(i,j,1)>1/D  )
            PH(i,j,1)= PH(i,j,1)*(1 - 1/D); % décroissance linéaire
            %PH(i,j,1)= PH(i,j,1) - 1/D; % décroissance linéaire
        end
        if (PH(i,j,2)>1/D  )
            PH(i,j,2)= PH(i,j,2)*(1 - 1/D); % décroissance linéaire
            %PH(i,j,2)= PH(i,j,2) - 1/D; % décroissance linéaire
        end
        if (PH(i,j,1)<=1/D  )
            PH(i,j,1)=0;
        end
        if (PH(i,j,2)<=1/D  )
            PH(i,j,2)=0;
        end
    end
end
return;

%%%%%%%%%%%%%%%%%%%%
function [Mapped] = Influence(FormerPosition,N_ph,PP,Thrs);

Mapped = zeros(4,1);
N_ph(FormerPosition,1) = 0;
N_ph(FormerPosition,2) = Thrs;
MaxP = 0;
MinP = 0;

% if    (N_ph(1,1) = (N_ph(4,1) + N_ph(2,1) + N_ph(3,1))/3 )
%     MaxP = 1;
% elseif(N_ph(2,1) >= (N_ph(1,1) + N_ph(4,1) + N_ph(3,1))/3 )
%     MaxP = 2;
% elseif(N_ph(3,1) >= (N_ph(2,1) + N_ph(1,1) + N_ph(4,1))/3 )
%     MaxP = 3;
% elseif(N_ph(4,1) >= (N_ph(1,1) + N_ph(2,1) + N_ph(3,1))/3 )
%     MaxP = 4;
% end

mm = max(N_ph(:,1));

if    (N_ph(FormerPosition,1) == mm )
    MaxP = FormerPosition;
elseif(N_ph(mod(FormerPosition - 1 + 3,4) + 1,1) == mm  )
    MaxP =  mod(FormerPosition - 1 + 3,4) + 1;
elseif(N_ph(mod(FormerPosition - 1 + 1,4) + 1,1) == mm  )
    MaxP =  mod(FormerPosition - 1 + 1,4) + 1;
elseif(N_ph(mod(FormerPosition - 1 + 2,4) + 1,1) == mm  )
    MaxP = 0;
end



if (MaxP ~= 0)
    Mapped( MaxP) = + PP;
    Mapped( mod(MaxP - 1 + 1,4) + 1) = - PP/2;
    Mapped( mod(MaxP - 1 + 2,4) + 1) = - PP/2;
    Mapped( mod(MaxP - 1 + 3,4) + 1) = - PP/2;

    Mapped(FormerPosition) = 0;

end

return;