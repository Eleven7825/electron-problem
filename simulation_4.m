% simulation_4
% Consider the size of the ball and consider the collisions
% between balls 
% add the side walls
% consider the damping force without video recorder

% % parameters
r = 0.005; 
L = 0.5;
m = 0.002;
T = 100;
dt = 0.001;
q0 = 1e-8;
a = 0.1;
N = 300;
% Cd = 0.5;
DEN = 1.2;
if ~exist('n','var'); n = 600;     end
if ~exist('Dr','var'); Dr = 0.05;  end
if ~exist('Plot','var'); Plot = 1; end
if ~exist('Xrn','var'); Xrn = L/2;  end
if ~exist('Yrn','var'); Yrn = L;  end
if ~exist('q1','var'); q1 = 0;     end
q0 
Dr
% initial state of the balls
if ~exist('balls','var')
    balls = [2*Xrn*rand(1,n)-Xrn;Yrn*rand(1,n)];
    [f_x,f_y] = ball2ball(balls',q0*ones(n,1));
    balls = [balls; zeros(2,n);f_x./m;f_y./m;q0*ones(1,n)];
end
t = 0;
Q = [];
E =[];
xs = [];
ys = [];

% plotting part
if Plot ==1
    n = size(balls,2);
    fig = scatter(balls(1,:),balls(2,:), [],balls(7,:)'.* [balls(7,:)'./q0^2,ones(n,2)], 'filled');
    xlim([-Xrn Xrn])
    ylim([0 Yrn])
end

while t<T 
    % Merge the electron of the balls
    Dists = dist(balls(1:2,:));
    Dists = triu(Dists);
    [I,J] = find(Dists<2*r&Dists~=0);
    xs = [xs, balls(1,1)];
    ys = [ys, balls(2,1)];   
    
    % consider the collision of the balls
    if ~isempty(I)
        for i = 1:length(I)
            % Move the balls because of the round off error
            Dist = Dists(I(i),J(i));
            xlen = 2*r*(balls(1,I(i))-balls(1,J(i)))/Dist;
            ylen = 2*r*(balls(2,I(i))-balls(2,J(i)))/Dist;
            x_mid = (balls(1,I(i))+balls(1,J(i)))/2;
            y_mid = (balls(2,I(i))+balls(2,J(i)))/2;
            
            if abs(Dist)>=r
                balls(1,I(i)) = x_mid+xlen/2;
                balls(2,I(i)) = y_mid+ylen/2;
                balls(1,J(i)) = x_mid-xlen/2;
                balls(2,J(i)) = y_mid-ylen/2;
            else
                balls(1,I(i)) = x_mid+xlen/2+r;
                balls(2,I(i)) = y_mid+ylen/2+r;
                balls(1,J(i)) = x_mid-xlen/2-r;
                balls(2,J(i)) = y_mid-ylen/2-r;                
            end          
            
            % electron merge
            q = (balls(7,I(i))+balls(7,J(i)))/2;
            balls(7,I(i)) = q;
            balls(7,J(i)) = q;
            
            % retrieve the speeds
            ux = balls(3,I(i));
            uy = balls(4,I(i));
            vx = balls(3,J(i));
            vy = balls(4,J(i));
            
            % recalculate the speeds
            balls(3,I(i)) = (ux+vx)/2-a*(ux-vx)/2;
            balls(4,I(i)) = (uy+vy)/2-a*(uy-vy)/2;
            balls(3,J(i)) = (ux+vx)/2+a*(ux-vx)/2;
            balls(4,J(i)) = (uy+vy)/2+a*(uy-vy)/2;
        end
    end
        
    % Merge the ball with the board
    I = find(max(balls(2,:)<r,balls(2,:)>Yrn-r));
    if ~isempty(I)
        q = sum([balls(7,I),q1])./(N+length(balls(7,I)));
        balls(7,I) = q;
        q1 = N*q;
        balls(4,I) = -a*balls(4,I);
    end
    
    J = find(balls(1,:)<-Xrn+r|balls(1,:)>Xrn-r);
    if ~isempty(J)
        q = sum([balls(7,J),q1])./(N+length(balls(7,J)));
        balls(7,J) = q;
        q1 = N*q;
        balls(3,J) = -a*balls(3,J);
    end
    [fx,fy] = force(balls,Xrn,Yrn,q1,m);
    balls(5,:) = fx;
    balls(6,:) = fy;
     
    % the damping force propotional to the speed
    balls(3,:) = balls(3,:) -Dr*balls(3,:).^2.*sign(balls(3,:));
    balls(4,:) = balls(4,:) - Dr*balls(4,:).^2.*sign(balls(4,:));
    
    %update the velocity
    balls(3,:) = balls(5,:).*dt+balls(3,:);
    balls(4,:) = balls(6,:).*dt+balls(4,:);
    
    % If touch the walls, do the reflection
    balls(1,balls(1,:)<r-Xrn) = 2*r-2*Xrn-balls(1,balls(1,:)<r-Xrn);
    balls(1,balls(1,:)>Xrn-r) = 2*Xrn-2*r-balls(1,balls(1,:)>Xrn-r);
    balls(2,balls(2,:)>Yrn-r) = 2*Yrn-2*r-balls(2,balls(2,:)>Yrn-r);
    balls(2,balls(2,:)<r) = 2*r-balls(2,balls(2,:)<r);
    
    % If the ball is too far, delete it
    q1 = q1+ sum(balls(7,balls(1,:)<-Xrn|balls(1,:)>Xrn|balls(2,:)<0|balls(2,:)>Yrn));
    balls(:,balls(1,:)<-Xrn|balls(1,:)>Xrn|balls(2,:)<0|balls(2,:)>Yrn)=[];
    disp(t/T)
    
    % update the position
    balls(1,:) = balls(3,:).*dt+balls(1,:);
    balls(2,:) = balls(4,:).*dt+balls(2,:);
    
    %plotting
   if Plot == 1
        clf
        n = size(balls,2);
        scatter(balls(1,:),balls(2,:), [], balls(7,:)'.*[balls(7,:)'./q0^2,ones(n,2)], 'filled');
        axis([-Xrn, Xrn,0,Yrn])
        axis equal
        drawnow
   end
   
    t = t+dt;
    Q = [Q sum([balls(7,:),q1])];
    E = [E sum(balls(3,:).^2+balls(4,:).^2)];
end