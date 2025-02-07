function [ si_barrier_certificate ] = create_si_connectivity_barrier_certificate_with_obstacles_stiff(varargin)

        
       parser = inputParser;
    parser.addParameter('BarrierGain', 1e2);
    parser.addParameter('MaxSeparation', 0.0015);
    parser.addParameter('SafetyDistance', 0);
    parser.addParameter('N', 60);
    parser.addParameter('tri_verices',cell(0,2));
    parse(parser, varargin{:})
    opts = optimoptions(@quadprog,'Display','off','ConstraintTolerance', 1e-15);

    gamma = parser.Results.BarrierGain;
    maxSeparation = parser.Results.MaxSeparation;
    safetyDistance = parser.Results.SafetyDistance;
    N = parser.Results.N;
    tri_verices = parser.Results.tri_verices;

    numObstacles = size(tri_verices,2);
 
    H = speye(2*(N+1));
    
    si_barrier_certificate = @barrier_certificate;

    function [ dx ,du, Bcon, Btis] = barrier_certificate(dxi, x, r, dr)

        
        assert(size(dxi,1)==2&&size(x,1)==2,'Expected inputs to be 2-by-N')

        N = size(dxi, 2);

%***********************************h_con*********************************%        
       stateDiff = (x(:,2:end)-x(:,1:end-1)).';
       normStateDiff = vecnorm(stateDiff,2,2);
       DeltaMinusNorm = 0.5*(maxSeparation.^2 - normStateDiff.^2)./maxSeparation.^2; 
       Bcon = zeros(N,1);
       Bcon(1:N-1) = DeltaMinusNorm;

       stateDiff0 = (x(:,1)-r(:)).';
       normStateDiff0 = vecnorm(stateDiff0,2,2);
       Bcon(end) =  0.5*(maxSeparation.^2-normStateDiff0.^2)./(maxSeparation.^2);

        dBdxcon = zeros(N,2*(N+1));

        for j = 1:(N-1)
            dBdxcon(j,2*j-1:2*j) = stateDiff(j,:)./maxSeparation.^2;
            dBdxcon(j,2*(j+1)-1:2*(j+1)) = -stateDiff(j,:)./maxSeparation.^2;
        end
            dBdxcon(end,1:2) = -stateDiff0./(maxSeparation.^2);

%*******************************h_con enhanced****************************%
       stateDiff_ext = (x(:,3:end)-x(:,1:end-2)).';
       normStateDiff_ext = vecnorm(stateDiff_ext,2,2);
       DeltaMinusNorm_ext = 0.5*(4*maxSeparation.^2 - normStateDiff_ext.^2)./(4*maxSeparation.^2); 
       Bcon_ext = zeros(N-1,1);
       Bcon_ext(1:N-2) = DeltaMinusNorm_ext;

       stateDiff_ext0 = (x(:,2)-r(:)).';
       normStateDiff_ext0 = vecnorm(stateDiff_ext0,2,2);
       Bcon_ext(end) =  0.5*(4*maxSeparation.^2-normStateDiff_ext0.^2)./(4*maxSeparation.^2);

        dBdxcon_ext = zeros(N-1,2*(N+1));

        for j = 1:(N-2)
            dBdxcon_ext(j,2*j-1:2*j) = stateDiff_ext(j,:)./(4*maxSeparation.^2);
            dBdxcon_ext(j,2*j+3:2*j+4) = -stateDiff_ext(j,:)./(4*maxSeparation.^2);
        end
            dBdxcon_ext(end,3:4) = -stateDiff_ext0./(4*maxSeparation.^2);

%***********************************h_stiff*********************************%
       DeltaMinusNorm_stiff = 0.5*(normStateDiff_ext.^2 - (2.0*maxSeparation).^2 ).^2; %0.5 * (x^2+y^2 - D^2)^2
       Bcon_stiff = zeros(N-1,1);
       Bcon_stiff(1:N-2) = DeltaMinusNorm_stiff;

       Bcon_stiff(end) =  0.5*(normStateDiff_ext0.^2 -(2.0*maxSeparation).^2 ).^2;

        dBdxcon_stiff = zeros(N-1,2*(N+1));

        for j = 1:(N-2)
            dBdxcon_stiff(j,2*j-1:2*j) = -2*(normStateDiff_ext(j,:).^2 - (2.0*maxSeparation).^2 ).* stateDiff_ext(j,:);  %(x^2+y^2 - D^2)*2*x
            dBdxcon_stiff(j,2*j+3:2*j+4) = 2*(normStateDiff_ext(j,:).^2 - (2.0*maxSeparation).^2 ).* stateDiff_ext(j,:);
        end
            dBdxcon_stiff(end,3:4) = 2*(normStateDiff_ext0.^2 - (2.0*maxSeparation).^2 ).* stateDiff_ext0;

%***********************************h_obs*********************************%

    Btis = zeros(numObstacles*(N+1),1);
    dBtisdx = zeros(numObstacles*(N+1),2*(N+1));
    X = [x,r];       
         for ii = 1:numObstacles 

                   [z_opt , dist] = closestPointOnTriangle_vectorized([x,r],tri_verices{ii} );
                   Btis((N+1)*(ii-1)+(1:(N+1))) = (dist.^2 - safetyDistance^2)./safetyDistance^2;

                for jj = 1:N+1
                    dBtisdx((N+1)*(ii-1)+jj,2*jj-1:2*jj) = 2*(X(:,jj) - z_opt(:,jj))'./safetyDistance^2;
                end
         end


%%
        %Eq 19: A matrix
        dBdx = [dBdxcon; dBdxcon_ext ; -dBdxcon_stiff; dBtisdx ];

        gammaB3=zeros(N,1);
        gammaB3(1:N-1) = gamma*Bcon(1:N-1);
        gammaB3(end) = 1*gamma*(Bcon(end))+ stateDiff0*dr(:)./maxSeparation.^2;

        gammaB_ext=zeros(N-1,1);
        gammaB_ext(1:N-2) = 1*gamma*Bcon_ext(1:end-1);
        gammaB_ext(end) = 1*gamma*(Bcon_ext(end))+ stateDiff_ext0*dr(:)./(4*maxSeparation.^2);

        gammaB_stiff=zeros(N-1,1);
        gammaB_stiff(1:N-2) = 1*gamma*Bcon_stiff(1:end-1);
        gammaB_stiff(end) = 1*gamma*(Bcon_stiff(end)) - 2*(normStateDiff_ext0.^2 - (2.0*maxSeparation).^2 ).* stateDiff_ext0 *dr(:);
        
        %Eq 19: B matrix
        gammaB = [gammaB3 ; gammaB_ext ;-gammaB_stiff+0.1 ; 1e1*Btis+1e-10 ];


%%
        %Solve QP program generated earlier
        vhat =[dxi(:) ; dr(:)];
        b = gammaB;
        f = -H*vhat;
        A = sparse(-dBdx);
    
%%
%slack for connectivity
% Number of constraints (same as the length of b)
num_constraints1 = length([gammaB3;gammaB_ext]);
num_constraints2 = length(gammaB_stiff);
num_constraints = num_constraints1 + num_constraints2;

% Extend H to include slack variables, penalty weight k
k1 = 1e5; %CBF
k2 = 1e0; %CLF
H_extended = blkdiag(H, k1 * speye(num_constraints1) , k2 * speye(num_constraints2));

% Extend f to include slack variables
f_extended = [f; zeros(num_constraints, 1)];

% Extend A to include slack variables
A_extended = [A, -[speye(num_constraints); zeros(length(Btis), num_constraints)]];

% Set bounds for slack variables
lb = [-inf(2*N+2, 1); zeros(num_constraints1, 1); -inf(num_constraints2, 1)];
ub = [inf(2*N+2, 1); inf(num_constraints1, 1); inf(num_constraints2, 1)];

%%
% Solve the quadratic programming problem with extended variables
opts = optimoptions('quadprog', 'Display', 'none'); % Adjust options as needed
[vnew_extended] = quadprog(H_extended, double(f_extended), A_extended, b, [], [], lb, ub, [], opts);
% Extract the solution for the original variables vnew
vnew = vnew_extended(1:2*(N+1));
slack = vnew_extended(2*(N+1)+1:end);

%%

        v = reshape(vnew, 2, N+1);
        dx = v(:,1:N);
        du = v(:,end);
      
        end
end