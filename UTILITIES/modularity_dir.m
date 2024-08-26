function [Ci,Q]=modularity_dir(A,gamma)
    % Modularity for Directed Networks
    % This function finds the optimal community structure and calculates modularity for a directed network
    % Inputs:
    %   A: directed weighted/binary connection matrix
    %   gamma: resolution parameter (optional, default = 1)
    % Outputs:
    %   Ci: optimal community structure
    %   Q: maximized modularity

% Set default gamma if not provided
if ~exist('gamma','var')
    gamma = 1;
end

% Initialize variables
N=length(A);                            %number of vertices
Ki=sum(A,1);                            %in-degree
Ko=sum(A,2);                            %out-degree
m=sum(Ki);   
aa=(Ko*Ki).'/m;

% Compute modularity matrix
b=A-gamma*(aa);
B=b+b.';                            	%directed modularity matrix

% Initialize community structure
Ci=ones(N,1);                           %community indices
cn=1;                                   %number of communities
U=[1 0];                                %array of unexamined communites

ind=1:N;
Bg=B;
Ng=N;

% Main community detection loop
while U(1)                              %examine community U(1)
    % Compute eigenvector corresponding to largest positive eigenvalue
    [V,D]=eig(Bg);
    [~,i1]=max(real(diag(D)));         %maximal positive (real part of) eigenvalue of Bg
    v1=V(:,i1);                         %corresponding eigenvector
    
    % Initial community division based on eigenvector sign
    S=ones(Ng,1);
    S(v1<0)=-1;
    q=S.'*Bg*S;                         %contribution to modularity
    
    % Community is divisible
    if q>1e-10                       	%contribution positive: U(1) is divisible
        qmax=q;                         %maximal contribution to modularity
        Bg(logical(eye(Ng)))=0;      	%Bg is modified, to enable fine-tuning
        indg=ones(Ng,1);                %array of unmoved indices
        Sit=S;
        % Fine-tuning loop
        while any(indg)                 %iterative fine-tuning
            Qit=qmax-4*Sit.*(Bg*Sit); 	%this line is equivalent to:
            [qmax,imax]=max(Qit.*indg); %for i=1:Ng
            Sit(imax)=-Sit(imax);       %	Sit(i)=-Sit(i);
            indg(imax)=nan;             %	Qit(i)=Sit.'*Bg*Sit;
            if qmax>q                   %	Sit(i)=-Sit(i);
                q=qmax;                 %end
                S=Sit;
            end
        end
         % Update community structure
        if abs(sum(S))==Ng              %unsuccessful splitting of U(1)
            U(1)=[];
        else
            cn=cn+1;
            Ci(ind(S==1))=U(1);         %split old U(1) into new U(1) and into cn
            Ci(ind(S==-1))=cn;
            U=[cn U];                   %#ok<AGROW>
        end
    else  % Community is indivisible                               %contribution nonpositive: U(1) is indivisible
        U(1)=[];
    end
     % Prepare for next iteration
    ind=find(Ci==U(1));                 %indices of unexamined community U(1)
    bg=B(ind,ind);
    Bg=bg-diag(sum(bg));                %modularity matrix for U(1)
    Ng=length(ind);                     %number of vertices in U(1)
end
% Compute final modularity
s=Ci(:,ones(1,N));                      %compute modularity
Q=~(s-s.').*B/(2*m);
Q=sum(Q(:));
