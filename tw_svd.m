function [ticktock,core,node,R,L] = tw_svd(X,Tol)
    tstart=tic;
    dims = size(X);
    n = length(dims);
    C = X;
    ep=Tol/sqrt(n);
    node=cell(1,n);
    core=cell(1,1);
    R = ones(n,1);
    L = ones(n,1);
    
    for i=1:n-1
        if i==1
            C = reshape(C,dims(1),numel(C)/dims(1));
            [U,S,V] = svd(C,'econ');
            S = diag(S);
            rc = my_chop2(S,sqrt(3)*ep*norm(S));
            temp = factor(rc);
            if length(temp) < 3
                R(1:length(temp)) = temp;     
            else
                R(1) = max(temp);
                L(1) = min(temp);
                R(2) = prod(temp)/(R(1)*L(1));
            end
            U = U(:,1:R(2)*L(1)*R(1));
            U = reshape(U,dims(1),R(2),L(1),R(1));
            node{1} = permute(U,[4 1 3 2]);
            S = S(1:R(2)*L(1)*R(1));
            V = V(:,1:R(2)*L(1)*R(1));
            V = V*diag(S);
            C = reshape(V,dims(2),prod(dims(3:end)),R(2),L(1),R(1));
            C = permute(C,[3 1 2 5 4]);
        else
            m=R(i)*dims(i); C=reshape(C,[m,numel(C)/m]);
            [U,S,V] = svd(C,'econ');
            S = diag(S);
            rc=my_chop2(S,sqrt(2)*ep*norm(S));
            temp = factor(rc);
            L(i)=min(temp); R(i+1)=rc/L(i);
            U=U(:,1:L(i)*R(i+1));
            node{i}=reshape(U,[R(i),dims(i),L(i),R(i+1)]);
            S = S(1:L(i)*R(i+1));
            V = V(:,1:L(i)*R(i+1));
            V = V*diag(S);
            C = reshape(V,dims(i+1),prod(dims((i+2):end)),R(1),prod(L(1:i)),R(i+1));
            C = permute(C,[5 1 2 3 4]);
        end  
    end
        m=R(n)*dims(n)*R(1); C=reshape(C,[m,numel(C)/m]);
        [U,S,V] = svd(C,'econ');
        S = diag(S);
        rc=my_chop2(S,ep*norm(S));
        L(n)=rc;
        U=U(:,1:L(n));
        node{n}=permute(reshape(U,[R(n),dims(n),R(1),L(n)]),[1 2 4 3]);
        S = S(1:L(n));
        V = V(:,1:L(n));
        V = V*diag(S);
        C = reshape(V,L(1:n)');
        core{1} = C;
        ticktock = toc(tstart);
        
end

