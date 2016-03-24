function z = maxprod(A, B, w, it)
nodeNumA = size(A,1);
nodeNumB = size(B,1);
%initialize
messages = ones([nodeNumA, nodeNumA, nodeNumB]); % m(i, j, xj) all is 1

%update messages
for t = 1:it % it(1000) times
    new_mes = zeros([nodeNumA, nodeNumA, nodeNumB]); % place new messages
    for i = 1:nodeNumA
        for j = 1:nodeNumA
            for xj = 1:nodeNumB
                for xi = 1:nodeNumB
                    new_mes(i,j,xj) = new_mes(i,j,xj) + exp(w(i,xi))*B(xi,xj)*N(A,messages,i,j,xi,nodeNumA);
                end
            end
        end
        messages = new_mes;
    end
end
%compute belief
bi = randn(nodeNumA, nodeNumB);
for i = 1:nodeNumA %i
    for xi = 1:nodeNumB % xi
        bi(i,xi) = exp(w(i,xi)) * N(A,messages,i,nodeNumA+1,xi,nodeNumA);
    end
end
%normalization
beliefSum = sum(bi,2);
for i = 1:nodeNumA %i
    for xi = 1:nodeNumB % xi
        bi(i,xi) = bi(i,xi)./beliefSum(i);
    end
end
%compute belief
bij = randn(nodeNumA,nodeNumA,nodeNumB,nodeNumB);
for i = 1: nodeNumA
    for j = 1:nodeNumA
        for xi = 1: nodeNumB
            for xj = 1:nodeNumB
                bij(i,j,xi,xj) = A(i,j) * B(xi,xj) * exp(w(i,xi)) * exp(w(j,xj)) * N(A,messages,i,j,xi,nodeNumA) * N(A,messages,j,i,xj,nodeNumA);
            end
        end
    end
end
%normalize
beliefSum = sum(sum(bij,4),3);
for i = 1: nodeNumA
    for j = 1:nodeNumA
        for xi = 1: nodeNumB
            for xj = 1:nodeNumB
                bij(i,j,xi,xj) = bij(i,j,xi,xj)./beliefSum(i,j);
            end
        end
    end
end
%compute max
for i = 1:nodeNumA
    maxBi = bi(i,1);
    for xi = 2:nodeNumB
        if bi(i,xi)>maxBi
            maxBi = bi(i,xi);
        elseif bi(i,xi) == maxBi
            z=0;
            return
        else
        end
    end
end
[~,I] = max(bi,[],2);
z=I;
end

function p = N(A,messages,i,j,xi,nodeNumA)
p = 1;
for Neighbor = 1:nodeNumA
    if Neighbor == j %N(i) = j
    elseif A(i,Neighbor) == 1
        p = p * messages(Neighbor,i,xi); %product all neighbor
    else
    end
end
end
