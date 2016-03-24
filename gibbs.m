function m = gibbs(A, B, w, burnin, its)
%declare
nodeNumA = size(A, 1);
nodeNumB = size(B, 1);
sample_burnin = randn(1, nodeNumA);

%initialize
for i = 1:nodeNumA
    sample_burnin(1, i) = 1;
end

%burnin
for i = 1:burnin % i is burnin round
    for j = 1:nodeNumA % j is node of A
        prob_candicate = randn(1, nodeNumB);
        for k = 1:nodeNumB % k is candicate of j
            edge_prod = 1;
            %count edge product
            for n = 1:nodeNumA
                if A(j, n) == 1 %has edge between j and n in A
                    edge_prod = edge_prod * B(k, sample_burnin(1, n));
                end
            end
            prob_candicate(1,k) = exp(w(j,k)) * edge_prod;
        end
        %random select candicate
        sumAll = sum(prob_candicate);
        prob_candicate = prob_candicate./sumAll; %normalize
        random = rand;
        for c = 1:nodeNumB
            if random <= prob_candicate(1, c)
                sample_burnin(1, j) = c;
                break
            else
                random = random - prob_candicate(1, c);
            end
        end
    end
end

%sampling
samples = randn(its, nodeNumA);
samples(1, :) = sample_burnin;
for i = 2:its % i is iterations
    for j = 1:nodeNumA % j is node of A
        prob_candicate = randn(1, nodeNumB);
        for k = 1:nodeNumB % k is candicate of j
            edge_prod = 1;
            %count edge product
            for n = 1:nodeNumA
                if A(j, n) == 1 %has edge between j and n
                    edge_prod = edge_prod * B(k, sample_burnin(1, n));
                end
            end
            prob_candicate(1,k) = exp(w(j,k)) * edge_prod;
        end
        %random select candicate
        sumAll = sum(prob_candicate);
        prob_candicate = prob_candicate./sumAll; %normalize
        random = rand;
        for c = 1:nodeNumB
            if random <= prob_candicate(1, c)
                sample_burnin(1, j) = c;
                break
            else
                random = random - prob_candicate(1, c);
            end
        end
    end
    samples(i, :) = sample_burnin;
end

%count samples
prob = zeros(nodeNumA, nodeNumB);
for its = 1:its
    for a = 1:nodeNumA
        prob(a, samples(its, a)) = prob(a, samples(its, a)) + 1;
    end
end

%compute probility
for a = 1:nodeNumA
    sumA = sum(prob(a,:));
    prob(a,:) = prob(a,:)./sumA;
end

m = prob;
end