function [Child] = Rank_Mutation(A, ...
                                Child, ...
                                IndicesInteractionProtein, NumInteractionProtein, ...
                                Pm)

    % Rank the proteins using a PageRank-like algorithm
    ProRankScores = RankProteins(A);

    %------------------  ProRank+ based Heuristic Mutation ----------------%

    N = length(A);
    for Node_i = 1:N
        if NumInteractionProtein(Node_i) > 0 && rand <= Pm
            % Get the current complex of Node_i
            CurrentComplex = Child.CmplxID(Node_i);

            % Find the best new complex for Node_i based on ProRank+ scores
            BestComplex = CurrentComplex;
            BestScore = ProRankScores(Node_i);
            MaxClusterID = max(Child.CmplxID);

            % Iterate over possible new clusters
            for NewCluster = 1:MaxClusterID
                if NewCluster ~= CurrentComplex
                    NewClusterNodes = find(Child.CmplxID == NewCluster);
                    if isempty(NewClusterNodes)
                        NewScoreSum = ProRankScores(Node_i);
                    else
                        NewScoreSum = sum(ProRankScores(NewClusterNodes)) + ProRankScores(Node_i);
                    end

                    % Prefer clusters with higher sum of ProRank+ scores
                    if NewScoreSum > BestScore
                        BestScore = NewScoreSum;
                        BestComplex = NewCluster;
                    end
                end
            end

            % Update the complex of Node_i if a better complex is found
            if BestComplex ~= CurrentComplex
                Child.CmplxID(Node_i) = BestComplex;

                % Update the Chromosome based on new complex connections
                ConnectedNode = Child.Chromosome(Node_i);
                for ConnectedNodeCounter = 1:NumInteractionProtein(Node_i)
                    if Child.CmplxID(IndicesInteractionProtein(Node_i, ConnectedNodeCounter)) == BestComplex
                        ConnectedNode = IndicesInteractionProtein(Node_i, ConnectedNodeCounter);
                    end
                end
                Child.Chromosome(Node_i) = ConnectedNode;
            end
        end
    end
end

% Protein Ranking: Using a PageRank-Like Algorithm
function [RankedProteins] = RankProteins(A)
    alpha = 0.85; % Damping factor
    N = size(A, 1);
    R = ones(N, 1) / N; % Initial rank
    M = A ./ sum(A, 1); % Column-stochastic matrix

    for iter = 1:100
        R = alpha * M * R + (1 - alpha) / N;
    end

    [~, RankedProteins] = sort(R,'descend');
end
