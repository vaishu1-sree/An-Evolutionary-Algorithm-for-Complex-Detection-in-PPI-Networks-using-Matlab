function [Population] = CreatePopulation(N, NumInteractionProtein, IndicesInteractionProtein, ...
                                         PopulationSize)


for IndividualCounter = 1 : PopulationSize
    for ProteinCounter = 1 : N
        if(NumInteractionProtein(ProteinCounter) > 0)
            random = randi(NumInteractionProtein(ProteinCounter));
            Population(IndividualCounter).Chromosome(ProteinCounter) =  IndicesInteractionProtein(ProteinCounter, random);
        else
            Population(IndividualCounter).Chromosome(ProteinCounter) =  0;
        end;
    end;
end;