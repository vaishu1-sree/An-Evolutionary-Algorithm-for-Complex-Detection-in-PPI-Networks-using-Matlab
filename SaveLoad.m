function [ComplexProteinLabel, ...
           NumberOfProteinsInComplexes, NumberOfKnownProteinsInComplexes, ...
           OverlapComplexesFlag,ProteinLabel,ProteinLabelPairs,A,N, ...
	       NumInteractionProtein,MaxNumInteractionProtein,IndicesInteractionProtein,NetworkNumber] = SaveLoad

DataInput = input('\n To view the datasets for Complex press 1: ');

if (DataInput == 1)
    if (length(dir('DataSets\Complex\*.mat')) == 0)
        InputData = input('\n No datasets available! \n To saving the datasets Press 1, Otherwise, Exit \n')
            if(InputData == 1)
                SelectedComplexProteinType = input('\n 1: Complex-D1-Yeast-D1 \n 2: Complex-D2-Yeast-D2 \n 3: Complex-Human-Protein \n')
                NetworkNumber = SelectedComplexProteinType;
                [ComplexType, ProteinType,ComplexName,ProteinName] = LoadComplexProteinDataSet(SelectedComplexProteinType)
                [ProteinLabel, ProteinLabelPairs, A, N, NumInteractionProtein,...
                    MaxNumInteractionProtein, IndicesInteractionProtein] = ComputeProteinDataSet(ProteinType)
                [ComplexProteinLabel,NumberOfProteinsInComplexes,NumberOfKnownProteinsInComplexes,OverlapComplexesFlag] = ...
                ComputeComplexDataSet(ComplexType,ProteinLabel,ProteinLabelPairs,A,N,NumInteractionProtein,MaxNumInteractionProtein,IndicesInteractionProtein)
                SaveDataSet(ComplexProteinLabel,NumberOfProteinsInComplexes,NumberOfKnownProteinsInComplexes,...
                    OverlapComplexesFlag,ProteinLabel,ProteinLabelPairs,A,N,NumInteractionProtein,MaxNumInteractionProtein,...
                    IndicesInteractionProtein,ComplexName,ProteinName)
            else
                return
            end
    else
    ComplexData = dir('DataSets\Complex\*.mat');
    ProteinData = dir('DataSets\Protein\*.mat');
    ComplexName = [ComplexData.name];
    ProteinName = ['Dataset for proteins are: ', ProteinData.name];
    fprintf('\n %s \n ', ComplexData.name); 
    ComplexDataSet = input('\n Select which dataset for complexes you need to load? \n Complex DataSet:');
    NetworkNumber = ComplexDataSet;
    fprintf('\n %s \n', ProteinData.name); 
    ProteinDataSet = input('\n Select which dataset for proteins you need to load? \n Protein DataSet:');
    SelectedComplexDataSet = convertStringsToChars(ComplexData(ComplexDataSet,1).name)
    SelectedProteinDataSet = convertStringsToChars(ProteinData(ProteinDataSet,1).name)
    SelectedComplexDataSet= append('DataSets/Complex/',SelectedComplexDataSet)
    SelectedProteinDataSet= append('DataSets/Protein/',SelectedProteinDataSet)
    load(SelectedComplexDataSet,'ComplexProteinLabel','NumberOfProteinsInComplexes', 'NumberOfKnownProteinsInComplexes', 'OverlapComplexesFlag');
    load(SelectedProteinDataSet,'ProteinLabel','ProteinLabelPairs','A','N','NumInteractionProtein','MaxNumInteractionProtein','IndicesInteractionProtein');
    end
end