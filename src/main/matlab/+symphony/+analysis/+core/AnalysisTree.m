classdef AnalysisTree < tree
    
    properties
        filename    % h5 file name
        filepath    % h5 file path
        name        % Descriptive name of analysis
    end
    
    methods
        
        function cellName = getCellName(obj, nodeInd)
            cellName = obj.getParameterValue('cellName', nodeInd);
        end
        
        function mode = getMode(obj, nodeInd)
            mode = obj.getParameterValue('ampModeParam', nodeInd);
        end
        
        function device = getDevice(obj, nodeInd)
            device = obj.getParameterValue('deviceName', nodeInd);
        end
        
        function className = getClassName(obj, nodeInd)
            className = obj.getParameterValue('class', nodeInd);
        end
        
        function value = getParameterValue(obj, parameter, nodeInd)
            nodeData = obj.get(nodeInd);
            value = [];
            while ~ isfield(nodeData, parameter);
                if nodeInd == 0
                    return;
                end
                
                nodeInd = obj.getparent(nodeInd);
                nodeData = obj.get(nodeInd);
            end
            value = nodeData.(parameter);
        end
        
        function nodeData = updateNodeDataInStructure(obj, nodeId, in, out)
            nodeData = obj.get(parent);
            ref = obj.get(nodeId);
            nodeData.(out) = [];
            inStructure = ref.(in);
            fnames = fieldnames(inStructure);
            
            for i = 1 : length(fnames)
                field = fnames{i};
                inValue = inStructure.(field);
                
                if strcmp(field, 'units') || stcmp(field, 'type')
                    nodeData.(out) = inValue;
                elseif length(inValue) < 2 && isfield(nodeData.(out), field)
                    % append to end of existing value vector
                    nodeData.(out).(field)(end + 1) = inValue;
                end
            end
        end
        
    end
end
