classdef App < handle

    methods (Static)
        
        function n = name()
            n = 'Symphony Analysis';
        end
        
        function d = description()
            d = 'Data Analysis System for Symphony data format';
        end

        function v = version()
            v = '2'; % i.e. 2.1-r1
        end
        
        function o = owner()
            o = 'ala-laurila-lab';
        end
        
        function r = repo()
            r = 'symphony-analysis-core';
        end

        function u = documentationUrl()
            u = symphony.analysis.app.App.getResource('docs', 'Home.html');
        end

        function p = getResource(varargin)
            resourcesPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))), 'resources');
            p = fullfile(resourcesPath, varargin{:});
        end

    end

end
