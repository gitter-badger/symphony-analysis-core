classdef Repository < appbox.Settings & mdepin.Bean
    
    properties
        startupFile
        searchPath
        analysisFolder
        rawDataFolder
        preferenceFolder
    end
    
    methods
        
         function obj = Repository(config)
            obj = obj@mdepin.Bean(config);
        end 
        
        function f = get.startupFile(obj)
            f = obj.get('startupFile', '');
        end
        
        function set.startupFile(obj, f)
            validateattributes(f, {'char', 'function_handle'}, {'2d'});
            obj.put('startupFile', f);
        end
        
        function set.searchPath(obj, p)
            validateattributes(p, {'char', 'function_handle'}, {'2d'});
            obj.put('searchPath', p);
        end
        
        function p = get.searchPath(obj)
            p = obj.get('searchPath', symphony.analysis.app.App.getResource('examples'));
        end
        
        function f = get.analysisFolder(obj)
            f = obj.get('analysisFolder', fullfile(char(java.lang.System.getProperty('user.home')), 'data', 'analysis'));
        end
        
        function set.analysisFolder(obj, f)
            validateattributes(f, {'char', 'function_handle'}, {'2d'});
            obj.put('analysisFolder', f);
        end
        
        function f = get.rawDataFolder(obj)
            f = obj.get('analysisFolder', fullfile(char(java.lang.System.getProperty('user.home')), 'data', 'analysis'));
        end
        
        function set.rawDataFolder(obj, f)
            validateattributes(f, {'char', 'function_handle'}, {'2d'});
            obj.put('analysisFolder', f);
        end
        
        function f = get.preferenceFolder(obj)
            f = obj.get('preferenceFolder', fullfile(char(java.lang.System.getProperty('user.home')), 'data', 'PreferenceFiles'));
        end
        
        function set.preferenceFolder(obj, f)
            validateattributes(f, {'char', 'function_handle'}, {'2d'});
            obj.put('preferenceFolder', f);
        end
    end
    
end

