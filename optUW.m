% optUW
%
%  holds all support functions for UW optimization functions
% project.
%
% functions can be called from outside with 'optUW.<function name>'
% written GMB
classdef optUW
    methods(Static)

        function [params,err] = fit(funName, params, freeList, varargin)
            % [params,err] = fit(funName, params, freeList, var1, var2,...)
            %
            % Helpful interface to MATLAB's 'fminsearch' function.
            %
            % Inputs:
            %   funName        Function to be optimized. Must have form
            %                  [err] = <funName>(params, var1, var2, ...)
            %
            %   params         A structure of with field names that correspond to
            %                  parameter values for fitted function. Params are allowed
            %                  to be matrices.
            %       options    A struc]ture with options for MATLAB's fminsearch program
            %                  (see OPTIMSET)
            %
            %   freeList       Cell array containing list of parameter names (strings)
            %                  to be free in fitting. Free strings can contain certain
            %                  values / ranges within the 'params' matrices. For
            %                  example, the following are valid.
            %
            %                  {'x(1)','y(3:4)', 'z(1:2,4:5)'}
            %
            %   var<n>         Extra variables to be sent into fitted function
            %                  'funName'
            %
            % Outputs:
            %   params         A structure with best fitting parameters as fieldnames
            %
            %   err            Error value at minimum, numeric
            %
            % Notes:
            % - Dependencies: params2var.m, var2params.m, fitFunc.m

            % Written by Geoffrey M. Boynton, Summer of '00
            % Edited by Kelly Chang, February 10, 2017

            %% Input Control

            options = []; % options for fminsearch (see OPTIMSET)
            if isfield(params, 'options')
                options = params.options;
            end

            if isempty(freeList)
                freeList = fieldnames(params);
            end

            %% Fit Function and Calculate Final Error
            options = optimset('MaxFunEvals', 2000, 'TolFun', 1.e-2);
            % turn free parameters in to 'var'
            vars = optUW.params2var(params, freeList);

            % calling fminsearch
            vars = fminsearch('fitfunc', vars, options, funName, params, freeList,varargin);

            % assign final parameters into 'params'
            params = optUW.var2params(vars, params, freeList);

            % evaluate the function 'funName' for error at minimum
            err = fitfunc(vars, funName, params, freeList, varargin);
        end

        function [params] = var2params(var, params, freeList)
            % [params] = var2params(var, params, freeList)
            %
            % Turns varues 'varues' into a field within 'params' with a field name
            % given in order from 'freeList'. Support function for 'fit.m' and
            % 'fitcon.m'.
            %
            % Inputs:
            %   var         New varues to be stored in the 'params' structure under
            %               field names (in order) from 'freeList'
            %
            %   params      A structure of parameter varues with field names that
            %               correspond with the parameter names in 'freeList'
            %
            %   freeList    Cell array containing list of parameter names (strings)
            %               that match the field names in 'params'
            %
            % Output:
            %   params      Same 'params' structure with parameter varues as field
            %               names that correspond with the parameter names in
            %               'freeList' with the varues from 'var'
            %
            % Notes:
            % - Dependencies: str2vec.m

            % Written by G.M. Boynton - Summer of '00
            % Edited by Kelly Chang - June 21, 2016

            %% Input Control

            if ischar(freeList)
                freeList = {freeList};
            end

            %% Assign 'var' into 'params' Structure

            count = 1;
            freeList = regexprep(freeList, '[= ]', '');
            varStr = regexprep(freeList, '(\(.*\))', '');
            numList = regexp(freeList, '(\(.*\))', 'match');
            for i = 1:length(varStr)
                indx = optUW.str2vec(params.(varStr{i}), char(numList{i}));
                params.(varStr{i})(indx) = var(count:(count + length(indx) - 1));
                count = count + length(indx);
            end
        end

        function [indx] = str2vec(param, str)
            % [indx] = str2vec(param, str)
            %
            % Extracts linearized 'indx' indices from 'param' parameter corresponding
            % with 'str'. Support function to 'params2var.m', 'params2varcon.m', and
            % 'var2params.m'.
            %
            % Inputs:
            %   param       The parameter that has values that are to be fitted
            %
            %   str         A string of the values within 'params' that to be fitted
            %               (i.e., if params were a vector of [1 5] size and only the
            %               4th and 5th values are to be fitted, the corresponding
            %               string would be '(4:5)').
            %
            % Output:
            %   indx        Linearized index into 'param' of the values to be fitted.
            %               If 'str' was empty, then returns all indices

            % Written by Kelly Chang - November 20, 2017

            %% Linearize Index for Value(s) to be Fitted in 'param'

            if isempty(str)
                indx = 1:numel(param);
            else
                if ~isempty(regexp(str, '(\(|,):', 'once'))
                    dims = discretize(regexp(str, '(\(|,):'), [0 regexp(str, ',') length(str)]);
                    for i = 1:length(dims)
                        indx = regexp(str, '(\(|,):', 'once');
                        str = [str(1:indx) sprintf('1:size(param,%d)',dims(i)) str((indx+2):end)];
                    end
                end
                indx = eval(sprintf('CombVec%s;', str));
                nIndx = arrayfun(@(x) sprintf('indx(%d,:)',x), 1:size(indx,1), 'UniformOutput', false);
                indx = eval(sprintf('sub2ind(size(param),%s);', strjoin(nIndx, ',')));
            end
        end

        function [var] = params2var(params, freeList)
            % [var] = params2var(params, freeList)
            %
            % Extracts 'freeList' fields from 'params' structure and outputs into a row
            % vector 'var'. Support function for 'fit.m' and 'fitcon.m'.
            %
            % Inputs:
            %   params      A structure of parameter values with field names that
            %               correspond with the parameter names in 'freeList'
            %
            %   freeList    Cell array containing list of parameter names (strings)
            %               that match the field names in 'params'
            %
            % Output:
            %   var         Values extracted from the 'params' structure with field
            %               names (in order) from 'freeList'
            %
            % Notes:
            % - Dependencies: str2vec.m

            % Written by G.M Boynton - Summer of '00
            % Edited by Kelly Chang - February 10, 2017

            %% Input Control

            if ischar(freeList)
                freeList = {freeList};
            end

            freeList = regexprep(freeList, '[= ]', '');
            varStr = regexprep(freeList, '(\(.*\))', '');
            numList = regexp(freeList, '(\(.*\))', 'match');
            if ~all(ismember(varStr, fieldnames(params)))
                errFlds = setdiff(fieldnames(params), varStr);
                error('Unknown ''freeList'' parameters: %s', strjoin(errFlds, ', '));
            end

            %% Extract 'freeList' Values from 'params' Structure

            var = cell(1, length(varStr));
            for i = 1:length(varStr)
                indx = optUW.str2vec(params.(varStr{i}), char(numList{i}));
                var{i} = params.(varStr{i})(indx);
            end
            var = cell2mat(var);
        end
    end
end