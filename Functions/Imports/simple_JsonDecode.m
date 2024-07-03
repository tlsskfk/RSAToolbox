
function jsonStruct = simple_JsonDecode(varargin)
% simple decoder of json files for older Matlab versions <2017
% please use 'jsondecode' in newer MATLAB versions
%
% Usage:
%   simple_JsonDecode(jsonFilePath)
%   simple_JsonDecode()                         % select json file via UI
%
%   input: jsonFilePath: full path of json file
%
% Assumed format of json-file:
% "FIELDNAME1": DOUBLE,
% "FIELDNAME2": "STRING",
% "FIELDNAME3": [ARRAYVALUE1,ARRAYVALUE2],
%
% Author: Lukas Pirpamer 
% lukas.pirpamer@medunigraz.at
%
% Acknowledgements:
% Many thanks to thrynae for the bug report, the varargin mistake is
% corrected in v.1.1
%
% Updates:
% v.1.1: Handling of varargin for GUI and command-line call
% v.1.2: Sub-level support implemented
%
% ------------------------------------------------------------------------------
%% BSD 3-Clause License
% Copyright (c) 2021, Lukas Pirpamer
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
%
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%
if isempty(varargin)
    [jFileName, jFolder] = uigetfile('*.json','Please select a json file');
    if jFileName==0
        disp('error while selecting file')
    else
        jsonFilePath=[jFolder,jFileName];
    end
elseif ~iscell(varargin)
    jsonFilePath=varargin;
elseif length(varargin)==1
    jsonFilePath=varargin{:};
else
    if length(varargin)>1
        disp('Attention: this function allows only one argument')
    end
    error(['ERROR: input could not be read.',char(10)...
        ' Please try: simple_JsonDecode (without any argument) ',char(10)...
        ' OR ',char(10)...
        ' simple_JsonDecode(''jsonfile_example.json'')',char(10)...
        char(10)...
        ' if the error still persists, please open an issue on the github page'])
end
jsonTxt = fileread(jsonFilePath);
jsonStruct = myJsonRead(jsonTxt);
%% Main function: parse JSON string
    function jsonStruct = myJsonRead(textRep)
        %% Reshape input string to be parsed
        textRep=regexprep(textRep,';','SEMICOLON_REPLACE'); % temporary replace of semicolons
        textRep=regexprep(textRep,[char(10),'\t"'],';');
        textRep=regexprep(textRep,'{"',';');
        textRep=regexprep(textRep,[char(10),'\s*"'],';');
        textRep=regexprep(textRep,'SEMICOLON_SUBLEVEL_REPLACE',';'); % temporary replace of semicolons       
        textRep=regexprep(textRep,'{;','{');
        textRep=regexprep(textRep,char(10),'');
        textRep=regexprep(textRep,'\t','');
        startIdcs=strfind(textRep,'{');
        endIdcs=strfind(textRep,'}');
        
        stairVec=zeros(length(textRep),1);
        for k = 1:length(startIdcs)
            stairVec(startIdcs(k):end)=stairVec(startIdcs(k):end)+1;
        end
        for k = 1:length(endIdcs)
            stairVec(endIdcs(k):end)=stairVec(endIdcs(k):end)-1;
        end
        
        % find start indices of outer level (OL)
        startIdcsOL=find(and(stairVec==1, [diff(stairVec);0]==1));
        endIdcsOS=find(and(stairVec==1, [0;diff(stairVec)]==-1));
        if ~isempty(startIdcsOL)
            textRepOS='';
            for k = 1:length(startIdcsOL)
                if k==1
                    textRepOS=[textRep(1:startIdcsOL(k)-1),regexprep(textRep(startIdcsOL(k):endIdcsOS(k)),';','SEMICOLON_SUBLEVEL_REPLACE')]; % temporary replace of semicolons
                else
                    textRepOS=[textRepOS,textRep(endIdcsOS(k-1)+1:startIdcsOL(k)-1),regexprep(textRep(startIdcsOL(k):endIdcsOS(k)),';','SEMICOLON_SUBLEVEL_REPLACE')]; % temporary replace of semicolons
                end
            end
            textRepOS=[textRepOS,textRep(endIdcsOS(length(startIdcsOL))+1:end)];
        else
            textRepOS=textRep;
        end
        if strcmp(textRepOS(1),'{')
            textRepOS=textRepOS(2:end);
        end
        if strcmp(textRepOS(end),'}')
            textRepOS=textRepOS(1:end-1);
        end
        jsonTextCell=textscan(textRepOS,'%s','delimiter',';');
        jsonTextCell=jsonTextCell{:};
        
        %% PARSE OUTER LAYER
        myJsonCell=cell(size(jsonTextCell,1),2);
        for k=1:size(jsonTextCell,1)
            currentline=jsonTextCell{k};
            currentline = regexprep(currentline,'SEMICOLON_REPLACE',';'); % re-substitute semicolons
           
            sepIdx=strfind(currentline,'": ');
             if isempty(sepIdx)
                 sepIdx=strfind(currentline,'":');
                 sepLength=2;
             else
                 sepLength=3;
             end
            if ~isempty(sepIdx)
                cellName=currentline(1:sepIdx-1);
                myJsonCell{k,1}=cellName;
                if strcmp(currentline(end),',')
                    valueRaw=currentline(sepIdx+sepLength:end-1);
                else
                    valueRaw=currentline(sepIdx+sepLength:end);
                end
                if strcmp(valueRaw(1),'{')
                    % multiple dimensions expected, recall main function
                     myJsonCell{k,2} = myJsonRead(valueRaw);
                else
                    if ~isempty(valueRaw)
                        myJsonCell{k,2} = readValue(valueRaw);
                    end
                end
            end
        end
        emptyCell=cellfun(@(x) isempty(x),myJsonCell(:,1));
        myJsonCell(emptyCell,:)=[];
        
        %% convert to struct
        for k=1:size(myJsonCell,1)
            if ~isvarname(myJsonCell{k,1})
                % Matlab struct does not allow special characters in the fieldname
                myJsonCell{k,1} = regexprep(myJsonCell{k,1},':','_');
                myJsonCell{k,1} = regexprep(myJsonCell{k,1},'\\','_');
                myJsonCell{k,1} = regexprep(myJsonCell{k,1},'//','_');
                myJsonCell{k,1} = regexprep(myJsonCell{k,1},'(','_');
                myJsonCell{k,1} = regexprep(myJsonCell{k,1},')','_');
                myJsonCell{k,1} = regexprep(myJsonCell{k,1},'+','_');
                myJsonCell{k,1} = regexprep(myJsonCell{k,1},'-','_');
                disp( ['invalid fieldname! fieldname changed to: ',myJsonCell{k,1}])
            end
            if k==1
                jsonStruct=struct(myJsonCell{k,1},myJsonCell{k,2});
            else
                jsonStruct.(myJsonCell{k,1})=myJsonCell{k,2};
            end
        end
    end
%% Read value of current field
    function returnVal = readValue(valueRaw)
        valueDouble=str2double(valueRaw);
        if ~isnan(valueDouble)
            returnVal=valueDouble;
        else
            if strcmp(valueRaw(1),'"')
                returnVal=valueRaw(2:strfind(valueRaw(2:end),'"'));
            elseif strcmp(valueRaw(1),'[')
                if strcmp(valueRaw(2),']')
                    returnVal=zeros(0);
                else
                    valueArray=regexprep(valueRaw,'\t','');
                    valueArray=valueArray(2:end-1); % take string between brackets
                    valueArray=regexprep(valueArray,'^ *',''); % remove empty characters at the beginning
                    valueArray=regexprep(valueArray,'} *, *{','},{'); % ensure no spaces are between the seperator
                    if strcmp(valueArray(1),'{')
                        seperatorIdcs= strfind(valueArray,'},{');
                        if ~isempty(seperatorIdcs)
                            %% merge structs to a struct array: 
                            % returnVal returns a stuct array
                            seperatorIdcs=[seperatorIdcs,length(valueArray)];
                            for k=1:length(seperatorIdcs)
                                if k==1
                                    currentStr=valueArray(1:seperatorIdcs(k));
                                    returnVal = myJsonRead(currentStr);
                                else
                                    currentStr=valueArray(seperatorIdcs(k-1)+2:seperatorIdcs(k));
                                    try
                                        currentReturnVal=myJsonRead(currentStr);
                                        % merge structs
                                        fcr = fieldnames(currentReturnVal);
                                        fr= fieldnames(returnVal);
                                        for i = 1:length(fcr)
                                            fieldIdx=find(cellfun(@(x) strcmp(x,fcr{i}),fr));
                                            if isfield(returnVal,fr{i})
                                               returnVal(k).(fr{fieldIdx}) = currentReturnVal.(fcr{i}); 
                                            else
                                                returnVal(k).(length(fieldnames(returnVal))+1) = currentReturnVal.(fcr{i});
                                            end
                                        end
                                    catch
                                        disp('DEBUG')
                                    end
                                end
                            end
                        else
                            returnVal = myJsonRead(valueArray);
                        end
                    else
                        if isempty(strfind(valueArray,'"')) % assume the value is a string when " are found
                            valueArray=regexprep(valueArray,'\ ',''); % a numeric array should not contain spaces 
                        end
                        valueArray=textscan(valueArray,'%s','delimiter',',');                        
                        valueArray=valueArray{:};
                        valueArrayDouble=cellfun(@(x) str2double(x),valueArray);
                        returnVal=valueArrayDouble;
                    end
                end
            elseif strcmp(valueRaw,'true')
                returnVal=true;
            elseif strcmp(valueRaw,'false')
                returnVal=false;
            elseif strcmp(valueRaw,'null')
                returnVal=[];
            end
        end
    end
end