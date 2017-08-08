%Outputs are written sequentially, and file column width is the same
%regardless of number of outputs, so you need to define the number:

%Type vector defines the way the hex number are read.  For position data
%from the encoders, just use the default [0 0].

function [data2] = hexinterp(pathname,filename,number_of_outputs,type_vector)

if isempty(number_of_outputs) == 1
    number_of_outputs = 2; %x and y position
else
    number_of_outputs = number_of_outputs;
end

if isempty(type_vector) == 1
    type_vector = [0 0];  %set conversion type for each output type
else
    type_vector = type_vector;
end

%[filename,pathname] = uigetfile;

%Read in as characters:
all_text = fileread(fullfile(pathname,filename));
all_text(13*16+1:13*16+1:end) = [];  %remove all end of lines

%Initialize:
index_time = 0;
data2 = zeros(floor(length(all_text)/floor((13*number_of_outputs))),number_of_outputs);

%Advance by thirteen characters * output to start data points for next time
%step:
for index_character = 1:13*number_of_outputs:length(all_text)-13*number_of_outputs
    index_time = index_time+1;

    %For each output desired, convert based on output type:
    for index_output = 1:number_of_outputs

        switch type_vector(index_output)
            case 0
                output_start = index_character + (index_output-1)*13;

                data_low = hex2dec(all_text(output_start+6:output_start+11));
                data_high= hex2dec(all_text(output_start:output_start+5));

                if data_high > 8388608
                    data_high = data_high - 16777216;
                end

                data2(index_time,index_output) = 16777216*data_high + data_low ;
              
        end %switch
        
    end
end