function surrogate = generate_constrained_surrogate(ts,order)

%original_ts = the original time series
%order = Markov order to be verified
%surrogates = a surrogate time series exactly preserving the Markov
%properties in the original time series

%Example:
%Create a surrogate that exactly preserves second-order Markov properties of
%the original time series
% >> ts = [1 2 1 1 2 1 2 1 1 2];
% >> surrogate = generate_constrained_surrogate(ts,2);

%converts sequence of order 'order' to a first order sequence
[orderOneSequence, dict] = convert_to_order_one(ts, order);

L = length(orderOneSequence);

%estimation of the number of swaps
nSwaps = log(2*L)/(log(3*L)-log(3*L-4));

%swap surrogate algorithm
surrogate_order1 = swap_algorithm(orderOneSequence', nSwaps);

%convert back first order surrogate to original order
surrogate = convert_to_higher_order(surrogate_order1, dict);

end

function timeseries = swap_algorithm(timeseries, n_swaps)

%this version of the code works with arrays of integers
%0s and non-sucessive numbers are not allowed. Examples:
%[0,1,2,1,2] -> [1,2,3,2,3]
%[1,3,5,3,2,1] -> [1,3,4,3,2,1]

%--------
unique_symbs = unique(timeseries);
n_symbs = length(unique_symbs);

%estimate the probability of each symbol to be use in the sampling part
symbolFrequency = histcounts(timeseries, [unique_symbs n_symbs+1]);

counts = (symbolFrequency - 1).*(symbolFrequency);
prob_symbs = counts/sum(counts);

vec_idx = 1:length(timeseries);
for n=1:n_swaps
    
    %sample two distinct symbols according to their probability
    rand_symbs = datasample(unique_symbs,2,'Weights',prob_symbs);
    
    s_1 = rand_symbs(1);
    s_2 = rand_symbs(2);
    
    %find locations where corresponding symbols occurr in the time series
    indices_s1_logical = timeseries==s_1;
    indices_s2_logical = timeseries==s_2;
    indices_s1 = vec_idx(indices_s1_logical);
    indices_s2 = vec_idx(indices_s2_logical);
    
    %if there are at least two symbols each case... (we always need two
    %pairs to perform a swap
    if((numel(indices_s1)>1)&&(numel(indices_s2)>1))
        
        %sort two distinct instances of each symbol
        A = sort(datasample(indices_s1,2, 'Replace', false));
        B = sort(datasample(indices_s2,2, 'Replace', false));
    
        %an adjust to facilitate the verification of an allowable swap
        %(make pairs in crescent order of indexes)
        if(A(1)>B(1))
            tmp = A;
            A = B;
            B = tmp;
        end

        %now let's check if a swap if possible
        unique_idx = unique([A B]);
        n_different = numel(unique_idx);

        % making them in the order X1 Yi X2 Yii
        indexes = [A(1) B(1) A(2) B(2)];
        indexes_sorted = sort(indexes);
        sum_flag = sum(abs(indexes-indexes_sorted));

        %if there is no overlap
        if((n_different==4)&&(sum_flag==0))
            string1 = timeseries((A(1)+1):(B(1)-1));
            string2 = timeseries((A(2)+1):(B(2)-1));
            between = timeseries(B(1):A(2));
            timeseries = [timeseries(1:A(1)) string2 between string1 timeseries(B(2):end)]; 
        %if there is a overlap and we have the situation A....A....A
        elseif(n_different==3)
            A(1) = unique_idx(1);
            A(2) = unique_idx(2);
            B(2) = unique_idx(3);
            string1 = timeseries((A(1)+1):(A(2)-1));
            string2 = timeseries((A(2)+1):(B(2)-1));
            between = timeseries(A(2));
            timeseries = [timeseries(1:A(1)) string2 between string1 timeseries(B(2):end)];
        end
      
    end         
end
       
end

function [orderOneSequence, dict] = convert_to_order_one(higherOrderSequence, order)
% %A function to convert a sequence hypothesised to be of Markov order
% %'order' to an equivalent first order Markov sequence.
% %
% %Example of commands to call this function:
% originalString = 'adadbcdabcadbadbcddcdbacbdabcdadadabdcbcbdadadbcbdbcdabdcbdadcbdbcbdadbabdadbcdbdbcabcdadbcdaacdabcdbacdcbadbcad';
% order = 2;
% [orderOneSequence, dict] = convert_to_order_one(originalString, order); 
L = length(higherOrderSequence);
numObs = L - (order - 1);
observs = NaN(numObs, order);
for ii = 1:order
    observs(:, ii) = higherOrderSequence(ii:(ii + (numObs - 1)));
end

[dict, ~, orderOneSequence] = unique(observs, 'rows');

end

function higherOrderSequence = convert_to_higher_order(orderOneSequence, dict)

% %A function to reverse the process performed by the function
% 'convert_to_order_one'.
% %
% %Example of commands to call this function:
% originalString = 'adadbcdabcadbadbcddcdbacbdabcdadadabdcbcbdadadbcbdbcdabdcbdadcbdbcbdadbabdadbcdbdbcabcdadbcdaacdabcdbacdcbadbcad';
% order = 2;
% [orderOneSequence, dict] = convert_to_order_one(originalString, order);
% higherOrderSequence = convert_to_higher_order(orderOneSequence, dict);
% higherOrderSequence = char(higherOrderSequence);
% isequal(higherOrderSequence, originalString)%Are the original and the converted then reverted strings the same? 

higherOrderSequence = [dict(orderOneSequence, 1)', dict(orderOneSequence(end), 2:end)];

end