function [p,options] = riffle_shuffle(p,options)

% riffe_shuffle(p) returns the outcome of applying one riffle shuffle
% with in-shuffle to vector p.
% k = nombre de shuffles
% For more information on riffle shuffles, see http://mathworld.wolfram.com/RiffleShuffle.html

% P.-A. Absil, 2015-07-30
% Update E.Massart September 2015.

% The implementation is a bit intricate but it avoids for loops. It is inspired
% from http://www.mathworks.com/matlabcentral/fileexchange/16919-interleave

if mod(length(p),2)==0      %if the number of elements to shuffle is even
    % Let's say that p(1) is the bottom of the deck and p(end) the top.
    % Bottom half deck (in the right hand):
    p_1 = p(1:floor(length(p)/2));
    % Top half deck (in the left hand):
    p_2 = p(floor(length(p)/2)+1:end);
    % Prepare a cell in which to put the two half decks "side by side":
    p_tmp = cell(2,max(length(p_1),length(p_2)));
    % Interleave the two half decks, starting with top half (the left hand):
    p_tmp(1,1:length(p_2)) = num2cell(p_2);
    p_tmp(2,1:length(p_1)) = num2cell(p_1);
    p = [p_tmp{:}];
else
    if mod(options.count,2)==1
        p_1 = p(1:floor(length(p)/2));
        p_2 = p(floor(length(p)/2)+1:end);
        
        p_tmp = cell(2,max(length(p_1),length(p_2)));
        p_tmp(1,1:length(p_2)) = num2cell(p_2);
        p_tmp(2,1:length(p_1)) = num2cell(p_1);
        p = [p_tmp{:}];
    else
        p_1 = p(1:ceil(length(p)/2));
        p_2 = p(ceil(length(p)/2)+1:end);
        
        p_tmp = cell(2,max(length(p_1),length(p_2)));
        p_tmp(1,1:length(p_1)) = num2cell(p_1);
        p_tmp(2,1:length(p_2)) = num2cell(p_2);
        p = [p_tmp{:}];
    end
    options.count = options.count + 1;
end
end

