% The minimum coin problem 
 
D = [1 2 5 10 20 50 100 200]';   % coin denominations
X = [2 2 2 1 3 2 2 1]'; 			% number of corresponding coins
p = 170; % item cost

n = length(X); % useful
 
 
 
% problem #1
 
% linear problem: 
% maximize sum(G) subject to 0 <= G <= X and D' * G = p
% where G = [G(0) ... G(n-1)]'

C = ones(n,1);
lb = zeros(n, 1); % lower bound
ub = X; % upper bound
ctype = ["S"]; % equality constraint on  D' * G = p
vartype = repmat("I", n, 1); % vars are integers
sense = -1; % maximization problem

param.msglev = 0; % verbosity

[g, fmin, status, extra] = glpk (C, D', p, lb, ub, ctype, vartype, sense, param);

% print results

total = 0;
printf("Problem #1:\nYou should use:\n");
for i = 1:n
	coins = g(i,1);
	if (coins == 0)
		continue;
	end
	total += coins;
	printf("	%d coins(s) of value %d\n", coins, D(i));
end
printf("Total number of coins = %d\n", total);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions for 2 & 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [maximal_weight coins_used change_coins] =  maximize(cost, denominations, highest_index_to_use, coins, weights, pi, bp_coins, bp_change_coins)
% whats the maximal weight (after change) that can be used to pay a cost.
% use only denomination up to highest_index_to_use
% weight - a vector with weight for each coin type (e.g. physical weight)
% pi is a dynamic programming table index by the cost, holding the maximal weight (or -inf if not calculated yet)
% bp_x is a back pointer tables which hold the coins & change coins to be used for each maximal weight

	n = length(denominations);
	change_coins = zeros(1, n);
	coins_used = zeros(1, n);
	
	if (cost <= 0) % give it in change
		[weight change_coins] =  get_change(-cost, denominations, weights);
		maximal_weight = - weight; % change has negative weight
	elseif (pi(cost) > -inf)
		% already calculated
		maximal_weight = pi(cost);
		coins = bp_coins(cost, :);
		change_coins = bp_change_coins(cost, :);
	else
		
		maximal_weight = -inf;
		% find the maximum by itterating through all denominations
		for i = highest_index_to_use:-1:1
			if (coins(i) == 0) % no coins in that denomination
				continue;
			end
			denom  = denominations(i);
			coins2 = coins; coins2(i) -= 1;
			[maximal_weight2 coins_used2 change_coins2] =  maximize(cost - denom, denominations, i, coins2, weights, pi, bp_coins, bp_change_coins);
			maximal_weight2 += weights(i);
			if (maximal_weight2 > maximal_weight)
				maximal_weight = maximal_weight2;
				coins_used = coins_used2; coins_used(i) += 1;
				change_coins = change_coins2;
			end
		end
		
		pi(cost) = maximal_weight;
		bp_coins(cost, :) = coins_used;
		bp_change_coins(cost, :) = change_coins;
		
	end

end

function [weight coins_used] =  get_change(cost, denominations, weights)
% we assumsume that denominations are increasing
% that change is given greedily
% and that the cost can be represented as a greedy sum of the denom-s

	weight = 0;
	coins_used = zeros(1, length(denominations));
	denom_index = length(denominations);

	while (cost > 0)
		if (cost - denominations(denom_index) >= 0)
			cost -= denominations(denom_index);
			weight += weights(denom_index);
			coins_used(denom_index) += 1;
		else 
			if (--denom_index < 1)
				printf("Error in change calculation\n");
				exit(1);
			end
		end
	end
end



% problem # 2

W = ones(1, n); % weights
pi = ones(p, 1) * (-inf); % dynamic programming table
bp_coins = bp_change_coins = zeros(p,n); % dynamic programming back pointers


[maximal_weight coins_used change_coins] =  maximize(p, D, n, X, W, pi, bp_coins, bp_change_coins);

total = 0;
printf("\n\nProblem #2:\nYou should use:\n");
for i = 1:n
	coins = coins_used(i);
	change = change_coins(i);
	if (coins == 0 && change == 0)
		continue;
	end
	total += coins - change;
	printf("	value %d: use %d coin(s), get %d coins in change\n",  D(i), coins, change);
end
printf("Total coins used = %d\n\n", maximal_weight);


% problem # 3

W = rand(1, n); % weights
pi = ones(p, 1) * (-inf); % dynamic programming table
bp_coins = bp_change_coins = zeros(p,n); % dynamic programming back pointers


[maximal_weight coins_used change_coins] =  maximize(p, D, n, X, W, pi, bp_coins, bp_change_coins);

total = 0;
printf("\n\nProblem #3:\n");
weights = W
printf("\nYou should use:\n");
for i = 1:n
	coins = coins_used(i);
	change = change_coins(i);
	if (coins == 0 && change == 0)
		continue;
	end
	total += coins - change;
	printf("	value %d: use %d coin(s), get %d coins in change\n",  D(i), coins, change);
end
printf("Total weight used = %f\n\n", maximal_weight);


