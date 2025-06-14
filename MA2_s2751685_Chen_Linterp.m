% This function returns a look-up table (as a matrix) suitable for 
% performing (N−1)th order Lagrange interpolation. 
% 
% The idea here is to pre-compute the Lagrange polynomials for a given N, 
% at Q sampled intervals between α=−1/2 and α= 1/2. Thus, we will create an
% interpolation look-up table.The advantage of this approach is the ability 
% to easily generalise N (as opposed to coding each case manually). 
% This approach also has the advantage of efficiency and speed within the 
% audio loop (at the cost of extra computer memory used), at least for large N.
% 
% Author: Zuyu Chen
% Date: 29/11/2024

function T = MA2_s2751685_Chen_Linterp(N, Q, fmode)
    
m = (0:N-1)';            % Indices of N samples
a_m_vec = -(N-1)/2 + m;  % vector for the N sample locations
T = zeros(Q, N);         % Initialize the table     

    switch fmode
        case 1
            % main mode that evaluates the polynomials for Q samples
            % in the range of [-1/2, 1/2)  

            % create the sample location vector for evaluating the
            % polynomials
            a = (-Q/2 + (0:Q-1))/Q ; 
            % loop along the columns
            for m = 0:N-1
                % the m-th sample location 
                a_m = -(N-1)/2 + m; 
                % sample location vector eliminating the current sample location
                a_q = a_m_vec(a_m_vec ~= a_m);
                % calculates the Lagrange polynomial and assign it to a column
                % in the table
                T(:,m+1) = prod(a - a_q)./prod(a_m - a_q);
            end
        case 2
            % plot each individual Lagrange polynomials functions in the
            % range of [-(N-1)/2, (N-1)/2)

            % create the sample location vector for evaluating the
            % polynomials 
            a = linspace(-(N-1)/2, (N-1)/2, Q+1);
            a = a(1:end-1); % exclude the right boundary
            % loop along the columns
            for m = 0:N-1
                % the m-th sample location 
                a_m = -(N-1)/2 + m;
                % select out the locations excluding the current one
                a_q = a_m_vec(a_m_vec ~= a_m);
                % calculate the polynomial and assign it to a column of the
                % table
                T(:,m+1) = prod(a - a_q)./prod(a_m - a_q);
                % plot each column
                plot(a, T(:,m+1));  
                hold on
            end
            xlabel('alpha')
            ylabel('value')
            title(sprintf("Lagrange Polynomials (N = %d, Q = %d)", N, Q))
            legend(num2str((-(N-1)/2:(N-1)/2)'))
    end
end