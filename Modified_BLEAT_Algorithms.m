%
% Modified BLEAT algorithm for generating hulls in the upper and lower half
% planes
%

% Receiving input:
%
% Driving function: real part
fprintf('Enter the real part of the driving function in t, beginning with "@(t)". \n');
real_function_input = input('Real part: ', 's');
f_real = str2func(real_function_input);

% Driving function: imaginary part
fprintf('Enter the imaginary part of the driving function in t, beginning with "@(t)". \n');
imag_function_input = input('Imaginary part: ', 's');
f_imag = str2func(imag_function_input);

% Upper bound of interval:
fprintf('Enter a value T > 0 such that the input function is continuous on the closed interval [0, T]. \n');
T = input('Your value: ');

% Number of subdivisions:
fprintf('Enter an initial number of subdivisions. \n');
N = input('Your number of subdivisions: ');

% Choice to refine:
fprintf('Do you want to refine the program at the tips? (N will increase)\n');
refine_choice = input('Y/N: ', 's');

% Choice to error check:
fprintf('Do you want to check for error points? (Some points will be removed from the hull)\n');
error_check = input('Y/N: ', 's');

% Defining d and t as arrays to store the time and value subdivisions:
d = zeros(1, N + 1);
t = zeros(1, N + 1);

% Filling in all values of t and d in the arrays
%
% These values correspond to the the time partitions and the corresponding
% values that the driving function takes on
%
% Fills in t with all subdivisions:
for r = 1:N+1
    t(r) = (r-1)*(T./N);
end

% If the user selects to refine at the tip, adds more points in this area:
if refine_choice == 'Y'
    % Declares a temporary variable to hold the initial number of subdivisions:
    J = N;
    
    i = floor(0.99*J);
    while i <= N+1
        while abs((f_real(t(i)) + 1j*f_imag(t(i)) - (f_real(t(i-1)) + 1j*f_imag(t(i-1))))) >= 1/(2*J)
            N = N + 1;
            for x = N+1:-1:i+1
                t(x) = t(x-1);
            end
            t(i) = (t(i) + t(i-1))/2;
        end
        i = i + 1;
    end
    % Prints the new number of sudvisions:
    fprintf('Final number of subdivisions: %d\n', N);
end

% Fills in d with values of the driving function at each subinterval:
for k = 1:N+1
    d(k) = f_real(t(k)) + 1j*f_imag(t(k));
end

% Creates lists for tips in upper and lower half planes
upper_list = zeros(1, N+1);
lower_list = zeros(1, N+1);

% Includes 0 as it is the starting point of the hulls in this algorithm
upper_list(N+1) = 0 + 0j;
lower_list(N+1) = 0 + 0j;

% Main loop: finds inverse comformal map for each subinterval and applies 
% them to all necessary points
for j = 1:N
    
    % First calculates R and s values
    R = d(N + 2 - j) - d(N + 1 - j);
    s = t(N + 2 - j) - t(N + 1 - j);
            
    % Finds formula for comformal map h
    h = @(z) (0.5)*(2*z + R + sqrt(R.^2 + 16*s))*((2*z + R - sqrt(R.^2 + 16*s))./(2*z + R + sqrt(R.^2 + 16*s))).^(0.5-((R)./(2*sqrt(R.^2 + 16*s))));
    
    %Finds c value
    c = c_value(R, s);
     
    %Finds a, b, alpha, and beta values
    a = a_value(c);
    b = b_value(c);
    alpha = alpha_value(c);
    beta = beta_value(c);
        
    % Checks if this is the first iteration, because the
    % procedure will be slightly different if so
    if j == 1
        % Calculates the first two tips first
        upper_tip = tip_value(a, b, alpha, beta, s, R);
        lower_tip = l_tip_value(a, b, alpha, beta, s, R);
       
        % Adds the tips to the list:
        upper_list(j) = upper_tip;
        lower_list(j) = lower_tip;
    else
        % Applies map for this iteration to all points in both lists
        %(except the origin, which is a special case)
        for k = 1:(j-1)
            u = h(upper_list(k));
            l = h(lower_list(k));
            upper_list(k) = u;
            lower_list(k) = l;
        end
        
        % Calculates new tips:
        upper_tip = tip_value(a, b, alpha, beta, s, R);
        lower_tip = l_tip_value(a, b, alpha, beta, s, R);
        
        % Adds new tips to the list:
        upper_list(j) = upper_tip;
        lower_list(j) = lower_tip;
        
    end
end

% Finds the value the driving function takes on at time 0
f_real_0 = f_real(0);
f_imag_0 = f_imag(0);

% Shifts the entire hull over by this value (translation property)
% (this is needed since the algorithm only approximates the change in value
% that the driving function takes on)
upper_list = upper_list + f_real_0 + 1j*(f_imag_0);
lower_list = lower_list + f_real_0 + 1j*(f_imag_0);

% Removes points that have been calculated incorrectly due to approximation
% error if prompted by the user:
if (error_check == 'Y')
    x = 20;
    while ((x > 0) && (abs(upper_list(x) - upper_list(x+1)) < 500/N))
        x = x - 1;
    end
    
    if (x ~= 1)
        for a = 1:x
            upper_list(1) = [];
        end
    end
end

%Creates plots, setting the window to have our desired dimensions:
plot(real(upper_list), imag(upper_list), real(lower_list), imag(lower_list));
xlim([-5 5]);
ylim([-5 5]);

%
% Functions:
%
function c = c_value(R, s)
    c = (-1)*(R./sqrt(s));
end

function a = a_value(c)
    a = c./(2*sqrt(c.^2 + 16)) - 0.5;
end

function b = b_value(c)
    b = (-1)*c./(2*sqrt(c.^2 + 16)) - 0.5;
end

function alpha = alpha_value(c)
    alpha = (c + sqrt(c.^2 + 16))./(2);
end

function beta = beta_value(c)
    beta = (c - sqrt(c.^2 + 16))./(2);
end

% Hard-coded formula to find the "upper" tips of hulls, modified for choice
% of log branch (in this case, (0, 2pi))
function tip = tip_value(a, b, alpha, beta, s, R)
    
    neg_beta = -beta;
    
    log_arg_beta = angle(neg_beta);
    
    log_beta_real = log(abs(neg_beta));
    
    beta_result_angle = (-1)*imag(a)*log_beta_real - real(a)*log_arg_beta;
    
    beta_result_real =(-1)*real(a)*log_beta_real + imag(a)*log_arg_beta;
    
    log_arg_alpha = angle(alpha);
    
    log_alpha_real = log(abs(alpha));
    
    alpha_result_angle = (-1)*imag(b)*log_alpha_real - real(b)*log_arg_alpha;

    alpha_result_real = (-1)*real(b)*log_alpha_real + imag(b)*log_arg_alpha;
    
    final_real = alpha_result_real + beta_result_real;
    
    final_angle = alpha_result_angle + beta_result_angle;
           
    pre_translate_tip = sqrt(s)*exp(final_real + 1j*(final_angle - (a*pi)));
    
    tip = pre_translate_tip + R;
    
end

% Hard-coded formula to find the "lower" tips of hulls, modified for choice
% of log branch (in this case, (-2pi, 0))
function lower_tip = l_tip_value(a, b, alpha, beta, s, R)
    
    neg_beta = -beta;
    
    log_arg_beta = angle(neg_beta);
    
    log_beta_real = log(abs(neg_beta));
    
    beta_result_angle = (-1)*imag(a)*log_beta_real - real(a)*log_arg_beta;
    
    beta_result_real =(-1)*real(a)*log_beta_real + imag(a)*log_arg_beta;
    
    log_arg_alpha = angle(alpha);
    
    log_alpha_real = log(abs(alpha));
    
    alpha_result_angle = (-1)*imag(b)*log_alpha_real - real(b)*log_arg_alpha;

    alpha_result_real = (-1)*real(b)*log_alpha_real + imag(b)*log_arg_alpha;
    
    final_real = alpha_result_real + beta_result_real;
    
    final_angle = alpha_result_angle + beta_result_angle;
           
    pre_translate_tip = sqrt(s)*exp(final_real + 1j*(final_angle + (a*pi)));
    
    lower_tip = pre_translate_tip + R;
    
end