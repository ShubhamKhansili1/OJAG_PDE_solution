function matrix = Linear_connection_mixed(N)

% This function file create a matrix which establishes a linear connection
% between the agents.

% Let i = 0, 1, 2,..., N be the agents. (Keep N even)

% Here, we consider total of three leaders. Two are at the boundary and one in middle of domain.
% That is, if we assume N to be even, then the third leader agent 
% will be i = (N/2).

% Define h 
h = 1/N; 

% The linear function for x in domain [a,b] is given by
% f(x) = f(a)(b-x/b-a) + f(b)(x-a/b-a).
% f(x) = f(a) L(x) + f(b)(1-L(x)), where L(x) = -(x-b/b-a).



% Lets define the lagrange linear function for x in [0,7*h ]
    l0= @(x) -(x-7*h)/(7*h);
    l1= @(x) ones(size(x))-l0(x);
% Lets define the linear function for x in [7*h,14*h ]
    l2= @(x) -(x-14*h)/(7*h);
    l3= @(x) ones(size(x))-l2(x);
% Lets define the linear function for x in [14*h,21*h ]
    l4= @(x) -(x-21*h)/(7*h);
    l5= @(x) ones(size(x))-l4(x);
% Lets define the linear function for x in [21*h,28*h ]
    l6= @(x) -(x-28*h)/(7*h);
    l7= @(x) ones(size(x))-l6(x);
% Lets define the linear function for x in [28*h,35*h ]
    l8= @(x) -(x-35*h)/(7*h);
    l9= @(x) ones(size(x))-l8(x);
% Lets define the linear function for x in [35*h,42*h ]
    l10= @(x) -(x-42*h)/(7*h);
    l11= @(x) ones(size(x))-l10(x);


% Creating array which contains these function values

l0val = l0((0: 7)*h) ; 
l1val = l1((0: 7)*h) ;  

l2val = l2((8: 14)*h) ;
l3val = l3((8: 14)*h) ;  

l4val = l4((15: 21)*h) ; 
l5val = l5((15: 21)*h) ; 

l6val = l6((22: 28)*h) ; 
l7val = l7((22: 28)*h) ;

l8val = l8((29: 35)*h) ; 
l9val = l9((29: 35)*h) ;  

l10val = l10((36: 42)*h) ;
l11val = l11((36: 42)*h) ;  




% creating the matrix
 matrix = zeros(43);
 matrix(:,1) = [l0val';zeros(35,1)]; % containing first leader values
 matrix(:, 8) = [l1val'; l2val'; zeros(28,1)];
 matrix(:, 15) = [zeros(8,1);l3val';l4val'; zeros(21,1)];
 matrix(:, 22) = [zeros(15,1);l5val';l6val'; zeros(14,1)];% containing mid leader values
 matrix(:, 29) = [zeros(22,1);l7val';l8val'; zeros(7,1)];
 matrix(:, 36) = [zeros(29,1);l9val';l10val'];
matrix(:, 43) = [zeros(36,1);l11val']; % containing last leader values







