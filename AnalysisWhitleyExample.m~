% This script computes the all marginal configurations of order up to 10
% for the additive decomposable function of order k=3, n=10 included in
% Whitley et al:2015

% Codomain vectors in index form
TableIndices = [1,2,5,8; 4,2,5,8; 2,3,6,8; 7,5,4,8; 2,5,7,8; 1,5,3,8; 2,3,6,8; 1,2,5,8; 1,3,4,8; 1,6,7,8];

% Indices of the 10 subfunctions
thefunctions = [1,2,3;
                2,3,4;
                3,4,5;
                4,5,6;
                5,6,7;
                6,7,8;
                7,8,9;
                8,9,10;
                9,10,1;
                10,1,2];
            

% We compute the Table containing the binary values
% for the 10 subfunctions. They will be
% contained in the table thevals 

for i=1:10,
   thevals(:,i) = zeros(8,1); 
   thevals(TableIndices(i,:),i) = 1;
end
    
% They look pretty much the same as the Codomain vectors included in the paper  Whittley et al 2015
% thevals'
% ans =
%     1     1     0     0     1     0     0     1
%     0     1     0     1     1     0     0     1
%     0     1     1     0     0     1     0     1
%     0     0     0     1     1     0     1     1
%     0     1     0     0     1     0     1     1
%     1     0     1     0     1     0     0     1
%     0     1     1     0     0     1     0     1
%     1     1     0     0     1     0     0     1
%     1     0     1     1     0     0     0     1
%     1     0     0     0     0     1     1     1

% We generate the complete space of 2^10 solutions
% and evaluate them using the additive decomposable function
    
nboolean = 10;
spacesize = 2^nboolean;

% We create all possible binary vectors of size nboolean and save them in  Pop
Card = 2*ones(1,nboolean);   % Cardinality (needed for generating the complete space of solutions)
AccCard = cumprod(Card)/2;
AccCard = AccCard([nboolean:-1:1]);
Pop = zeros(spacesize,nboolean);   % All space of solutions (binary vectors) are saved in Pop

% We fill the population with all possible binary vectors
for j=1:spacesize,
 num = IndexconvertCard(j-1,nboolean,AccCard); % Convert ordinal number to binary
 Pop(j,:) = num;                               % We save in the population
end 

% We evaluate all solutions using the additive function composed by the 10 subfunctions
% and save the results  in sumval
sumval = zeros(spacesize,1);
for j=1:1024, 
  for k=1:10,
    vars = thefunctions(k,:);               % First we get the 3 variables involved
    vals = Pop(j,vars);                     % We get their configuration in solution j
    theind = 4*vals(1)+2*vals(2)+vals(3)+1; % We compute the index for the subfunction table     
    sumval(j) =  sumval(j) + thevals(theind,k); % The local contribution is added to value of solution j    
    %[j,k,vars,vals,theind,sumval(j)]
  end 
end  


% We will compute the statistics for building blocks of size 
% between 1 and 6. The variables defining the blocks
% are saved in Blocks
% The Frequencies for all blocks are saved in Frequencies

largevals = [1:10,1:10];

for BB=1:10,
   for i=1:10,
     Blocks{BB}(i,:) = largevals(i:i+BB-1);     
   end
   Frequencies{BB} = zeros(2^BB,10);
end    


% Finally, we compute the frequencies  for blocks of up to 10 variables
% (relevant here those with 5 variables)

for BB=1:10,
 Card = 2*ones(1,BB);   % Cardinality (needed for generating the complete space of solutions)
 AccCard = cumprod(Card)/2;
 AccCard = AccCard([BB:-1:1]);   
 for k=1:10,
  for j=1:1024, 
    vars = Blocks{BB}(k,:);                                              % For block BB, the kth factor of variables
    vals = Pop(j,vars);                                                  % The configuration of these variables in the solution j
    theind = NumconvertCard(vals,BB,AccCard) + 1;                        % What is the index of the configuration  
    Frequencies{BB}(theind,k) = Frequencies{BB}(theind,k) +  sumval(j);  % The contribution to the hyperplane
    %[j,k,vars,vals,theind,sumval(j),Frequencies{BB}(theind,k)]
  end 
 end 
end  

% We visualize for each of the factors the highest marginal values
% Hyperplanes with the highest contributions

BB = 5;
for k=1:10,
  maxBB = max(Frequencies{BB}(:,k));
  [k,maxBB,Blocks{BB}(k,:)]  
  [[find(Frequencies{BB}(:,k) == maxBB)-1]']  
end


% We save the Population and the function values 
csvwrite('solutions_and_values.csv',[Pop,sumval]);

% We print all the tables for building blocks of size 5
mytable = [];

for i=1:16,
  sprintf('& %d ', Frequencies{4}(i,[10,1:9]))    
end  

for i=1:32,
  sprintf('& %d ', Frequencies{5}(i,[10,1:9]))    
end  
  

% Now we compute a factorization of the original problem
% using the triangulations of the original interaction graph

% The following matrix encodes the problem structure

matrix10 = [1 1 1 0 0 0 0 0 1 1;
            1 1 1 1 0 0 0 0 0 1;
            1 1 1 1 1 0 0 0 0 0; 
            0 1 1 1 1 1 0 0 0 0;
            0 0 1 1 1 1 1 0 0 0;
            0 0 0 1 1 1 1 1 0 0; 
            0 0 0 0 1 1 1 1 1 0; 
            0 0 0 0 0 1 1 1 1 1;
            1 0 0 0 0 0 1 1 1 1; 
            1 1 0 0 0 0 0 1 1 1];
        
    

% We triangulate the graph. The relevant factors are saved in "cliques"
order = [1:10];
[G, cliques, fill_ins] = triangulate(matrix10, order)
 
% Now we compute the frequencies  for the five cliques of (5 variables each)
  
 FactorizationFrequencies = zeros(32,6); 
 Card = 2*ones(1,5);   % Cardinality (needed for generating the complete space of solutions)
 AccCard = cumprod(Card)/2;
 AccCard = AccCard([5:-1:1]);   
 for k=1:6,
  for j=1:1024, 
    vars = cliques{k};                                              % For block BB, the kth factor of variables
    vals = Pop(j,vars);                                                  % The configuration of these variables in the solution j
    theind = NumconvertCard(vals,5,AccCard) + 1;                        % What is the index of the configuration  
    FactorizationFrequencies(theind,k) = FactorizationFrequencies(theind,k) +  sumval(j);  % The contribution to the hyperplane
   
  end   
 end  

 % The tables are saved
fid = fopen('Factorization_Frequencies.txt','w');
%fprintf(fid,'psize %d  gen %d \n',psizes(npsize),ngen(npsize)); 
   for i=1:32,  
      A = [sprintf('  %s  ',dec2bin(i-1,5)),sprintf('&  %d  ', FactorizationFrequencies(i,1:6)), '\\ \hline'];     
      fprintf(fid,'%s \n',A);
    end   
  fprintf(fid,'\n \n');
fclose(fid);

  


  
