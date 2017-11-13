function [B, Edges] = BGRMI_Uneven(mRNA, Times, Replicates, Transcription_Factors, Prior)

% Some function with 2 required inputs, 3 optional.

% Check number of inputs.
if nargin > 6
    error('myfuns:BADR:TooManyInputs', ...
        'requires at most 2 optional inputs');
end

% Fill in unset optional values.
switch nargin
    
    case 3
        Transcription_Factors = 0;
        Prior = 0;
    case 4
       Prior = 0;
end

No_of_Genes = size(mRNA,2);

Number_of_Replicates = unique(Replicates);
No_Rows = 0;
for i = 1:length(Number_of_Replicates)
    
No_Rows = No_Rows+length(find(Replicates==i))-1;

end

Size = [No_Rows,No_of_Genes];

dmRNA_dt = zeros(Size);

mRNA_to = zeros(Size);



%% Preprocessing:




%% Preprocess prior information about the nature of the genes is available

if  sum(sum(Prior)) == 0 && sum(sum(Transcription_Factors)) ==0
  
  
   Prior = 0.5*(ones(No_of_Genes)-eye(No_of_Genes));
  
    elseif sum(sum(Transcription_Factors)) ~= 0 && sum(sum(Prior)) ==0
    
    No_TFs = length(Transcription_Factors);
   Prior = 0.5*(ones(No_of_Genes, No_TFs)-eye(No_of_Genes,No_TFs));
    
    
    
end
   
    j = 1;

for i = 1:length(Number_of_Replicates)
    
   Idx = find(Replicates==Number_of_Replicates(i));
   
   X = mRNA(Idx,:);
   
   Ts = Times(Idx);
   
   [Ts, Order] = sortrows(Ts);
    
   X = X(Order,:);
   
   
   dmRNA_dt(j:j+length(Ts)-2,:) = (X(2:end,:)-X(1:end-1,:))./repmat(Ts(2:end)- Ts(1:end-1),1,No_of_Genes);
   
   mRNA_to(j:j+length(Ts)-2,:) = X(1:end-1,:);
   
    j = j+length(Ts)-1;
    
end



if Transcription_Factors == 0

K =  1:No_of_Genes;

else 
    
    K = Transcription_Factors;
    

end


%% Substract the mean of the Gene expression, this ensures that the Intercept is orthogonal to the variables  

 Mean_mRNA_to = mean(mRNA_to);
 Mean_mRNA_to = repmat(Mean_mRNA_to, size(mRNA_to,1),1);
 mRNA_to = mRNA_to-Mean_mRNA_to;
 
 
 [~, B, Edges]=  Occams_Model_Proposal_with_F_and_Steel_Sort_Berno_Predictors(mRNA_to, dmRNA_dt, K, Prior);


