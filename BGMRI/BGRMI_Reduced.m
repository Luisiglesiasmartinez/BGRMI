function [B, Edges] = BGRMI_Reduced(mRNA, dt, No_Time_Points, No_Replicates, Transcription_Factors, Prior)

if nargin > 6
    error('myfuns:BADR:TooManyInputs', ...
        'requires at most 2 optional inputs');
end

% Fill in unset optional values.
switch nargin
    
    case 4
        Transcription_Factors = 0;
        Prior = 0;
    case 5
       Prior = 0;
end



%% Preprocessing:
No_of_Genes = size(mRNA,2);

%% Preprocess prior information about the nature of the genes is available

if  sum(sum(Prior)) == 0 && sum(sum(Transcription_Factors)) ==0
  
  
   Prior = 0.5*(ones(No_of_Genes)-eye(No_of_Genes));
  
elseif sum(sum(Transcription_Factors)) ~= 0 && sum(sum(Prior)) ==0
    
    No_TFs = length(Transcription_Factors);
   Prior = 0.5*(ones(No_of_Genes, No_TFs)-eye(No_of_Genes,No_TFs));
    
    
else 
    

    
end

%% Check that the dt is a vector of the correct size

if length(dt) ~= No_Time_Points -1
    
    error('dt must be a vector of length equal to the number of time points minus one')
    
end

%% Preprocessing:
No_of_Genes = size(mRNA,2);

Size = [No_Replicates*(length(dt)),No_of_Genes];

dmRNA_dt = zeros(Size);

mRNA_to = zeros(Size);

%% Resize the dt vector

dt = repmat(dt, 1 ,No_of_Genes);

%% Calculate the Derivative and Scale it and centre it to zero


    j = 1;
k = 1;
for i = 1:No_Replicates
    
   Idx = j:j+No_Time_Points-1;
   
   X = mRNA(Idx,:);  
   
   dmRNA_dt(k:k+No_Time_Points-2,:) = (X(2:end,:)-X(1:end-1,:))./dt;
   
   
   mRNA_to(k:k+No_Time_Points-2,:) = X(1:end-1,:);
   
    j = j+No_Time_Points;
    
    k = k+No_Time_Points-1;
    
end





 
%% Create a List of Genes which are Candidate Regulators

if Transcription_Factors == 0

K = 1:No_of_Genes;
else K = Transcription_Factors;

end


%% Substract the mean of the Gene expression, this ensures that the Intercept is orthogonal to the variables  

 Mean_mRNA_to = mean(mRNA_to);
 Mean_mRNA_to = repmat(Mean_mRNA_to, size(mRNA_to,1),1);
 mRNA_to = mRNA_to-Mean_mRNA_to;
 
 [~, B, Edges]=  Occams_Model_Proposal_with_F_and_Steel_Sort_Berno_Predictors(mRNA_to, dmRNA_dt, K, Prior);
 