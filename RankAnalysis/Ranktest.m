%---------------------------
%C Melian July 2020
%Malena's restauration test

%Wilcoxon signed rank test (Proxy for Restauration ranking test)
%good example -- https://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test


%Rank the importance of each species for "functional" restauration with input data containing species level continuous and discrete values -- each of the 18 species has one value for variable for a total of 10 variables.One option is to use the Wilcoxon signed-rank test (a non-parametric statistical hypothesis test) -- the test is used to compare two related samples whether their population median ranks differ. The null hypothesis is that data in x and y are samples from continuous distributions with equal medians, against the alternative that they are not. The test assumes that the two samples are independent. If one species has the largest median value and significantly differ from all the others for a specific variable, then that species contributes the most for that specific variable to the functional restauration. We can then rank each species for each variable and the total species rank is the rank for all the species for all the variables. The challenge for this specific dataset is that we need to resample each variable with an underlying (normal) distribution with mean given by the data and std chosen for a specific value to compare two distributions from two species -- pval = wilcoxon_test(x,y) size of x and size of y must match, it is tested whether the difference x-y is significantly different to m=0. The test assumes that the data in x and y come from a continuous distribution symmetric about its median.
%--------------------------

pkg load dataframe
pkg load statistics

%https://stackoverflow.com/questions/28407344/reading-text-number-mixed-csv-files-as-tables-in-octave
%data = dataframe ("Rankdata.csv");
Rankdata;

N = 10;%#variables
S = 18;%#Species

%Rank test
%for b = 0.01:0.05:10;
%Transient matrix
P = zeros(S,S,N);%p-value matrix
M = zeros(S,S,N);%median matrix
Me = zeros(S,S,N);%mean matrix

%Output matrix
R = zeros(S,N);
Rmedian = zeros(S,N+1);
Rmean = zeros(S,N+1);


for i = 1:S;
    for k = i+1:S;
        for j = 1:N;
        b = 5;%std
        x = b.*randn(100,1) + D(i,j);%Input data given by D(i,j) is assumed to be the mean
        y = b.*randn(100,1) + D(i+1,j);
        
        %subplot(2,2,1)
        %hist(x)%Plot distributions
        
        %subplot(2,2,2)
        %hist(y)%Plot distributions
        
        [pval, z] = wilcoxon_test (x, y);
        
        P(i,k,j) = pval;
        P(k,i,j) = pval;
        %pause
        
        if pval <=0.05;
           Md = median(x) - median(y);
           M(i,k,j) = Md; 
           M(k,i,j) = Md;
           %pause
        
           Mn = mean(x) - mean(y);
           Me(i,k,j) = Mn;
           Me(k,i,j) = Mn;
        end
        end        
            
    end
    %Fill in Output matrix
    for n = 1:N;   
        Ptest = P(i,:,n);
        Mtest = M(i,:,n);
        Total = find(P(i,:,n) <= 0.05); 
        Rmedian(i,n) = sum(M(i,Total,n)); 
        Rmean(i,n) = sum(Me(i,Total,n)); 
        %pause
    end
end


%The importance of each species for each variable 
%M2 = zeros(S,S,N);
%M1 = sum(abs(M),2);

%for i = 1:S;
%    for j = 1:N;
%        for k = 1:S;
%            if M1(i,1,j) > 0;
%               M2(i,k,j) = M(i,k,j)/M1(i,1,j);
%            end  
%        end
%        M(i,:,j);
%        M1(i,1,j);
%        M2(i,:,j); 
%        M3(i,:,j) = sum(abs(M2(i,:,:)));
        %pause 
%    end
%end

%The importance of each species for all variables
%TODO



for sp = 1:S;
    Rmedian(sp,11) = sum(Rmedian(sp,1:10));
    Rmean(sp,11) = sum(Rmean(sp,1:10));
end
hold on
red=unifrnd(0,1,1); green=unifrnd(0,1,1); blue=unifrnd(0,1,1);
subplot(2,2,1)
plot(Rmedian(:,11),Rmean(:,11),'.','color',[red green blue],'Markersize',20)
%set(gca,'XTick',[],'YTick',[])
xlabel('Median','fontsize', 16)
ylabel('Mean','fontsize', 16)
%end%b std

%[R, P, LCI, HCI] = corrcoef (Rmedian(:,11),Rmean(:,11))


%Compare original rank with WT rank
for u = 1:S;
    D(u,11) = sum(D(u,1:10));
end

subplot(2,2,2)
plot(D(:,11),Rmedian(:,11),'.','Markersize',20)
%set(gca,'XTick',[],'YTick',[])
xlabel('Raw Rank','fontsize', 16)
ylabel('WT Rank','fontsize', 16)

[Rc, Pc, LCI, HCI] = corrcoef (D(:,11),Rmean(:,11))



save WTRankMatrix.csv Rmean



