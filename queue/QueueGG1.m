function [MeanCustTime] = QueueGG1(lambda, x, runlength)


%% -- INPUTS:
% x is a theta value chosen, within limits, represent mean #arrivals
% lambda is parameter for service times
% runlength is the number of replications to simulate
% seed is the index of the substreams to use (integer >= 1)
% other is not used
%Note: RandStream.setGlobalStream(stream) can only be used for Matlab
%versions 2011 and later
%For earlier versions, use the method RandStream.setDefaultStream(stream)
%

%% -- OUTPUTS:
% RETURNS fn and FnVar, throughput are variance in customers/minute


    %% *********************PARAMETERS*********************

    warmup=500;      %people warmup before start count  1000
    people=36;      %measure avg sojourn time over next fifty people.   36
    total=warmup+people;         
    nRuns=runlength;%number of repetitions

    
    %% Generate random data
    arrival=exprnd(1/x,total,nRuns);                      % Generate Interarrival time
    service=exprnd(1/lambda,total,nRuns);                       % Generate Service time
    
    
    %% RUN SIMULATIONS
    %Tracks performance times for each run
    MeanCustTime=zeros(nRuns,1);
    for k=1:nRuns                           
                                            % Customers:: 
        Customer=zeros(total,4);            % Column 1-time of arrival to queue, Column 2-service time, Column
        Customer(:,1)=cumsum(arrival(:,k)); % 3-time service complete, Column 4-sojourn time

        Customer(:,2)=service(:,k);                     % Insert service times into matrix
       
        Customer(1,3)=Customer(1,1)+Customer(1,2);      % Put first customer through empty system
        Customer(1,4)=Customer(1,2);
        
        for i=2:total
            Customer(i,3)=max([Customer(i,1),Customer(i-1,3)])+Customer(i,2);
            Customer(i,4)=Customer(i,3)-Customer(i,1);  % remaining customers through system
        end
        MeanCustTime(k)=mean(Customer((warmup+1):total,4));     % Calculate mean sojourn time for last 50 customers
    end
%     fn=mean(MeanCustTime)+ lambda^2;                      % Mean and Variance of Alpha Function=Mean sojourn time + theta^2
%     FnVar=var(MeanCustTime)/nRuns;
    
end

