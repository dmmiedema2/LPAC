%% L-PAC: MATLAB function for determining absolute copy number alterations (CNAs) from multi-region data
% L-PAC leverages information from high-purity samples to improve the
% inference of CNAs from low-purity samples from the same tumor.

%% PROCEDURE
% 1) Single-sample inference of absolute CNAs by the CNH.m algorithm
% (available at https://github.com/dmmiedema/CNH)
% 2) Filtering of samples for which single-sample inference in step (1)
% failed based on the output purity == 1 and ploidy < 2.5.
% 3) Identification of the best sample from the multi-region data of a
% tumor as the sample with the lowest intratumor heterogeneity as defined by CNH in step (1), and for which inference was succesfull according to the
% criteria in step (2).
% 4) Alignment of relative CNAs to the absolute CNAs in a minimization
% procedure analogous to that described in Van Dijk et al., Nature
% Communications 2021, for CNH. But with minimization of candidate absolute
% CNA values to the best sample identified in step (3), instead of to
% integer values as in step (1).

%% DEPENDENCY

% CNH.m         available at https://github.com/dmmiedema/CNH. Description in Van
%               Dijk et al., Nature Communications 2021.

%% INPUT
% bin_val       Nb x Ns size array containing relative CNA values of bins
%               of equal size. With Nb the number of bins, and Ns the
%               number of multi-region samples. The number of input samples should be 2 or more (Ns > 1).
%               NOTE: the input bin-values should be non-log values, and
%               should be normalized, i.e. have an average of 1 per sample.

%% OUTPUT
% LPAC_success  Binary variable with 1 if LPAC inference was succesfull and
%               0 if unsuccesfull
% ploidies      Ns x 1 size vector containing the tumor ploidies
% purities      Ns x 1 size vector containing the sample purities
% CNH           Ns x 1 size vector containing the copy number
%               heterogeneity
% best_sample   Ns x 1 size binary vector in which the best sample used for
%               multi-region inference is identified.



%% MAIN FUNCTION

function [LPAC_success, ploidies, purities, CNHout, best_sample] = LPAC(bin_val)
    % determine sample size and number of bins
    [Nb,Ns] = size(bin_val);
    % check if data consists of at least two samples
    if Ns < 2
        error('Input should consist of multiple samples (=columns)');
    end
    % check if data is normalized
    for i = 1:Ns
        if abs(mean(bin_val(:,i)) - 1) > 0.02
            tekst = sprintf( 'input bin values (column %d) are not-normalized, i.e. have an average value of 1',i);
            if mean(bin_val(:,i)) < 0.5 &&  abs( mean(2.^bin_val(:,i)) - 1) < 0.1
                tekst = sprintf('%s\n input data appears to consist of log transformed data. Please provide non-log transformed data',tekst);
            end
            error(tekst);
        end
    end
    
    % STEP 1: single-sample inference of absolute CNAs
    CNHout = zeros(Ns,1);
    ploidies = zeros(Ns,1);
    purities = zeros(Ns,1);
    for i = 1:Ns
        [CNHout(i), ploidies(i), purities(i)] = CNH(bin_val(:,i), ones(Nb,1), [], []);
    end

    % STEP 2: criterion for success of single-sample inference
    step1_success = ~(ploidies < 2.5 & purities == 1);
    if sum(step1_success) == 0
        LPAC_success = false;
        best_sample = false(Ns,1);
    else
        LPAC_success = true;
    end
    % Continue to STEP 3 if for at least one sample the inference of STEP 2
    % was successful
    
    % STEP 3: identification of best single-sample
    if LPAC_success
        id = find( CNHout == min(CNHout(step1_success)),1);
        best_sample = false(Ns,1);
        best_sample(id) = true;
        
    % STEP 4: infer CNAs by alignment to absolute CNAs of best samples of
    % relative CNAs of other samples.
        % define reference absolute CNA
        a1 = (purities(id) * ploidies(id) + (1-purities(id))*2) / purities(id);
        a2 = -2*(1-purities(id))/purities(id);
        cna_abs_ref = a1 * bin_val(:,id) + a2;
        
        % define purities and ploidies for grid search
        alphas = 0.1:0.01:1;
        taus = 1.5:0.05:5;
        % loop over non-best samples
        ids = find(~best_sample);
        for i = 1:numel(ids)
            Dmin = 1;
            bin_val_test = bin_val(:,ids(i));
            % loop over purities
            for i2 = i:numel(alphas)
                alpha = alphas(i2);
                % loop over ploidies
                for i3 = 1:numel(taus)
                    tau = taus(i3);
                    a1 = (alpha*tau+(1-alpha)*2)/alpha;
                    a2 = -2*(1-alpha)/alpha;
                    bin_val_test_abs = a1*bin_val_test+a2;
                    % determine distance between test and reference CNA
                    Dtest = mean(abs(bin_val_test_abs - cna_abs_ref));
                    % select sample with minimum distance to reference as
                    % best sample
                    if Dtest < Dmin
                        Dmin = Dtest;
                        purities(ids(i)) = alpha;
                        ploidies(ids(i)) = tau;
                    end                
                end
            end
            % Determine CNH from L-PAC inferred absolute CNAs
            a1 = (purities(ids(i))*ploidies(ids(i))+(1-purities(ids(i)))*2)/purities(ids(i));
            a2 = -2*(1-purities(ids(i)))/purities(ids(i));
            cna_abs_test = a1*bin_val_test + a2;
            CNHout(ids(i)) = mean( min( mod(cna_abs_test,1), 1-mod(cna_abs_test,1)) );
        end
    end
end
        
    
    
    
    
        
   
    


