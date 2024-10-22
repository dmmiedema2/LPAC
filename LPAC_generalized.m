%% L-PAC GENERALIZED: MATLAB function for determining absolute copy number alterations (CNAs) from multi-region data
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

% CNH.m                 available at https://github.com/dmmiedema/CNH. Description in Van
%                       Dijk et al., Nature Communications 2021.

% cnv_expand_segment.m  available at https://github.com/dmmiedema2/LPAC. 

%% INPUT
% seg_val                   Nb x Ns size array containing relative segmented CNA values. With Nb the number of segments, and Ns the
%                           number of multi-region samples. The number of input samples should be 2 or more (Ns > 1).
%                           NOTE: the input bin-values should be non-log values, and
%                           should be normalized, i.e. have an average of 1 per sample.
% seg_len                   Nb x 1 size array containing the genomic length of each segment
% candidate_best_samples    Nb x 1 size binary array containing potential
%                           best samples (=1). Leave empty ("[]") in case
%                           of no pre-selection of best samples.
% purity_ref                single scalar value of purity of pre-selected
%                           sample.  Leave empty ("[]") in case
%                           of no pre-selection of best samples.
% ploidy_ref                single scalar value of ploidyy of pre-selected
%                           sample.  Leave empty ("[]") in case
%                           of no pre-selection of best samples.

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

function [LPAC_success, ploidies, purities, CNHout, best_sample] = LPAC_generalized(seg_val, seg_len, candidate_best_samples, purity_ref, ploidy_ref)
    % parameters:
    bin_size = 10^5; %size of underlying bins in segmented CNA data. Used only in case segmented data is given is input.

    % determine sample size and number of bins
    [Nseg,Ns] = size(seg_val);    
    CNHout = zeros(Ns,1);
    ploidies = zeros(Ns,1);
    purities = zeros(Ns,1);
    % check if data consists of at least two samples
    if Ns < 2
        error('Input should consist of multiple samples (=columns)');
    end
    % check if input consist of segments (seg_len > 1) or bins (seg_len ==
    % 1): expand segments to bins
    if sum( seg_len > 1) > 0
        Nbin = ceil(sum(seg_len) / bin_size);
        bin_val = zeros(Nbin,Nbin);
        for i = 1:Ns
            cc = 1;
            for i2 = 1:Nseg
                db = round(seg_len(i2) / bin_size);
                bin_val(cc:cc+db-1,i) = seg_val(i2,i);
                cc = cc+db;
            end
            bin_val = bin_val(1:Nbin,:);
        end
    else
        bin_val = seg_val;
    end        
    % check if there is a candidate best sample
    if isempty(candidate_best_samples)
        candidate_best_samples = ones(Ns,1);
    else        
        if ~(size(candidate_best_samples,1) == Ns && size(candidate_best_samples,2)==1)
            error('candidate_best_samples input should be a vector of size: Nsamples x 1');
        elseif sum(candidate_best_samples == 1) < 1
            error('at least one candidate samples should be identified (element =1) or field should be left empty: []');
        end
    end
    % check if reference input ploidy and purity are specified: if so, skip
    % steps 1-3 of LPAC.
    if isempty(purity_ref)
        % do steps 1,2 and 3 of LPAC
        do_steps123 = 1;
        if ~isempty(ploidy_ref)
            warning('reference ploidy is specified, but reference purity is not. LPAC proceeds without reference sample');
        end
    elseif  numel(purity_ref)~=1 || numel(ploidy_ref)~=1 || ~(isnumeric(purity_ref) || isnumeric(ploidy_ref))
        error('input reference ploidy and purity should be single scalar value'); 
    elseif sum(candidate_best_samples==1) ~= 1
        error('one and only one reference sample is expected in input candidate_best_sample (vector with one element == 1'); 
    else
        % skip steps 1,2 and 3 of LPAC and proceed with input reference
        % ploidies and purities from reference sample
        do_steps123 = 0;
        LPAC_success = 1;
        ploidies(candidate_best_samples==1) = ploidy_ref;
        purities(candidate_best_samples==1) = purity_ref;
    end
    
    % check if data is normalized
    for i = 1:Ns
        if abs(mean(bin_val(:,i)) - 1) > 0.05
            tekst = sprintf( 'input bin values (column %d) are not-normalized, i.e. have an average value of 1',i);
            if mean(bin_val(:,i)) < 0.5 &&  abs( mean(2.^seg_val(:,i)) - 1) < 0.1
                tekst = sprintf('%s\n input data appears to consist of log transformed data. Please provide non-log transformed data',tekst);
            end
            error(tekst);
        end
    end
    
    if do_steps123
        % STEP 1: single-sample inference of absolute CNAs
        for i = 1:Ns
            [CNHout(i), ploidies(i), purities(i)] = CNH(bin_val(:,i), ones(size(bin_val,1),1), [], []);
        end

        % STEP 2: criterion for success of single-sample inference
        step1_success = ~(ploidies < 2.5 & purities == 1) & candidate_best_samples;
        %step1_success = ones(Ns,1);
        if sum(step1_success) == 0
            LPAC_success = false;
            best_sample = false(Ns,1);
        else
            LPAC_success = true;
        end
    end
    % Continue to STEP 3 if for at least one sample the inference of STEP 2
    % was successful
    
    % STEP 3: identification of best single-sample
    if LPAC_success
        if do_steps123
            id = find( CNHout == min(CNHout(step1_success)),1);
            best_sample = false(Ns,1);
            best_sample(id) = true;
        else
            id = find(candidate_best_samples == 1);
            best_sample = candidate_best_samples == 1;
        end
        % STEP 4: infer CNAs by alignment to absolute CNAs of best samples of
        % relative CNAs of other samples.

        % define reference absolute CNA
        a1 = (purities(id) * ploidies(id) + (1-purities(id))*2) / purities(id);
        a2 = -2*(1-purities(id))/purities(id);
        cna_abs_ref = a1 * bin_val(:,id) + a2;
        
        % define purities and ploidies for grid search
        alphas = 0.02:0.01:1;
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
        
    
    
    
    
        
   
    


