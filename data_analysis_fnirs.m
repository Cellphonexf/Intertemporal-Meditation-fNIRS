%% fNIRS data analysis
%% Required to run in MATLAB version R2017b
%% Required to load Homer3 processed data files
%% Programmed by Feng Xiao (2024.6.25)
clear all,
clc,
%% Parameter settings
samplingrate = 5.1;
Channels = 1:36;
MediCond = [2 3]; %2 for mindfulness; 3 for intertemporal meditation
Mindfulness_first_subjs = [6 10 11 15 16 32 34 35 36 38 42 43 44 46 47 48 53 54 56 58 60 62 64 66 68 70 72];
Intertemporal_first_subjs = [17 18 19 20 21 23 24 26 27 29 33 37 40 41 50 51 52 55 57 59 61 65 67 69 71];
subj = [Mindfulness_first_subjs, Intertemporal_first_subjs];
sclConc = 1e6; %convert Conc from Molar to uMolar
%% Session 1
cd session1\derivatives\homer\
interval1 = [0 300]; %in second, for calculating means
ses1_mindfulness_hbo = [];
ses1_intertemporal_hbo = [];
for i_subj = subj
load(num2str(i_subj), '-mat')
channelnum=size(output.dc.dataTimeSeries,3)/36;
setlength=channelnum*36; 
t=output.dcAvg.time;
temp_chPrune = cell2mat(output.misc.mlActAuto);
chPrune = temp_chPrune(1:36,3); %1:valid channels; 0: bad channels
subj_mindfulness_hbo = [];
subj_intertemporal_hbo = [];
HbO_mindfulness = [];
HbO_intertemporal = [];
    for j_Cond = MediCond
        for k_Ch = Channels
            SigHbO = output.dcAvg.dataTimeSeries(:,(j_Cond-1)*setlength*3+(k_Ch-1)*3+1)*sclConc;
            intind = [find(t >= interval1(1),1,'first') find(t <= interval1(end),1,'last')]; %extract data from 15s before stimulus onset to 300s after it
            HbO_bas = mean(SigHbO(1:intind(1)));
            HbO_act = mean(SigHbO(intind(1):intind(end)));
            HbO_act = SigHbO(intind(1):intind(end));
            HbO = HbO_act - HbO_bas; %calculate the relative activation based on the baseline
            if chPrune(k_Ch,1) == 1
                if ismember(i_subj,Mindfulness_first_subjs)==1&&j_Cond==2 || ismember(i_subj,Intertemporal_first_subjs)==1&&j_Cond==2
                    HbO_mindfulness = [HbO_mindfulness, HbO];
                elseif ismember(i_subj,Intertemporal_first_subjs)==1&&j_Cond==3 || ismember(i_subj,Mindfulness_first_subjs)==1&&j_Cond==3
                    HbO_intertemporal = [HbO_intertemporal, HbO];
                end
            else
                if ismember(i_subj,Mindfulness_first_subjs)==1&&j_Cond==2 || ismember(i_subj,Intertemporal_first_subjs)==1&&j_Cond==2
                    HbO_mindfulness = [HbO_mindfulness, zeros(size(HbO,1),1)];
                elseif ismember(i_subj,Intertemporal_first_subjs)==1&&j_Cond==3 || ismember(i_subj,Mindfulness_first_subjs)==1&&j_Cond==3
                    HbO_intertemporal = [HbO_intertemporal, zeros(size(HbO,1),1)];
                end
            end
        end
    end
    subj_mindfulness_hbo = [subj_mindfulness_hbo, HbO_mindfulness];
    subj_intertemporal_hbo = [subj_intertemporal_hbo, HbO_intertemporal];
    
    ses1_mindfulness_hbo = [ses1_mindfulness_hbo, subj_mindfulness_hbo];
    ses1_intertemporal_hbo = [ses1_intertemporal_hbo, subj_intertemporal_hbo]; %col: channels * subj
end

clear subj_mindfulness_hbo subj_intertemporal_hbo HbO_mindfulness HbO_intertemporal

Ses1_mindfulness_hbo = [];
Ses1_intertemporal_hbo = [];
Ses1_mindfulness_hbo_se = [];
Ses1_intertemporal_hbo_se = [];
Ses1_channel_hbo_t = ones(size(ses1_mindfulness_hbo, 1), 36);
Ses1_channel_hbo_p = ones(size(ses1_mindfulness_hbo, 1), 36);
for i = Channels
    i_temp_mindfulness_hbo = [];
    i_temp_intertemporal_hbo = [];
    i_mean_mindfulness_hbo = [];
    i_mean_intertemporal_hbo = [];
    i_se_mindfulness_hbo = [];
    i_se_intertemporal_hbo = [];
    ch_temp_t = ones(size(ses1_mindfulness_hbo, 1), 1);
    ch_temp_p = ones(size(ses1_mindfulness_hbo, 1), 1);
    for j = 1:size(subj, 2)
        temp_mindfulness_hbo = ses1_mindfulness_hbo(:, i+size(Channels, 2)*(j-1));
        temp_intertemporal_hbo = ses1_intertemporal_hbo(:, i+size(Channels, 2)*(j-1));
        
        i_temp_mindfulness_hbo = [i_temp_mindfulness_hbo, temp_mindfulness_hbo];
        i_temp_intertemporal_hbo = [i_temp_intertemporal_hbo, temp_intertemporal_hbo];
    end
    i_temp_mindfulness_hbo = i_temp_mindfulness_hbo(:, any(i_temp_mindfulness_hbo)); 
    i_temp_intertemporal_hbo = i_temp_intertemporal_hbo(:, any(i_temp_intertemporal_hbo)); %delete column with all zeros
    
    i_mean_mindfulness_hbo = mean(i_temp_mindfulness_hbo, 2);%average subjects' data for each channel
    i_mean_intertemporal_hbo = mean(i_temp_intertemporal_hbo, 2); %average subjects' data for each channel
    
    i_se_mindfulness_hbo = std(i_temp_mindfulness_hbo, 0, 2)./sqrt(size(i_temp_mindfulness_hbo, 2)); %calculate standard error
    i_se_intertemporal_hbo = std(i_temp_intertemporal_hbo, 0, 2)./sqrt(size(i_temp_intertemporal_hbo, 2)); %calculate standard error   
    for k = 1:size(ses1_mindfulness_hbo,1)
        [h,p_hbo,ci,stats_hbo] = ttest(i_temp_mindfulness_hbo(k,:), i_temp_intertemporal_hbo(k,:)); %paired ttest between two meditation for each time point in each channel
        ch_temp_hbo_t(k) = stats_hbo.tstat; %t>0 denote mindfulness > intertemporal
        ch_temp_hbo_p(k) = p_hbo;
    end
    Ses1_mindfulness_hbo = [Ses1_mindfulness_hbo, i_mean_mindfulness_hbo];
    Ses1_intertemporal_hbo = [Ses1_intertemporal_hbo, i_mean_intertemporal_hbo];
    Ses1_mindfulness_hbo_se = [Ses1_mindfulness_hbo_se, i_se_mindfulness_hbo];
    Ses1_intertemporal_hbo_se = [Ses1_intertemporal_hbo_se, i_se_intertemporal_hbo];
    Ses1_channel_hbo_t(:, i) = ch_temp_hbo_t;
    Ses1_channel_hbo_p(:, i) = ch_temp_hbo_p;
end

clear i_temp_mindfulness_hbo
clear i_temp_intertemporal_hbo
clear i_mean_mindfulness_hbo 
clear i_mean_intertemporal_hbo 
clear i_se_mindfulness_hbo 
clear i_se_intertemporal_hbo
clear temp_mindfulness_hbo 
clear temp_intertemporal_hbo 
clear ch_temp_hbo_t 
clear ch_temp_hbo_p 
clear h ci p_hbo stats_hbo 

tTest_Ses1_mindfulness_hbo = zeros(5,36);
tTest_Ses1_intertemporal_hbo = zeros(5,36);
tTest_Ses1_contrast_hbo = zeros(5,36);
for i = 1:36
    temp_col_m = Ses1_mindfulness_hbo(:, i);
    [h,p_m,ci_m,stats_m] = ttest(temp_col_m); %one-sample ttest
    tTest_Ses1_mindfulness_hbo(1, i) = stats_m.tstat;
    tTest_Ses1_mindfulness_hbo(2, i) = p_m;
    tTest_Ses1_mindfulness_hbo(3, i) = ci_m(1,1);
    tTest_Ses1_mindfulness_hbo(4, i) = ci_m(2,1);
    if p_m < 0.001 %alpha level of .001
       tTest_Ses1_mindfulness_hbo(5, i) = 1; %accept h1
    else
       tTest_Ses1_mindfulness_hbo(5, i) = 0; %accept h0 
    end
    
    temp_col_i = Ses1_intertemporal_hbo(:, i);
    [h,p_i,ci_i,stats_i] = ttest(temp_col_i); %one-sample ttest
    tTest_Ses1_intertemporal_hbo(1, i) = stats_i.tstat;
    tTest_Ses1_intertemporal_hbo(2, i) = p_i;
    tTest_Ses1_intertemporal_hbo(3, i) = ci_i(1,1);
    tTest_Ses1_intertemporal_hbo(4, i) = ci_i(2,1);
    if p_i < 0.001 %alpha level of .001
       tTest_Ses1_intertemporal_hbo(5, i) = 1; %accept h1
    else
       tTest_Ses1_intertemporal_hbo(5, i) = 0; %accept h0 
    end

    [h,p_c,ci_c,stats_c] = ttest(temp_col_i, temp_col_m); %intertemporal - mindfulness; paired ttest
    tTest_Ses1_contrast_hbo(1,i) = stats_c.tstat; 
    tTest_Ses1_contrast_hbo(2,i) = p_c;
    tTest_Ses1_contrast_hbo(3,i) = ci_c(1,1);
    tTest_Ses1_contrast_hbo(4,i) = ci_c(2,1); 
    if p_c < 0.001 %alpha level of .001
       tTest_Ses1_contrast_hbo(5, i) = 1; %accept h1
    else
       tTest_Ses1_contrast_hbo(5, i) = 0; %accept h0 
    end
end

clear temp_col_m temp_col_i
clear h p_m p_i p_c ci_m ci_i ci_c stats_m stats_i stats_c
%% Session 2
cd ..\..\..
cd session2\derivatives\homer\
interval2 = [0 60]; %in second, for calculating means
ses2_mindfulness_hbo = [];
ses2_intertemporal_hbo = [];
for i_subj = subj
load(num2str(i_subj), '-mat')
channelnum=size(output.dc.dataTimeSeries,3)/36;
setlength=channelnum*36; 
t=output.dcAvg.time;
temp_chPrune = cell2mat(output.misc.mlActAuto);
chPrune = temp_chPrune(1:36,3); %1:valid channels; 0: bad channels
subj_mindfulness_hbo = [];
subj_intertemporal_hbo = [];
HbO_mindfulness = [];
HbO_intertemporal = [];
    for j_Cond = MediCond
        for k_Ch = Channels
            SigHbO = output.dcAvg.dataTimeSeries(:,(j_Cond-1)*setlength*3+(k_Ch-1)*3+1)*sclConc;
            intind = [find(t >= interval2(1),1,'first') find(t <= interval2(end),1,'last')]; %extract data from 15s before stimulus onset to 60s after it
            HbO_bas = mean(SigHbO(1:intind(1))); %baseline
            HbO_act = SigHbO(intind(1):intind(end));
            HbO = HbO_act - HbO_bas; %calculate the relative activation based on the baseline
            if chPrune(k_Ch,1) == 1
                if ismember(i_subj,Mindfulness_first_subjs)==1&&j_Cond==2 || ismember(i_subj,Intertemporal_first_subjs)==1&&j_Cond==2
                    HbO_mindfulness = [HbO_mindfulness, HbO];
                elseif ismember(i_subj,Intertemporal_first_subjs)==1&&j_Cond==3 || ismember(i_subj,Mindfulness_first_subjs)==1&&j_Cond==3
                    HbO_intertemporal = [HbO_intertemporal, HbO];
                end
            else
                if ismember(i_subj,Mindfulness_first_subjs)==1&&j_Cond==2 || ismember(i_subj,Intertemporal_first_subjs)==1&&j_Cond==2
                    HbO_mindfulness = [HbO_mindfulness, zeros(size(HbO,1),1)];
                elseif ismember(i_subj,Intertemporal_first_subjs)==1&&j_Cond==3 || ismember(i_subj,Mindfulness_first_subjs)==1&&j_Cond==3
                    HbO_intertemporal = [HbO_intertemporal, zeros(size(HbO,1),1)];
                end
            end
        end
    end
    subj_mindfulness_hbo = [subj_mindfulness_hbo, HbO_mindfulness];
    subj_intertemporal_hbo = [subj_intertemporal_hbo, HbO_intertemporal];
    
    ses2_mindfulness_hbo = [ses2_mindfulness_hbo, subj_mindfulness_hbo];
    ses2_intertemporal_hbo = [ses2_intertemporal_hbo, subj_intertemporal_hbo]; %col: channels * subj
end

clear subj_mindfulness_hbo subj_intertemporal_hbo HbO_mindfulness HbO_intertemporal

Ses2_mindfulness_hbo = [];
Ses2_intertemporal_hbo = [];
Ses2_mindfulness_hbo_se = [];
Ses2_intertemporal_hbo_se = [];
Ses2_channel_hbo_t = ones(size(ses2_mindfulness_hbo, 1), 36);
Ses2_channel_hbo_p = ones(size(ses2_mindfulness_hbo, 1), 36);
for i = Channels
    i_temp_mindfulness_hbo = [];
    i_temp_intertemporal_hbo = [];
    i_mean_mindfulness_hbo = [];
    i_mean_intertemporal_hbo = [];
    i_se_mindfulness_hbo = [];
    i_se_intertemporal_hbo = [];
    ch_temp_t = ones(size(ses2_mindfulness_hbo, 1), 1);
    ch_temp_p = ones(size(ses2_mindfulness_hbo, 1), 1);
    for j = 1:size(subj, 2)
        temp_mindfulness_hbo = ses2_mindfulness_hbo(:, i+size(Channels, 2)*(j-1));
        temp_intertemporal_hbo = ses2_intertemporal_hbo(:, i+size(Channels, 2)*(j-1));
        
        i_temp_mindfulness_hbo = [i_temp_mindfulness_hbo, temp_mindfulness_hbo];
        i_temp_intertemporal_hbo = [i_temp_intertemporal_hbo, temp_intertemporal_hbo];
    end
    i_temp_mindfulness_hbo = i_temp_mindfulness_hbo(:, any(i_temp_mindfulness_hbo)); %delete column with all zeros
    i_temp_intertemporal_hbo = i_temp_intertemporal_hbo(:, any(i_temp_intertemporal_hbo)); %delete column with all zeros
    
    i_mean_mindfulness_hbo = mean(i_temp_mindfulness_hbo, 2); %average subjects' data for each channel
    i_mean_intertemporal_hbo = mean(i_temp_intertemporal_hbo, 2); %average subjects' data for each channel
    
    i_se_mindfulness_hbo = std(i_temp_mindfulness_hbo, 0, 2)./sqrt(size(i_temp_mindfulness_hbo, 2)); %calculate standard error
    i_se_intertemporal_hbo = std(i_temp_intertemporal_hbo, 0, 2)./sqrt(size(i_temp_intertemporal_hbo, 2)); %calculate standard error   
    for k = 1:size(ses2_mindfulness_hbo,1)
        [h,p_hbo,ci,stats_hbo] = ttest(i_temp_mindfulness_hbo(k,:), i_temp_intertemporal_hbo(k,:)); %paired ttest between two meditation for each time point in each channel
        ch_temp_hbo_t(k) = stats_hbo.tstat; %t>0 denote mindfulness > intertemporal
        ch_temp_hbo_p(k) = p_hbo;
    end
    Ses2_mindfulness_hbo = [Ses2_mindfulness_hbo, i_mean_mindfulness_hbo];
    Ses2_intertemporal_hbo = [Ses2_intertemporal_hbo, i_mean_intertemporal_hbo];
    Ses2_mindfulness_hbo_se = [Ses2_mindfulness_hbo_se, i_se_mindfulness_hbo];
    Ses2_intertemporal_hbo_se = [Ses2_intertemporal_hbo_se, i_se_intertemporal_hbo];
    Ses2_channel_hbo_t(:, i) = ch_temp_hbo_t;
    Ses2_channel_hbo_p(:, i) = ch_temp_hbo_p;
end

clear i_temp_mindfulness_hbo
clear i_temp_intertemporal_hbo 
clear i_mean_mindfulness_hbo 
clear i_mean_intertemporal_hbo
clear i_se_mindfulness_hbo 
clear i_se_intertemporal_hbo 
clear temp_mindfulness_hbo 
clear temp_intertemporal_hbo 
clear ch_temp_hbo_t
clear ch_temp_hbo_p 
clear h ci p_hbo stats_hbo

tTest_Ses2_mindfulness_hbo = zeros(5,36);
tTest_Ses2_intertemporal_hbo = zeros(5,36);
tTest_Ses2_contrast_hbo = zeros(5,36);
for i = 1:36
    temp_col_m = Ses2_mindfulness_hbo(:, i);
    [h,p_m,ci_m,stats_m] = ttest(temp_col_m); %one-sample ttest
    tTest_Ses2_mindfulness_hbo(1, i) = stats_m.tstat;
    tTest_Ses2_mindfulness_hbo(2, i) = p_m;
    tTest_Ses2_mindfulness_hbo(3, i) = ci_m(1,1);
    tTest_Ses2_mindfulness_hbo(4, i) = ci_m(2,1);
    if p_m < 0.001 %alpha level of .001
       tTest_Ses2_mindfulness_hbo(5, i) = 1; %accept h1
    else
       tTest_Ses2_mindfulness_hbo(5, i) = 0; %accept h0 
    end
    
    temp_col_i = Ses2_intertemporal_hbo(:, i);
    [h,p_i,ci_i,stats_i] = ttest(temp_col_i); %one-sample ttest
    tTest_Ses2_intertemporal_hbo(1, i) = stats_i.tstat;
    tTest_Ses2_intertemporal_hbo(2, i) = p_i;
    tTest_Ses2_intertemporal_hbo(3, i) = ci_i(1,1);
    tTest_Ses2_intertemporal_hbo(4, i) = ci_i(2,1);
    if p_i < 0.001 %alpha level of .001
       tTest_Ses2_intertemporal_hbo(5, i) = 1; %accept h1
    else
       tTest_Ses2_intertemporal_hbo(5, i) = 0; %accept h0 
    end

    [h,p_c,ci_c,stats_c] = ttest(temp_col_i, temp_col_m); %intertemporal - mindfulness; paired ttest
    tTest_Ses2_contrast_hbo(1,i) = stats_c.tstat; 
    tTest_Ses2_contrast_hbo(2,i) = p_c;
    tTest_Ses2_contrast_hbo(3,i) = ci_c(1,1);
    tTest_Ses2_contrast_hbo(4,i) = ci_c(2,1);
    if p_c < 0.001 %alpha level of .001
       tTest_Ses2_contrast_hbo(5, i) = 1; %accept h1
    else
       tTest_Ses2_contrast_hbo(5, i) = 0; %accept h0 
    end
end

clear temp_col_m temp_col_i
clear h p_m p_i p_c ci_m ci_i ci_c stats_m stats_i stats_c
%% Session 3
cd ..\..\..
cd session3\derivatives\homer\
interval3 = [0 60]; %in second, for calculating means
ses3_mindfulness_hbo = [];
ses3_intertemporal_hbo = [];
for i_subj = subj
load(num2str(i_subj), '-mat')
channelnum=size(output.dc.dataTimeSeries,3)/36;
setlength=channelnum*36; 
t=output.dcAvg.time;
temp_chPrune = cell2mat(output.misc.mlActAuto);
chPrune = temp_chPrune(1:36,3); %1:valid channels; 0: bad channels
subj_mindfulness_hbo = [];
subj_intertemporal_hbo = [];
HbO_mindfulness = [];
HbO_intertemporal = [];
    for j_Cond = MediCond
        for k_Ch = Channels
            SigHbO = output.dcAvg.dataTimeSeries(:,(j_Cond-1)*setlength*3+(k_Ch-1)*3+1)*sclConc;
            intind = [find(t >= interval3(1),1,'first') find(t <= interval3(end),1,'last')]; %extract data from 15s before stimulus onset to 60s after it
            HbO_bas = mean(SigHbO(1:intind(1))); %baseline
            HbO_act = SigHbO(intind(1):intind(end));
            HbO = HbO_act - HbO_bas; %calculate the relative activation based on the baseline
            if chPrune(k_Ch,1) == 1
                if ismember(i_subj,Mindfulness_first_subjs)==1&&j_Cond==2 || ismember(i_subj,Intertemporal_first_subjs)==1&&j_Cond==2
                    HbO_mindfulness = [HbO_mindfulness, HbO];
                elseif ismember(i_subj,Intertemporal_first_subjs)==1&&j_Cond==3 || ismember(i_subj,Mindfulness_first_subjs)==1&&j_Cond==3
                    HbO_intertemporal = [HbO_intertemporal, HbO];
                end
            else
                if ismember(i_subj,Mindfulness_first_subjs)==1&&j_Cond==2 || ismember(i_subj,Intertemporal_first_subjs)==1&&j_Cond==2
                    HbO_mindfulness = [HbO_mindfulness, zeros(size(HbO,1),1)];
                elseif ismember(i_subj,Intertemporal_first_subjs)==1&&j_Cond==3 || ismember(i_subj,Mindfulness_first_subjs)==1&&j_Cond==3
                    HbO_intertemporal = [HbO_intertemporal, zeros(size(HbO,1),1)];
                end
            end
        end
    end
    subj_mindfulness_hbo = [subj_mindfulness_hbo, HbO_mindfulness];
    subj_intertemporal_hbo = [subj_intertemporal_hbo, HbO_intertemporal];
    
    ses3_mindfulness_hbo = [ses3_mindfulness_hbo, subj_mindfulness_hbo];
    ses3_intertemporal_hbo = [ses3_intertemporal_hbo, subj_intertemporal_hbo]; %col: channels * subj
end

clear subj_mindfulness_hbo subj_intertemporal_hbo HbO_mindfulness HbO_intertemporal 

Ses3_mindfulness_hbo = [];
Ses3_intertemporal_hbo = [];
Ses3_mindfulness_hbo_se = [];
Ses3_intertemporal_hbo_se = [];
Ses3_channel_hbo_t = ones(size(ses3_mindfulness_hbo, 1), 36);
Ses3_channel_hbo_p = ones(size(ses3_mindfulness_hbo, 1), 36);
for i = Channels
    i_temp_mindfulness_hbo = [];
    i_temp_intertemporal_hbo = [];
    i_mean_mindfulness_hbo = [];
    i_mean_intertemporal_hbo = [];
    i_se_mindfulness_hbo = [];
    i_se_intertemporal_hbo = [];
    ch_temp_t = ones(size(ses3_mindfulness_hbo, 1), 1);
    ch_temp_p = ones(size(ses3_mindfulness_hbo, 1), 1);
    for j = 1:size(subj, 2)
        temp_mindfulness_hbo = ses3_mindfulness_hbo(:, i+size(Channels, 2)*(j-1));
        temp_intertemporal_hbo = ses3_intertemporal_hbo(:, i+size(Channels, 2)*(j-1));
        
        i_temp_mindfulness_hbo = [i_temp_mindfulness_hbo, temp_mindfulness_hbo];
        i_temp_intertemporal_hbo = [i_temp_intertemporal_hbo, temp_intertemporal_hbo];
    end
    i_temp_mindfulness_hbo = i_temp_mindfulness_hbo(:, any(i_temp_mindfulness_hbo)); %delete column with all zeros
    i_temp_intertemporal_hbo = i_temp_intertemporal_hbo(:, any(i_temp_intertemporal_hbo)); %delete column with all zeros
    
    i_mean_mindfulness_hbo = mean(i_temp_mindfulness_hbo, 2); %average subjects' data for each channel
    i_mean_intertemporal_hbo = mean(i_temp_intertemporal_hbo, 2); %average subjects' data for each channel
    
    i_se_mindfulness_hbo = std(i_temp_mindfulness_hbo, 0, 2)./sqrt(size(i_temp_mindfulness_hbo, 2)); %calculate standard error
    i_se_intertemporal_hbo = std(i_temp_intertemporal_hbo, 0, 2)./sqrt(size(i_temp_intertemporal_hbo, 2)); %calculate standard error   
    for k = 1:size(ses3_mindfulness_hbo,1)
        [h,p_hbo,ci,stats_hbo] = ttest(i_temp_mindfulness_hbo(k,:), i_temp_intertemporal_hbo(k,:)); %paired ttest between two meditation for each time point in each channel
        ch_temp_hbo_t(k) = stats_hbo.tstat; %t>0 denote mindfulness > intertemporal
        ch_temp_hbo_p(k) = p_hbo;
    end
    Ses3_mindfulness_hbo = [Ses3_mindfulness_hbo, i_mean_mindfulness_hbo];
    Ses3_intertemporal_hbo = [Ses3_intertemporal_hbo, i_mean_intertemporal_hbo];
    Ses3_mindfulness_hbo_se = [Ses3_mindfulness_hbo_se, i_se_mindfulness_hbo];
    Ses3_intertemporal_hbo_se = [Ses3_intertemporal_hbo_se, i_se_intertemporal_hbo];
    Ses3_channel_hbo_t(:, i) = ch_temp_hbo_t;
    Ses3_channel_hbo_p(:, i) = ch_temp_hbo_p;
end

clear i_temp_mindfulness_hbo
clear i_temp_intertemporal_hbo 
clear i_mean_mindfulness_hbo 
clear i_mean_intertemporal_hbo
clear i_se_mindfulness_hbo 
clear i_se_intertemporal_hbo
clear temp_mindfulness_hbo 
clear temp_intertemporal_hbo 
clear ch_temp_hbo_t 
clear ch_temp_hbo_p 
clear h ci p_hbo stats_hbo

tTest_Ses3_mindfulness_hbo = zeros(5,36);
tTest_Ses3_intertemporal_hbo = zeros(5,36);
tTest_Ses3_contrast_hbo = zeros(5,36);
for i = 1:36
    temp_col_m = Ses3_mindfulness_hbo(:, i);
    [h,p_m,ci_m,stats_m] = ttest(temp_col_m); %one-sample ttest
    tTest_Ses3_mindfulness_hbo(1, i) = stats_m.tstat;
    tTest_Ses3_mindfulness_hbo(2, i) = p_m;
    tTest_Ses3_mindfulness_hbo(3, i) = ci_m(1,1);
    tTest_Ses3_mindfulness_hbo(4, i) = ci_m(2,1);
    if p_m < 0.001 %alpha level of .001
       tTest_Ses3_mindfulness_hbo(5, i) = 1; %accept h1
    else
       tTest_Ses3_mindfulness_hbo(5, i) = 0; %accept h0 
    end
    
    temp_col_i = Ses3_intertemporal_hbo(:, i);
    [h,p_i,ci_i,stats_i] = ttest(temp_col_i); %one-sample ttest
    tTest_Ses3_intertemporal_hbo(1, i) = stats_i.tstat;
    tTest_Ses3_intertemporal_hbo(2, i) = p_i;
    tTest_Ses3_intertemporal_hbo(3, i) = ci_i(1,1);
    tTest_Ses3_intertemporal_hbo(4, i) = ci_i(2,1);
    if p_i < 0.001 %alpha level of .001
       tTest_Ses3_intertemporal_hbo(5, i) = 1; %accept h1
    else
       tTest_Ses3_intertemporal_hbo(5, i) = 0; %accept h0 
    end

    [h,p_c,ci_c,stats_c] = ttest(temp_col_i, temp_col_m); %intertemporal - mindfulness; paired ttest
    tTest_Ses3_contrast_hbo(1,i) = stats_c.tstat; 
    tTest_Ses3_contrast_hbo(2,i) = p_c;
    tTest_Ses3_contrast_hbo(3,i) = ci_c(1,1);
    tTest_Ses3_contrast_hbo(4,i) = ci_c(2,1);
    if p_c < 0.001 %alpha level of .001
       tTest_Ses3_contrast_hbo(5, i) = 1; %accept h1
    else
       tTest_Ses3_contrast_hbo(5, i) = 0; %accept h0 
    end
end

clear temp_col_m temp_col_i
clear h p_m p_i p_c ci_m ci_i ci_c stats_m stats_i stats_c
%% Total meditation process (session 1+2)
Total12_mindfulness_hbo = [Ses1_mindfulness_hbo; Ses2_mindfulness_hbo];
Total12_intertemporal_hbo = [Ses1_intertemporal_hbo; Ses2_intertemporal_hbo];
Total12_mindfulness_hbo_se = [Ses1_mindfulness_hbo_se; Ses2_mindfulness_hbo_se];
Total12_mindfulness_hbo_upper = Total12_mindfulness_hbo + Total12_mindfulness_hbo_se;
Total12_mindfulness_hbo_lower = Total12_mindfulness_hbo - Total12_mindfulness_hbo_se;
Total12_intertemporal_hbo_se = [Ses1_intertemporal_hbo_se; Ses2_intertemporal_hbo_se];
Total12_intertemporal_hbo_upper = Total12_intertemporal_hbo + Total12_intertemporal_hbo_se;
Total12_intertemporal_hbo_lower = Total12_intertemporal_hbo - Total12_intertemporal_hbo_se;
Total12_channel_hbo_t = [Ses1_channel_hbo_t; Ses2_channel_hbo_t];
Total12_channel_hbo_p = [Ses1_channel_hbo_p; Ses2_channel_hbo_p];

tTest_Total12_mindfulness_hbo = zeros(5,36);
tTest_Total12_intertemporal_hbo = zeros(5,36);
tTest_Total12_contrast_hbo = zeros(5,36);
for i = 1:36
    temp_col_m = Total12_mindfulness_hbo(:, i);
    [h,p_m,ci_m,stats_m] = ttest(temp_col_m); %one-sample ttest
    tTest_Total12_mindfulness_hbo(1, i) = stats_m.tstat;
    tTest_Total12_mindfulness_hbo(2, i) = p_m;
    tTest_Total12_mindfulness_hbo(3, i) = ci_m(1,1);
    tTest_Total12_mindfulness_hbo(4, i) = ci_m(2,1);
    if p_m < 0.001 %alpha level of .001
       tTest_Total12_mindfulness_hbo(5, i) = 1; %accept h1
    else
       tTest_Total12_mindfulness_hbo(5, i) = 0; %accept h0 
    end
    
    temp_col_i = Total12_intertemporal_hbo(:, i);
    [h,p_i,ci_i,stats_i] = ttest(temp_col_i); %one-sample ttest
    tTest_Total12_intertemporal_hbo(1, i) = stats_i.tstat;
    tTest_Total12_intertemporal_hbo(2, i) = p_i;
    tTest_Total12_intertemporal_hbo(3, i) = ci_i(1,1);
    tTest_Total12_intertemporal_hbo(4, i) = ci_i(2,1);
    if p_i < 0.001 %alpha level of .001
       tTest_Total12_intertemporal_hbo(5, i) = 1; %accept h1
    else
       tTest_Total12_intertemporal_hbo(5, i) = 0; %accept h0 
    end

    [h,p_c,ci_c,stats_c] = ttest(temp_col_i, temp_col_m); %intertemporal - mindfulness; paired ttest
    tTest_Total12_contrast_hbo(1,i) = stats_c.tstat; 
    tTest_Total12_contrast_hbo(2,i) = p_c;
    tTest_Total12_contrast_hbo(3,i) = ci_c(1,1);
    tTest_Total12_contrast_hbo(4,i) = ci_c(2,1);
    if p_c < 0.001 %alpha level of .001
       tTest_Total12_contrast_hbo(5, i) = 1; %accept h1
    else
       tTest_Total12_contrast_hbo(5, i) = 0; %accept h0 
    end
end

clear temp_col_m temp_col_i
clear h p_m p_i p_c ci_m ci_i ci_c stats_m stats_i stats_c
%% Total meditation process (session 1+2+3)
Total_mindfulness_hbo = [Ses1_mindfulness_hbo; Ses2_mindfulness_hbo; Ses3_mindfulness_hbo];
Total_intertemporal_hbo = [Ses1_intertemporal_hbo; Ses2_intertemporal_hbo; Ses3_intertemporal_hbo];
Total_mindfulness_hbo_se = [Ses1_mindfulness_hbo_se; Ses2_mindfulness_hbo_se; Ses3_mindfulness_hbo_se];
Total_mindfulness_hbo_upper = Total_mindfulness_hbo + Total_mindfulness_hbo_se;
Total_mindfulness_hbo_lower = Total_mindfulness_hbo - Total_mindfulness_hbo_se;
Total_intertemporal_hbo_se = [Ses1_intertemporal_hbo_se; Ses2_intertemporal_hbo_se; Ses3_intertemporal_hbo_se];
Total_intertemporal_hbo_upper = Total_intertemporal_hbo + Total_intertemporal_hbo_se;
Total_intertemporal_hbo_lower = Total_intertemporal_hbo - Total_intertemporal_hbo_se;
Total_channel_hbo_t = [Ses1_channel_hbo_t; Ses2_channel_hbo_t; Ses3_channel_hbo_t];
Total_channel_hbo_p = [Ses1_channel_hbo_p; Ses2_channel_hbo_p; Ses3_channel_hbo_p];

tTest_Total_mindfulness_hbo = zeros(5,36);
tTest_Total_intertemporal_hbo = zeros(5,36);
tTest_Total_contrast_hbo = zeros(5,36);
for i = 1:36
    temp_col_m = Total_mindfulness_hbo(:, i);
    [h,p_m,ci_m,stats_m] = ttest(temp_col_m); %one-sample ttest
    tTest_Total_mindfulness_hbo(1, i) = stats_m.tstat;
    tTest_Total_mindfulness_hbo(2, i) = p_m;
    tTest_Total_mindfulness_hbo(3, i) = ci_m(1,1);
    tTest_Total_mindfulness_hbo(4, i) = ci_m(2,1);
    if p_m < 0.001 %alpha level of .001
       tTest_Total_mindfulness_hbo(5, i) = 1; %accept h1
    else
       tTest_Total_mindfulness_hbo(5, i) = 0; %accept h0 
    end
    
    temp_col_i = Total_intertemporal_hbo(:, i);
    [h,p_i,ci_i,stats_i] = ttest(temp_col_i); %one-sample ttest
    tTest_Total_intertemporal_hbo(1, i) = stats_i.tstat;
    tTest_Total_intertemporal_hbo(2, i) = p_i;
    tTest_Total_intertemporal_hbo(3, i) = ci_i(1,1);
    tTest_Total_intertemporal_hbo(4, i) = ci_i(2,1);
    if p_i < 0.001 %alpha level of .001
       tTest_Total_intertemporal_hbo(5, i) = 1; %accept h1
    else
       tTest_Total_intertemporal_hbo(5, i) = 0; %accept h0 
    end

    [h,p_c,ci_c,stats_c] = ttest(temp_col_i, temp_col_m); %intertemporal - mindfulness; paired ttest
    tTest_Total_contrast_hbo(1,i) = stats_c.tstat; 
    tTest_Total_contrast_hbo(2,i) = p_c;
    tTest_Total_contrast_hbo(3,i) = ci_c(1,1);
    tTest_Total_contrast_hbo(4,i) = ci_c(2,1);
    if p_c < 0.001 %alpha level of .001
       tTest_Total_contrast_hbo(5, i) = 1; %accept h1
    else
       tTest_Total_contrast_hbo(5, i) = 0; %accept h0 
    end
end

clear temp_col_m temp_col_i
clear h p_m p_i p_c ci_m ci_i ci_c stats_m stats_i stats_c
%% Effect size calculation (Cohen's d)
hbo_mindfulness_ses1 = zeros(3,36);
hbo_intertemporal_ses1 = zeros(3,36);

hbo_mindfulness_ses12 = zeros(3,36);
hbo_intertemporal_ses12 = zeros(3,36);

hbo_mindfulness = zeros(3,36);
hbo_intertemporal = zeros(3,36);

hbo_mindfulness_ses2 = zeros(3,36);
hbo_intertemporal_ses2 = zeros(3,36);

hbo_mindfulness_ses3 = zeros(3,36);
hbo_intertemporal_ses3 = zeros(3,36);

for i = 1:36
    temp_hbo_mean = mean(Total_mindfulness_hbo(:,i));
    temp_hbo_se = std(Total_mindfulness_hbo(:,i))/sqrt(size(Total_mindfulness_hbo(:,i),1));
    hbo_mindfulness(1,i) = temp_hbo_mean;
    hbo_mindfulness(2,i) = temp_hbo_se;
    hbo_mindfulness(3,i) = std(Total_mindfulness_hbo(:,i));
 
    clear temp_hbo_mean temp_hbo_se
    
    temp_hbo_mean = mean(Total_intertemporal_hbo(:,i));
    temp_hbo_se = std(Total_intertemporal_hbo(:,i))/sqrt(size(Total_intertemporal_hbo(:,i),1));
    hbo_intertemporal(1,i) = temp_hbo_mean;
    hbo_intertemporal(2,i) = temp_hbo_se;
    hbo_intertemporal(3,i) = std(Total_intertemporal_hbo(:,i));
    
    clear temp_hbo_mean temp_hbo_se
    
    temp_hbo_mean = mean(Total12_mindfulness_hbo(:,i));
    temp_hbo_se = std(Total12_mindfulness_hbo(:,i))/sqrt(size(Total12_mindfulness_hbo(:,i),1));
    hbo_mindfulness_ses12(1,i) = temp_hbo_mean;
    hbo_mindfulness_ses12(2,i) = temp_hbo_se;
    hbo_mindfulness_ses12(3,i) = std(Total12_mindfulness_hbo(:,i));
 
    clear temp_hbo_mean temp_hbo_se
    
    temp_hbo_mean = mean(Total12_intertemporal_hbo(:,i));
    temp_hbo_se = std(Total12_intertemporal_hbo(:,i))/sqrt(size(Total12_intertemporal_hbo(:,i),1));
    hbo_intertemporal_ses12(1,i) = temp_hbo_mean;
    hbo_intertemporal_ses12(2,i) = temp_hbo_se;
    hbo_intertemporal_ses12(3,i) = std(Total12_intertemporal_hbo(:,i));
    
    clear temp_hbo_mean temp_hbo_se
    
    temp_hbo_mean = mean(Ses1_mindfulness_hbo(:,i));
    temp_hbo_se = std(Ses1_mindfulness_hbo(:,i))/sqrt(size(Ses1_mindfulness_hbo(:,i),1));
    hbo_mindfulness_ses1(1,i) = temp_hbo_mean;
    hbo_mindfulness_ses1(2,i) = temp_hbo_se;
    hbo_mindfulness_ses1(3,i) = std(Ses1_mindfulness_hbo(:,i));
    
    clear temp_hbo_mean temp_hbo_se
    
    temp_hbo_mean = mean(Ses1_intertemporal_hbo(:,i));
    temp_hbo_se = std(Ses1_intertemporal_hbo(:,i))/sqrt(size(Ses1_intertemporal_hbo(:,i),1));
    hbo_intertemporal_ses1(1,i) = temp_hbo_mean;
    hbo_intertemporal_ses1(2,i) = temp_hbo_se;
    hbo_intertemporal_ses1(3,i) = std(Ses1_intertemporal_hbo(:,i));
    
    clear temp_hbo_mean temp_hbo_se
    
    temp_hbo_mean = mean(Ses2_mindfulness_hbo(:,i));
    temp_hbo_se = std(Ses2_mindfulness_hbo(:,i))/sqrt(size(Ses2_mindfulness_hbo(:,i),1));
    hbo_mindfulness_ses2(1,i) = temp_hbo_mean;
    hbo_mindfulness_ses2(2,i) = temp_hbo_se;
    hbo_mindfulness_ses2(3,i) = std(Ses2_mindfulness_hbo(:,i));
    
    clear temp_hbo_mean temp_hbo_se
    
    temp_hbo_mean = mean(Ses2_intertemporal_hbo(:,i));
    temp_hbo_se = std(Ses2_intertemporal_hbo(:,i))/sqrt(size(Ses2_intertemporal_hbo(:,i),1));
    hbo_intertemporal_ses2(1,i) = temp_hbo_mean;
    hbo_intertemporal_ses2(2,i) = temp_hbo_se;
    hbo_intertemporal_ses2(3,i) = std(Ses2_intertemporal_hbo(:,i));
    
    clear temp_hbo_mean temp_hbo_se
    
    temp_hbo_mean = mean(Ses3_mindfulness_hbo(:,i));
    temp_hbo_se = std(Ses3_mindfulness_hbo(:,i))/sqrt(size(Ses3_mindfulness_hbo(:,i),1));
    hbo_mindfulness_ses3(1,i) = temp_hbo_mean;
    hbo_mindfulness_ses3(2,i) = temp_hbo_se;
    hbo_mindfulness_ses3(3,i) = std(Ses3_mindfulness_hbo(:,i));
    
    clear temp_hbo_mean temp_hbo_se
    
    temp_hbo_mean = mean(Ses3_intertemporal_hbo(:,i));
    temp_hbo_se = std(Ses3_intertemporal_hbo(:,i))/sqrt(size(Ses3_intertemporal_hbo(:,i),1));
    hbo_intertemporal_ses3(1,i) = temp_hbo_mean;
    hbo_intertemporal_ses3(2,i) = temp_hbo_se;
    hbo_intertemporal_ses3(3,i) = std(Ses3_intertemporal_hbo(:,i));
    
    clear temp_hbo_mean temp_hbo_se
end

efs_hbo_contrast = zeros(5,36); %overall,ses1,ses2,ses3,ses1+2
n_overall = size(Total_intertemporal_hbo,1);
n_overall12 = size(Total12_intertemporal_hbo,1);
n_ses1 = size(Ses1_intertemporal_hbo,1);
n_ses2 = size(Ses2_intertemporal_hbo,1);
n_ses3 = size(Ses3_intertemporal_hbo,1);
for i = 1:36
    temp_mean_overall = hbo_intertemporal(1,i) - hbo_mindfulness(1,i); %Intertemporal - Mindfulness
    temp_mean_overall12 = hbo_intertemporal_ses12(1,i) - hbo_mindfulness_ses12(1,i); %Intertemporal - Mindfulness
    temp_mean_ses1 = hbo_intertemporal_ses1(1,i) - hbo_mindfulness_ses1(1,i);
    temp_mean_ses2 = hbo_intertemporal_ses2(1,i) - hbo_mindfulness_ses2(1,i);
    temp_mean_ses3 = hbo_intertemporal_ses3(1,i) - hbo_mindfulness_ses3(1,i);
    
    temp_poolstd_overall = sqrt((hbo_intertemporal(3,i)^2 + hbo_mindfulness(3,i)^2) / 2);
    temp_poolstd_overall12 = sqrt((hbo_intertemporal_ses12(3,i)^2 + hbo_mindfulness_ses12(3,i)^2) / 2);
    temp_poolstd_ses1 = sqrt((hbo_intertemporal_ses1(3,i)^2 + hbo_mindfulness_ses1(3,i)^2) / 2);
    temp_poolstd_ses2 = sqrt((hbo_intertemporal_ses2(3,i)^2 + hbo_mindfulness_ses2(3,i)^2) / 2);
    temp_poolstd_ses3 = sqrt((hbo_intertemporal_ses3(3,i)^2 + hbo_mindfulness_ses3(3,i)^2) / 2);
    
    efs_hbo_contrast(1,i) = temp_mean_overall/temp_poolstd_overall;
    efs_hbo_contrast(2,i) = temp_mean_ses1/temp_poolstd_ses1;
    efs_hbo_contrast(3,i) = temp_mean_ses2/temp_poolstd_ses2;
    efs_hbo_contrast(4,i) = temp_mean_ses3/temp_poolstd_ses3;
    efs_hbo_contrast(5,i) = temp_mean_overall12/temp_poolstd_overall12;
end

clear temp_mean_overall temp_mean_ses1 temp_mean_ses2 temp_mean_ses3
clear temp_poolstd_overall temp_poolstd_ses1 temp_poolstd_ses2 temp_poolstd_ses3

efs_hbo_intertemporal = zeros(5,36);
efs_hbo_mindfulness = zeros(5,36); %overall,ses1,ses2,ses3,ses1+2
for i = 1:36
    efs_hbo_intertemporal(1,i) = hbo_intertemporal(1,i)/hbo_intertemporal(3,i);
    efs_hbo_intertemporal(2,i) = hbo_intertemporal_ses1(1,i)/hbo_intertemporal_ses1(3,i);
    efs_hbo_intertemporal(3,i) = hbo_intertemporal_ses2(1,i)/hbo_intertemporal_ses2(3,i);
    efs_hbo_intertemporal(4,i) = hbo_intertemporal_ses3(1,i)/hbo_intertemporal_ses3(3,i);
    efs_hbo_intertemporal(5,i) = hbo_intertemporal_ses12(1,i)/hbo_intertemporal_ses12(3,i);
    
    efs_hbo_mindfulness(1,i) = hbo_mindfulness(1,i)/hbo_mindfulness(3,i);
    efs_hbo_mindfulness(2,i) = hbo_mindfulness_ses1(1,i)/hbo_mindfulness_ses1(3,i);
    efs_hbo_mindfulness(3,i) = hbo_mindfulness_ses2(1,i)/hbo_mindfulness_ses2(3,i);
    efs_hbo_mindfulness(4,i) = hbo_mindfulness_ses3(1,i)/hbo_mindfulness_ses3(3,i);
    efs_hbo_mindfulness(5,i) = hbo_mindfulness_ses12(1,i)/hbo_mindfulness_ses12(3,i);
end
%% Bilateral activation difference comparisons (contrast)
Ses1_contrast_hbo = Ses1_intertemporal_hbo - Ses1_mindfulness_hbo;
Ses12_contrast_hbo = Total12_intertemporal_hbo - Total12_mindfulness_hbo;
Ses123_contrast_hbo = Total_intertemporal_hbo - Total_mindfulness_hbo;
left_ch = [16 17 11 12 1 6 7 18 22 29 35];
right_ch = [15 19 8 9 2 4 5 20 24 32 36];
label_ba = [9 9 10 10 11 11 46 45 21 21 17];
L_contrast_ses1 = [];
L_contrast_ses12 = [];
L_contrast_ses123 = [];
R_contrast_ses1 = [];
R_contrast_ses12 = [];
R_contrast_ses123 = [];

for i = left_ch
    temp_left = Ses1_contrast_hbo(:, i);
    L_contrast_ses1 = [L_contrast_ses1, temp_left];
    clear temp_left
    
    temp_left = Ses12_contrast_hbo(:, i);
    L_contrast_ses12 = [L_contrast_ses12, temp_left];
    clear temp_left
    
    temp_left = Ses123_contrast_hbo(:, i);
    L_contrast_ses123 = [L_contrast_ses123, temp_left];
    clear temp_left
end

for i = right_ch
    temp_right = Ses1_contrast_hbo(:, i);
    R_contrast_ses1 = [R_contrast_ses1, temp_right];
    clear temp_right
    
    temp_right = Ses12_contrast_hbo(:, i);
    R_contrast_ses12 = [R_contrast_ses12, temp_right];
    clear temp_right
    
    temp_right = Ses123_contrast_hbo(:, i);
    R_contrast_ses123 = [R_contrast_ses123, temp_right];
    clear temp_right
end

bi_contrast_ses1 = [label_ba; zeros(1,11)]; %col2:p; col3:t-value; col4: efs; 
bi_contrast_ses12 = [label_ba; zeros(1,11)];
bi_contrast_ses123 = [label_ba; zeros(1,11)];

for i = 1:size(label_ba,2)
    [h,p,ci,stats] = ttest(L_contrast_ses1(:,i), R_contrast_ses1(:,i)); %left - right; paired ttest
    temp_mean_df = mean(L_contrast_ses1(:,i)) - mean(R_contrast_ses1(:,i));
    std_left = std(L_contrast_ses1(:,i));
    std_right = std(R_contrast_ses1(:,i));
    temp_poolstd = sqrt((std_left^2 + std_right^2) / 2);
    bi_contrast_ses1(2,i) = p;
    bi_contrast_ses1(3,i) = stats.tstat;
    bi_contrast_ses1(4,i) = temp_mean_df / temp_poolstd; %left - right; Cohen's d
    clear h p ci stats temp_mean_df std_left std_right temp_poolstd
    
    [h,p,ci,stats] = ttest(L_contrast_ses12(:,i), R_contrast_ses12(:,i)); %left - right; paired ttest
    temp_mean_df = mean(L_contrast_ses12(:,i)) - mean(R_contrast_ses12(:,i));
    std_left = std(L_contrast_ses12(:,i));
    std_right = std(R_contrast_ses12(:,i));
    temp_poolstd = sqrt((std_left^2 + std_right^2) / 2);
    bi_contrast_ses12(2,i) = p;
    bi_contrast_ses12(3,i) = stats.tstat;
    bi_contrast_ses12(4,i) = temp_mean_df / temp_poolstd; %left - right; Cohen's d
    clear h p ci stats temp_mean_df std_left std_right temp_poolstd
    
    [h,p,ci,stats] = ttest(L_contrast_ses123(:,i), R_contrast_ses123(:,i)); %left - right; paired ttest
    temp_mean_df = mean(L_contrast_ses123(:,i)) - mean(R_contrast_ses123(:,i));
    std_left = std(L_contrast_ses123(:,i));
    std_right = std(R_contrast_ses123(:,i));
    temp_poolstd = sqrt((std_left^2 + std_right^2) / 2);
    bi_contrast_ses123(2,i) = p;
    bi_contrast_ses123(3,i) = stats.tstat;
    bi_contrast_ses123(4,i) = temp_mean_df / temp_poolstd; %left - right; Cohen's d
    clear h p ci stats temp_mean_df std_left std_right temp_poolstd
end
%% Figure 5 Plotting (Hemispheric activation contrast between two meditations)
cd ..\..\..
cd Figure5\
mediPeriodLabels = {'1','1 & 2','1 & 2 & 3'};
for i = [2 4 6 8 9 10] %the pairs of lateral channels with significant differences
    main_plot = figure('color', 'w');
    set(main_plot, 'Units', 'inches', 'Position', [0, 0, 2, 1.5]);
    x = 1:3;
    x_offsets = [-0.15, 0.15];
    temp_hbo_L = [mean(L_contrast_ses1(:,i)),mean(L_contrast_ses12(:,i)),mean(L_contrast_ses123(:,i))];
    temp_hbo_L_se = [std(L_contrast_ses1(:,i))/sqrt(length(L_contrast_ses1(:,i))),std(L_contrast_ses12(:,i))/sqrt(length(L_contrast_ses12(:,i))),std(L_contrast_ses123(:,i))/sqrt(length(L_contrast_ses123(:,i)))];
    temp_hbo_R = [mean(R_contrast_ses1(:,i)),mean(R_contrast_ses12(:,i)),mean(R_contrast_ses123(:,i))];
    temp_hbo_R_se = [std(R_contrast_ses1(:,i))/sqrt(length(R_contrast_ses1(:,i))),std(R_contrast_ses12(:,i))/sqrt(length(R_contrast_ses12(:,i))),std(R_contrast_ses123(:,i))/sqrt(length(R_contrast_ses123(:,i)))];
    bar(x + x_offsets(1), temp_hbo_L, 0.3, 'FaceColor', [0.7, 0.1, 0.1], 'EdgeColor', 'none', 'DisplayName', 'Left hemisphere');
    hold on;
    errorbar(x + x_offsets(1), temp_hbo_L, temp_hbo_L_se, 'Color', [0, 0, 0], 'linewidth', 0.3, 'LineStyle', 'none');
    bar(x + x_offsets(2), temp_hbo_R, 0.3, 'FaceColor', [0.3, 0.4, 0.9], 'EdgeColor', 'none', 'DisplayName', 'Right hemisphere');
    errorbar(x + x_offsets(2), temp_hbo_R, temp_hbo_R_se, 'Color', [0, 0, 0], 'linewidth', 0.3, 'LineStyle', 'none');
    title(['Bilateral pair',num2str(i)]);
    xlabel('Meditation session', 'Color', 'k');
    ylabel('HbO amplitude (uM)', 'Color', 'k');
    xlim([0 4]);
    xticks(1:3);
    xticklabels(mediPeriodLabels);
    set(gca,'fontsize',7,'FontName','Times New Roman', 'XColor', 'k', 'YColor', 'k');
    set(findall(gcf,'type','line'),'LineWidth',0.5, 'Color', 'k');
    set(get(gca,'xlabel'),'fontsize',7);
    set(get(gca,'ylabel'),'fontsize',7);
    set(get(gca,'title'),'fontsize',7);
    grid off
    box off
    print(main_plot,['Pair',num2str(i)],'-dpdf','-r600')
    close all
end