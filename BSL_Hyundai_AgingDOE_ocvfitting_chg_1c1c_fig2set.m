clc;clear;close all


addpath = 'C:\Users\GSPARK\Desktop\Hcode';
data_folder = "G:\공유 드라이브\GSP_Data\1C1C";
cd(data_folder);
[save_folder,save_name] = fileparts(data_folder);


    n_points = 100; % 95 for coincells 400 for NE cell data
    w_ocv_scale = 1;
    w_dvdq_scale = 5;
    N_iter = 100;
    N_multistart = 24;
    t_pause_plot = 0;

Crate_Q = [56.2549, 54.9090, 54.3720, 53.7353, 52.9119];

% load OCP datas
    load ("G:\공유 드라이브\BSL-Data\Processed_data\Hyundai_dataset\Coin_cell\OCV\AHC_(5)_OCV_C20.mat")
    ocpn_raw = OCV_golden.OCVchg;
    x_raw = ocpn_raw(:,1);
    x = linspace(min(x_raw),max(x_raw),n_points)';
    ocpn_raw = ocpn_raw(:,2);
    ocpn_mva = movmean(ocpn_raw,round(length(ocpn_raw)/n_points));
    ocpn = interp1(x_raw,ocpn_mva,x);
    ocpn = [x ocpn];
    clear OCV_golden OCV_all Q_cell x_raw x ocpn_raw ocpn_mva

    load ("G:\공유 드라이브\BSL-Data\Processed_data\Hyundai_dataset\Coin_cell\OCV\CHC_(5)_OCV_C20.mat")
    ocpp_raw = OCV_golden.OCVchg;
    y_raw = ocpp_raw(:,1);
    y = linspace(min(y_raw),max(y_raw),n_points)';
    ocpp_raw = ocpp_raw(:,2);
    ocpp_mva = movmean(ocpp_raw,round(length(ocpp_raw)/n_points));
    ocpp = interp1(y_raw,ocpp_mva,y);
    ocpp = [y ocpp];
    clear OCV_golden OCV_all Q_cell y_raw y ocpp_raw ocpp_mva

    % folder = 'G:\공유 드라이브\BSL-Data\Processed_data\Hyundai_dataset\Coin_cell\OCV';
    % filename_ocpn = 'AHC_(5)_OCV_C20.mat';
    % filename_ocpp = 'CHC_(5)_OCV_C20.mat';
    % ocpn1 = load([folder filesep filename_ocpn]);
    % ocpp1 = load([folder filesep filename_ocpp]);
    % 
    % ocpn1 = ocpn1.OCV_golden.OCVchg;
    % ocpp1 = ocpp1.OCV_golden.OCVdis;
    % figure(1)
    % plot(ocpn1(:,1),ocpn1(:,2))
    % figure(2)
    % plot(ocpp1(:,1),ocpp1(:,2))

  %   SOC = -0.1:0.01:1.2;
  % 
  %  x_vec = x0 + (x1-x0)*SOC;
  %  y_vec = y0 + (y1-y0)*SOC;
  % 
  % ocpp_vec = interp1(ocpp(:,1),ocpp(:,2),y_vec);
  % ocpn_vec = interp1(ocpn(:,1),ocpn(:,2),x_vec);
  % 
  % ocv_vec = ocpp_vec - ocpn_vec;
  % 
  % figure(3)
  % yyaxis left
  % plot(SOC, ocv_vec,'k'); hold on
  % plot(SOC, ocpp_vec,'b')
  % yyaxis right
  % plot(SOC, ocpn_vec,'r'); 
  % xline(1)
  % xline(0)

% load Merged data
    % see BSL_hyundai_agingDOE_merge.m
    load("G:\공유 드라이브\BSL-Data\Processed_data\Hyundai_dataset\현대차파우치셀 (rOCV,Crate)\1C1C\NE_OCV_Merged.mat")
 

% check data
    % figure(1)
    % yyaxis left
    % plot(ocpn(:,1),ocpn(:,2)); hold on7
    % yyaxis right
    % plot(ocpp(:,1),ocpp(:,2)); hold on
% -----------------------------------
j_count = 0; % to count OCV steps
tic;

for i = 1:size(data_merged,1) % loop over steps

    data_merged(i).step = i; % add step number field


    % detect OCV step (**charging)
    if data_merged(i).OCVflag == 1

    j_count = j_count+1; % number of OCV now


    ocv_raw = data_merged(i).V; % OCV
    i_cc2cv = find(ocv_raw >= max(ocv_raw),1,'first');
    ocv_cc = ocv_raw(1:i_cc2cv);
    ocv_mva = movmean(ocv_cc,round(length(ocv_cc)/n_points));

    q_raw = data_merged(i).cumQ; % cum-capacity **charging
    q_cc = q_raw(1:i_cc2cv);
    q = linspace(min(q_cc),max(q_cc),n_points)';

    ocv = interp1(q_cc,ocv_mva,q);

    data_merged(i).q_ocv = [q ocv]; % add a field for later use


    %% Initial guess

            if j_count == 1 % first RPT
                Q_cell = abs(data_merged(i).Q); 
                x_guess = [0,Q_cell,1,Q_cell]; % x0, Qn, y0, Qp
                x_lb = [0,  Q_cell*0.5, 0.8,  Q_cell*0.5];
                x_ub = [0.2,Q_cell*2,   1,  Q_cell*2];

            else % non-first RPT

                % detect the first and the last OCV results
                i_first = find([data_merged(1:i-1).OCVflag] ==1,1,'first'); % firt OCV result -> for upper bounds
                i_last = find([data_merged(1:i-1).OCVflag] ==1,1,'last'); % last OCV result -> for initial guess

                Qp_first = data_merged(i_first).ocv_para_hat(4);
                Qn_first = data_merged(i_first).ocv_para_hat(2);

                Q_cell = abs(data_merged(i).Q); % 1**
                x_guess =   data_merged(i_last).ocv_para_hat; % initial guess from the last RPT
                x_lb =      [0,     Q_cell*0.5,          0.5,     Q_cell*0.5];
                x_ub =      [0.2,   Qn_first,    1,      Qp_first];

            end

    %% Weighting

        % OCV weighting
        w_ocv = w_ocv_scale*ones(size(q)); 

        % dvdq weighting
        w_dvdq = w_dvdq_scale*ones(size(q)); %



        %% Fitting 

        options = optimoptions(@fmincon,'MaxIterations',N_iter, 'StepTolerance',1e-10,'ConstraintTolerance', 1e-10, 'OptimalityTolerance', 1e-10);

        problem = createOptimProblem('fmincon', 'objective', @(x)func_ocvdvdq_cost(x,ocpn,ocpp,[q ocv],w_dvdq,w_ocv), ...
            'x0', x_guess, 'lb', x_lb, 'ub', x_ub, 'options', options);
        ms = MultiStart('Display','iter','UseParallel',true,'FunctionTolerance',1e-100,'XTolerance',1e-100);

        [x_hat, f_val, exitflag, output] = run(ms,problem,N_multistart);
        [cost_hat, ocv_hat, dvdq_mov, dvdq_sim_mov] = func_ocvdvdq_cost(x_hat,ocpn,ocpp,[q ocv],w_dvdq,w_ocv);


%         figure()
%         plot(soc,ocv); hold on
%         plot(soc,ocv_hat); hold on
% 
%         figure()
%         plot(soc,dvdq_mov); hold on
%         plot(soc,dvdq_sim_mov); hold on
        % Set Y lim
%         ylim_top = 2*max(dvdq_mov((soc > 0.2) & (soc < 0.8)));
%         ylim([0 ylim_top])

        % save the result to the struct
        data_merged(i).ocv_para_hat = x_hat;
        data_merged(i).ocv_hat = ocv_hat;
        data_merged(i).dvdq_mov = dvdq_mov;
        data_merged(i).dvdq_sim_mov = dvdq_sim_mov;


    end
end
toc


%% LOOP over OCV Steps 

data_ocv = data_merged([data_merged.OCVflag]' == 1);
J = size(data_ocv,1);
c_mat = lines(J);

figure_handle = figure();
subplot_idx = 0;

for j = [1,3,5]
    subplot_idx = subplot_idx + 1;

    % Plot OCV fitting results: OCV and dVdQ plots
    if data_ocv(j).q_ocv(1,2) < data_ocv(j).q_ocv(end,2) % charging ocv
        soc_now =  data_ocv(j).q_ocv(:,1)/data_ocv(j).q_ocv(end,1);
    else
        soc_now =  (data_ocv(j).q_ocv(end,1)-data_ocv(j).q_ocv(:,1))/data_ocv(j).q_ocv(end,1);
    end
    ocv_now = data_ocv(j).q_ocv(:,2);
    ocv_sim_now = data_ocv(j).ocv_hat;
    dvdq_now = data_ocv(j).dvdq_mov;
    dvdq_sim_now = data_ocv(j).dvdq_sim_mov;


%   figure(3)
%     SOC = -0.2:0.01:1.2;
% 
%     x_vec = data_ocv(1).ocv_para_hat(1) + ((data_ocv(1).ocv_para_hat(1)+(Q_cell/data_ocv(1).ocv_para_hat(2)*SOC)));
%     y_vec = data_ocv(1).ocv_para_hat(3) + ((data_ocv(1).ocv_para_hat(3)+(Q_cell/data_ocv(1).ocv_para_hat(4)*SOC)));
% 
% ocpp1_vec = interp1(ocpp1(:,1), ocpp1(:,2), y_vec, 'spline', 'extrap');
% ocpn1_vec = interp1(ocpn1(:,1), ocpn1(:,2), x_vec, 'spline', 'extrap');
% 
% 
%     ocv1_vec = ocpp1_vec - ocpn1_vec;
% 
%     figure(3)
%     % yyaxis left
%     plot(SOC, ocv1_vec,'-','Color','k'); hold on;
%     plot(SOC, ocpp1_vec,'o')
%     yyaxis right
%     plot(SOC, ocpn1_vec,'o'); 
%     xline(1)
%     xline(0)




    % set(gcf,'position',[100,100,1600,800])

   subplot(2,3,subplot_idx);

    plot(soc_now,ocv_now,'color',[0.455, 0.463, 0.471]); hold on
    plot(soc_now,ocv_sim_now(:,3),'color',[0.027, 0.451, 0.761]);
    plot(soc_now,ocv_sim_now(:,1),'--','color',[0.937, 0.753, 0.000]);
    plot(soc_now,ocv_sim_now(:,2),'-r');

    ylabel('PE, FC Voltage [V]');
    xticks(0:0.2:1);
    yticks(3:0.2:4.4);
    ylim([3 4.4]);

    yyaxis right
    % xlabel('SOC');  % x축 레이블
    ylabel('NE Voltage [V]');  % y축 레이블
    plot(soc_now,ocv_sim_now(:,2),'color',[ 0.804, 0.325, 0.298]);
    h =legend('FC data','PE fit','FC fit', 'NE fit', 'Location', 'north','Box','off','NumColumns', 2);
    h.ItemTokenSize(1) = 15;
    h.FontSize = 4;


     yticks(0:0.1:0.6);
     ylim([0 0.6]);
     set(gca, 'XColor', 'k', 'YColor', 'k');
     % figuresettings8('0cyc_OCV', 1200);


    subplot(2,3,subplot_idx+3);

    plot(soc_now,dvdq_now(:,1),'color',[0.455, 0.463, 0.471]); hold on
    plot(soc_now,dvdq_sim_now(:,3),'color',[0.027, 0.451, 0.761]);
    plot(soc_now,dvdq_sim_now(:,1),'--','color',[0.937, 0.753, 0.000]);
    plot(soc_now,-dvdq_sim_now(:,2),'color',[ 0.804, 0.325, 0.298]); hold on;

    xlabel('SOC');  % x축 레이블
    ylabel('dV/dQ [V/mAh]\times10');
    h =legend('FC data','PE fit','FC fit', 'NE fit', 'Location', 'north','Box','off','NumColumns', 2);
    h.ItemTokenSize(1) = 15;
    h.FontSize = 4;



        % Set Y lim
        ylim_top = 2*max(dvdq_now((soc_now > 0.2) & (soc_now < 0.8)));
        ylim([0 ylim_top]);

    ylim([0 0.035]);
    xticks(0:0.2:1);

    ylim([0 0.04]);
    yticks(0:0.005:0.035);
    xticks(0:0.2:1);

    yticklabels(arrayfun(@(y) sprintf('%.2f', y*10), 0:0.005:0.035, 'UniformOutput', false));
    set(gca, 'XColor', 'k', 'YColor', 'k');
end

    filenames = sprintf('1c1cset');
    figuresettingsSET(filenames, 1200);




% ------------------------------------------------------fitting2
 % Crate_Q = [56.5790, 55.9170 55.5610 55.0090 54.1250 51.4230 25.1460];


    % folder = 'G:\공유 드라이브\BSL-Data\Processed_data\Hyundai_dataset\Coin_cell\OCV';
    % filename_ocpn = 'AHC_(5)_OCV_C20.mat';
    % filename_ocpp = 'CHC_(5)_OCV_C20.mat';
    % ocpn1 = load([folder filesep filename_ocpn]);
    % ocpp1 = load([folder filesep filename_ocpp]);
    % 
    % ocpn1 = ocpn1.OCV_golden.OCVchg;
    % ocpp1 = ocpp1.OCV_golden.OCVdis;
    % figure(1)
    % plot(ocpn1(:,1),ocpn1(:,2))
    % figure(2)
    % plot(ocpp1(:,1),ocpp1(:,2))

  %   SOC = -0.1:0.01:1.2;
  % 
  %  x_vec = x0 + (x1-x0)*SOC;
  %  y_vec = y0 + (y1-y0)*SOC;
  % 
  % ocpp_vec = interp1(ocpp(:,1),ocpp(:,2),y_vec);
  % ocpn_vec = interp1(ocpn(:,1),ocpn(:,2),x_vec);
  % 
  % ocv_vec = ocpp_vec - ocpn_vec;
  % 
  % figure(3)
  % yyaxis left
  % plot(SOC, ocv_vec,'k'); hold on
  % plot(SOC, ocpp_vec,'b')
  % yyaxis right
  % plot(SOC, ocpn_vec,'r'); 
  % xline(1)
  % xline(0)

% load Merged data
    % see BSL_hyundai_agingDOE_merge.m
%     load("C:\Users\GSPARK\Desktop\QC1C\QC1C사이클(C10)\NE_OCV_Merged.mat")
%     %load('G:\Shared drives\BSL_Data2\HNE_AgingDOE_Processed\HNE_FCC\4CPD 1C (25-42)\25degC\HNE_FCC_4CPD 1C (25-42)_25degC_s01_3_6_Merged.mat')
% j_count = 0; % to count OCV steps
% tic;
% 
% for i = 1:size(data_merged,1) % loop over steps
% 
%     data_merged(i).step = i; % add step number field
% 
% 
%     % detect OCV step (**charging)
%     if data_merged(i).OCVflag == 1
% 
%     j_count = j_count+1; % number of OCV now
% 
% 
%     ocv_raw = data_merged(i).V; % OCV
%     i_cc2cv = find(ocv_raw >= max(ocv_raw),1,'first');
%     ocv_cc = ocv_raw(1:i_cc2cv);
%     ocv_mva = movmean(ocv_cc,round(length(ocv_cc)/n_points));
% 
%     q_raw = data_merged(i).cumQ; % cum-capacity **charging
%     q_cc = q_raw(1:i_cc2cv);
%     q = linspace(min(q_cc),max(q_cc),n_points)';
% 
%     ocv = interp1(q_cc,ocv_mva,q);
% 
%     data_merged(i).q_ocv = [q ocv]; % add a field for later use
% 
% 
%     %% Initial guess
% 
%             if j_count == 1 % first RPT
%                 Q_cell = abs(data_merged(i).Q); 
%                 x_guess = [0,Q_cell,1,Q_cell]; % x0, Qn, y0, Qp
%                 x_lb = [0,  Q_cell*0.5, 0.8,  Q_cell*0.5];
%                 x_ub = [0.2,Q_cell*2,   1,  Q_cell*2];
% 
%             else % non-first RPT
% 
%                 % detect the first and the last OCV results
%                 i_first = find([data_merged(1:i-1).OCVflag] ==1,1,'first'); % firt OCV result -> for upper bounds
%                 i_last = find([data_merged(1:i-1).OCVflag] ==1,1,'last'); % last OCV result -> for initial guess
% 
%                 Qp_first = data_merged(i_first).ocv_para_hat(4);
%                 Qn_first = data_merged(i_first).ocv_para_hat(2);
% 
%                 Q_cell = abs(data_merged(i).Q); % 1**
%                 x_guess =   data_merged(i_last).ocv_para_hat; % initial guess from the last RPT
%                 x_lb =      [0,     Q_cell*0.5,          0.5,     Q_cell*0.5];
%                 x_ub =      [0.2,   Qn_first,    1,      Qp_first];
% 
%             end
% 
%     %% Weighting
% 
%         % OCV weighting
%         w_ocv = w_ocv_scale*ones(size(q)); 
% 
%         % dvdq weighting
%         w_dvdq = w_dvdq_scale*ones(size(q)); %
% 
% 
%         %% Fitting 
% 
%         options = optimoptions(@fmincon,'MaxIterations',N_iter, 'StepTolerance',1e-10,'ConstraintTolerance', 1e-10, 'OptimalityTolerance', 1e-10);
% 
%         problem = createOptimProblem('fmincon', 'objective', @(x)func_ocvdvdq_cost(x,ocpn,ocpp,[q ocv],w_dvdq,w_ocv), ...
%             'x0', x_guess, 'lb', x_lb, 'ub', x_ub, 'options', options);
%         ms = MultiStart('Display','iter','UseParallel',true,'FunctionTolerance',1e-100,'XTolerance',1e-100);
% 
%         [x_hat, f_val, exitflag, output] = run(ms,problem,N_multistart);
%         [cost_hat, ocv_hat, dvdq_mov, dvdq_sim_mov] = func_ocvdvdq_cost(x_hat,ocpn,ocpp,[q ocv],w_dvdq,w_ocv);
% 
% 
% %         figure()
% %         plot(soc,ocv); hold on
% %         plot(soc,ocv_hat); hold on
% % 
% %         figure()
% %         plot(soc,dvdq_mov); hold on
% %         plot(soc,dvdq_sim_mov); hold on
%         % Set Y lim
% %         ylim_top = 2*max(dvdq_mov((soc > 0.2) & (soc < 0.8)));
% %         ylim([0 ylim_top])
% 
%         % save the result to the struct
%         data_merged(i).ocv_para_hat = x_hat;
%         data_merged(i).ocv_hat = ocv_hat;
%         data_merged(i).dvdq_mov = dvdq_mov;
%         data_merged(i).dvdq_sim_mov = dvdq_sim_mov;
% 
% 
%     end
% end
% toc


%% LOOP over OCV Steps 

% data_ocv = data_merged([data_merged.OCVflag]' == 1);
% J = size(data_ocv,1);
% c_mat = lines(J);
% figure_handle = figure();
% subplot_idx = 0;
% 
% for j = [1 4 6]
% 
%     subplot_idx = subplot_idx + 1;
%     % Plot OCV fitting results: OCV and dVdQ plots
%     if data_ocv(j).q_ocv(1,2) < data_ocv(j).q_ocv(end,2) % charging ocv
%         soc_now =  data_ocv(j).q_ocv(:,1)/data_ocv(j).q_ocv(end,1);
%     else
%         soc_now =  (data_ocv(j).q_ocv(end,1)-data_ocv(j).q_ocv(:,1))/data_ocv(j).q_ocv(end,1);
%     end
%     ocv_now = data_ocv(j).q_ocv(:,2);
%     ocv_sim_now = data_ocv(j).ocv_hat;
%     dvdq_now = data_ocv(j).dvdq_mov;
%     dvdq_sim_now = data_ocv(j).dvdq_sim_mov;


%   figure(3)
%     SOC = -0.2:0.01:1.2;
% 
%     x_vec = data_ocv(1).ocv_para_hat(1) + ((data_ocv(1).ocv_para_hat(1)+(Q_cell/data_ocv(1).ocv_para_hat(2)*SOC)));
%     y_vec = data_ocv(1).ocv_para_hat(3) + ((data_ocv(1).ocv_para_hat(3)+(Q_cell/data_ocv(1).ocv_para_hat(4)*SOC)));
% 
% ocpp1_vec = interp1(ocpp1(:,1), ocpp1(:,2), y_vec, 'spline', 'extrap');
% ocpn1_vec = interp1(ocpn1(:,1), ocpn1(:,2), x_vec, 'spline', 'extrap');
% 
% 
%     ocv1_vec = ocpp1_vec - ocpn1_vec;
% 
%     figure(3)
%     % yyaxis left
%     plot(SOC, ocv1_vec,'-','Color','k'); hold on;
%     plot(SOC, ocpp1_vec,'o')
%     yyaxis right
%     plot(SOC, ocpn1_vec,'o'); 
%     xline(1)
%     xline(0)

    % set(gcf,'position',[100,100,1600,800])

%     subplot(2,3,subplot_idx);
%     plot(soc_now,ocv_now,'color',[0.455, 0.463, 0.471]); hold on
%     plot(soc_now,ocv_sim_now(:,3),'color',[0.027, 0.451, 0.761]);
%     plot(soc_now,ocv_sim_now(:,1),'--','color',[0.937, 0.753, 0.000]);
%     plot(soc_now,ocv_sim_now(:,2),'-r');
% 
%     ylabel('PE, FC Voltage [V]');
%     xticks(0:0.2:1);
%     yticks(3:0.2:4.4);
%     ylim([3 4.4]);
% 
%     yyaxis right
%     xlabel('SOC');  
%     ylabel('NE Voltage [V]');  
%     plot(soc_now,ocv_sim_now(:,2),'color',[ 0.804, 0.325, 0.298]);
%     h =legend('FC data','PE fit','FC fit', 'NE fit', 'Location', 'north','Box','off','NumColumns', 2);
%     h.ItemTokenSize(1) = 15;
%     h.FontSize = 4;
%      yticks(0:0.1:0.6);
%      ylim([0 0.6]);
%      set(gca, 'XColor', 'k', 'YColor', 'k');
%      % figuresettings8('0cyc_OCV', 1200);
% 
%     subplot(2,3,subplot_idx+3);
%     plot(soc_now,dvdq_now(:,1),'color',[0.455, 0.463, 0.471]); hold on
%     plot(soc_now,dvdq_sim_now(:,3),'color',[0.027, 0.451, 0.761]);
%     plot(soc_now,dvdq_sim_now(:,1),'--','color',[0.937, 0.753, 0.000]);
%     plot(soc_now,-dvdq_sim_now(:,2),'color',[ 0.804, 0.325, 0.298]); hold off;
% 
%     xlabel('SOC');  % x축 레이블
%     ylabel('dV/dQ [V/mAh] \times10');
%     h =legend('FC data','PE fit','FC fit', 'NE fit', 'Location', 'north','Box','off','NumColumns', 2);
%     h.ItemTokenSize(1) = 15;
%     h.FontSize = 4;
% 
%         % Set Y lim
%         ylim_top = 2*max(dvdq_now((soc_now > 0.2) & (soc_now < 0.8)));
%         ylim([0 ylim_top]);
%     % 
%     ylim([0 0.04]);
%     yticks(0:0.005:0.035);
%     xticks(0:0.2:1);
% 
%     yticklabels(arrayfun(@(y) sprintf('%.2f', y*10), 0:0.005:0.035, 'UniformOutput', false));
%     set(gca, 'XColor', 'k', 'YColor', 'k');
% 
% end
%     filenames = sprintf('qc1cset');
%     figuresettingsSET(filenames, 1200);


    % filenames = sprintf('fig2_%d', j);
    % figuresettings14(filenames, 1200);

   % %  % fig1 = sprintf('G:\\공유 드라이브\\GSP_Data\\1C1C\\%d.fig', j);
   % %  % savefig(gcf, fig1);
   % % % ax = gca;
   % % % ax.FontWeight = 'bold';
   
   % pause(t_pause_plot)

   % save_fullpath = [save_folder filesep save_name '.mat'];
   % save(save_fullpath)
   % OCV_golden = OCV_all;
   % save(save_fullpath,'OCV_golden')


%     %% LLI, LAMp
%     data_ocv(j).LAMp = data_ocv(1).ocv_para_hat(4)...
%                         -data_ocv(j).ocv_para_hat(4);
% 
%     data_ocv(j).LAMn = data_ocv(1).ocv_para_hat(2)...
%                         -data_ocv(j).ocv_para_hat(2); 
% 
%     data_ocv(j).LLI = (data_ocv(1).ocv_para_hat(4)*data_ocv(1).ocv_para_hat(3) + data_ocv(1).ocv_para_hat(2)*data_ocv(1).ocv_para_hat(1))...
%                         -(data_ocv(j).ocv_para_hat(4)*data_ocv(j).ocv_para_hat(3) + data_ocv(j).ocv_para_hat(2)*data_ocv(j).ocv_para_hat(1)); 
% 
%     data_ocv(j).dQ_LLI = (data_ocv(1).ocv_para_hat(4)*(data_ocv(1).ocv_para_hat(3)-1))...
%                         -(data_ocv(j).ocv_para_hat(4)*(data_ocv(j).ocv_para_hat(3)-1));
% 
% 
%     data_ocv(j).dQ_LAMp = (data_ocv(1).ocv_para_hat(4) - data_ocv(j).ocv_para_hat(4))...
%                             *(1-data_ocv(1).ocv_para_hat(3)+data_ocv(1).Q/data_ocv(1).ocv_para_hat(4));
% 
% 
%     data_ocv(j).dQ_data = data_ocv(1).Q - data_ocv(j).Q;
% 
%     dQ_data_now = data_ocv(j).dQ_data;
% 
%     dQ_total_now = data_ocv(j).dQ_LLI + data_ocv(j).dQ_LAMp;
% 
%     data_ocv(j).R = [abs(data_merged(4).Q)] - [Crate_Q(j)] - [data_ocv(j).dQ_data];
% 
%     data_ocv(j).Crate_Q = [Crate_Q(j)];
% 
%     % manipulate loss scale to be consistent with the data Q
%     scale_now = dQ_data_now/dQ_total_now;
%     data_ocv(j).dQ_LLI = data_ocv(j).dQ_LLI*scale_now;
%     data_ocv(j).dQ_LAMp = data_ocv(j).dQ_LAMp*scale_now;
% 
% end


% data_ocv(1).dQ_LLI = 0;
% data_ocv(1).dQ_LAMp = 0;
% 
% figure()
% bar([data_ocv.cycle],[[data_ocv.dQ_LLI];[data_ocv.dQ_LAMp];[abs(data_merged(4).Q)] - [Crate_Q] - [data_ocv.dQ_data]],'stacked')
% hold on
% plot([data_ocv.cycle],[data_ocv.dQ_data], '-sc', 'LineWidth', 2);
% plot([data_ocv.cycle],[data_ocv.dQ_data] + [abs(data_merged(4).Q)] - [Crate_Q]- [data_ocv.dQ_data],'-sm', 'LineWidth', 2)
% legend({'Loss by LLI', 'Loss by LAMp', 'Loss by R', 'Loss data (C/10)', 'Loss data (C/3)'},'Location', 'northwest'); 
% xlabel('Cycle')
% ylabel('$\Delta Q$ [Ah]', 'Interpreter', 'latex','FontName', 'Arial','FontWeight', 'bold')
% ylim([0 4.5])
% title('완속층전(1C1C) 대형셀 열화 인자');
% fig2 = sprintf('G:\\공유 드라이브\\BSL-Data\\Processed_data\\Hyundai_dataset\\현대차파우치셀 (rOCV,Crate)\\1C1C\\bar_plot.fig');
% saveas(gcf, fig2);
% 
% save_path = 'G:\공유 드라이브\BSL-Data\Processed_data\Hyundai_dataset\현대차파우치셀 (rOCV,Crate)\1C1C\NE_data_ocv';
% data_merged = data_ocv; 
% save(save_path,'data_merged');
% xlabel('Cycle')
% ylabel('$\Delta Q$ [Ah]', 'Interpreter', 'latex','FontName', 'Arial','FontWeight', 'bold')
% 
% save_folder = data_folder; 
% 
% % 파일 전체 경로
% fig_path = fullfile(save_folder, fig_name);
% 
% % Figure 저장
% savefig(fig_path)
% print(fig_path, '-dpdf', '-r2400')
