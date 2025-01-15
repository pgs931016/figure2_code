function figuresettingsSET(filename, dpi, width, height)

    alw = 0.5;    
    fsz = 8; 
    lw = 1; 

    if nargin < 4, height = 8; end
    if nargin < 3, width = 17; end
    if nargin < 2, dpi = 1200; end
    if nargin < 1, filename = ['unnamed_', datestr(datetime, 'yyyymmdd_HHMMSS')]; end

    fig = gcf; 
    fig.Units = 'centimeters';
    fig.Position = [0, 0, width, height];

    % Subplot 핸들 생성
    axHandles = gobjects(2, 3); 
    for col = 1:3
        % 첫 번째 행
        axHandles(1, col) = subplot(2, 3, col);
    
        % 두 번째 행
        axHandles(2, col) = subplot(2, 3, col + 3);
    
    end

    % 각 subplot 설정
    for col = 1:3
        % 첫 번째 행
        set(axHandles(1, col), 'LineWidth', lw, 'FontSize', fsz);
    
        set(axHandles(1, col), 'XTickLabel', []); 

        % 두 번째 행
        set(axHandles(2, col), 'LineWidth', lw, 'FontSize', fsz);
    end

   
  for col = 1:3
        linkaxes([axHandles(1, col), axHandles(2, col)], 'x');
  end


    % Subplot 간 간격 조정 (행 간격 제거)
    for row = 1:2
        for col = 1:3
            xPos = 0.06 + (col - 1) * 0.34; 
            height = 0.43;
             if row == 1
             yPos = 0.53; 
             else
             yPos = 0.53-height ; 
             end
            axHandles(row, col).Position = [xPos, yPos, 0.2, height];
       
        end
    end


    % Figure 저장 설정
    set(gcf, 'InvertHardcopy', 'on');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperSize', [width, height]);
    % set(gcf, 'PaperPosition', [0, 0, width, height]);
    set(gcf, 'PaperPositionMode', 'auto');

    % 파일 저장
    savefig(gcf, [filename, '.fig']);
    print(gcf, [filename, '.tiff'], '-dtiff', ['-r', num2str(dpi)]);
end



% print('figure','-dpdf','-r2400');  
% 
% date = datetime('now','Format','yyyy-mm-dd');


 % idx = 1;
 %    while true
 %        filename = sprintf('%s_%s_%d', customname, date, idx);
 %        if ~isfile([filename '.pdf']) 
 %            break;
 %        end
 %        idx = idx + 1; 
 %    end
 % 
 %    % Save the figure as a PDF
 %    print(filename, '-dpdf', '-r2400');  


% if ispc % Use Windows ghostscript call
%   system('gswin64c -o -q -sDEVICE=png256 -dEPSCrop -r2400 -oimprovedExample_eps.png improvedExample.eps');
% else % Use Unix/OSX ghostscript call
%   system('gs -o -q -sDEVICE=png256 -dEPSCrop -r2400 -oimprovedExample_eps.png improvedExample.eps');
% end
% end
