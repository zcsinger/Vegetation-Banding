% Plot angle phi vs. front speed s in two dimensional planar system

% Zachary Singer, Univerisy of Minnesota - Twin Cities, 8/22/17

%%%%%%%%%%%%%%%%%%%% Specifications and Guide %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Graphs angle vs. front speed s in two dimensional planar system. 

% Given data (s, theta') want to transform to (phi, s):

% We undo scaling as follows:

% theta' = sigma * theta = 1 / sqrt(s + c*cos(phi)) * (s + c*cos(phi)) w_+

%%%%%%%%%%%%%%%%%%% Model and Parameters: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ansatz in 2D : 
% (b, w)( cos(phi)*x + sin(phi)*y - s*t)

% sigma^2 = 1 / (s + c cos(phi))


%%%%%%%%%%%%%%%%%%% Initializing Constants: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clf
clear all
set(gcf, 'Position', [100, 800, 800, 800])
set(0, 'defaulttextinterpreter', 'latex')

mode = input('Mode (1 = c only, 2 = w_plus only, 3 = both): ');

fontsize = 24;
tick_fontsize = 20;

% Upper edge heteroclinic (het1) data: (s, theta) values from continuation

cont_file = 'st_h1m_81717_std_01.txt';
data_ct = textread(cont_file, '', 'delimiter', ',', 'emptyvalue', NaN);
s_res = data_ct(1, :); % s data
t_res = data_ct(2, :); % theta data

% Important Constants

c_start = 2;
d_c = .1;
c_stop = 10;
w_plus_start = 6;
d_w = .1;
w_plus_stop = 10;

phi_val = @(c, w_plus) acos( ((t_res / w_plus).^2 - s_res) / c);

c = c_start;
w_plus = w_plus_start;

c_range = size(c_start : d_c : c_stop);
c_dirnum = c_range(2);
cm = colormap(jet(c_dirnum));

w_range = size(w_plus_start : d_w : w_plus_stop);
w_dirnum = w_range(2);
wm = colormap(jet(w_dirnum));

hold on

switch mode 
    
    case 1
        
    %%%%% JUST CHANGE C WITH FIXED W_PLUS %%%%%
    
    % w_plus > 6?
    
    for i = 1 : c_dirnum

        phi = phi_val (c, w_plus);

        plot(phi, s_res, 'Color', cm(i, :), 'LineWidth', 2)
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', tick_fontsize)

        ylabel('Front speed $$s$$', 'FontSize', fontsize)
        xlabel('Angle $$\phi$$', 'FontSize', fontsize)
        title(['Angle $$\phi$$ vs. Front speed $$s$$' ...
                ' with ($$w_+$$ fixed, $$c$$ varied)'], ...
                'FontSize', fontsize)
            
        drawnow

        c = c + d_c;

    end
    
    case 2
        
    %%%%% JUST CHANGE W_PLUS WITH FIXED C %%%%%
    
    % c > 3 ?
    
    for j = 1 : w_dirnum

        phi = phi_val (c, w_plus);

        plot(phi, s_res, 'Color', wm(j, :), 'LineWidth', 2)
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', tick_fontsize)

        ylabel('Front speed $$s$$', 'FontSize', fontsize)
        xlabel('Angle $$\phi$$', 'FontSize', fontsize)
        title(['Angle $$\phi$$ vs. Front speed $$s$$' ...
                ' with ($$c$$ fixed, $$w_+$$ varied)'], ...
                'FontSize', fontsize)

        drawnow

        w_plus = w_plus + d_w;

    end
    
    otherwise
        
    warning(['Unexpected mode: Type 1 for c only, 2 for w_+ only, 3 for both']);
        
end

