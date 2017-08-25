% Plot angle phi vs. front speed s

% Zachary Singer, University of Minnesota Twin Cities, 8/23/17

% Set up window, interpreter, figure for plotting

clc
clf
clear all
set(gcf, 'Position', [100, 800, 800, 800])
set(0, 'defaulttextinterpreter', 'latex')

%%% AESTHETIC PLOTTING CONSTANTS %%%

fontsize = 36;
tick_fontsize = 20;
linewidth = 4;

% Read files from continuation_h1.m

cont_file_c = 'c_results_82317_01.txt';
cont_file_w = 'w_results_82317_01.txt';
cont_file_m = 'mat_results_82317_01.txt';
c_data = textread(cont_file_c, '', 'delimiter', ',', 'emptyvalue', NaN);
w_data = textread(cont_file_w, '', 'delimiter', ',', 'emptyvalue', NaN);
m_data = textread(cont_file_m, '', 'delimiter', ',', 'emptyvalue', NaN);


