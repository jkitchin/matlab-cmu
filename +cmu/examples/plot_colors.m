%% using rgbcolors in plots
% John Kitchin
%
% The standard colors available in plotting are not that pretty:           
% b     blue         
% g     green        
% r     red          
% c     cyan         
% m     magenta      
% y     yellow       
% k     black         
% w     white 
%
% yellow is particularly hard to see.

clear all; close all; clc;

x = linspace(0,1);
y = sin(x);

figure
plot(x,y,'y-')

%%
plot(x,y,'c-')

%% using cmu.colors
% the cmu.colors provides a function to get a lot of colors
% deep carrot orange is much easier to see than the default yellow. Note
% this feature is only available for +cmu.version >= 1.7.
c = @cmu.colors;
figure
plot(x,y,'color',c('deep carrot orange'),'linewidth',3)

%% List of available colors
c('list')