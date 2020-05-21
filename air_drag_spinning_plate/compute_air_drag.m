function [drag_force_x,drag_force_y]=compute_air_drag(theta,velocity)

%the purpose of this script is to model the drag force on a metal plate
%which can change its orientation relative to the direction of air passing
%over it. Here we will have the direction of wind blowing upwards (+Y)
%direction and the plate angle will change relative to the vertical (+Y)
%axis as well. Note this is a very simple model, desinged for a
%reinforcement learning agent testing. It can be improved i'm sure. 

%% equation for force exerted
% here is the equation 

%begin by declaring our constants
coefficient_of_drag=0.8; %this can change
density_of_air=1.225;
area_of_plate=0.3; %[m^2]

air_force=(0.5)*coefficient_of_drag*density_of_air*area_of_plate*sind(theta)*velocity.^2;

%% components of this force

drag_force_y=(sind(theta))*air_force; %one which will slow vertical 
drag_force_x=cosd(theta)*air_force; %one which will induce the rotation 

