function [sigma, pres] = tool_hydrostatic_pressure_2D(stress, axis, varargin)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
    sigma=zeros(size(stress));
    pres=zeros(size(stress,1),1);
    if isempty(varargin)
        identity = [1,0,0,1];
    elseif varargin{1}==3
        identity = [1, 0, 1];
    end
    for i = 1:size(stress,1)
        pres(i,1)  = stress(i,axis);
        sigma(i,:) = stress(i,:) - pres(i,1)*identity;
    end
end

