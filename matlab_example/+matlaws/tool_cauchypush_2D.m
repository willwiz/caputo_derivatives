function sigma = tool_cauchypush_2D(stress,strain)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
    sigma=zeros(size(stress));
    for i = 1:size(stress,1)
        sigma(i,:) = strain(i,:).*stress(i,:).*strain(i,:);
    end
end

