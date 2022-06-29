%%% Function created by Nawara T. (Mattheyses lab - 05/20/2022) compatible
%%% with MATLAB R2020b / Function calculates gamma correction factor
%%% based on microscope setup x would need to be adjusted

%This work is licensed under the Creative Commons Attribution 4.0
%International License. To view a copy of this license, visit
%http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative
%Commons, PO Box 1866, Mountain View, CA 94042, USA.

function [gamma] = gamma_calculator(y, yourpath)
lambda_short = 488;     %shorter wavlegth used for imaging
lambda_long = 647;      %longer wavlegth used for imaging
n0 = 1;                 %RI of glass
n1 = 1.512;             %RI of prism
n2 = 1.37;              %RI of cell
x = 114;                %distnas from objective to the wall

z = sqrt((x^2) + (y^2));
sin_alpha = sin(y/z);

theta = 90 - rad2deg(sin((sin_alpha * n0) / n1));

d1 = lambda_short / ((4 * pi) * (sqrt((n1^2) * (sin(deg2rad(theta))^2) - ((n2)^2))));
d2 = lambda_long / ((4 * pi) * (sqrt((n1^2) * (sin(deg2rad(theta))^2) - ((n2)^2))));

gamma = (d2-d1)/(d2*d1);

end