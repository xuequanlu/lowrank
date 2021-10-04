function [ V,Vn ] = readPLY2( fileName )
%READPLY2 此处显示有关此函数的摘要
%   此处显示详细说明

data = readPLY(fileName,'');
V = [data.vertex.x data.vertex.y data.vertex.z];
Vn = [data.vertex.nx data.vertex.ny data.vertex.nz];

end

