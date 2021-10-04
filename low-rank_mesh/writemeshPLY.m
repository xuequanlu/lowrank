function [succ] = writemeshPLY( pos, F, Normal, filename )
%WRITEPLY Summary of this function goes here
%   Detailed explanation goes here

% write ply file
[numPoints, dimPoints] = size(pos);
[numF, dim] = size(F);

fid = fopen(filename, 'w');
fprintf(fid, 'ply\n');
fprintf(fid, 'format ascii 1.0\n');
fprintf(fid, 'comment VCGLIB generated\n');
fprintf(fid, 'element vertex %d\n', numPoints);
fprintf(fid, 'property float x\n');
fprintf(fid, 'property float y\n');
fprintf(fid, 'property float z\n');
fprintf(fid, 'property float nx\n');
fprintf(fid, 'property float ny\n');
fprintf(fid, 'property float nz\n');
fprintf(fid, 'element face %d\n', numF);
fprintf(fid, 'property list uchar int vertex_indices\n');
fprintf(fid, 'end_header\n');

%fprintf(fid, '%0.6f %0.6f %0.6f\n', [ V1(:,1)'; V1(:,2)'; V1(:,3)']);
fprintf(fid, '%0.6f %0.6f %0.6f %d %d %d\n', [ pos(:,1)';pos(:,2)';pos(:,3)';... 
    Normal(:,1)';Normal(:,2)';Normal(:,3)']);

fprintf(fid, '3 %d %d %d\n', [F(:,1)'; F(:,2)'; F(:,3)';]);
fclose(fid);

succ = 1;

end

