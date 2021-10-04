function [succ] = writePLY( pos, normal, color, filename )
%WRITEPLY Summary of this function goes here
%   Detailed explanation goes here

% write ply file
[numPoints, dimPoints] = size(pos);
numF = 0;

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
fprintf(fid, 'property uchar red\n');
fprintf(fid, 'property uchar green\n');
fprintf(fid, 'property uchar blue\n');
fprintf(fid, 'property uchar alpha\n');
fprintf(fid, 'element face %d\n', numF);
fprintf(fid, 'property list uchar int vertex_indices\n');
fprintf(fid, 'end_header\n');

%fprintf(fid, '%0.6f %0.6f %0.6f\n', [ V1(:,1)'; V1(:,2)'; V1(:,3)']);
fprintf(fid, '%0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %d %d %d %d\n', [ pos(:,1)'; pos(:,2)'; pos(:,3)'; ...
                normal(:,1)'; normal(:,2)'; normal(:,3)';  color(:,1)';color(:,2)';color(:,3)';color(:,4)']);

% fprintf(fid, '%0.6f %0.6f %0.6f %0.6f %0.6f %0.6f\n', [ Point.pos(:,1)'; Point.pos(:,2)'; Point.pos(:,3)'; ...
%     Point.localstru(:,1)'; Point.localstru(:,2)'; Point.localstru(:,3)']);
%fprintf(fid, 'f %d %d %d\n', [F(:,1)'; F(:,2)'; F(:,3)';]);
fclose(fid);

succ = 1;

end

