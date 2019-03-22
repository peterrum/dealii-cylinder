function [] = trans(P3,name)

%% Reference element (not refined)
P1 = [
 0.0 0.0 0.0;  
 1.0 0.0 0.0;
 0.0 1.0 0.0;
 1.0 1.0 0.0;
 0.0 0.0 1.0;  
 1.0 0.0 1.0;
 0.0 1.0 1.0;
 1.0 1.0 1.0
];

%% Reference element (refined)
P2 = [
 % layer 0
 0.00 0.00 0.00; %  0
 0.50 0.00 0.00; %  1
 1.00 0.00 0.00; %  2
 0.00 0.50 0.00; %  3  
 0.50 0.50 0.00; %  4
 1.00 0.50 0.00; %  5
 0.00 1.00 0.00; %  6  
 0.50 1.00 0.00; %  7
 1.00 1.00 0.00; %  8
 
 % layer 1
 0.50 0.00 0.25; %  9
 1.00 0.00 0.25; % 10
 0.50 0.50 0.25; % 11
 1.00 0.50 0.25; % 12
 
 % layer 2
 0.00 0.00 0.50; % 13
 0.50 0.00 0.50; % 14 
 1.00 0.00 0.50; % 15
 0.00 0.50 0.50; % 16
 0.50 1.00 0.50; % 17
 1.00 1.00 0.50; % 18
 
 % layer 3
 0.00 0.00 1.00; % 19
 1.00 0.00 1.00; % 20
 0.00 1.00 1.00; % 21
 1.00 1.00 1.00  % 22
];

%% Create Interpolation matrix
N = zeros(size(P1,1), size(P2,1));
dim = size(P1,2);
for i=1:size(P1,1)
    for j=1:size(P2,1)
        acc = 1.0;
        for d=1:dim
            if P1(i,d)==0.0
                acc = acc * (1-P2(j,d));
            else
                acc = acc * P2(j,d);
            end
        end
        N(i,j) = acc;
    end
end

%% Isoparametric concept
P3 = N'*P3;

%% Move points onto circle (correction step)

for i=1:size(P3, 1)
   if abs(norm(P3(i,1:2)) - 1.00*cos(22.5/360*2*pi)) < 1e-6
       P3(i,1:2) = P3(i,1:2) / norm(P3(i,1:2));
   end
end

%% Output vtk-file
fileID = fopen(name, 'w');

fprintf(fileID,'# vtk DataFile Version 2.0\n');
fprintf(fileID,'Unstructured Grid Example\n');
fprintf(fileID,'ASCII\n');
fprintf(fileID,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fileID,'POINTS 23 float\n');
fprintf(fileID,'%12.8f %12.8f %12.8f\n', P3');

fprintf(fileID,'CELLS 7 63\n');
PP = [
    8  0  1  3  4 13  9 16 11;
    8  1  2  4  5  9 10 11 12;
    8  3  4  6  7 16 11 21 17;
    8  4  5  7  8 11 12 17 18;
    8  9 10 11 12 14 15 17 18;
    8 13  9 16 11 19 14 21 17;
    8 14 15 17 18 19 20 21 22
    ];


fprintf(fileID,'%3d %3d %3d %3d %3d %3d %3d %3d %3d\n', PP');

fprintf(fileID,'CELL_TYPES 7\n');
PP = [
    11;
    11;
    11;
    11;
    11;
    11;
    11
    ];

fprintf(fileID,'%4d\n', PP);

fclose(fileID);

end

