%% Cell 1
ycor  = 0.55 * cos(pi/12)/(sin(pi/12)+cos(pi/12));
P1 = [
 0.00 0.55 0.00;  
 0.00 0.00 0.00;
 ycor ycor 0.00;
 0.55 0.00 0.00;
 0.00 0.55 0.50;  
 0.00 0.00 0.50;
 ycor ycor 0.50;
 0.55 0.00 0.50;
];
trans(P1, 'a.vtk')

%% Cell 2
ycos  = 0.5^0.5;
P1 = [
 +ycos ycos 0.00;  
 +1.00 0.00 0.00;
 +ycor ycor 0.00;
 +0.55 0.00 0.00;
 +ycos ycos 0.50;  
 +1.00 0.00 0.50;
 +ycor ycor 0.50;
 +0.55 0.00 0.50;
];
trans(P1, 'b.vtk')

%% Cell 3
ycos  = 0.5^0.5;
P1 = [
 +0.00 0.55 0.00;  
 +0.00 1.00 0.00;
 +ycor ycor 0.00;
 +ycos ycos 0.00;
 +0.00 0.55 0.50;  
 +0.00 1.00 0.50;
 +ycor ycor 0.50;
 +ycos ycos 0.50;
];
trans(P1, 'c.vtk')