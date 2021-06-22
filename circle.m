function h = circle(x,y,r)

%function plotting the unitary circle, to check te position of the
%eigenvalues easily

hold on;
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit, 'k', 'LineWidth' , 0.7);

hold off;

end

