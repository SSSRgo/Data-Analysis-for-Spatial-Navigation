function A=circ_rmoutliers(A)

if ~isempty(A)

    Mean=nancirc_mean(A*pi/180)*180/pi;
    [s,std] = circ_std(A*pi/180);
    std=std*180/pi;
    A=A(find(abs(mod(A-Mean+180,360)-180)<(3*std))); % 3 times standard deviation 

else 

    A=A;
    
end

end