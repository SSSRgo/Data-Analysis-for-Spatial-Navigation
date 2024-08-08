function DiffDeg=DiffAng(aa,bb)
% aa to bb, clockwise is positive 

adaf = @(a,b) diff(unwrap([a,b]/180*pi)*180/pi);
DiffDeg=adaf(aa,bb);
end