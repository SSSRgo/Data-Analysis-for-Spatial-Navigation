function D=Viewdistance(ego)

pmaxX=50;
nmaxX=-50;
pmaxY=50;
nmaxY=-50;

ego.x(find(ego.x>=pmaxX))=pmaxX*0.999;
ego.y(find(ego.y>=pmaxY))=pmaxY*0.999;
ego.x(find(ego.x<=nmaxX))=nmaxX*0.999;
ego.y(find(ego.y<=nmaxY))=nmaxY*0.999;


y1=tan(ego.hd).*(50-ego.x)+ego.y; %(50,y1);
x2=(50-ego.y)./tan(ego.hd)+ego.x;%(x2,50);
y3=tan(ego.hd).*(-50-ego.x)+ego.y;%(-50,y3);
x4=(-50-ego.y)./tan(ego.hd)+ego.x;%(x4,-50);

y1(find(abs(y1)>50))=nan;
x2(find(abs(x2)>50))=nan;
y3(find(abs(y3)>50))=nan;
x4(find(abs(x4)>50))=nan;

theta=nan(length(ego.x),4);
theta(:,1)=mod(atan2((y1-ego.y),(50-ego.x))+2*pi,2*pi);
theta(:,2)=mod(atan2((50-ego.y),(x2-ego.x))+2*pi,2*pi);
theta(:,3)=mod(atan2((y3-ego.y),(-50-ego.x))+2*pi,2*pi);
theta(:,4)=mod(atan2((-50-ego.y),(x4-ego.x))+2*pi,2*pi);

Distance=nan(length(ego.x),4);
Distance(:,1)=sqrt((y1-ego.y).^2+(50-ego.x).^2);
Distance(:,2)=sqrt((50-ego.y).^2+(x2-ego.x).^2);
Distance(:,3)=sqrt((y3-ego.y).^2+(-50-ego.x).^2);
Distance(:,4)=sqrt((-50-ego.y).^2+(x4-ego.x).^2);

for i=    1:length(Distance)

    [row,col,v] = find(abs(theta(i,:)-ego.hd(i))<0.01);
    if ~isempty(row)
        col=min(col); % avoid the same distance 
        D(i,1)=Distance(i,col);
    else
        D(i,1)=nan;
    end

end

end