function [output] = bin(series,bins)
%BIN  series: Number,Channel

 [l,c]=size(series);
 n=floor(l/bins);
 output=zeros(n+1,c);
 for i=1:n
     output(i,:)=mean(series(bins*(i-1)+1:i*bins,:));
 end

output(n+1,:)=mean(series(bins*(n-1)+1:end,:));
end

