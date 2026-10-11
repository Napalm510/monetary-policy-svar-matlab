function [vd vdk] = TopicsVD2(sirfp,cirf, hor)
n=size(sirfp,2);

expvar=sum(sirfp.^2,3); % (numerator)
totvar = sum(sum(cirf.^2,3),2); % total variance of forecast error (denominator)

vd=expvar./totvar*100; % quota percentuale di varianza spiegata da ciascuno shock 

totvark = squeeze(sum(cumsum(cirf.^2,3),2)); % 
expvark = cumsum(sirfp.^2,3);

for j=1:n
    tempvd(:,:,j)=squeeze(expvark(:,j,:))./totvark;
end

vdk=tempvd(:,hor,:);