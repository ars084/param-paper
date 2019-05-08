function BlandAltman(y,x,titlestr)
S(:,1) = (x + y)/2;
S(:,2) = 100*(x - y)./x;
avg = mean(S(:,2));
stdev = std(S(:,2));
plot(S(:,1),S(:,2),'k.','markersize',20)
hold on
plot([min(S(:,1))*0.9-1, max(S(:,1)*1.1+1)],[avg avg],'k-')
plot([min(S(:,1))*0.9-1, max(S(:,1)*1.1+1)],[avg-1.95*stdev,avg-1.95*stdev],'k--')
plot([min(S(:,1))*0.9-1, max(S(:,1)*1.1+1)],[avg+1.95*stdev,avg+1.95*stdev],'k--')

xlim([min(S(:,1))*0.9-1, max(S(:,1)*1.1+1)])
%truemax = max(abs(avg+2.1*stdev),abs(avg-2.1*stdev));
limit = max(max(abs(S(:,2))+10),abs(avg+1.95*stdev));
ylim([-limit,limit])
if nargin == 3
    title(titlestr,'fontsize',24)
end
ylabel('% Difference','fontsize',24)
xlabel('Mean','fontsize',24)
end