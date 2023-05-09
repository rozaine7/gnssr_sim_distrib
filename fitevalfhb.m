function output=fitevalfhb(xdata)%,ydata)
xdata=xdata(:);


inds1=find(xdata<=.022061);
inds2=find(xdata>.022061);
xd1=xdata(inds1);
xd2=xdata(inds2);

yd1=107.1*xd1+72.9;
yd2=577.4*exp(-319.2*xd2)+64.6*exp(6.5*xd2);

%
% These were values in the latest version of code that I got from JVoo (JLG
% 12/2010)
%
%yd1=1.005768785523987e+02*xd1+72.884254625333028;
%yd2=5.774239737798026e+02*exp(-3.191556944841277e+02*xd2)+64.639311451748767*exp(6.491948570146627*xd2);


%yd1 = 27.944770197190334*xd1+74.128355948087261;
%yd2 = 1.070181545854494e+002*exp(-2.020581772992435e+002*xd2)+65.649411575136384*exp(5.122974466333582*xd2);

output(inds1)=yd1;
output(inds2)=yd2;