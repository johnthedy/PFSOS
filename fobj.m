function y=fobj(x,MdlC)

y1=1/predict(MdlC,x);

%Price per cubic in USD
pr1=0.9;
pr2=0.8;
pr3=0.7;
pr=pr2;
Cprice=150;
NA1price=30;
NA2price=100;
RA1price=pr*NA1price;
RA2price=pr*NA2price;

NAprice=-(NA2price-NA1price)*(x(5)-10)/20+NA2price;
RAprice=-(RA2price-RA1price)*(x(4)-10)/20+RA2price;

Cquantity=1/(x(1)+x(2)+1);
Aquantity=x(2)/(x(1)+x(2)+1);
NAquantity=(100-x(3))*Aquantity/100;
RAquantity=x(3)*Aquantity/100;
y2=Cquantity*Cprice+NAquantity*NAprice+RAquantity*RAprice;

y=[y1.*100 y2];
end

