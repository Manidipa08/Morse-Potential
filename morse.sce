clc
clear
n=1001
rmin=1D-15
rmax=1D-9
r=linspace(rmin,rmax,n)
new=r(2:n-1)
rnew=new'
dr=(rmax-rmin)/(n-1)
h=6.626D-34
h_bar=(h/%pi)
e=1.6D-19
p=1.67D-27
/*[HF HCL HBr HI CO NO ICl ]
xe=[0.0218,0.0174,0.0171,0.0172,0.0061,0.073,0.0038]
k=[966,516,412,314,1902,1595,238]
req=[0.0927,0.1274,0.1414,0.16609,0.1131,0.1151,0.2321]*1d-9
m1=[1,1,1,1,12,14,127]
m2=[19,35.5,80,127,16,16,35.5]
*/
xe=input("Enter the value of xe : ")
k=input("Enter the value of k : : ") //FORCE CONSTANT
req=input("Enter req : ")//internuclear distance
m1=input("Enter m1 : ")
m2=input("Enter m2 : ")
m_reduced=((m1*m2)/(m1+m2))*p;
we=sqrt(k/m_reduced); // rad per sec
Deq=h_bar*we/(4*xe) // in J
De=Deq/e // in eV
a=sqrt(k/(2*Deq)) // m^-1
Vm=round(1/(2*xe)-(1/2)) //max vibrational quantum no 
for i=1:n
 V0(i)=De*((1-exp(a*(req-r(i))))^2)
end
vnew=diag(V0(2:n-1));
C=-(h_bar^2)/(2*m_reduced*dr*dr*e)
KE=diag(-2*C*ones(1,n-2))+diag(C*ones(1,n-3),1)+diag(C*ones(1,n-3),-1)
H=KE+vnew;
[a1,b1]=spec(H);
Z=spec(H);
counter=0; //for no of Bound States 
for i=1:n-2
 if Z(i)<De
 counter=counter+1
 end
end
plot(r'*1D+9,De*ones(1,n),'o-c')
plot(r'*1D+9,V0','g','linewidth',3)
for i=1:counter
 plot(rnew*1D+9,Z(i)*ones(n-2,1),':r','linewidth',2)
 ss=gca()
 legend('D-equlibrium','Potential curve','Bound States')
end
at= gca()
at.font_size = 4; //tics label font size
at.data_bounds=[0,0;1,De+1]
title("< MORSE POTENTIAL VS R PLOT >",'color','red','Fontsize',5)
xlabel("R in nm >>>>",'color','blue','Fontsize',4)
ylabel("Morse Potential (V) in eV >>>>",'color','blue','Fontsize',4)
xgrid()
