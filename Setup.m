HP = HestonModel();

clear B
clear ans
B = RostBarrierClass(HP,200,1000);

B.price(@(t) t)

T = B.H*ones(1,length(B.time_vec))+B.G;
T2 = T - ones(length(B.x),1)*(B.time_vec.^2/2);

surf(T2(1:50:2000,1:200:6001))


pause()

clear all

HP = HestonModel();

B = RostBarrierClass(HP,2000,2000);
f = @(t) 1.0*(t>0.05);

B.price(f)

[Pyff,h_f,iv_f] = B.SimulateHedgeFinal(f,2000,10);

hist(Pyff)
pause;
hist(h_f-Pyff)


