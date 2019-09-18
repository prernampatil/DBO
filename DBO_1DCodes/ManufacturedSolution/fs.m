function f=fs(t)

global NModes
global nu
global Nr
global x
t = double(t);
ue = u_exact(t);       ude  = ud_exact(t);
Ue = UDO_exact(t);     Ude  = UDOd_exact(t);
Ye = YDO_exact(t);     Yde  = YDOd_exact(t);

ue  = ue(x);
ude = ude(x);
Ue  = Ue(x);
Ude = Ude(x);

due     = fourdifft(ue,1);
ddue    = fourdifft(ue,2);
du      = Diff_rk4(Ue,1);
ddu     = Diff_rk4(Ue,2);


%% ---------------- Driving Force for Burgers Equation ------------------
LUe = nu*ddu  - repmat(ue,1,NModes).*du - Ue.*repmat(due,1,NModes);
BUe    = [];
Qe     = [];
Le     = [];

for j=1:NModes
    for k=1:NModes
        BUe= [BUe Ue(:,j).*du(:,k)];
        Qe  = [Qe Ye(:,j).*Ye(:,k)];
    end
end

f = repmat(ude +  ue.*due - nu*ddue,1,Nr);
f = f +  Ue*Yde' + Ude*Ye' + BUe*Qe' - LUe*Ye';

%% ----------------- Driving Force for Static Functions -----------
% f = repmat(ude,1,Nr)+ Ue*Yde' + Ude*Ye';

