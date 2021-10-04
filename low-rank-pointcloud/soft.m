% define the soft threshold function, which is used above.
function y = soft(x,tau)

y = sign(x).*max(abs(x)-tau,0);

%y = max(x-tau,0); % by xqLu
