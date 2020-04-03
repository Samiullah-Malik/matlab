function y = periodicTopHat(x, alpha, beta)  
    xp = mod(x,2*pi);
    y=zeros(size(xp));
    i=xp > alpha & xp < beta;
    y(i==1) = 1;
    j=xp == alpha | xp == beta;
    y(j==1) = 1/2.; 
end