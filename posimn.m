% Code for posimn (equation (19)
function value=posimn(p,q,x,y,r)
value=besselk(p,x*r).*besseli(q,y*r)+(-1)^(p-q+1)*besselk(q,y*r).*besseli(p,x*r);