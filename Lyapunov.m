% Calculate Lyapunov's function for linear system of differential equations

if isOctave
    pkg load symbolic;
end

syms x y
syms p11 p12 p22

X = [x y].';

A = [2 -6; 2 6] * X;
W = -(x^2 + y^2);

nonlinear = [ (x - 3 * y)^2, (x + 3 * y)^2 ].';

P = [p11 p12; p12 p22];
V = X.' * P * X;

dV = [diff(V,x) diff(V,y)] * A;
[cV, tV] = coeffs(dV, X);
[cW, tW] = coeffs(-W, X);

S = cV;
found = 0;
for i = 1:length(cV)
  idx = find(tW == tV(i));
  if length(idx) > 0
    rhs = cW(idx);
    found += 1;
  else
    rhs = 0;
  end
  
  S(i) = S(i) == rhs;
end

if found < length(cW)
  printf("ERROR: rhs has terms that lhs does not")
  return
end

[r11 r12 r22] = solve(S);
V = expand(X.' * [r11 r12; r12 r22] * X);
printf("Result:\n")
V
W

printf("Verification...");

dV = [diff(V,x) diff(V,y)] * A;
[cV, tV] = coeffs(dV, X);
passed = all(isAlways(cV == cW));

if passed
  printf(" passed\n");
else
  printf(" not passed:\n%s != %s\n", char(cV), char(cW));
end

printf("\nNon-linear part:\n");
dV = [diff(V,x) diff(V,y)] * nonlinear;
[cV, tV] = coeffs(dV, X);

expand(dV)

