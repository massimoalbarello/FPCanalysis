A = [ 1 2 3 4;
    5 6 7 8;
    9 10 11 12];

c=[10,21,35,7];

P=[0 0 1 0;
  0 1 0 0;
  1 0 0 0;
  0 0 0 1];

A_permuted = P*A*P'
A_perm = permute_matrix(c,2,A)