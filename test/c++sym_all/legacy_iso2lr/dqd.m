def2ch[2];

Hc = gammaPolCh[1] hop[f[0], d[]] + gammaPolCh[2] hop[f[1], a[]];
Had = tpp hop[a[], d[]];

lrchain = {f[0], d[], a[], f[1]};

H = H0 + Hc + (H1 + Ha + Had);

nnop[f[0]] = 0;
nnop[d[]] = 1;
nnop[a[]] = 0;
nnop[f[1]] = 1;
