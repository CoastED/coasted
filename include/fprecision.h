#ifdef DOUBLEP
#define dp_to_nat(X) X
#define nat_to_dp(X) X
#else
#define dp_to_nat(X) sngl(X)
#define nat_to_dp(X) dble(X)
#endif
