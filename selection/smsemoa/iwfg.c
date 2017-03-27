#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "iwfg.h"
int i_fr = 0;   // current depth

int i_sorter(const void *a, const void *b)
//sort point indexes into groups of last objective
{
	int i = *(int*)a;
	int j = *(int*)b;
	if (torder[i][i_n-1]==torder[j][i_n-1]) {
		return tcompare[j][torder[j][i_n-1]] - tcompare[i][torder[i][i_n-1]];
	}
	else {
		return torder[i][i_n-1] - torder[j][i_n-1];
	}
}

int i_greater(const void *v1, const void *v2)
// this sorts points worsening in the last objective
{
	int i;
	POINT p = *(POINT*)v1;
	POINT q = *(POINT*)v2;
	for (i = i_n - 1; i >= 0; i--) {
		if BEATS(p.objectives[i],q.objectives[i]) {
			return -1;
		}
		else if BEATS(q.objectives[i],p.objectives[i]) {
			return  1;
		}
	}
	return 0;
}

int i_same(const void *v1, const void *v2)
// this sorts points worsening in the last objective for a certain objective ordering
{
	int i;
	POINT p = *(POINT*)v1;
	POINT q = *(POINT*)v2;
	for ( i = i_n - 1; i >= 0; i--) {
		if (p.objectives[i]!=q.objectives[i]) {
			return 0;
		}
	}
	return 1;
}
int i_greaterorder(const void *v1, const void *v2)
// this sorts points worsening in the last objective for a certain objective ordering
{
	int i;
	POINT p = *(POINT*)v1;
	POINT q = *(POINT*)v2;
	for ( i = i_n - 1; i >= 0; i--) {
		if BEATS(p.objectives[gorder[i]],q.objectives[gorder[i]]) {
			return -1;
		}
		else if BEATS(q.objectives[gorder[i]],p.objectives[gorder[i]]) {
			return  1;
		}
	}
	return 0;
}

int i_greaterabbrev(const void *v1, const void *v2)
// this sorts points worsening in the penultimate objective
{
	int i;
	POINT p = *(POINT*)v1;
	POINT q = *(POINT*)v2;
	for ( i = i_n - 2; i >= 0; i--) {
		if BEATS(p.objectives[i],q.objectives[i]) {
			return -1;
		}
		else if BEATS(q.objectives[i],p.objectives[i]) {
			return  1;
		}
	}
	return 0;
}

int i_greaterabbrevorder(const void *v1, const void *v2)
// this sorts points worsening in the penultimate objective for a certain objective ordering
{
	int i;
	POINT p = *(POINT*)v1;
	POINT q = *(POINT*)v2;
	for ( i = i_n - 2; i >= 0; i--) {
		if BEATS(p.objectives[gorder[i]],q.objectives[gorder[i]]) {
			return -1;
		}
		else if BEATS(q.objectives[gorder[i]],p.objectives[gorder[i]]) {
			return  1;
		}
	}
	return 0;
}

int i_dominates2way(POINT p, POINT q, int k)
// returns -1 if p dominates q, 1 if q dominates p, 2 if p == q, 0 otherwise
// k is the highest index inspected
{
	int i,j;
	for (i = k; i >= 0; i--) {
		if BEATS(p.objectives[i],q.objectives[i]) {
			for ( j = i - 1; j >= 0; j--) {
	 			if BEATS(q.objectives[j],p.objectives[j]) {
	 				return 0;
	 			}
			}
			return -1;
		}
		else if BEATS(q.objectives[i],p.objectives[i]) {
			for ( j = i - 1; j >= 0; j--) {
				if BEATS(p.objectives[j],q.objectives[j]) {
					return 0;
				}
			}
			return  1;
		}
	}
	return 2;
}

bool i_dominates1way(POINT p, POINT q, int k)
// returns true if p dominates q or p == q, false otherwise
// the assumption is that q doesn't dominate p
// k is the highest index inspected
{
	int i;
	for ( i = k; i >= 0; i--) {
		if BEATS(q.objectives[i],p.objectives[i]) {
			return false;
		}
	}
	return true;
}

bool i_dominates1wayOrder(POINT p, POINT q, int k, int* order)
// returns true if p dominates q or p == q, false otherwise
// the assumption is that q doesn't dominate p
// k is the highest index inspected
{
	int i;
	for ( i = k; i >= 0; i--) {
		if BEATS(q.objectives[order[i]],p.objectives[order[i]]) {
			return false;
		}
	}
	return true;
}

void i_removeDominated(int l, int limit)
{
	POINT t;
	int i,j;
  // points below l are all equal in the last objective; points above l are all worse
  // points below l can dominate each other, and we don't need to compare the last objective
  // points above l cannot dominate points that start below l, and we don't need to compare the last objective
	i_fs[i_fr].nPoints = 1;
	for ( i = 1; i < l; i++) {
		j = 0;
		while (j < i_fs[i_fr].nPoints) {
			switch (i_dominates2way(i_fs[i_fr].points[i], i_fs[i_fr].points[j], i_n-2)) {
				case  0:
					j++;
					break;
				case -1: // AT THIS POINT WE KNOW THAT i CANNOT BE DOMINATED BY ANY OTHER PROMOTED POINT j
					// SWAP i INTO j, AND 1-WAY DOM FOR THE REST OF THE js
					t = i_fs[i_fr].points[j];
					i_fs[i_fr].points[j] = i_fs[i_fr].points[i];
					i_fs[i_fr].points[i] = t;
					while(j < i_fs[i_fr].nPoints - 1 && i_dominates1way(i_fs[i_fr].points[j], i_fs[i_fr].points[i_fs[i_fr].nPoints - 1], i_n-1)) {
						i_fs[i_fr].nPoints--;
					}
					int k = j+1;
					while (k < i_fs[i_fr].nPoints) {
						if(i_dominates1way(i_fs[i_fr].points[j], i_fs[i_fr].points[k], i_n-2)) {
							t = i_fs[i_fr].points[k];
							i_fs[i_fr].nPoints--;
							i_fs[i_fr].points[k] = i_fs[i_fr].points[i_fs[i_fr].nPoints];
							i_fs[i_fr].points[i_fs[i_fr].nPoints] = t;
						}
						else {
							k++;
						}
					}
				default:
					j = i_fs[i_fr].nPoints + 1;
			}
		}
		if (j == i_fs[i_fr].nPoints) {
			t = i_fs[i_fr].points[i_fs[i_fr].nPoints];
			i_fs[i_fr].points[i_fs[i_fr].nPoints] = i_fs[i_fr].points[i];
			i_fs[i_fr].points[i] = t;
			i_fs[i_fr].nPoints++;
		}
	}
	i_safe = WORSE(l,i_fs[i_fr].nPoints);
	for ( i = l; i < limit; i++) {
		j = 0;
		while (j < i_safe) {
			if(i_dominates1way(i_fs[i_fr].points[j], i_fs[i_fr].points[i], i_n-2)) {
				j = i_fs[i_fr].nPoints + 1;
			}
			else {
				j++;
			}
		}
		while (j < i_fs[i_fr].nPoints) {
			switch (i_dominates2way(i_fs[i_fr].points[i], i_fs[i_fr].points[j], i_n-1)) {
				case  0:
					j++;
					break;
				case -1: // AT THIS POINT WE KNOW THAT i CANNOT BE DOMINATED BY ANY OTHER PROMOTED POINT j
					// SWAP i INTO j, AND 1-WAY DOM FOR THE REST OF THE js
					t = i_fs[i_fr].points[j];
					i_fs[i_fr].points[j] = i_fs[i_fr].points[i];
					i_fs[i_fr].points[i] = t;
					while(j < i_fs[i_fr].nPoints - 1 && i_dominates1way(i_fs[i_fr].points[j], i_fs[i_fr].points[i_fs[i_fr].nPoints - 1], i_n-1)) {
						i_fs[i_fr].nPoints--;
					}
					int k = j+1;
					while (k < i_fs[i_fr].nPoints) {
						if(i_dominates1way(i_fs[i_fr].points[j], i_fs[i_fr].points[k], i_n-1)) {
							t = i_fs[i_fr].points[k];
							i_fs[i_fr].nPoints--;
							i_fs[i_fr].points[k] = i_fs[i_fr].points[i_fs[i_fr].nPoints];
							i_fs[i_fr].points[i_fs[i_fr].nPoints] = t;
						}
						else {
							k++;
						}
					}
				default:
					j = i_fs[i_fr].nPoints + 1;
			}
		}
		if (j == i_fs[i_fr].nPoints) {
			t = i_fs[i_fr].points[i_fs[i_fr].nPoints];
			i_fs[i_fr].points[i_fs[i_fr].nPoints] = i_fs[i_fr].points[i];
			i_fs[i_fr].points[i] = t;
			i_fs[i_fr].nPoints++;
		}
	}
	i_fr++;
}

void i_makeDominatedBit(FRONT ps, int p)
// creates the front ps[0 .. p-1] in i_fs[i_fr], with each point bounded by ps[p] and dominated points removed
{
	int i,j;
	int l = 0;
	int u = p - 1;
	for (i = p - 1; i >= 0; i--) {
		if (BEATS(ps.points[p].objectives[i_n - 1],ps.points[i].objectives[i_n - 1])) {
			i_fs[i_fr].points[u].objectives[i_n - 1] = ps.points[i].objectives[i_n - 1];
			for ( j = 0; j < i_n - 1; j++) {
				i_fs[i_fr].points[u].objectives[j] = WORSE(ps.points[p].objectives[j],ps.points[i].objectives[j]);
			}
			u--;
		}
		else {
			i_fs[i_fr].points[l].objectives[i_n - 1] = ps.points[p].objectives[i_n - 1];
			for ( j = 0; j < i_n - 1; j++) {
				i_fs[i_fr].points[l].objectives[j] = WORSE(ps.points[p].objectives[j],ps.points[i].objectives[j]);
			}
			l++;
		}
	}
	i_removeDominated(l,p);
}


double i_hv2(FRONT ps, int k)
// returns the hypervolume of ps[0 .. k-1] in 2D
// assumes that ps is sorted improving
{
	int i;
	double volume = ps.points[0].objectives[0] * ps.points[0].objectives[1];
	for ( i = 1; i < k; i++) {
		volume += ps.points[i].objectives[1] * (ps.points[i].objectives[0] - ps.points[i - 1].objectives[0]);
	}
	return volume;
}


double i_inclhv(POINT p)
// returns the inclusive hypervolume of p
{
	int i;
	double volume = 1;
	for ( i = 0; i < i_n; i++) {
		volume *= p.objectives[i];
	}
	return volume;
}

double i_inclhvOrder(POINT p, int* order)
// returns the inclusive hypervolume of p
{
	int i;
	double volume = 1;
	for ( i = 0; i < i_n; i++) {
		volume *= p.objectives[order[i]];
	}
	return volume;
}

double i_inclhv2(POINT p, POINT q)
// returns the hypervolume of {p, q}
{
	int i;
	double vp  = 1;
	double vq  = 1;
	double vpq = 1;
	for ( i = 0; i < i_n; i++) {
		vp  *= p.objectives[i];
		vq  *= q.objectives[i];
		vpq *= WORSE(p.objectives[i],q.objectives[i]);
	}
	return vp + vq - vpq;
}


double i_inclhv3(POINT p, POINT q, POINT r)
// returns the hypervolume of {p, q, r}
{
	int i;
	double vp = 1;
	double vq = 1;
	double vr = 1;
	double vpq = 1;
	double vpr = 1;
	double vqr = 1;
	double vpqr = 1;
	for ( i = 0; i < i_n; i++) {
		vp *= p.objectives[i];
		vq *= q.objectives[i];
		vr *= r.objectives[i];
		if (BEATS(p.objectives[i],q.objectives[i])) {
			if (BEATS(q.objectives[i],r.objectives[i])) {
				vpq  *= q.objectives[i];
				vpr  *= r.objectives[i];
				vqr  *= r.objectives[i];
				vpqr *= r.objectives[i];
			}
			else {
				vpq  *= q.objectives[i];
				vpr  *= WORSE(p.objectives[i],r.objectives[i]);
				vqr  *= q.objectives[i];
				vpqr *= q.objectives[i];
			}
		}
		else if (BEATS(p.objectives[i],r.objectives[i])) {
			vpq  *= p.objectives[i];
			vpr  *= r.objectives[i];
			vqr  *= r.objectives[i];
			vpqr *= r.objectives[i];
		}
		else {
			vpq  *= p.objectives[i];
			vpr  *= p.objectives[i];
			vqr  *= WORSE(q.objectives[i],r.objectives[i]);
			vpqr *= p.objectives[i];
		}

	}
	return vp + vq + vr - vpq - vpr - vqr + vpqr;
}


double i_inclhv4(POINT p, POINT q, POINT r, POINT s)
// returns the hypervolume of {p, q, r, s}
{
	double vp = 1;
	double vq = 1;
	double vr = 1;
	double vs = 1;
	double vpq = 1;
	double vpr = 1;
	double vps = 1;
	double vqr = 1;
	double vqs = 1;
	double vrs = 1;
	double vpqr = 1;
	double vpqs = 1;
	double vprs = 1;
	double vqrs = 1;
	double vpqrs = 1;
	int i;
	for ( i = 0; i < i_n; i++) {
		vp *= p.objectives[i];
		vq *= q.objectives[i];
		vr *= r.objectives[i];
		vs *= s.objectives[i];
		if (BEATS(p.objectives[i],q.objectives[i])) {
			if (BEATS(q.objectives[i],r.objectives[i])) {
				if (BEATS(r.objectives[i],s.objectives[i])) {
					vpq *= q.objectives[i];
					vpr *= r.objectives[i];
					vps *= s.objectives[i];
					vqr *= r.objectives[i];
					vqs *= s.objectives[i];
					vrs *= s.objectives[i];
					vpqr *= r.objectives[i];
					vpqs *= s.objectives[i];
					vprs *= s.objectives[i];
					vqrs *= s.objectives[i];
					vpqrs *= s.objectives[i];
				}
				else {
					OBJECTIVE z1 = WORSE(q.objectives[i],s.objectives[i]);
					vpq *= q.objectives[i];
					vpr *= r.objectives[i];
					vps *= WORSE(p.objectives[i],s.objectives[i]);
					vqr *= r.objectives[i];
					vqs *= z1;
					vrs *= r.objectives[i];
					vpqr *= r.objectives[i];
					vpqs *= z1;
					vprs *= r.objectives[i];
					vqrs *= r.objectives[i];
					vpqrs *= r.objectives[i];
				}
			}
			else if (BEATS(q.objectives[i],s.objectives[i])) {
				vpq *= q.objectives[i];
				vpr *= WORSE(p.objectives[i],r.objectives[i]);
				vps *= s.objectives[i];
				vqr *= q.objectives[i];
				vqs *= s.objectives[i];
				vrs *= s.objectives[i];
				vpqr *= q.objectives[i];
				vpqs *= s.objectives[i];
				vprs *= s.objectives[i];
				vqrs *= s.objectives[i];
				vpqrs *= s.objectives[i];
			}
			else {
				OBJECTIVE z1 = WORSE(p.objectives[i],r.objectives[i]);
				vpq *= q.objectives[i];
				vpr *= z1;
				vps *= WORSE(p.objectives[i],s.objectives[i]);
				vqr *= q.objectives[i];
				vqs *= q.objectives[i];
				vrs *= WORSE(r.objectives[i],s.objectives[i]);
				vpqr *= q.objectives[i];
				vpqs *= q.objectives[i];
				vprs *= WORSE(z1,s.objectives[i]);
				vqrs *= q.objectives[i];
				vpqrs *= q.objectives[i];
			}
		}
		else if (BEATS(q.objectives[i],r.objectives[i])) {
			if (BEATS(p.objectives[i],s.objectives[i])) {
				OBJECTIVE z1 = WORSE(p.objectives[i],r.objectives[i]);
				OBJECTIVE z2 = WORSE(r.objectives[i],s.objectives[i]);
				vpq *= p.objectives[i];
				vpr *= z1;
				vps *= s.objectives[i];
				vqr *= r.objectives[i];
				vqs *= s.objectives[i];
				vrs *= z2;
				vpqr *= z1;
				vpqs *= s.objectives[i];
				vprs *= z2;
				vqrs *= z2;
				vpqrs *= z2;
			}
			else {
				OBJECTIVE z1 = WORSE(p.objectives[i],r.objectives[i]);
				OBJECTIVE z2 = WORSE(r.objectives[i],s.objectives[i]);
				vpq *= p.objectives[i];
				vpr *= z1;
				vps *= p.objectives[i];
				vqr *= r.objectives[i];
				vqs *= WORSE(q.objectives[i],s.objectives[i]);
				vrs *= z2;
				vpqr *= z1;
				vpqs *= p.objectives[i];
				vprs *= z1;
				vqrs *= z2;
				vpqrs *= z1;
			}
		}
		else if (BEATS(p.objectives[i],s.objectives[i])) {
				vpq *= p.objectives[i];
				vpr *= p.objectives[i];
				vps *= s.objectives[i];
				vqr *= q.objectives[i];
				vqs *= s.objectives[i];
				vrs *= s.objectives[i];
				vpqr *= p.objectives[i];
				vpqs *= s.objectives[i];
				vprs *= s.objectives[i];
				vqrs *= s.objectives[i];
				vpqrs *= s.objectives[i];
			}
		else {
			OBJECTIVE z1 = WORSE(q.objectives[i],s.objectives[i]);
			vpq *= p.objectives[i];
			vpr *= p.objectives[i];
			vps *= p.objectives[i];
			vqr *= q.objectives[i];
			vqs *= z1;
			vrs *= WORSE(r.objectives[i],s.objectives[i]);
			vpqr *= p.objectives[i];
			vpqs *= p.objectives[i];
			vprs *= p.objectives[i];
			vqrs *= z1;
			vpqrs *= p.objectives[i];
		}
	}
	return vp + vq + vr + vs - vpq - vpr - vps - vqr - vqs - vrs + vpqr + vpqs + vprs + vqrs - vpqrs;
}

double i_inclhv5(POINT p, POINT q, POINT r, POINT s, POINT t)
// returns the hypervolume of {p, q, r, s, t}
{
	double vp = 1;
	double vq = 1;
	double vr = 1;
	double vs = 1;
	double vt = 1;

	double vpq = 1;
	double vpr = 1;
	double vps = 1;
	double vpt = 1;
	double vqr = 1;
	double vqs = 1;
	double vqt = 1;
	double vrs = 1;
	double vrt = 1;
	double vst = 1;

	double vpqr = 1;
	double vpqs = 1;
	double vpqt = 1;
	double vprs = 1;
	double vprt = 1;
	double vpst = 1;
	double vqrs = 1;
	double vqrt = 1;
	double vqst = 1;
	double vrst = 1;

	double vpqrs = 1;
	double vpqrt = 1;
	double vpqst = 1;
	double vprst = 1;
	double vqrst = 1;

	double vpqrst = 1;
	int i;

	for ( i = 0; i < i_n; i++) {
		vp *= p.objectives[i];
		vq *= q.objectives[i];
		vr *= r.objectives[i];
		vs *= s.objectives[i];
		vt *= t.objectives[i];
		vpq *= WORSE(p.objectives[i],q.objectives[i]);
		vpr *= WORSE(p.objectives[i],r.objectives[i]);
		vps *= WORSE(p.objectives[i],s.objectives[i]);
		vpt *= WORSE(p.objectives[i],t.objectives[i]);
		vqr *= WORSE(q.objectives[i],r.objectives[i]);
		vqs *= WORSE(q.objectives[i],s.objectives[i]);
		vqt *= WORSE(q.objectives[i],t.objectives[i]);
		vrs *= WORSE(r.objectives[i],s.objectives[i]);
		vrt *= WORSE(r.objectives[i],t.objectives[i]);
		vst *= WORSE(s.objectives[i],t.objectives[i]);
		vpqr *= WORSE(p.objectives[i],
                      WORSE(q.objectives[i],r.objectives[i]));
		vpqs *= WORSE(p.objectives[i],
                      WORSE(q.objectives[i],s.objectives[i]));
		vpqt *= WORSE(p.objectives[i],
                      WORSE(q.objectives[i],t.objectives[i]));
		vprs *= WORSE(p.objectives[i],
                      WORSE(r.objectives[i],s.objectives[i]));
		vprt *= WORSE(p.objectives[i],
                      WORSE(r.objectives[i],t.objectives[i]));
		vpst *= WORSE(p.objectives[i],
                      WORSE(s.objectives[i],t.objectives[i]));
		vqrs *= WORSE(q.objectives[i],
                      WORSE(r.objectives[i],s.objectives[i]));
		vqrt *= WORSE(q.objectives[i],
                      WORSE(r.objectives[i],t.objectives[i]));
		vqst *= WORSE(q.objectives[i],
                      WORSE(s.objectives[i],t.objectives[i]));
		vrst *= WORSE(r.objectives[i],
                      WORSE(s.objectives[i],t.objectives[i]));
		vpqrs *= WORSE(WORSE(p.objectives[i],q.objectives[i]),
                         WORSE(r.objectives[i],s.objectives[i]));
		vpqrt *= WORSE(WORSE(p.objectives[i],q.objectives[i]),
                         WORSE(r.objectives[i],t.objectives[i]));
		vpqst *= WORSE(WORSE(p.objectives[i],q.objectives[i]),
                         WORSE(s.objectives[i],t.objectives[i]));
		vprst *= WORSE(WORSE(p.objectives[i],r.objectives[i]),
                         WORSE(s.objectives[i],t.objectives[i]));
		vqrst *= WORSE(WORSE(q.objectives[i],r.objectives[i]),
                         WORSE(s.objectives[i],t.objectives[i]));
		vpqrst *= WORSE(WORSE(p.objectives[i],
                         WORSE(q.objectives[i],r.objectives[i])),
                         WORSE(s.objectives[i],t.objectives[i]));
	}
	return vp + vq + vr + vs + vt
- vpq  - vpr  - vps  - vpt  - vqr  - vqs  - vqt  - vrs  - vrt  - vst
+ vpqr + vpqs + vpqt + vprs + vprt + vpst + vqrs + vqrt + vqst + vrst
- vpqrs - vpqrt - vpqst - vprst - vqrst
+ vpqrst;
}


double i_exclhv(FRONT ps, int p)
// returns the exclusive hypervolume of ps[p] relative to ps[0 .. p-1]
{
	i_makeDominatedBit(ps, p);
	double volume = i_inclhv(ps.points[p]) - i_hv(i_fs[i_fr - 1]);
	i_fr--;
	return volume;
}

double i_hv_contribution(FRONT ps, int id,double whole_hv)
{
    int i,j;
/*
	printf("\n------%d-------\n",id);
	for(i=0;i<ps.nPoints;i++)
	{
		for(j=0;j<ps.n;j++)
			printf("%lf ",ps.points[i].objectives[j]);
		printf("\n");
	}
*/
	for(i=0;i<ps.n;i++)
		ps.points[id].objectives[i]=0;

/*
	printf("\n-------%d--------\n",id);
	for(i=0;i<ps.nPoints;i++)
	{
		for(j=0;j<ps.n;j++)
			printf("%lf ",ps.points[i].objectives[j]);
		printf("\n");
	}
	printf("\n--------------\n");
*/


    whole_hv = whole_hv - i_hv(ps);


	/*
	for(i=0;i<ps.nPoints;i++)
	{
		for(j=0;j<ps.n;j++)
			printf("%lf ",ps.points[i].objectives[j]);
		printf("\n");
	}
	printf("\n-----++++------\n");
    */
    return  whole_hv;
}
double i_hv(FRONT ps)
// returns the hypervolume of ps[0 ..]
{
	int i;
	// process small fronts with the IEA
	switch (ps.nPoints) {
		case 1:
			return i_inclhv (ps.points[0]);
		case 2:
			return i_inclhv2(ps.points[0], ps.points[1]);
		case 3:
			return i_inclhv3(ps.points[0], ps.points[1], ps.points[2]);
		case 4:
			return i_inclhv4(ps.points[0], ps.points[1], ps.points[2], ps.points[3]);
	}

	// these points need sorting
	qsort(&ps.points[i_safe], ps.nPoints - i_safe, sizeof(POINT), i_greater);
	// n = 2 implies that safe = 0
	if (i_n == 2) {
		return i_hv2(ps, ps.nPoints);
	}
	// these points don't NEED sorting, but it helps
	qsort(ps.points, i_safe, sizeof(POINT), i_greaterabbrev);

	if (i_n == 3 && i_safe > 0) {
		double volume = ps.points[0].objectives[2] * i_hv2(ps, i_safe);
		i_n--;
		for ( i = i_safe; i < ps.nPoints; i++) {
			// we can ditch dominated points here, but they will be ditched anyway in makeDominatedBit
			volume += ps.points[i].objectives[i_n] * i_exclhv(ps, i);
		}
		i_n++;
		return volume;
	}
	else {
		double volume = i_inclhv4(ps.points[0], ps.points[1], ps.points[2], ps.points[3]);
		i_n--;
		for ( i = 4; i < ps.nPoints; i++) {
			// we can ditch dominated points here, but they will be ditched anyway in makeDominatedBit
			volume += ps.points[i].objectives[i_n] * i_exclhv(ps, i);
		}
		i_n++;
		return volume;
	}
}

void i_makeDominatedBitPoint(FRONT ps, POINT p, int* order)
// creates the front ps in i_fs[i_fr], with each point bounded by p and dominated points removed
{
	int i,j;
	int l = 0;
	int u = ps.nPoints - 1;
	for ( i = ps.nPoints - 1; i >= 0; i--) {
		if (BEATS(p.objectives[order[i_n - 1]],ps.points[i].objectives[order[i_n - 1]])) {
			i_fs[i_fr].points[u].objectives[i_n - 1] = ps.points[i].objectives[order[i_n - 1]];
			for ( j = 0; j < i_n - 1; j++) {
				i_fs[i_fr].points[u].objectives[j] = WORSE(p.objectives[order[j]],ps.points[i].objectives[order[j]]);
			}
			u--;
		}
		else {
			i_fs[i_fr].points[l].objectives[i_n - 1] = p.objectives[order[i_n - 1]];
			for ( j = 0; j < i_n - 1; j++) {
				i_fs[i_fr].points[l].objectives[j] = WORSE(p.objectives[order[j]],ps.points[i].objectives[order[j]]); 
			}
			l++;
		}
	}
	i_removeDominated(l,ps.nPoints);
}

double i_exclhvPoint(FRONT ps, POINT p, int* order)
// returns the exclusive hypervolume of p relative to ps
{
	i_makeDominatedBitPoint(ps, p, order);
	double volume = i_inclhvOrder(p,order) - i_hv(i_fs[i_fr - 1]);
	i_fr--;
	return volume;
}

void i_heapify(int location, int index)
// restores heap property starting at location and working downwards to place index in heap
{
	while (2*location+2<heapsize) {
		bool left = false;
		bool right = false;
		if (partial[heap[2*location+1]] < partial[index]) {
			left = true;
		}
		if (partial[heap[2*location+2]] < partial[index]) {
			right = true;
		}
		if (left) {
			if (right && partial[heap[2*location+2]] < partial[heap[2*location+1]]) {
				heap[location] = heap[2*location+2];
				location = 2*location+2;
			}
			else {
				heap[location] = heap[2*location+1];
				location = 2*location+1;
			}
		}
		else if (right) {
			heap[location] = heap[2*location+2];
			location = 2*location+2;
		}
		else {
			break;
		}
	}
	if (2*location+1<heapsize && partial[heap[2*location+1]] < partial[index]) {
		heap[location] = heap[2*location+1];
		location = 2*location+1;
	}
	heap[location] = index;
}

int i_peekFromHeap(void)
{
	return heap[0];
}

void i_initialiseHeap(int capacity)
// creates the heap with the indexes 0..(capacity-1) 
{
	int i;
	heapsize = capacity;
	for ( i=heapsize-1; i>=0; i--) {
		i_heapify(i,i);
	}
}

void i_insert(POINT p, int k, FRONT pl, int i, int j, int* order)
// inserts p into pl with the result in stacks[i][j]
{
	int place = 0;
	while (place < pl.nPoints && pl.points[place].objectives[order[k]] > p.objectives[order[k]]) {
		stacks[i][j].front.points[place] = pl.points[place];
		place++;
	}
	POINT pp = pl.points[place];
	stacks[i][j].front.points[place] = p;
	int placeNext = place + 1;
	POINT ppn = pl.points[place+1];
	while (place < pl.nPoints) {
		if (!i_dominates1wayOrder(p,pp,k,order)) {
			stacks[i][j].front.points[placeNext] = pp;
			placeNext++;
		}
		place++;
		pp = ppn;
		ppn = pl.points[place+1];
	}
	stacks[i][j].front.nPoints = placeNext;
}

void i_sliceOrder(int nPoints)
// slice using a separate objective ordering per point
{
	int i,j;
	int sorder[nPoints];
	for ( i=0; i<nPoints; i++) {
		sorder[i] = i;
	}
	qsort(sorder,nPoints,sizeof(int),i_sorter);
	int seen = 0;
	int p;
	for ( p=0; p<nPoints; p++) {
		i = sorder[p];
		if (p==0 || torder[i][i_n-1]!=torder[sorder[p-1]][i_n-1]) {
			seen = 0;
			stacks[i][1].front.nPoints = 0;
		}
		else {
			for (j=0; j<stacks[sorder[p-1]][1].front.nPoints; j++) {
				stacks[i][1].front.points[j] = stacks[sorder[p-1]][1].front.points[j];
			}
			stacks[i][1].front.nPoints = stacks[sorder[p-1]][1].front.nPoints;
		}
		/////////
		/*
		int ii,kk;
		printf("\ntcompare:\n");
		for (ii=0; ii<nPoints; ii++) {
			for (kk=0; kk<i_n; kk++) {
				printf("%d ",tcompare[ii][kk]);
			}
			printf("\n");
		}
		printf("\ntorder:\n");
		for (ii=0; ii<nPoints; ii++) {
			for (kk=0; kk<i_n; kk++) {
				printf("%d ",torder[ii][kk]);
			}
			printf("\n");
		}
		 */
		/////////
		int pos = nPoints-1-tcompare[i][torder[i][i_n-1]];
		for (j=seen; j<pos; j++) {
			stacks[i][1].front.points[stacks[i][1].front.nPoints+j-seen] = stacks[i][0].front.points[j];
		}
		int start = stacks[i][1].front.nPoints;
		int end = stacks[i][1].front.nPoints+pos-seen;
		seen = pos;
		POINT temp;
		for (j=start; j<end; j++) {
			int k = 0;
			while (k<stacks[i][1].front.nPoints) {
				if (i_dominates1wayOrder(stacks[i][1].front.points[j],stacks[i][1].front.points[k],i_n-2,torder[i])) {
					temp = stacks[i][1].front.points[k];
					stacks[i][1].front.points[k] = stacks[i][1].front.points[j];
					stacks[i][1].front.points[j] = temp;
					while(k<stacks[i][1].front.nPoints-1 &&
						  i_dominates1wayOrder(stacks[i][1].front.points[k],stacks[i][1].front.points[stacks[i][1].front.nPoints-1],i_n-2,torder[i])) {
						stacks[i][1].front.nPoints--;
					}
					int l = k+1;
					while (l < stacks[i][1].front.nPoints) {
						if(i_dominates1wayOrder(stacks[i][1].front.points[k],stacks[i][1].front.points[l],i_n-2,torder[i])) {
							temp = stacks[i][1].front.points[l];
							stacks[i][1].front.nPoints--;
							stacks[i][1].front.points[l] = stacks[i][1].front.points[stacks[i][1].front.nPoints];
							stacks[i][1].front.points[stacks[i][1].front.nPoints] = temp;
						}
						else {
							l++;
						}
					}
					k = stacks[i][1].front.nPoints + 1;
				}
				else {
					k++;
				}
			}
			if (k==stacks[i][1].front.nPoints) {
				temp = stacks[i][1].front.points[stacks[i][1].front.nPoints];
				stacks[i][1].front.points[stacks[i][1].front.nPoints] = stacks[i][1].front.points[j];
				stacks[i][1].front.points[j] = temp;
				stacks[i][1].front.nPoints++;
			}
		}
		stacks[i][1].index = pos+1;
		if (pos<nPoints-1) {
			stacks[i][1].width = fabs(stacks[i][0].front.points[pos].objectives[torder[i][i_n-1]]-
									  stacks[i][0].front.points[pos+1].objectives[torder[i][i_n-1]]);
		}
		else {
			stacks[i][1].width = stacks[i][0].front.points[pos].objectives[torder[i][i_n-1]];
		}
	}
	for (i=0; i<nPoints; i++) {
		gorder = torder[i];
		qsort(stacks[i][1].front.points,stacks[i][1].front.nPoints,sizeof(POINT),i_greaterabbrevorder);
	}
}

void i_slice(FRONT pl)
// slice in the last objective
{
	int i;

	stacks[0][1].front.nPoints = 0;

	for ( i=0; i<pl.nPoints-1; i++) {
		stacks[i][1].width = fabs(pl.points[i].objectives[i_n-1]-pl.points[i+1].objectives[i_n-1]);
		stacks[i][1].index = i+1;
		i_insert(pl.points[i],i_n-2,stacks[i][1].front,i+1,1,torder[i]);
	}
	stacks[pl.nPoints-1][1].width = pl.points[pl.nPoints-1].objectives[i_n-1];
	stacks[pl.nPoints-1][1].index = pl.nPoints;
}

int i_binarySearch(POINT p, int d)
{
	int min = 0;
	int max = fsorted[d].nPoints-1;
	gorder = torder[d];
	while (min<=max) {
		int mid = (max+min)/2;
		if (i_same(&p,&(fsorted[d].points[mid]))) {
			return mid;
		}
		else if (i_greaterorder(&p,&fsorted[d].points[mid])==-1) {
			max = mid-1;
		}
		else {
			min = mid+1;
		}
	}
	print_error(1,1,"EE: Cannot find the point in binary search");
	int i,j;
	printf("p:");
	for(i = 0; i< number_objective; i++)
		printf("%lf ",p.objectives[i]);
	printf("\n");

	printf("ps:");
	for(j = 0 ; j< fsorted[d].nPoints;j++) {
		for (i = 0; i < number_objective; i++)
			printf("%lf ", fsorted[d].points[j].objectives[i]);
		printf("\n");
	}

	return -1;
}

void i_runHeuristic(FRONT ps)
{
	int i,j,k;
	for ( i=0; i<i_n-1; i++) {
		torder[i][i_n-1] = i;
		torder[i][i] = i_n-1;
	}
	for ( i=i_n-1; i>=0; i--) {
		for ( j=0; j<ps.nPoints; j++) {
			fsorted[i].points[j] = ps.points[j];
			tcompare[j][i] = 0;
		}
		fsorted[i].nPoints = ps.nPoints;
		gorder = torder[i];
		qsort(fsorted[i].points,ps.nPoints,sizeof(POINT),i_greaterorder);
	}
	for ( i=0; i<ps.nPoints; i++) {
		for ( k=0; k<i_n; k++) {
			tcompare[i][k] = ps.nPoints - 1 - i_binarySearch(ps.points[i],k);

		}
	}
	for ( i=0; i<ps.nPoints; i++) {
		for ( j=1; j<i_n; j++) {
			int x = torder[i][j];
			k = j;
			while (k>0 && tcompare[i][x] < tcompare[i][torder[i][k-1]]) {
				torder[i][k] = torder[i][k-1];
				k--;
			}
			torder[i][k] = x;
		}
	}
}

int i_slicingDepth(int d)
{
	if (d <=  5) return 1;
	if (d <=  7) return 2;
	if (d <= 12) return 3;
	return 4;
}

void i_ihv2(FRONT ps, double *min)
{
  // returns the minimum exclusive hypervolume of points in ps for the 2D case
  qsort(ps.points,ps.nPoints,sizeof(POINT),i_greater);
  double vol = ps.points[0].objectives[0] * (ps.points[0].objectives[1] - ps.points[1].objectives[1]);
  double kvol;
  int sm = 0;
	int k;
  for ( k = 1; k < ps.nPoints - 1; k++)
    {
      kvol = (ps.points[k].objectives[0] - ps.points[k-1].objectives[0]) * (ps.points[k].objectives[1] - ps.points[k+1].objectives[1]);
      if (kvol <= 0)
	{
	  min[0] = ps.points[k].objectives[0];
	  min[1] = ps.points[k].objectives[1];
	  min[2] = 0;

	  return;
	}
      else
      if (kvol < vol)
	{
	  vol = kvol;
	  sm = k;
	}
    }
  kvol = (ps.points[ps.nPoints - 1].objectives[0] - ps.points[ps.nPoints - 2].objectives[0]) * ps.points[ps.nPoints - 1].objectives[1];
  if (kvol < vol)
    {
      vol = kvol;
      sm = ps.nPoints - 1;
    }
  min[0] = ps.points[sm].objectives[0];
  min[1] = ps.points[sm].objectives[1];
  min[2] = vol;

}



void i_ihv(FRONT ps, double *min)
// returns the minimum exclusive hypervolume of points in ps
{
	int i,j,k,z;
    int maxStackSize = MIN(i_slicingDepth(i_n),i_n-2)+1;   // ****
    for ( i=0; i<MAX(ps.nPoints,i_n); i++)
        for ( j=0; j<i_n; j++)
            torder[i][j] = j;
    i_runHeuristic(ps);  // ****
    for ( i=0; i<ps.nPoints; i++) {
        stacks[i][0].front = fsorted[torder[i][i_n-1]];
        stacks[i][0].width = 1;
        stacksize[i] = 2;
    }
    i_sliceOrder(ps.nPoints); // ****
    i_n--;
    for ( i=0; i<ps.nPoints; i++) {
        SLICE top = stacks[i][stacksize[i]-1];
        while (stacksize[i] < maxStackSize && top.front.nPoints > SLICELIMIT) {
            stacks[i][stacksize[i]].front.nPoints = 0;
            int index = 0;
            while (index < top.front.nPoints && ps.points[i].objectives[torder[i][i_n-1]] < top.front.points[index].objectives[torder[i][i_n-1]]) {
                i_insert(top.front.points[index],i_n-2,stacks[i][stacksize[i]].front,i,stacksize[i],torder[i]);
                index++;
            }
            if (index < top.front.nPoints) {
                stacks[i][stacksize[i]].width = ps.points[i].objectives[torder[i][i_n-1]]-top.front.points[index].objectives[torder[i][i_n-1]];
            }
            else {
                stacks[i][stacksize[i]].width = ps.points[i].objectives[torder[i][i_n-1]];
            }
            stacks[i][stacksize[i]].index = index;
            top = stacks[i][stacksize[i]];
            stacksize[i]++;
            i_n--;
        }
        double width = 1;
        for ( j=0; j<stacksize[i]; j++) {
            width *= stacks[i][j].width;
        }
        if (top.front.nPoints == 0) {
            partial[i] = width * i_inclhvOrder(ps.points[i],torder[i]); // ****
        }
        else {
            partial[i] = width * i_exclhvPoint(top.front,ps.points[i],torder[i]);
        }
        i_n += stacksize[i]-2;
        while (stacksize[i]>1 && (top.index==stacks[i][stacksize[i]-2].front.nPoints ||
                                  i_dominates1wayOrder(stacks[i][stacksize[i]-2].front.points[top.index],ps.points[i],i_n-stacksize[i]+1,torder[i]))) {
            stacksize[i]--;
            top = stacks[i][stacksize[i]-1];
        }
    }
    i_initialiseHeap(ps.nPoints);
    maxStackSize = 2;
    while (true) {
        // int i = removeFromHeap();
		i = i_peekFromHeap();
        if (stacksize[i]<=1) {
            for ( z=0; z<ps.n; z++) {
                min[z] = ps.points[i].objectives[z];
            }
            min[ps.n] = partial[i];
            break;
        }
        i_n -= stacksize[i]-2;
        j = stacks[i][stacksize[i]-1].index;
        if (j<stacks[i][stacksize[i]-2].front.nPoints-1) {
            stacks[i][stacksize[i]-1].width = stacks[i][stacksize[i]-2].front.points[j].objectives[torder[i][i_n]] -
                                              stacks[i][stacksize[i]-2].front.points[j+1].objectives[torder[i][i_n]];
        }
        else {
            stacks[i][stacksize[i]-1].width = stacks[i][stacksize[i]-2].front.points[j].objectives[torder[i][i_n]];
        }
        i_insert(stacks[i][stacksize[i]-2].front.points[j],i_n-1,stacks[i][stacksize[i]-1].front,i,stacksize[i]-1,torder[i]);
        stacks[i][stacksize[i]-1].index = j+1;
        SLICE top = stacks[i][stacksize[i]-1];

        double width = 1;
        for (k=0; k<stacksize[i]; k++) {
            width *= stacks[i][k].width;
        }
        if (top.front.nPoints == 0) {
            partial[i] += width * i_inclhvOrder(ps.points[i],torder[i]);
        }
        else {
            partial[i] += width * i_exclhvPoint(top.front,ps.points[i],torder[i]);
        }
        i_n += stacksize[i]-2;
        while (stacksize[i]>1 && (top.index==stacks[i][stacksize[i]-2].front.nPoints ||
                                  i_dominates1wayOrder(stacks[i][stacksize[i]-2].front.points[top.index],ps.points[i],i_n-stacksize[i]+1,torder[i]))) {
            stacksize[i]--;
            top = stacks[i][stacksize[i]-1];
        }

        i_heapify(0,i);
    }
    i_n++;
}


/*
int main(int argc, char *argv[])
// processes each front from the file
{
  if (argc == 1)
    {
      printf("There are two modes of use:\n");
      printf("- exe filename r1 ... rn uses the reference point (r1, ..., rn)\n");
      printf("- exe filename           uses the reference point (0,  ..., 0)\n");
      printf("No checks are performed on the data or the reference point.\n");
    }
  else
    {
	FILECONTENTS *f = readFile(argv[1]);

	// find the biggest fronts
	maxm = 0;
	maxn = 0;
	for (int i = 0; i < f->nFronts; i++) {
		if (f->fronts[i].nPoints > maxm) {
			maxm = f->fronts[i].nPoints;
		}
		if (f->fronts[i].n       > maxn) {
			maxn = f->fronts[i].n;
		}
	}

	// allocate memory
	int maxdepth = maxn - 2;
	fs = malloc(sizeof(FRONT) * maxdepth);
	for (int i = 0; i < maxdepth; i++) {
		fs[i].points = malloc(sizeof(POINT) * maxm);
		for (int j = 0; j < maxm; j++) {
			fs[i].points[j].objectives = malloc(sizeof(OBJECTIVE) * (maxn - i - 1));
		}
	}
	partial = malloc(sizeof(double) * maxm);
	heap = malloc(sizeof(int) * maxm);
	stacksize = malloc(sizeof(int) * maxm);
	stacks = malloc(sizeof(SLICE*) * maxm);
	int maxStackSize = MIN(maxn-2,i_slicingDepth(maxn))+1;
	for (int i=0; i<maxm; i++) {
		stacks[i] = malloc(sizeof(SLICE) * maxStackSize);
		for (int j=1; j<maxStackSize; j++) {
			stacks[i][j].front.points = malloc(sizeof(POINT) * maxm);
		}
	}

	fsorted = malloc(sizeof(FRONT) * maxn);
	for (int i=0; i<maxn; i++) {
		fsorted[i].points = malloc(sizeof(POINT) * maxm);
	}
	torder = malloc(sizeof(int*) * MAX(maxm,maxn));
	tcompare = malloc(sizeof(int*) * maxm);
	for (int i=0; i<MAX(maxn,maxm); i++) {
		torder[i] = malloc(sizeof(int) * maxn);
	}
	for (int i=0; i<maxm; i++) {
		tcompare[i] = malloc(sizeof(int) * maxn);
	}

	// initialise the reference point
	ref.objectives = malloc(sizeof(OBJECTIVE) * maxn);
	if (argc == 2) {
		printf("No reference point provided: using the origin\n");
		for (int i = 0; i < maxn; i++) {
			ref.objectives[i] = 0;
		}
	}
	else if (argc - 2 != maxn) {
		printf("Your reference point should have %d values\n", maxn);
		return 0;
	}
	else {
		for (int i = 2; i < argc; i++) {
			ref.objectives[i - 2] = atof(argv[i]);
		}
	}

	// record the directions of the objectives, so we can restore the smallest point later
	// TRUE indicates maximisation
	// ASSUMES that the first point in the first front has no zero objectives
	dirs.objectives = malloc(sizeof(OBJECTIVE) * maxn);
	for (int k = 0; k < maxn; k++)
	  dirs.objectives[k] = f->fronts[0].points[0].objectives[k] > ref.objectives[k];
	// modify the objective values relative to the reference point
	for (int i = 0; i < f->nFronts; i++)
		for(int j = 0; j < f->fronts[i].nPoints; j++)
			for(int k = 0; k < f->fronts[i].n; k++)
				f->fronts[i].points[j].objectives[k] = fabs(f->fronts[i].points[j].objectives[k] - ref.objectives[k]);

	// printContents(f);

	totaltime = 0;
	for (int i = 0; i < f->nFronts; i++) {
		struct timeval tv1, tv2;
		struct rusage ru_before, ru_after;
		getrusage (RUSAGE_SELF, &ru_before);

		n = f->fronts[i].n;
		double eh[n+1];
		if(n == 2)
		  ihv2(f->fronts[i],eh);
		else
		  ihv(f->fronts[i],eh);
		printf("mehv(%d) = %1.16f\n", i+1, eh[n]);
		printf("Smallest: ");
		printPoint(eh);

		getrusage (RUSAGE_SELF, &ru_after);
		tv1 = ru_before.ru_utime;
		tv2 = ru_after.ru_utime;
		// printf("Time: %f (s)\n", tv2.tv_sec + tv2.tv_usec * 1e-6 - tv1.tv_sec - tv1.tv_usec * 1e-6);
		totaltime += tv2.tv_sec + tv2.tv_usec * 1e-6 - tv1.tv_sec - tv1.tv_usec * 1e-6;
	}
	printf("Total time = %f (s)\n", totaltime);
    }
	return 0;
}
*/