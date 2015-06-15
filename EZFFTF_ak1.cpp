// EZFFTF_ak1.cpp - Program for computing the Fourier coefficients of a real periodic sequence (Fourier analysis.)
// Written in Microsoft Visual Studio Express 2013 for Windows Desktop
// 15 June 2015
//
// The sub-routines listed below are translations of FORTRAN routines included in FFTPACK, posted off the NETLIB site:
//
// http://www.netlib.org/fftpack/
//
// Specifically, EZFFTF.FOR and its dependencies, written by Paul N. Swarztrauber, National Center for Atmospheric Research, Boulder, CO.
//
// To distinguish the routines posted below from others, an _ak1 suffix has been appended to them.
//
// Following is a list of the major changes made in the course of translating the FFTPACK routines
// to the C++ versions posted below:
// 1) All global variables have been eliminated.
// 2) All "GO TO" statements have been eliminated.
//
// A small main program is included also, to provide an example of how to use EZFFTF_ak1. In this 
// example, data is input from a file to eliminate the need for a user to type data in via
// the console.

#include <iostream>
#include <fstream>
#include <cctype>
#include <cmath>
#include <vector>
#include <cfloat>

using namespace std;

#define PI 3.141592653589793238462643383279502884197169399375158209749445923
#define FACARSIZE 15

typedef vector<double> C1DArray;

void EZFFTI_ak1(int N, C1DArray& WA, int ifac_array[FACARSIZE]);
void RFFTF_ak1(int N, C1DArray& WA1, C1DArray& WA2, C1DArray& r2Ar, int ifac_array[FACARSIZE]);
void RADF2_ak1(int ID0, int L1, C1DArray CC, C1DArray& CH, double* WA1);
void RADF3_ak1(int ID0, int L1, C1DArray CC, C1DArray& CH, double* WA1, double* WA2);
void RADF4_ak1(int ID0, int L1, C1DArray CC, C1DArray& CH, double* WA1, double* WA2, double* WA3);
void RADF5_ak1(int ID0, int L1, C1DArray CC, C1DArray& CH, double* WA1, double* WA2, double* WA3, double* WA4);
void RADFG_ak1(int ID0, int IP, int L1, int IDL1, C1DArray& CC, C1DArray& C1, C1DArray& C2, C1DArray& CH, C1DArray& CH2, double* WA1);

void EZFFTI_ak1(int N, C1DArray& WA, int ifac_array[FACARSIZE]){
	// EZFFTI initializes the WA and IFAC arrays. 
	// The tabulation of the trigonometric functions are computed and stored in WA.
	// The prime factorization of N is stored in ifac_array.

	if (N > 1) {
		double arg1, ch1, ch1h, dch1, dsh1, sh1;
		double argh = 2.0 * PI / (double) N;
		int i, id0, ii, ip, ipm, is = 0, j = 0, k1 = 0, l1 = 1, l2, nf = 0, nfm1, nl = N, nq, nr, ntry, ntryh[4] = { 4, 2, 3, 5 };

		do {
			if (j < 4)  ntry = ntryh[j++];
			else  ntry += 2;
			do {
				nq = nl / ntry;
				nr = nl - ntry*nq;
				if (nr != 0) break;
				++nf;
				ifac_array[nf + 1] = ntry;
				nl = nq;
				if ((ntry == 2) && (nf != 1)) {
					for (i = nf; i > 1; --i) ifac_array[i + 1] = ifac_array[i];
					ifac_array[2] = 2;
				}  // End if ((ntry == 2) && (nf != 1))
			} while (nl != 1);
	    } while (nr != 0);

		ifac_array[0] = N;
		ifac_array[1] = nf;
		nfm1 = nf - 1;

		while (k1 < nfm1){
			ip = ifac_array[2 + k1++];
			l2 = l1*ip;
			id0 = N / l2;
			ipm = ip - 1;
			arg1 = (double)l1 * argh;
			ch1 = 1.0;
			sh1 = 0.0;
			dch1 = cos(arg1);
			dsh1 = sin(arg1);
			for (j = 0; j < ipm; ++j){
				ch1h = dch1*ch1 - dsh1*sh1;
				sh1 = dch1*sh1 + dsh1*ch1;
				ch1 = ch1h;
				i = is;
				WA[i++] = ch1;
				WA[i++] = sh1;
				ii = 4;
				while (ii < id0){
					WA[i] = ch1*WA[i - 2] - sh1*WA[i - 1];
					++i;
					WA[i] = ch1*WA[i - 2] + sh1*WA[i - 3];
					++i;
					ii += 2;
				} // End while (ii < id0)
				is += id0;
			}  // End for j
			l1 = l2;
		}  // End while (k1 < nfm1){

	} // End if (N > 1)
	return;
} // End EZFFTI_ak1

void RFFTF_ak1(int N, C1DArray& WA1, C1DArray& WA2, C1DArray& r2Ar, int ifac_array[FACARSIZE]){
	// RFFTF computes the Fourier coefficients of a real periodic sequence.

	if (N > 1) {

		bool na = 1;
		int iw = N - 1, l2 = N, nf = ifac_array[1], idl1, id0, ip, ix2, ix3, ix4, kh = nf + 1, l1;

		for (int k1 = 0; k1 < nf; ++k1){
			ip = ifac_array[kh--];
			l1 = l2 / ip;
			id0 = N / l2;
			iw -= (ip - 1)*id0;
			na = !na;

			switch (ip){
				case 4:	ix2 = iw + id0;
						ix3 = ix2 + id0;
						if (!na) RADF4_ak1(id0, l1, r2Ar, WA2, &WA1[iw], &WA1[ix2], &WA1[ix3]);
						else RADF4_ak1(id0, l1, WA2, r2Ar, &WA1[iw], &WA1[ix2], &WA1[ix3]);
						break;
				case 2:	if (na) RADF2_ak1(id0, l1, WA2, r2Ar, &WA1[iw]);
						else RADF2_ak1(id0, l1, r2Ar, WA2, &WA1[iw]);
						break;
				case 3:	ix2 = iw + id0;
						if (na) RADF3_ak1(id0, l1, WA2, r2Ar, &WA1[iw], &WA1[ix2]);
						else RADF3_ak1(id0, l1, r2Ar, WA2, &WA1[iw], &WA1[ix2]);
						break;
				case 5:	ix2 = iw + id0;
						ix3 = ix2 + id0;
						ix4 = ix3 + id0;
						if (na) RADF5_ak1(id0, l1, WA2, r2Ar, &WA1[iw], &WA1[ix2], &WA1[ix3], &WA1[ix4]);
						else RADF5_ak1(id0, l1, r2Ar, WA2, &WA1[iw], &WA1[ix2], &WA1[ix3], &WA1[ix4]);
						break;
				default:idl1 = id0*l1;
					    if (id0 == 1) na = !na;
						if (na) RADFG_ak1(id0, ip, l1, idl1, WA2, WA2, WA2, r2Ar, r2Ar, &WA1[iw]);
						else RADFG_ak1(id0, ip, l1, idl1, r2Ar, r2Ar, r2Ar, WA2, WA2, &WA1[iw]);
						na = !na;
			} // End switch
			l2 = l1;
		}  // End for k1

		if (!na) for (int i = 0; i < N; ++i) r2Ar[i] = WA2[i];

	} // End if (N > 1)
	return;
}  // End RFFTF_ak1

void RADF2_ak1(int ID0, int L1, C1DArray CC, C1DArray& CH, double* WA1){

	int cc_l_index, cc_u_index, ccZIncr = ID0*L1, twoido = 2*ID0, i, idij, idkk, idxz0, idz5, k = L1, l3_index, z_index = -1, zi0_index, z5_index;
	double dum1, tr2, ti2;

	cc_u_index = cc_l_index = -ID0;
	cc_u_index += ccZIncr;
	l3_index = -twoido;

	while (k){
		cc_l_index += ID0;
		cc_u_index += ID0;
		l3_index += twoido;
		//z_index = 2*ID0 - 1 + k*ccZIncr;
		z_index += twoido;
		CH[l3_index] = CC[cc_l_index] + CC[cc_u_index];
		CH[z_index] = CC[cc_l_index] - CC[cc_u_index];
		--k;
	} // End while

	if (ID0 >= 2){
		if (ID0 > 2) {

			idij = -ID0 + ccZIncr;
			idkk = -1;
			idxz0 = -ID0;
			idz5 = -twoido + 1;

			for (k = 0; k < L1; ++k){

				cc_l_index = -1;
				cc_u_index = idij = idij + ID0;
				z_index = idkk = idkk + twoido;
				zi0_index = idxz0 = idxz0 + ID0;
				z5_index = idz5 = idz5 + twoido;
				
				//cc_u_index = (k + L1)*ID0 + 1;
				//z_index = (1 + k) * twoido - 2;
				//zi0_index = 1 + k*ID0;
				//z5_index = 2 + k*twoido;

				i = (ID0 - 1) / 2;
				while (i){
					++cc_l_index;
					++cc_u_index;
					++zi0_index;
					++z5_index;
					--z_index;
					dum1 = CC[cc_u_index + 1];
					ti2 = WA1[cc_l_index + 1];
					tr2 = WA1[cc_l_index] * CC[cc_u_index] + ti2 * dum1;
					ti2 = WA1[cc_l_index++] * dum1 - ti2 * CC[cc_u_index++];
					dum1 = CC[zi0_index++];
					CH[z5_index] = CC[zi0_index] + ti2;
					CH[z_index--] = ti2 - CC[zi0_index];
					CH[z5_index++ - 1] = dum1 + tr2;
					CH[z_index] = dum1 - tr2;
					--i;
				} // End while
				
			}  // End for k
		} // End if (ID0 > 2)

		if ((ID0 % 2) == 0){

			cc_l_index = -1;
			cc_u_index = ccZIncr - 1;
			z_index = -(1 + ID0);

			//z_index = cc_l_index = ID0 - 1;
			//cc_u_index = cc_l_index + ccZIncr;

			k = L1;
			while (k){
				cc_l_index += ID0;
				cc_u_index += ID0;
				z_index += twoido;
				CH[z_index + 1] = -CC[cc_u_index];
				CH[z_index] = CC[cc_l_index];
				--k;
			}  // End while

		}  // End if ((ID0 % 2) == 0)

	} // End if (ID0 >= 2)

	return;
}  // End RADF2_ak1

void RADF3_ak1(int ID0, int L1, C1DArray CC, C1DArray& CH, double* WA1, double* WA2){

	int cc_l_index, cc_u_index, ccZIncr = ID0*L1, chzincr = 3*ID0, i, idkk, idm2, incr3row, izero2, k, l3_index, twoido = 2*ID0, z_index = 0, zi0_index, z4_index, z5_index;

	double ci2, cr2, di2, dr2, di3, dr3, ti2, tr2, ti3, tr3;

	static double taui = 0.8660254037844386467637231707529361834710262690519031402790348975;
	static double taur = -0.5;

	l3_index = cc_l_index = -ID0;
	l3_index += ccZIncr;
	cc_u_index = l3_index + ccZIncr;
	//cc_u_index = l3_index*2;

	z4_index = zi0_index = -ID0;
	z_index = -chzincr;
	--z4_index;

	//z4_index = zi0_index - 1;

	for (k = 0; k < L1; k++){

		cc_l_index += ID0;
		l3_index += ID0;
		cc_u_index += ID0;
		zi0_index += chzincr;
		z_index += chzincr;
		z4_index += chzincr;

		cr2 = CC[l3_index] + CC[cc_u_index];
		CH[z_index] = CC[cc_l_index] + cr2;
		CH[zi0_index] = taui*(CC[cc_u_index] - CC[l3_index]);
		CH[z4_index] = CC[cc_l_index] + taur*cr2;

	}  // End for k

	if (ID0 != 1) {

		idm2 = ID0 - 2;
		izero2 = -ID0;
		idkk = -chzincr;

		for (k = 0; k < L1; k++){

			z4_index = izero2 = izero2 + ID0;

			//z4_index = 1 + k*ID0;

			cc_u_index = cc_l_index = z4_index + ccZIncr;
			++cc_u_index;
			z_index = l3_index = cc_l_index + ccZIncr;
			++z_index;

			zi0_index = idkk = idkk + chzincr;

			//zi0_index = 1 + k * 3 * ID0;

			z5_index = incr3row = zi0_index + twoido;
			--z5_index;

			//z5_index = 2 * (ID0 - 1) + k * 3 * ID0;

			for (i = 0; i < idm2; i += 2){

				++cc_l_index;
				++cc_u_index;
				++z4_index;
				++l3_index;
				++z_index;
				++zi0_index;
				++incr3row;
				--z5_index;

				cr2 = WA1[i + 1];
				ci2 = WA2[i + 1];
				dr2 = WA1[i] * CC[cc_l_index] + cr2 * CC[cc_u_index];
				di2 = WA1[i] * CC[cc_u_index++] - cr2 * CC[cc_l_index++];
				dr3 = WA2[i] * CC[l3_index] + ci2 * CC[z_index];
				di3 = WA2[i] * CC[z_index++] - ci2 * CC[l3_index++];
				cr2 = dr2 + dr3;
				ci2 = di2 + di3;

				tr2 = CC[z4_index++];
				CH[zi0_index++] = tr2 + cr2;
				CH[zi0_index] = CC[z4_index] + ci2;

				tr2 += taur*cr2;
				ti2 = CC[z4_index] + taur*ci2;
				tr3 = taui*(di2 - di3);
				ti3 = taui*(dr3 - dr2);

				CH[incr3row++] = tr2 + tr3;
				CH[z5_index--] = ti3 - ti2;
				CH[z5_index] = tr2 - tr3;
				CH[incr3row] = ti2 + ti3;

			}  // End for i
		}  // End for k
	} // End if (ID0 != 1)

	return;
}  // End RADF3_ak1

void RADF4_ak1(int ID0, int L1, C1DArray CC, C1DArray& CH, double* WA1, double* WA2, double* WA3){

	int cc_l_index, cc_u_index, i, ccZIncr = L1*ID0, twoido = 2*ID0, idm2, z4_index, incr3row, izero2, k, l3_index, l4_index, chzincr = 2*twoido, z_index, zi0_index, z3_index;

	double cr2, ci2, cr3, ci3, cr4, ci4, tr1, tr2, tr4, ti1, ti4, ti2, ti3, tr3;
	static double HSQT2 = 0.70710678118654752440084436210484903928483593768474036588339869;

	cc_l_index = -ID0;
	l3_index = cc_l_index + ccZIncr;
	l4_index = l3_index + ccZIncr;
	cc_u_index = l4_index + ccZIncr;

	z_index = -chzincr;
	z3_index = z4_index = z_index + twoido;
	--z3_index;
	zi0_index = -1;

	for (k = 0; k < L1; k++){

		cc_l_index += ID0;
		cc_u_index += ID0;
		l3_index += ID0;
		l4_index += ID0;

		z_index += chzincr;
		z4_index += chzincr;
		zi0_index += chzincr;
		z3_index += chzincr;

		tr1 = CC[l3_index] + CC[cc_u_index];
		tr2 = CC[cc_l_index] + CC[l4_index];
		CH[z_index] = tr1 + tr2;
		CH[zi0_index] = tr2 - tr1;
		CH[z3_index] = CC[cc_l_index] - CC[l4_index];
		CH[z4_index] = CC[cc_u_index] - CC[l3_index];

	} // End for k
	
	if (ID0 >= 2){
		if (ID0 > 2) {

			idm2 = ID0 - 2;
			izero2 = -ID0;
			incr3row = -chzincr;

			for (k = 0; k < L1; k++){

				//cc_l_index = 1 + k*ID0;
				cc_l_index = izero2 = izero2 + ID0;

				l3_index = cc_l_index + ccZIncr;
				l4_index = l3_index + ccZIncr;
				cc_u_index = l4_index + ccZIncr;

				zi0_index = incr3row = incr3row + chzincr;

				//z_index = zi0_index = 1 + k*chzincr;
				z_index = zi0_index + twoido;

				//z3_index = zi0_index + twoido - 1;
				z3_index = z_index - 1;
				z4_index = z3_index + twoido;

				for (i = 0; i < idm2; i += 2){

					++cc_l_index;
					++l3_index;
					++l4_index;
					++cc_u_index;

					++zi0_index;
					++z_index;
					--z4_index;
					--z3_index;

					tr2 = WA1[i + 1];
					tr3 = WA2[i + 1];
					tr4 = WA3[i + 1];

					ci2 = CC[l3_index++];
					cr2 = WA1[i] * ci2 + tr2*CC[l3_index];
					ci2 = WA1[i] * CC[l3_index] - tr2*ci2;

					ci3 = CC[l4_index++];
					cr3 = WA2[i] * ci3 + tr3*CC[l4_index];
					ci3 = WA2[i] * CC[l4_index] - tr3*ci3;

					ci4 = CC[cc_u_index++];
					cr4 = WA3[i] * ci4 + tr4*CC[cc_u_index];
					ci4 = WA3[i] * CC[cc_u_index] - tr4*ci4;

					tr1 = cr2 + cr4;
					tr4 = cr4 - cr2;

					ti1 = ci2 + ci4;
					ti4 = ci2 - ci4;

					tr2 = CC[cc_l_index] + cr3;
					tr3 = CC[cc_l_index++] - cr3;

					ti2 = CC[cc_l_index] + ci3;
					ti3 = CC[cc_l_index] - ci3;

					CH[zi0_index++] = tr1 + tr2;
					CH[zi0_index] = ti1 + ti2;
					CH[z4_index--] = ti1 - ti2;
					CH[z4_index] = tr2 - tr1;
					CH[z_index++] = ti4 + tr3;
					CH[z_index] = tr4 + ti3;
					CH[z3_index--] = tr4 - ti3;
					CH[z3_index] = tr3 - ti4;
					//if (i == 24){
					//	cout << "\ni = 24. \n";
					//}
				} // End for i

			}  // End for k

		} // End if (ID0 > 2)
		
		if ((ID0 % 2) == 0){

			//cc_u_index = z3_index = zi0_index = cc_l_index = ID0 - 1;
			cc_l_index = -1;

			cc_u_index = cc_l_index + ccZIncr;
			l3_index = cc_u_index + ccZIncr;
			l4_index = l3_index + ccZIncr;

			zi0_index = z4_index = -chzincr + ID0;
			--zi0_index;
			z3_index = z_index = z4_index + twoido;
			--z3_index;

			//zi0_index = ID0 - 1 - chzincr;
			//z_index = twoido + ID0 - chzincr;
			//z3_index += z_index;
			//z4_index = ID0 - chzincr;
			//z_index += ID0;

			for (k = 0; k < L1; k++){

				cc_l_index += ID0;
				cc_u_index += ID0;
				l3_index += ID0;
				l4_index += ID0;

				zi0_index += chzincr;
				z_index += chzincr;
				z4_index += chzincr;
				z3_index += chzincr;

				ti1 = -HSQT2*(CC[cc_u_index] + CC[l4_index]);
				tr1 = HSQT2*(CC[cc_u_index] - CC[l4_index]);
				CH[zi0_index] = tr1 + CC[cc_l_index];
				CH[z3_index] = CC[cc_l_index] - tr1;
				CH[z4_index] = ti1 - CC[l3_index];
				CH[z_index] = ti1 + CC[l3_index];

			} // End for k

		}  // End if ((ID0 % 2) == 0)

	} // End if (ID0 >= 2)

	return;
}  // End RADF4_ak1

void RADF5_ak1(int ID0, int L1, C1DArray CC, C1DArray& CH, double* WA1, double* WA2, double* WA3, double* WA4){

	int cc_l_index, cc_u_index, ccZIncr = ID0*L1, chzincr = ID0 * 5, idm2, izero2, l2_index, l3_index, l4_index, twoido, z_index, zi0_index, zi2, z4, z5, z7, i, k;

	double cr2, ci2, ci3, ci5, cr3, ci4, cr4, cr5, dr2, di2, dr3, di3, dr4, di4, dr5, di5, tr2, ti2, tr3, ti3, tr5, ti5, tr4, ti4;
	
	static double tr11 = 0.3090169943749474241022934171828195886015458990288143106772431137;
	static double ti11 = 0.9510565162951535721164393337938214340569863412575022244730564442;
	static double tr12 = -0.8090169943749474241022934171828190588601545899028814310677431135;
	static double ti12 = 0.5877852522924731291687059546390727685976524376431459107227248076;

	z7 = z4 = z5 = cc_l_index = -ID0;
	--z4;
	z7 += ccZIncr;
	l3_index = z7 + ccZIncr;
	l4_index = l3_index + ccZIncr;
	cc_u_index = l4_index + ccZIncr;

	z_index = -chzincr;
	zi0_index = zi2 = -(3 * ID0);
	--zi0_index;

	for (k = 0; k < L1; k++){

		cc_l_index += ID0;
		z7 += ID0;
		l3_index += ID0;
		l4_index += ID0;
		cc_u_index += ID0;

		z_index += chzincr;
		zi2 += chzincr;
		zi0_index += chzincr;
		z4 += chzincr;
		z5 += chzincr;

		cr2 = CC[cc_u_index] + CC[z7];
		ci5 = CC[cc_u_index] - CC[z7];
		cr3 = CC[l4_index] + CC[l3_index];
		ci4 = CC[l4_index] - CC[l3_index];
		CH[z_index] = CC[cc_l_index] + cr2 + cr3;
		CH[zi0_index] = CC[cc_l_index] + tr11*cr2 + tr12*cr3;
		CH[zi2] = ti11*ci5 + ti12*ci4;
		CH[z4] = CC[cc_l_index] + tr12*cr2 + tr11*cr3;
		CH[z5] = ti12*ci5 - ti11*ci4;

	}  // End for k

	if (ID0 != 1){

		idm2 = ID0 - 2;
		twoido = 2 * ID0;
		izero2 = -ID0;
		z7 = -chzincr;

		for (k = 0; k < L1; k++){

			cc_l_index = izero2 = izero2 + ID0;
			l2_index = cc_l_index + ccZIncr;
			l3_index = l2_index + ccZIncr;
			l4_index = l3_index + ccZIncr;
			cc_u_index = l4_index + ccZIncr;

			zi0_index = z7 = z7 + chzincr;
			z5 = z_index = zi0_index + twoido;
			zi2 = z_index + twoido;
			--z5;
			z4 = z5 + twoido;

			for (i = 0; i < idm2; i += 2){

				++cc_l_index;
				++l2_index;
				++l3_index;
				++l4_index;
				++cc_u_index;
				++zi0_index;
				++z_index;
				++zi2;
				--z4;
				--z5;

				tr2 = WA1[i + 1];
				tr3 = WA2[i + 1];
				tr4 = WA3[i + 1];
				tr5 = WA4[i + 1];

				di2 = CC[l2_index++];
				dr2 = WA1[i] * di2 + tr2*CC[l2_index];
				di2 = WA1[i] * CC[l2_index] - tr2*di2;

				di3 = CC[l3_index++];
				dr3 = WA2[i] * di3 + tr3*CC[l3_index];
				di3 = WA2[i] * CC[l3_index] - tr3*di3;

				di4 = CC[l4_index++];
				dr4 = WA3[i] * di4 + tr4*CC[l4_index];
				di4 = WA3[i] * CC[l4_index] - tr4*di4;

				di5 = CC[cc_u_index++];
				dr5 = WA4[i] * di5 + tr5*CC[cc_u_index];
				di5 = WA4[i] * CC[cc_u_index] - tr5*di5;

				cr2 = dr2 + dr5;
				ci5 = dr5 - dr2;
				cr5 = di2 - di5;
				ci2 = di2 + di5;
				cr3 = dr3 + dr4;
				ci4 = dr4 - dr3;
				cr4 = di3 - di4;
				ci3 = di3 + di4;

				tr5 = CC[cc_l_index++];
				CH[zi0_index++] = tr5 + cr2 + cr3;
				CH[zi0_index] = CC[cc_l_index] + ci2 + ci3;
				tr2 = tr5 + tr11 * cr2 + tr12*cr3;
				ti2 = CC[cc_l_index] + tr11*ci2 + tr12*ci3;
				tr3 = tr5 + tr12*cr2 + tr11*cr3;
				ti3 = CC[cc_l_index] + tr12*ci2 + tr11*ci3;
				tr5 = ti11*cr5 + ti12*cr4;
				ti5 = ti11*ci5 + ti12*ci4;
				tr4 = ti12*cr5 - ti11*cr4;
				ti4 = ti12*ci5 - ti11*ci4;
				
				CH[z_index++] = tr2 + tr5;
				CH[z_index] = ti2 + ti5;

				CH[zi2++] = tr3 + tr4;
				CH[zi2] = ti3 + ti4;

				CH[z4--] = ti4 - ti3;
				CH[z4] = tr3 - tr4;

				CH[z5--] = ti5 - ti2;
				CH[z5] = tr2 - tr5;

			} // End for i

		}  // End for k

	} // End if (ID0 != 1)

	return;
}  // End RADF5_ak1

void RADFG_ak1(int ID0, int IP, int L1, int IDL1, C1DArray& CC, C1DArray& C1, C1DArray& C2, C1DArray& CH, C1DArray& CH2, double* WA1){

	int cc_l_index, cc_u_index, i, idij, idkk, j, k, ipm1, idp2, idxj, incr3row = IP*ID0, ipph, ipp2, is, izero2, l3_index, nbd, twoido = 2*ID0, z_index;

	double AR1, ar2, AI1, ai2, ar1h, ar2h, ARG_G, dcp, dc2, ds2, dsp;

	static double TWO_PI = 2.0*PI;

	ARG_G = TWO_PI / (double) IP;
	dcp = cos(ARG_G);
	dsp = sin(ARG_G);
	//idm2 = ID0 - 2;
	ipm1 = IP - 1;
	ipph = ((IP + 1) / 2) - 1;
	nbd = (ID0 - 1) / 2;

	if (ID0 != 1){

		for (i = 0; i < IDL1; i++)  CH2[i] = C2[i];

		is = -ID0;
		for (j = 0; j < ipm1; j++){  // LOOP 103 IN FORTRAN CODE

			cc_l_index = is = is + IDL1;

			for (k = 0; k < L1; k++){
				cc_l_index += ID0;
				CH[cc_l_index] = C1[cc_l_index];
			}  // End for k
		}  //  End for j

		is = -(ID0 + 1);

		if (nbd <= L1){

			ipp2 = is;

			for (j = 0; j < ipm1; j++){

				idij = is = is + ID0;
				idp2 = ipp2 = ipp2 + IDL1;

				//for (i = 0; i < idm2; i += 2){
				for (i = 0; i < nbd; i++){

					++idij;
					ai2 = WA1[idij++];
					ds2 = WA1[idij];

					++idp2;
					++idp2;
					cc_l_index = idp2;
					//cc_l_index = 1 + (1 + j)*IDL1 + i;

					for (k = 0; k < L1; k++){

						cc_l_index += ID0;

						dc2 = C1[cc_l_index];
						ar2 = C1[cc_l_index + 1];
						CH[cc_l_index] = ai2 * dc2 + ds2 * ar2;
						CH[cc_l_index + 1] = ai2 * ar2 - ds2 * dc2;

					}  // End for k
				}  // End for i
			}  // End for j

		}  // End if (nbd <= L1)
		else {  // else nbd > L1

			ipp2 = 0;

			for (j = 0; j < ipm1; j++){

				is += ID0;

				ipp2 = ipp2 + IDL1;
				idp2 = ipp2 - ID0;

				for (k = 0; k < L1; k++){

					idij = is;

					cc_l_index = idp2 = idp2 + ID0;
					//cc_l_index = 1 + (1 + j)*IDL1 + k*ID0;

					//for (i = 0; i < idm2; i += 2){
					for (i = 0; i < nbd; i++){

						++idij;
						++cc_l_index;

						ai2 = WA1[idij++];
						dc2 = C1[cc_l_index];
						ar2 = C1[cc_l_index + 1];

						CH[cc_l_index++] = ai2 * dc2 + WA1[idij] * ar2;
						CH[cc_l_index] = ai2 * ar2 - WA1[idij] * dc2;

					}  // End for i
				}  // End for k
			}  // End for j

		}  // End else nbd > L1

		ipp2 = IP*IDL1;

		if (nbd >= L1){  // LINE 111 IN FORTRAN CODE

			is = 0;

			for (j = 0; j < ipph; j++){  // LOOP 114 IN FORTRAN CODE

				is = is + IDL1;
				izero2 = is - ID0;

				ipp2 = ipp2 - IDL1;
				idp2 = ipp2 - ID0;

				for (k = 0; k < L1; k++){

					cc_l_index = izero2 = izero2 + ID0;
					cc_u_index = idp2 = idp2 + ID0;

					//cc_l_index = 1 + (1 + j)*IDL1 + k*ID0;
					//cc_u_index = 1 + (IP - 1 - j)*IDL1 + k*ID0;

					//for (i = 0; i < idm2; i += 2){
					for (i = 0; i < nbd; i++){

						++cc_l_index;
						++cc_u_index;

						ai2 = CH[cc_l_index];
						ar2 = CH[cc_u_index];
						dc2 = CH[cc_u_index + 1];

						C1[cc_l_index++] = ai2 + ar2;
						C1[cc_l_index] = CH[cc_l_index] + dc2;

						C1[cc_u_index++] = CH[cc_l_index] - dc2;
						C1[cc_u_index] = ar2 - ai2;
						
					}  // End for i
				}  // End for k
			}  // End for j

		}  // End if (nbd >= L1)
		else {  // else nbd < L1

			is = -(ID0 + 1);
			ipp2 += is;
			//--ipp2;
			//ipp2 -= ID0;

			for (j = 0; j < ipph; j++){  // LOOP 118 IN FORTRAN CODE

				izero2 = is = is + IDL1;
				idp2 = ipp2 = ipp2 - IDL1;

				//for (i = 0; i < idm2; i += 2){
				for (i = 0; i < nbd; i++){

					++izero2;
					++izero2;
					cc_l_index = izero2;
					++idp2;
					++idp2;
					cc_u_index = idp2;

					//cc_l_index = 1 + (1 + j)*IDL1 + i;
					//cc_u_index = 1 + (IP - 1 - j)*IDL1 + i;

					for (k = 0; k < L1; k++){

						cc_l_index += ID0;
						cc_u_index += ID0;

						ai2 = CH[cc_l_index];
						ar2 = CH[cc_u_index];
						dc2 = CH[cc_u_index + 1];

						C1[cc_l_index] = ai2 + ar2;
						C1[cc_l_index + 1] = CH[cc_l_index + 1] + dc2;

						C1[cc_u_index] = CH[cc_l_index + 1] - dc2;
						C1[cc_u_index + 1] = ar2 - ai2;

					}  // End for k
				}  // End for i
			}  // End for j

		}  // End else nbd < L1

	} // End if (ID0 != 1){
	else { // else ID0 == 1	// 119 in FORTRAN CODE
		for (i = 0; i < IDL1; i++) C2[i] = CH2[i];
	}  // End else ID0 == 1

	// 121 in FORTRAN CODE

	is = -ID0;
	izero2 = IP*IDL1 - ID0;

	for (j = 0; j < ipph; j++){

		cc_l_index = is = is + IDL1;
		cc_u_index = izero2 = izero2 - IDL1;

		//cc_l_index = (1 + j)*IDL1;
		//cc_u_index = (IP - 1 - j)*IDL1;

		for (k = 0; k < L1; k++){
			cc_l_index += ID0;
			cc_u_index += ID0;
			C1[cc_l_index] = CH[cc_l_index] + CH[cc_u_index];
			C1[cc_u_index] = CH[cc_u_index] - CH[cc_l_index];
		}  // End for k
	}  // End for j

	AR1 = 1.0;
	AI1 = 0.0;

	is = idxj = -1;
	is += ipm1*IDL1;
	idkk = is + IDL1;

	//idkk = (ipm1 + 1) * IDL1 - 1;

	for (i = 0; i < ipph; i++){
		ar1h = dcp*AR1 - dsp*AI1;
		AI1 = dcp*AI1 + dsp*AR1;
		AR1 = ar1h;

		cc_l_index = IDL1;
		cc_u_index = is + 1;
		l3_index = idxj = idxj + IDL1;
		z_index = idkk = idkk - IDL1;
		//l3_index = (1 + i)*IDL1;
		//z_index = (ipm1 - i)*IDL1;

		for (k = 0; k < IDL1; k++){
			++l3_index;
			++z_index;
			CH2[l3_index] = C2[k] + AR1*C2[cc_l_index++]; 
			CH2[z_index] = AI1*C2[cc_u_index++];
		}  // End for k

		ar2 = dc2 = AR1;
		ai2 = ds2 = AI1;

		idij = IDL1 - 1;
		izero2 = is;

		for (j = 0; j < (ipph - 1) ; j++){

			ar2h = dc2*ar2 - ds2*ai2;
			ai2 = dc2*ai2 + ds2*ar2;
			ar2 = ar2h;

			cc_l_index = idij = idij + IDL1;
			cc_u_index = izero2 = izero2 - IDL1;
			l3_index = idxj;
			z_index = idkk;

			//cc_l_index = (2 + j) * IDL1;
			//cc_u_index = (ipm1 - 1 - j)*IDL1;
			//l3_index = (1 + i)*IDL1;
			//z_index = (ipm1 - i)*IDL1;

			for (k = 0; k < IDL1; k++){

				++cc_l_index;
				++cc_u_index;
				++l3_index;
				++z_index;

				CH2[l3_index] += ar2*C2[cc_l_index];
				CH2[z_index] += ai2*C2[cc_u_index];

			}  // End for k
		}  // End for j
	}  // End for i

	is = 0;

	for (j = 0; j < ipph; j++){  // LOOP 129 IN FORTRAN CODE

		cc_l_index = is = is + IDL1;

		//cc_l_index = (1 + j) * IDL1;

		for (k = 0; k < IDL1; k++)  CH2[k] += C2[cc_l_index++];

	}  // End for j

	if (ID0 < L1){

		for (i = 0; i < ID0; i++){

			cc_u_index = cc_l_index = i;

			for (k = 0; k < L1; k++){
				CC[cc_l_index] = CH[cc_u_index];

				cc_l_index += incr3row;
				cc_u_index += ID0;
			}  // End for k
		}  // End for i

	}  // End if (ID0 < L1)
	else {  // else (ID0 >= L1)
		is = -ID0;
		for (k = 0; k < L1; k++){
			cc_u_index = is = is + ID0;
			//cc_u_index = k*ID0;
			cc_l_index = cc_u_index*IP;

			for (i = 0; i < ID0; i++)  CC[cc_l_index++] = CH[cc_u_index++];

		}  // End for k

	}  // End else ID0 >= L1

	//  135 IN FORTRAN CODE

	idxj = is = 0;
	idkk = IP*IDL1;

	for (j = 0; j < ipph; j++){ // LOOP 137 IN FORTRAN CODE

		cc_l_index = cc_u_index = is = is + twoido;
		//cc_l_index = cc_u_index = (1 + j) * twoido;
		--cc_l_index;

		l3_index = idxj = idxj + IDL1;
		//l3_index = (1 + j) * IDL1;

		z_index = idkk = idkk - IDL1;
		//z_index = (IP - 1 - j) * IDL1;

		for (k = 0; k < L1; k++){
			CC[cc_l_index] = CH[l3_index];
			CC[cc_u_index] = CH[z_index];

			cc_l_index += incr3row;
			cc_u_index += incr3row;
			l3_index += ID0;
			z_index += ID0;
		}  // End for k
	}  // End for j

	if (ID0 != 1) {

		ipp2 = IP*IDL1;

		if (nbd < L1){  // GO TO LINE 141

			is = -(ID0 + 1);
			ipp2 += is;

			//--ipp2;
			//ipp2 -= ID0;

			izero2 = -(incr3row + 1);

			for (j = 0; j < ipph; j++){ // LINE 141, LOOP 144 IN FORTRAN CODE

				idij = is = is + IDL1;

				idp2 = ipp2 = ipp2 - IDL1;

				idxj = idkk = izero2 = izero2 + twoido;

				//for (i = 0; i < idm2; i += 2){
				for (i = 0; i < nbd; i++){  // Try to replace +2 increment with +1 increment

					++idij;
					++idij;
					cc_l_index = idij;

					++idp2;
					++idp2;
					cc_u_index = idp2;

					++idxj;
					++idxj;
					l3_index = idxj;

					--idkk;
					--idkk;
					z_index = idkk;

					//cc_l_index = 1 + (1 + j) * IDL1 + i;
					//cc_u_index = 1 + (IP - 1 - j) * IDL1 + i;
					//l3_index = 1 + (1 + j) * twoido + i;
					//z_index = (1 + j) * twoido - 3 - i;

					for (k = 0; k < L1; k++){

						cc_l_index += ID0;
						cc_u_index += ID0;
						l3_index += incr3row;
						z_index += incr3row;

						ai2 = CH[cc_l_index];
						ar2 = CH[cc_u_index];

						CC[l3_index] = ai2 + ar2;
						CC[l3_index + 1] = CH[cc_l_index + 1] + CH[cc_u_index + 1];
						CC[z_index] = ai2 - ar2;
						CC[z_index + 1] = CH[cc_u_index + 1] - CH[cc_l_index + 1];

					}  // End for k
				}  // End for i	
			}  // End for j				

		}  // end if nbd < L1
		else { // else nbd >= L1

			izero2 = is = 0;
			idxj = -1;

			for (j = 0; j < ipph; j++){ // LOOP 140 IN FORTRAN CODE

				is = is + IDL1;
				idij = is - ID0;

				ipp2 = ipp2 - IDL1;
				idp2 = ipp2 - ID0;

				izero2 = izero2 + twoido;
				ipm1 = izero2 - incr3row;

				idxj = idxj + twoido;
				idkk = idxj - incr3row;

				for (k = 0; k < L1; k++){

					cc_l_index = idij = idij + ID0;
					cc_u_index = idp2 = idp2 + ID0;
					l3_index = ipm1 = ipm1 + incr3row;
					z_index = idkk = idkk + incr3row;

					//cc_l_index = 1 + (1 + j) * IDL1 + k * ID0;
					//cc_u_index = 1 + (IP - 1 - j) * IDL1 + k * ID0;
					//l3_index = 1 + (1 + j) * 2 * ID0 + k * IP*ID0;
					//z_index = (1 + j) * 2 * ID0 + k * incr3row - 2;

					//for (i = 0; i < idm2; i += 2){
					for (i = 0; i < nbd; i++){

						++cc_l_index;
						++cc_u_index;
						++l3_index;
						--z_index;

						ai2 = CH[cc_l_index++];
						ar2 = CH[cc_u_index++];

						CC[l3_index++] = ai2 + ar2;
						CC[l3_index] = CH[cc_l_index] + CH[cc_u_index];

						CC[z_index--] = CH[cc_u_index] - CH[cc_l_index];
						CC[z_index] = ai2 - ar2;

					}  // End for i
				}  // End for k
			}  // End for j
		}  // end else nbd >= L1

	}  // End if (ID0 != 1)

	return;
}  // End RADFG_ak1

int main()
{
	char rflag = 0; //Readiness flag 

	cout << "                                   EZFFTF_ak1 (15 June 2015)\n";
	cout << "=========================================================================== \n";
	cout << "This program computes the Fourier coefficients of \n";
	cout << "a real periodic sequence (Fourier analysis.)\n";
	cout << "\n--------------------------------------------------------------------------- \n";
	cout << "\nThe (real) data entries should have\n";
	cout << "been saved beforehand in a file named Fourier_Data_In.txt.\n";
	cout << "Fourier_Data_In.txt should be in the same folder as the EZFFTF_ak1 executable. \n";
	cout << "--------------------------------------------------------------------------- \n";
	cout << "\nThe first entry of this file must be the number of data entries, N.\n";
	cout << "The actual data entries should follow.\n";
	cout << "\nThe data is assumed to be of type double. Variables used within this program\n";
	cout << "are type double. The data is also assumed to be UNIFORMLY-spaced.\n";
	cout << "\n--------------------------------------------------------------------------- \n";
	cout << "\nThe output is written to the file EZFFTF_ak1out.txt.\n";
	cout << "\n--------------------------------------------------------------------------- \n";
	cout << "\nAdditional information is posted at the following URL:\n";
	cout << "http://www.akiti.ca/FourierTransform.html\n";
	cout << "--------------------------------------------------------------------------- \n";
	cout << "\nIs everything ready (are you ready to continue)? If yes, Enter y. \n";
	cout << "Otherwise Enter any other key. \n";
	cin >> rflag;

	if (toupper(rflag) == 'Y') {

		int N; // The number of R Array entries.

		cout << "Appear to be ready. \n";

		ifstream in("Fourier_Data_In.txt", ios::in);

		if (!in) {
			cout << "Cannot open the input file.\n";
			return 0;
		}

		in >> N; //Input the number of data entries saved in the file
		if (N < 1) {
			cout << "Invalid number of data entries input. Program terminated. \n";
			in.close(); //Close the input file before terminating
			return 0;
		}

		ofstream out("EZFFTF_ak1out.txt", ios::out);
		if (!out) {
			cout << "Cannot open the output file. Program terminated.\n";
			in.close(); //Close the input file before terminating
			return 0;
		}

		C1DArray rArray;					// A real array of length N which contains the sequence to be transformed, R. rArray is not destroyed
		C1DArray A_Array, B_Array;			// Arrays containing the resultant A and B coefficients
		C1DArray WORK1_Array, WORK2_Array;	// Arrays containing information for the factorization of the input data
		C1DArray r_copy_Array;				// Initially, a copy of the R Array
		double azero;						// The sum from I = 1 to I = N of R(I)/N
		int KMAX;							// The size of the A and B arrays.

		try { // Resize the R Array to its required size
			rArray.resize(N);
		} // End of try block

		catch (bad_alloc& xa) { // Catch block, for exceptions
			in.close();
			out.close();
			cerr << "In catch block for resizing R Array: " << xa.what() << "\n";
			cout << "\nEnter any key to continue. \n";
			cin >> rflag;
			return 0;
		} // End of catch block

		//Input the data entries from the file and put them in the R Array
		for (int i = 0; i < N; i++) in >> rArray[i];

		in.close(); //Close the input file

		out.precision(DBL_DIG);

		out << "N = " << N << "\n\n";

		if (N == 1){
			azero = rArray[0];
			out << "azero = " << azero << "\n\n";
			out.close(); // Close the output file
			return 0;
		}  // End if (N == 1)

		KMAX = N / 2;

		try { // Resize the A Array to its required size
			A_Array.resize(KMAX);
		} // End of try block

		catch (bad_alloc& xa) { // Catch block, for exceptions
			out.close();
			cerr << "In catch block for resizing A Array: " << xa.what() << "\n";
			cout << "\nEnter any key to continue. \n";
			cin >> rflag;
			return 0;
		} // End of catch block

		if (N == 2){
			azero = 0.5*(rArray[0] + rArray[1]);
			A_Array[0] = 0.5*(rArray[0] - rArray[1]);
			out << "azero = " << azero << "\n\n";
			out << "A(1) = " << A_Array[0] << "\n";
			out.close(); // Close the output file
			return 0;
		}  // End if (N == 2)

		try { // Resize the other arrays to their required sizes
			WORK1_Array.resize(N);
			WORK2_Array.resize(N);
			r_copy_Array.resize(N);
			B_Array.resize(KMAX);
		} // End of try block

		catch (bad_alloc& xa) { // Catch block, for exceptions
			out.close();
			cerr << "In catch block, so an exception occurred: " << xa.what() << "\n";
			cout << "\nEnter any key to continue. \n";
			cin >> rflag;
			return 0;
		} // End of catch block

		double cf = 2.0/N, cfm;
		int ns2 = 1, ns2m = (N + 1) / 2 - 1;
		int IFAC[FACARSIZE];				// An array that holds factorization info about the number of data inputs

		EZFFTI_ak1(N, WORK1_Array, IFAC);

		for (int i = 0; i < N; i++) r_copy_Array[i] = rArray[i];

		RFFTF_ak1(N, WORK1_Array, WORK2_Array, r_copy_Array, IFAC);

		cfm = -cf;
		azero = 0.5*cf*r_copy_Array[0];

		for (int i = 0; i < ns2m; i++){
			A_Array[i] = cf*r_copy_Array[ns2++];
			B_Array[i] = cfm*r_copy_Array[ns2++];
		} // End for i

		if ((N%2) == 0){
			A_Array[ns2m] = 0.5*cf*r_copy_Array[N - 1];
			B_Array[ns2m] = 0.0;
		} // End if ((N%2) == 0)

		out << "The program has completed successfully.\n\n";
		out << "The A Array entries follow:\n\n";
		out << "azero = " << azero << "\n";
		for (int i = 0; i < KMAX; i++) out << A_Array[i] << "\n";
		out << "\n";
		out << "The B Array entries follow:\n\n";
		for (int i = 0; i < KMAX; i++) out << B_Array[i] << "\n";
		out.close(); // Close the output file
	} //End if rflag = 'Y'
	else cout << "\nNot ready. Try again when ready with information. \n";
	cout << "\nEnter any key to continue. \n";
	cin >> rflag;
	return 0;
} // End main program.