#include<bits/stdc++.h>
#include <stdint.h>
#include <random>  
#include <stdio.h>


using namespace std;

const double pi = acos(-1);
const double eps = 1e-10;
typedef long long LL;
typedef unsigned long long ULL;
typedef pair<int, int>pii;
typedef pair<LL, LL>pLL;

const int max_I = 2e4 + 5;
const int maxn = 2000 + 5;


struct PlainData {
	unsigned char key[16], plain[16], res[16];
};

struct Parameters {
	double A1, A2, A3, C1, C2, C3;
};

struct Parameters_second {
	double A, C;
};

struct Paras_duplication {
	double A_duplication, C_duplication;
};

//种群个体node
struct node {
	unsigned char key[16];
	double corr;
	bool operator <(const node& oth)const {
		return corr < oth.corr;
	}
};

//S盒
const unsigned char S_Table[16][16] =
{
	0x63, 0x7C, 0x77, 0x7B, 0xF2, 0x6B, 0x6F, 0xC5, 0x30, 0x01, 0x67, 0x2B, 0xFE, 0xD7, 0xAB, 0x76,
	0xCA, 0x82, 0xC9, 0x7D, 0xFA, 0x59, 0x47, 0xF0, 0xAD, 0xD4, 0xA2, 0xAF, 0x9C, 0xA4, 0x72, 0xC0,
	0xB7, 0xFD, 0x93, 0x26, 0x36, 0x3F, 0xF7, 0xCC, 0x34, 0xA5, 0xE5, 0xF1, 0x71, 0xD8, 0x31, 0x15,
	0x04, 0xC7, 0x23, 0xC3, 0x18, 0x96, 0x05, 0x9A, 0x07, 0x12, 0x80, 0xE2, 0xEB, 0x27, 0xB2, 0x75,
	0x09, 0x83, 0x2C, 0x1A, 0x1B, 0x6E, 0x5A, 0xA0, 0x52, 0x3B, 0xD6, 0xB3, 0x29, 0xE3, 0x2F, 0x84,
	0x53, 0xD1, 0x00, 0xED, 0x20, 0xFC, 0xB1, 0x5B, 0x6A, 0xCB, 0xBE, 0x39, 0x4A, 0x4C, 0x58, 0xCF,
	0xD0, 0xEF, 0xAA, 0xFB, 0x43, 0x4D, 0x33, 0x85, 0x45, 0xF9, 0x02, 0x7F, 0x50, 0x3C, 0x9F, 0xA8,
	0x51, 0xA3, 0x40, 0x8F, 0x92, 0x9D, 0x38, 0xF5, 0xBC, 0xB6, 0xDA, 0x21, 0x10, 0xFF, 0xF3, 0xD2,
	0xCD, 0x0C, 0x13, 0xEC, 0x5F, 0x97, 0x44, 0x17, 0xC4, 0xA7, 0x7E, 0x3D, 0x64, 0x5D, 0x19, 0x73,
	0x60, 0x81, 0x4F, 0xDC, 0x22, 0x2A, 0x90, 0x88, 0x46, 0xEE, 0xB8, 0x14, 0xDE, 0x5E, 0x0B, 0xDB,
	0xE0, 0x32, 0x3A, 0x0A, 0x49, 0x06, 0x24, 0x5C, 0xC2, 0xD3, 0xAC, 0x62, 0x91, 0x95, 0xE4, 0x79,
	0xE7, 0xC8, 0x37, 0x6D, 0x8D, 0xD5, 0x4E, 0xA9, 0x6C, 0x56, 0xF4, 0xEA, 0x65, 0x7A, 0xAE, 0x08,
	0xBA, 0x78, 0x25, 0x2E, 0x1C, 0xA6, 0xB4, 0xC6, 0xE8, 0xDD, 0x74, 0x1F, 0x4B, 0xBD, 0x8B, 0x8A,
	0x70, 0x3E, 0xB5, 0x66, 0x48, 0x03, 0xF6, 0x0E, 0x61, 0x35, 0x57, 0xB9, 0x86, 0xC1, 0x1D, 0x9E,
	0xE1, 0xF8, 0x98, 0x11, 0x69, 0xD9, 0x8E, 0x94, 0x9B, 0x1E, 0x87, 0xE9, 0xCE, 0x55, 0x28, 0xDF,
	0x8C, 0xA1, 0x89, 0x0D, 0xBF, 0xE6, 0x42, 0x68, 0x41, 0x99, 0x2D, 0x0F, 0xB0, 0x54, 0xBB, 0x16
};



inline int HextoDec(char ch) {
	if (isdigit(ch))return ch - '0';
	return ch - 'A' + 10;
}


void strtobyte(const char* s, unsigned char res[], int n) {
	for (int i = 0; i < n; i++) {
		res[i] = (HextoDec(s[i * 2]) << 4) + HextoDec(s[i * 2 + 1]);
	}
}

//计算相关性
double corr(int x[], double y[], int n) {
	double ex = 0, ey = 0;

	for (int i = 1; i <= n; i++) {
		ex += x[i];
		ey += y[i];
	}

	ex /= n;
	ey /= n;

	double ans = 0;
	for (int i = 1; i <= n; i++) {
		ans += (x[i] - ex) * (y[i] - ey);
	}
	ans /= n - 1;


	double t;
	t = 0;
	for (int i = 1; i <= n; i++) {
		t += (x[i] - ex) * (x[i] - ex);
	}
	t /= n - 1;
	ans /= sqrt(t);

	t = 0;
	for (int i = 1; i <= n; i++) {
		t += (y[i] - ey) * (y[i] - ey);
	}
	t /= n - 1;
	ans /= sqrt(t);

	return ans;
}

//计算汉明距离
int f[256];
int calc_HD(int x, int y) {
	x ^= y;
	if (f[x])return f[x];
	int t = x;
	//计算x二进制下多少个1
	while (t) {
		f[x]++;
		t -= t & -t;
	}
	return f[x];
}


int dis(const unsigned char* s1, const unsigned char* s2) { 
	int cnt = 0;
	for (int i = 0; i < 16; i++)
		cnt += calc_HD(s1[i], s2[i]);
	return cnt;
}

int cmp_equal(int x, int y) {//计算x和y是否相等
	if (x == y)
		return 0;
	else
		return 1;
}

int dis_byte(const unsigned char* s1, const unsigned char* s2) { 
	int cnt = 0;
	for (int i = 0; i < 16; i++)
		cnt += cmp_equal(s1[i], s2[i]);
	return cnt;
}



//计算汉明重量
int calc_HW(int x) {
	int ans = 0;
	while (x) {
		ans++;
		x -= x & -x; 
	}
	return ans;
}

//创建文件夹
void CreateFolder(const string& directoryPath)
{
	string cmd_mkdir = "mkdir -p ";
	cmd_mkdir += directoryPath;
	system(cmd_mkdir.c_str());
}

//输出当前时间
string currentTimetoStr() {
	char tmp[64];
	time_t t = time(NULL);
	tm* _tm = localtime(&t);
	int year = _tm->tm_year + 1900;
	int month = _tm->tm_mon + 1;
	int date = _tm->tm_mday;
	int hh = _tm->tm_hour;
	int mm = _tm->tm_min;
	int ss = _tm->tm_sec;
	sprintf(tmp, "%04d-%02d-%02d-%02d-%02d-%02d", year, month, date, hh, mm, ss);
	return string(tmp);
}


//Note: you should use your own data here
void InvokeData(PlainData Data[], double trace[][maxn], int trace_num, int start_trace) {
	'''
	Use your own data
	'''
}

void PrintData(int i, string dir, int n, PlainData Data[], double trace[][maxn], int start_trace) {
	InvokeData(Data, trace, n, start_trace);
	char fileName[100];

	sprintf(fileName, (dir + "BODIO-Real-Trace_%d.csv").c_str(), i);
	FILE* fp = fopen(fileName, "w");

	for (int i = 1; i <= n; i++) {
		fprintf(fp, "%lf,", trace[0][i]);

		for (int j = 0; j < 15; j++) {
			fprintf(fp, "%02X ", Data[i].plain[j]);//plaintext
		}
		fprintf(fp, "%02X,", Data[i].plain[15]);
		for (int j = 0; j < 15; j++) {
			fprintf(fp, "%02X ", Data[i].key[j]);//key
		}
		fprintf(fp, "%02X,", Data[i].key[15]);
		for (int j = 0; j < 15; j++) {
			fprintf(fp, "%02X ", Data[i].res[j]);//ciphertext
		}
		fprintf(fp, "%02X,", Data[i].res[15]);
		fprintf(fp, "\n");
	}
	fclose(fp);
}

double calc_corr_SboxPre(unsigned char key[], int trace_num, PlainData Data[], int HD_key[], double trace[][maxn]) {   //fitness function
	double corr1;
	for (int i = 2; i <= trace_num; i++) {
		HD_key[i] = 0;
		for (int j = 0; j < 16; j++) {
			int tmp = Data[i].plain[j] ^ key[j];//加密
			HD_key[i] += calc_HD(tmp, Data[i - 1].res[j]);
		}
	}
	corr1 = corr(HD_key + 1, trace[0] + 1, trace_num - 1);
	return corr1;
}

void creat_node(node w[], int pn, int trace_num, PlainData Data[], int HD_key[], double trace[][maxn]) {

	for (int i = 1; i <= pn; i++) {
		for (int j = 0; j < 16; j++) {
			w[i].key[j] = rand() % 256;
		}
		w[i].corr = calc_corr_SboxPre(w[i].key, trace_num, Data, HD_key, trace);  //fitness function
	}

	sort(w + 1, w + 1 + pn);  //调用node结构体里的operator < 进行排序

}

double C_method() {
	//srand(time(0));
	double r2, C;

	r2 = ((double)rand() / (double)RAND_MAX);  // rand()产生[0, RAND_MAX]之间的随机数
	//printf("r2: %f\n", r2);
	C = r2 + r2;
	return C;
}

double A_method(int evolution_max, int iter_num) {
	//srand(time(0));
	double r1, a, A;

	//a = 0.05 * (1 - ((double)iter_num / (double)evolution_max));
	double tmp_a = (double)iter_num / (double)evolution_max;
	a = tmp_a + tmp_a;
	r1 = ((double)rand() / (double)RAND_MAX); // rand()产生[0, RAND_MAX]之间的随机数
	//printf("rand: %d   r1: %f\n", rand(), r1);
	double tmp = a * r1;
	A = tmp + tmp - a;
	return A;
}

Parameters update_parameters(Parameters para, int evolution_max, int iter_num) {
	double C1, C2, C3, C_duplication, A1, A2, A3, A_duplication;

	A1 = A_method(evolution_max, iter_num);
	C1 = C_method();
	para.A1 = A1;
	para.C1 = C1;
	//printf("A1: %f  C1: %f\n", para.A1, para.C1);

	A2 = A_method(evolution_max, iter_num);
	C2 = C_method();
	para.A2 = A2;
	para.C2 = C2;
	//printf("A2: %f  C2: %f\n", para.A1, para.C1);

	A3 = A_method(evolution_max, iter_num);
	C3 = C_method();
	para.A3 = A3;
	para.C3 = C3;
	//printf("A3: %f  C3: %f\n", para.A1, para.C1);

	return para;
}

Paras_duplication update_paras_duplication(Paras_duplication para, int evolution_max, int iter_num) {
	double C_duplication, A_duplication, r1, r2, a;

	//C-method
	r2 = (double)rand() / (double)RAND_MAX;  // rand()产生[0, RAND_MAX]之间的随机数
	para.C_duplication = r2 + r2;

	//A-method
	//a = 2 * (1 - ((double)iter_num / (double)evolution_max));
	double tmp_a = (double)iter_num / (double)evolution_max;
	a = tmp_a + tmp_a;
	r1 = (double)rand() / (double)RAND_MAX; // rand()产生[0, RAND_MAX]之间的随机数
	double tmp = a * r1;
	para.A_duplication = tmp + tmp - a;

	return para;
}

Parameters_second update_parameters_second(Parameters_second para, int evolution_max, int iter_num) {
	double C_duplication, A_duplication, r1, r2, a;

	//C-method
	r2 = (double)rand() / (double)RAND_MAX;  // rand()产生[0, RAND_MAX]之间的随机数
	para.C = r2 + r2;

	//A-method
	//a = 2 * (1 - ((double)iter_num / (double)evolution_max));
	double tmp_a = (double)iter_num / (double)evolution_max;
	a = tmp_a + tmp_a;
	r1 = (double)rand() / (double)RAND_MAX; // rand()产生[0, RAND_MAX]之间的随机数
	double tmp = a * r1;
	para.A = tmp + tmp - a;

	return para;
}

void reduce_duplication(node w[], int pn, int evolution_max, int iter_num, int trace_num, PlainData Data[], int HD_key[], double trace[][maxn], Paras_duplication para_dup) {
	bitset<128> alpha, beta, delta, tmp;
	double d[128], X[128];
	double cmp_value_12 = memcmp(w[1].key, w[2].key, 16);
	double cmp_value_23 = memcmp(w[2].key, w[3].key, 16);

	if (cmp_value_12 == 0 || cmp_value_23 == 0) {
		for (int i = 0; i < 16; i++) {
			bitset<8> tmp(w[1].key[i]);
			//cout << tmp << endl;
			for (int j = i * 8; j < (i + 1) * 8; j++) {
				if (tmp[j % 8] == 1)
					alpha.set(j);
			}
		}

		for (int i = 0; i < 16; i++) {
			bitset<8> tmp(w[2].key[i]);
			//cout << tmp << endl;
			for (int j = i * 8; j < (i + 1) * 8; j++) {
				if (tmp[j % 8] == 1)
					beta.set(j);
			}
		}

		for (int i = 0; i < 16; i++) {
			bitset<8> tmp(w[3].key[i]);
			//cout << tmp << endl;
			for (int j = i * 8; j < (i + 1) * 8; j++) {
				if (tmp[j % 8] == 1)
					delta.set(j);
			}
		}

		if (cmp_value_12 == 0) {
			for (int i = 2; memcmp(w[1].key, w[i].key, 16) == 0; i++) {
				for (int j = 0; j < 128; j++) {
					para_dup = update_paras_duplication(para_dup, evolution_max, iter_num);
					if (rand() % 10 >= 5) {
						int tmp = delta[j];
						delta[j] = beta[j];
						beta[j] = tmp;
					}

					d[j] = fabs(para_dup.C_duplication * (double)beta[j] - (double)delta[j]);
					X[j] = fabs((double)alpha[j] - para_dup.A_duplication * d[j]);

					double r_sigmoid = (double)rand() / (double)RAND_MAX;
					double result = 1 / (1 + exp(-10 * X[j] + 5));
					if (result >= r_sigmoid) {
						tmp.set(j);
					}
				}

				for (int j = 0; j < 16; j++) {
					bitset<8> tmp_byte;
					for (int k = j * 8; k < (j + 1) * 8; k++) {
						if (tmp[k] == 1)
							tmp_byte.set(k % 8);
					}
					w[i].key[j] = tmp_byte.to_ulong();
					//w[i].key[j] = tmp.to_ulong();
				}

				w[i].corr = calc_corr_SboxPre(w[i].key, trace_num, Data, HD_key, trace);  //fitness function
			}
		}
		else if (cmp_value_23 == 0) {
			for (int i = 3; memcmp(w[2].key, w[i].key, 16) == 0; i++) {
				for (int j = 0; j < 128; j++) {
					para_dup = update_paras_duplication(para_dup, evolution_max, iter_num);
					if (rand() % 10 >= 5) {
						int tmp = delta[j];
						delta[j] = beta[j];
						beta[j] = tmp;
					}

					d[j] = fabs(para_dup.C_duplication * (double)beta[j] - (double)delta[j]);
					X[j] = fabs((double)alpha[j] - para_dup.A_duplication * d[j]);

					double r_sigmoid = (double)rand() / (double)RAND_MAX;
					double result = 1 / (1 + exp(-10 * X[j] + 5));
					if (result >= r_sigmoid) {
						tmp.set(j);
					}
				}

				for (int j = 0; j < 16; j++) {
					bitset<8> tmp_byte;
					for (int k = j * 8; k < (j + 1) * 8; k++) {
						if (tmp[k] == 1)
							tmp_byte.set(k % 8);
					}
					w[i].key[j] = tmp_byte.to_ulong();
					//w[i].key[j] = tmp.to_ulong();
				}
				w[i].corr = calc_corr_SboxPre(w[i].key, trace_num, Data, HD_key, trace);  //fitness function
			}
		}
		sort(w + 1, w + 1 + pn);
	}
}

void update_omega(node w[], int pn, Parameters para, int evolution_max, int iter_num, int trace_num, PlainData Data[], int HD_key[], double trace[][maxn], Paras_duplication para_dup) {
	bitset<128> alpha, beta, delta, omega, omega_tmp;
	double d_alpha[128], d_beta[128], d_delta[128], X1[128], X2[128], X3[128], X[128];

	for (int i = 0; i < 16; i++) {
		bitset<8> tmp(w[1].key[i]);
		//cout << tmp << endl;
		for (int j = i * 8; j < (i + 1) * 8; j++) {
			if (tmp[j % 8] == 1)
				alpha.set(j);
		}
	}

	for (int i = 0; i < 16; i++) {
		bitset<8> tmp(w[2].key[i]);
		//cout << tmp << endl;
		for (int j = i * 8; j < (i + 1) * 8; j++) {
			if (tmp[j % 8] == 1)
				beta.set(j);
		}
	}

	for (int i = 0; i < 16; i++) {
		bitset<8> tmp(w[3].key[i]);
		//cout << tmp << endl;
		for (int j = i * 8; j < (i + 1) * 8; j++) {
			if (tmp[j % 8] == 1)
				delta.set(j);
		}
	}

	for (int i = 4; i <= pn; i++) {
		for (int j = 0; j < 16; j++) {
			bitset<8> tmp(w[i].key[j]);
			//cout << tmp << endl;
			for (int k = j * 8; k < (j + 1) * 8; k++) {
				if (tmp[k % 8] == 1)
					omega.set(j);
			}
		}

		for (int j = 0; j < 128; j++) {

			para = update_parameters(para, evolution_max, iter_num);

			d_alpha[j] = fabs(para.C1 * (double)alpha[j] - (double)omega[j]);
			X1[j] = fabs((double)alpha[j] - para.A1 * d_alpha[j]);

			d_beta[j] = fabs(para.C2 * (double)beta[j] - (double)omega[j]);
			X2[j] = fabs((double)beta[j] - para.A2 * d_beta[j]);

			d_delta[j] = fabs(para.C3 * (double)delta[j] - (double)omega[j]);
			X3[j] = fabs((double)delta[j] - para.A3 * d_delta[j]);

			double r_sigmoid = (double)rand() / (double)RAND_MAX;
			X[j] = (X1[j] * 0.6 + X2[j] * 0.3 + X3[j] * 0.1);  
			double result = 1 / (1 + exp(-10 * X[j] + 5));
			if (result >= r_sigmoid) {
				omega_tmp.set(j);
			}
		}

		node w_tmp;
		w_tmp.corr = 0;

		for (int j = 0; j < 16; j++) {
			bitset<8> tmp;
			for (int k = j * 8; k < (j + 1) * 8; k++) {
				if (omega_tmp[k] == 1)
					tmp.set(k % 8);
			}
			w_tmp.key[j] = tmp.to_ulong();
			//w[i].key[j] = tmp.to_ulong();
		}
		w_tmp.corr = calc_corr_SboxPre(w_tmp.key, trace_num, Data, HD_key, trace);

		memcpy(w[i].key, w_tmp.key, 16);
		w[i].corr = w_tmp.corr;
	}
	sort(w + 1, w + 1 + pn);
	reduce_duplication(w, pn, evolution_max, iter_num, trace_num, Data, HD_key, trace, para_dup);
}


void precise_search(node best, node w_precise[], int pn, int trace_num, PlainData Data[], int HD_key[], double trace[][maxn]) {
	bitset<128> best_bit;
	for (int i = 0; i < 16; i++) {
		bitset<8> tmp(best.key[i]);
		//cout << tmp << endl;
		for (int j = i * 8; j < (i + 1) * 8; j++) {
			if (tmp[j % 8] == 1)
				best_bit.set(j);
		}
	}
	for (int i = 0; i < 128; i++) {
		bitset<128> omega;
		for (int j = 0; j < 128; j++) {
			if (i == j) {
				if(best_bit[i] == 0)
					omega.set(j);
			}
			else
			{
				omega[j] = best_bit[j];
			}
		}

		for (int j = 0; j < 16; j++) {
			bitset<8> tmp_byte;
			for (int k = j * 8; k < (j + 1) * 8; k++) {
				if (omega[k] == 1)
					tmp_byte.set(k % 8);
			}
			w_precise[i + 1].key[j] = tmp_byte.to_ulong();
		}
		w_precise[i + 1].corr = calc_corr_SboxPre(w_precise[i + 1].key, trace_num, Data, HD_key, trace);  //fitness function
	}
	sort(w_precise + 1, w_precise + 1 + 128);
}

void work_bGWO(int trace_num, node w[], PlainData Data[], int HD_key[], double trace[][maxn], int pn, int evolution_max, int iter_num) {
	Parameters para;
	Paras_duplication para_dup;
	Parameters_second para_2;
	para.A1 = 0;
	para.A2 = 0;
	para.A3 = 0;
	para.C1 = 0;
	para.C2 = 0;
	para.C3 = 0;

	para_dup.A_duplication = 0;
	para_dup.C_duplication = 0;

	para_2.A = 0;
	para_2.C = 0;

	update_omega(w, pn, para, evolution_max, iter_num, trace_num, Data, HD_key, trace, para_dup);
}
