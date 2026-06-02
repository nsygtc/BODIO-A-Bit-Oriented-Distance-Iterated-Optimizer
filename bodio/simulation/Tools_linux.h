#include<bits/stdc++.h>
#include <random>  
#include <stdio.h>
#include <stdint.h>
#include <io.h>
#include <thread>

using namespace std;

const double pi = acos(-1);
const double eps = 1e-10;
typedef long long LL;
typedef unsigned long long ULL;
typedef pair<int, int>pii;
typedef pair<LL, LL>pLL;

const int max_I = 2e4 + 5;
const int maxn = 1000 + 5;

double sigma = 5;//标准差 

struct PlainData {//plain文件里的内容  res是密文
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
		return corr > oth.corr;
	}
};


//16进制转换为10进制
inline int HextoDec(char ch) {
	if (isdigit(ch))return ch - '0';
	return ch - 'A' + 10;
}

//16进制字符串转 16个int  每个int表示一个字节
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


int f[256];
int calc_HD(int x, int y) {//计算x  y汉明距离
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


int dis(const unsigned char* s1, const unsigned char* s2) { //密钥与正确密钥差多少比特
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

int dis_byte(const unsigned char* s1, const unsigned char* s2) { //密钥与正确密钥差几个字节
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
		x -= x & -x;  //把x最低位变为0   二进制x的负数 -x = 原码 x 求反 得反码，反码+1为补码，即为-x  正负相与 
	}
	return ans;
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



//创建文件夹
void CreateFolder(const string& directoryPath)
{
	string cmd_mkdir = "mkdir -p ";
	cmd_mkdir += directoryPath;
	system(cmd_mkdir.c_str());
}


//生成随机明文和随机密钥
void RandomPlainAndKey(int n, PlainData Data[]) {
	int tmp;
	//srand(time(0));  //把当前的时间作为随机数种子
	int key[16];//随机生成的16字节密钥

	for (int j = 0; j < 16; j++) {  //生成16字节随机密钥
		key[j] = rand() % 256;
	}
	/*printf("random key is ");
	for (int i = 0; i < 16; i++) {
		printf("%02X ", key[i]);
	}
	putchar(10);*/

	for (int i = 1; i <= n; i++) {  //生成n条16字节随机明文
		for (int j = 0; j < 16; j++) {
			tmp = rand() % 256;
			Data[i].plain[j] = tmp;
			Data[i].key[j] = key[j];
		}
	}
}

//生成中间值
void GenerateIntermediateValue_SBoxPre(int n, PlainData Data[]) {
	for (int i = 1; i <= n; i++) {
		for (int j = 0; j < 16; j++) {
			int tmp = Data[i].key[j] ^ Data[i].plain[j];//加密
			Data[i].res[j] = tmp;
			//Data[i].res[j] = S_Table[tmp >> 4][tmp & ((1 << 4) - 1)];
		}
	}
}

//生成中间值
void GenerateIntermediateValue_SBoxAfter(int n, PlainData Data[]) {
	for (int i = 1; i <= n; i++) {
		for (int j = 0; j < 16; j++) {
			int tmp = Data[i].key[j] ^ Data[i].plain[j];//加密
			Data[i].res[j] = S_Table[tmp >> 4][tmp & ((1 << 4) - 1)];
		}
	}
}

//生成波形
void GenerateTrace(int n, PlainData Data[], int HW[], double trace[][maxn]) {
	default_random_engine generator;  //一个随机数引擎类,它定义的调用运算符返回一个随机的 unsigned 类型的值  

	for (int i = 1; i <= n; i++) {
		HW[i] = 0;
		for (int j = 0; j < 16; j++) {
			HW[i] += calc_HW(Data[i].res[j]);
		}
		normal_distribution<double> dist(HW[i], sigma);
		trace[0][i] = dist(generator);
	}
}

void PrintData(int i, string dir, int n, PlainData Data[], int HW[], double trace[][maxn]) {

	//srand(time(0));  //把当前的时间作为随机数种子

	RandomPlainAndKey(n, Data);
	GenerateIntermediateValue_SBoxPre(n, Data);
	//GenerateIntermediateValue_SBoxAfter(n, Data);
	GenerateTrace(n, Data, HW, trace);
	char fileName[100];
	sprintf(fileName, (dir + "bGWO-CPA-Trace_%d.csv").c_str(), i);
	FILE* fp = fopen(fileName, "w");

	for (int i = 1; i <= n; i++) {
		fprintf(fp, "%lf,", trace[0][i]);

		for (int j = 0; j < 15; j++) {
			fprintf(fp, "%02X ", Data[i].plain[j]);//明文
		}
		fprintf(fp, "%02X,", Data[i].plain[15]);
		for (int j = 0; j < 15; j++) {
			fprintf(fp, "%02X ", Data[i].key[j]);//密钥
		}
		fprintf(fp, "%02X,", Data[i].key[15]);
		for (int j = 0; j < 15; j++) {
			fprintf(fp, "%02X ", Data[i].res[j]);//密文
		}
		fprintf(fp, "%02X,", Data[i].res[15]);
		fprintf(fp, "\n");
	}
	fclose(fp);
}

double calc_corr_SboxPre(unsigned char key[], int trace_num, PlainData Data[], int HW_key[], double trace[][maxn]) {   //fitness function
	double corr1;
	for (int i = 1; i <= trace_num; i++) {
		HW_key[i] = 0;
		for (int j = 0; j < 16; j++) {
			int tmp = key[j] ^ Data[i].plain[j];//加密
			//int tmp2 = S_Table[tmp >> 4][tmp & ((1 << 4) - 1)];
			HW_key[i] += calc_HW(tmp);
		}
	}
	corr1 = corr(HW_key, trace[0], trace_num);
	return corr1;
}

double calc_corr_SboxAfter(unsigned char key[], int trace_num, PlainData Data[], int HW_key[], double trace[][maxn]) {   //fitness function
	double corr1;
	for (int i = 1; i <= trace_num; i++) {
		HW_key[i] = 0;
		for (int j = 0; j < 16; j++) {
			int tmp = key[j] ^ Data[i].plain[j];//加密
			int tmp2 = S_Table[tmp >> 4][tmp & ((1 << 4) - 1)];
			HW_key[i] += calc_HW(tmp2);
		}
	}
	corr1 = corr(HW_key, trace[0], trace_num);
	return corr1;
}

void creat_node(node w[], int pn, int trace_num, PlainData Data[], int HW_key[], double trace[][maxn]) {
	//char s[] = "0123456799ABCDEFFEDCBA9876543210";
	//strtobyte(s, w[1].key, 16);
	/*for (int i = 0; i < 16; i++) {
		w[1].key[i] = Data[1].key[i];
	}
	w[1].corr = calc_corr_SboxPre(w[1].key, trace_num, Data, HW_key, trace);*/
	//w[1].corr = calc_corr2(w[1].key);
	//printf("正确密钥的相关系数为 %lf \n", w[1].corr);
	for (int i = 1; i <= pn; i++) {
		for (int j = 0; j < 16; j++) {
			w[i].key[j] = rand() % 256;
		}
		w[i].corr = calc_corr_SboxPre(w[i].key, trace_num, Data, HW_key, trace);  //fitness function
	}

	sort(w + 1, w + 1 + pn);  //调用node结构体里的operator < 进行排序
	//printf("排序后第一个密钥的相关系数为 %lf  %lf  %lf\n", w[1].corr, w[2].corr, w[3].corr);
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

void reduce_duplication(node w[], int pn, int evolution_max, int iter_num, int trace_num, PlainData Data[], int HW_key[], double trace[][maxn], Paras_duplication para_dup) {
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

				w[i].corr = calc_corr_SboxPre(w[i].key, trace_num, Data, HW_key, trace);  //fitness function
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
				w[i].corr = calc_corr_SboxPre(w[i].key, trace_num, Data, HW_key, trace);  //fitness function
			}
		}
		sort(w + 1, w + 1 + pn);
	}
}

void update_omega(node w[], int pn, Parameters para, int evolution_max, int iter_num, int trace_num, PlainData Data[], int HW_key[], double trace[][maxn], Paras_duplication para_dup) {
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
	//cout << alpha << endl;

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
			X[j] = (X1[j] * 0.6 + X2[j] * 0.3 + X3[j] * 0.1);  // 平均 /  5 3 2
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
		w_tmp.corr = calc_corr_SboxPre(w_tmp.key, trace_num, Data, HW_key, trace);
		//w[i].corr = calc_corr_SboxPre(w[i].key, trace_num, Data, HW_key, trace);

		/*double r_compare = (double)rand() / (double)RAND_MAX;

		if (((double)iter_num / (double)evolution_max) <= 0.8) {   //进化前期
			if (w[i].corr < w_tmp.corr && r_compare < 0.8) {
				memcpy(w[i].key, w_tmp.key, 16);
				w[i].corr = w_tmp.corr;
			}
		}
		else {   //进化后期
			if (r_compare < 0.8) {
				memcpy(w[i].key, w_tmp.key, 16);
				w[i].corr = w_tmp.corr;
			}
		}*/
		//if (w[i].corr < w_tmp.corr) {// && r_compare < 0.8
			memcpy(w[i].key, w_tmp.key, 16);
			w[i].corr = w_tmp.corr;
		//}
	}
	sort(w + 1, w + 1 + pn);
	reduce_duplication(w, pn, evolution_max, iter_num, trace_num, Data, HW_key, trace, para_dup);
}

void precise_search(node best, node w_precise[], int pn, int trace_num, PlainData Data[], int HW_key[], double trace[][maxn]) {
	bitset<128> best_bit;
	//printf("----------------------------------------test precise_search-----------------------------------------\n");
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
			//w[i].key[j] = tmp.to_ulong();
		}
		w_precise[i + 1].corr = calc_corr_SboxPre(w_precise[i + 1].key, trace_num, Data, HW_key, trace);  //fitness function
	}
	sort(w_precise + 1, w_precise + 1 + 128);
}

void work_bGWO(int trace_num, node w[], PlainData Data[], int HW_key[], double trace[][maxn], int pn, int evolution_max, int iter_num) {
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

	update_omega(w, pn, para, evolution_max, iter_num, trace_num, Data, HW_key, trace, para_dup);
}
