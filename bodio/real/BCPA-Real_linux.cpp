#include<bits/stdc++.h>
#include <random>  
#include <iostream>
#include <stdint.h>
#include <sys/time.h>

using namespace std;

const double pi = acos(-1);
const double eps = 1e-10;
typedef long long LL;
typedef unsigned long long ULL;
typedef pair<int, int>pii;
typedef pair<LL, LL>pLL;


const int maxn = 2000 + 5;
const int max_I = 2e4 + 5;
const double Pc = 0.9, Pm = 0.07;   //


double trace[1][maxn];//trace[x][y]在第y个波形里  第x个点的值
int HD_key[maxn];

struct PlainData {//plain文件里的内容  res是密文
	unsigned char key[16], plain[16], res[16];
}Data[maxn];

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

//计算汉明重量
int calc_HW(int x) {
	int ans = 0;
	while (x) {
		ans++;
		x -= x & -x;  
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


void CreateFolder(const string& directoryPath)
{
	string cmd_mkdir = "mkdir -p ";
	cmd_mkdir += directoryPath;
	system(cmd_mkdir.c_str());
}

//invoke data (note: you should use your own real data)
void InvokeData(PlainData Data[], double trace[][maxn], int trace_num, int start_trace) {
	'''
	Use your own data
	'''
}

void PrintData(int i, string dir, int n, PlainData Data[], double trace[][maxn], int start_trace) {

	InvokeData(Data, trace, n, start_trace);
	char fileName[100];

	sprintf(fileName, (dir + "BCPA-Real-Trace_%d.csv").c_str(), i);
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

double calc_corr1(unsigned char key[], int trace_num, PlainData Data[], int HD_key[], double trace[][maxn]) {   //fitness function
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


struct node {
	unsigned char key[16];
	double corr;
	bool operator <(const node& oth)const {
		return corr < oth.corr;
	}
};

void creat_node(node w[], int n, int trace_num, PlainData Data[], int HD_key[], double trace[][maxn]) {


	for (int i = 1; i <= n; i++) {
		for (int j = 0; j < 16; j++) {
			w[i].key[j] = rand() % 256;
		}
		w[i].corr = calc_corr1(w[i].key, trace_num, Data, HD_key, trace);   //fitness function
	}

	sort(w + 1, w + 1 + n);  

}



void crossover_singlePoint_tournament(node w[]) {
	//if (rand() % 10 >= 5)
	//	swap(a, b);

	node res_1 = w[1];
	node res_2 = w[2];

	int tmp_bit = rand() % 128;//rand() % 128

	//printf("%d\n", tmp_bit);

	for (int i = (tmp_bit / 8) + 1; i < 16; i++) {
		res_1.key[i] = w[2].key[i];
	}
	res_1.key[tmp_bit / 8] = 0;
	for (int i = 0; i < tmp_bit % 8; i++) {
		unsigned char tmp;
		res_1.key[tmp_bit / 8] += (w[1].key[tmp_bit / 8] & (128 >> (i % 8)));
	}
	for (int i = tmp_bit % 8; i < 8; i++) {
		unsigned char tmp;
		res_1.key[tmp_bit / 8] += (w[2].key[tmp_bit / 8] & (128 >> (i % 8)));
	}

	for (int i = (tmp_bit / 8) + 1; i < 16; i++) {
		res_2.key[i] = w[1].key[i];
	}
	res_2.key[tmp_bit / 8] = 0;
	for (int i = 0; i < tmp_bit % 8; i++) {
		unsigned char tmp;
		res_2.key[tmp_bit / 8] += (w[2].key[tmp_bit / 8] & (128 >> (i % 8)));
	}
	for (int i = tmp_bit % 8; i < 8; i++) {
		unsigned char tmp;
		res_2.key[tmp_bit / 8] += (w[1].key[tmp_bit / 8] & (128 >> (i % 8)));
	}

	w[1] = res_1;
	w[2] = res_2;
}


//单字节单比特变异
void mutation_singleBit_tournament(node& x) {

	for (int i = 0; i < 16; i++) {
		if (rand() < Pm * RAND_MAX) {
			int index = rand() % 8;
			x.key[i] ^= (1 << index);   
		}
	}
}

void work_GA_tournament(int trace_num, int pn, int evolution_max, node w[], node w_tmp[], PlainData Data[], int HD_key[], double trace[][maxn], int pn_tournament) {
	for (int i = 1; i <= pn; i++) {
		node tmp;
		tmp.corr = 0;
		int index;
		for (int j = 1; j < pn_tournament; j++) {
			index = rand() % pn + 1;
			if (tmp.corr > w[index].corr) {
				memcpy(tmp.key, w[index].key, 16);
				tmp.corr = w[index].corr;
			}
		}
		memcpy(w_tmp[i].key, tmp.key, 16);
		w_tmp[i].corr = tmp.corr;
	}
	for (int i = 1; i <= pn / 2; i++) {
		node w_crossover[2 + 1];

		memcpy(w_crossover[1].key, w_tmp[i * 2 - 1].key, 16);
		memcpy(w_crossover[2].key, w_tmp[i * 2].key, 16);

		if (rand() < Pc * RAND_MAX) {
			crossover_singlePoint_tournament(w_crossover);
		}

		mutation_singleBit_tournament(w_crossover[1]);
		mutation_singleBit_tournament(w_crossover[2]);

		memcpy(w[i * 2 - 1].key, w_crossover[1].key, 16);
		w[i * 2 - 1].corr = calc_corr1(w[i * 2 - 1].key, trace_num, Data, HD_key, trace);
		memcpy(w[i * 2].key, w_crossover[2].key, 16);
		w[i * 2].corr = calc_corr1(w[i * 2].key, trace_num, Data, HD_key, trace);
	}
	sort(w + 1, w + 1 + pn);
}

void GA(unsigned char* key, int trace_num, int pn, int evolution_max, int pn_tournament, int experiment_max, FILE* fp1, string dir) {
	
	vector<double> runtimes(experiment_max);

	struct timeval tpstart, tpend;
	double timeuse;
	double total_runtime = 0.0;
	
	int count = 0, count_2 = 0, count_fail = 0;
	int fail_num = 0;
	int k, k_2 = 0;

	int start_trace_set[200 + 5];
	FILE* fp_crt; //ChooseRandomTraces
	fp_crt = fopen("./ChooseRandomTraces.csv", "r");
	int skip_round = trace_num / 50 - 1;
	for (int i = 1; i <= skip_round; i++) {
		for (int j = 1; j <= experiment_max; j++) {
			int tmp;
			fscanf(fp_crt, "%d", &tmp);
			fscanf(fp_crt, ",");
		}
	}
	for (int i = 1; i <= experiment_max; i++) {
		fscanf(fp_crt, "%d", &start_trace_set[i]);
		fscanf(fp_crt, ",");
	}
	fclose(fp_crt);


	node w[max_I], w_tmp[max_I];
	PlainData Data[maxn];
	double trace[1][maxn];
	int HD_key[maxn];


	for (k = 1; k <= experiment_max; k++) {  
		//srand(time(0));  

		int start_trace = start_trace_set[k];

		PrintData(k, dir, trace_num, Data, trace, start_trace);
		int i;

		gettimeofday(&tpstart, NULL);

		creat_node(w, pn, trace_num, Data, HD_key, trace);
		/*for (int m = 1; m < pn; m++) {
			printf("%lf   ", w[m].corr);
		}*/
		node best;
		best.corr = 0;

		//node correct;
		//memcpy(correct.key, Data[1].key, 16);
		//correct.corr = calc_corr1(correct.key, trace_num, Data, HD_key, trace);

		for (i = 1; memcmp(w[1].key, Data[1].key, 16) && i <= evolution_max; i++) {

			work_GA_tournament(trace_num, pn, evolution_max, w, w_tmp, Data, HD_key, trace, pn_tournament);

			//printf("%d ", i);
			int HD = 0;

			best = min(best, w[1]);
			best.corr = w[1].corr;


			//printf("%lf  %lf  %lf  %lf  ", correct.corr, best.corr, w[1].corr, w[2].corr);
			//for (int j = 0; j < 16; j++) {
			//	printf("%02X", best.key[j]);  
			//	HD += calc_HD(Data[1].key[j], best.key[j]);
			//}
			//printf(" %d\n", dis(best.key, Data[1].key));

		}
		memcpy(key, best.key, 16);

		if (i == evolution_max + 1) {
			//fprintf(fp, "fail\n");
			fail_num++;
		}
		else {
			count += i - 1;
			//fprintf(fp, "success\n");
		}

		gettimeofday(&tpend, NULL);

		timeuse = 1000000 * (tpend.tv_sec - tpstart.tv_sec) + tpend.tv_usec - tpstart.tv_usec;
		timeuse /= 1000000.0;
		runtimes[k - 1] = timeuse;
		total_runtime += timeuse;

		printf("第%d次实验的代数为：%d\n", k, i - 1);

	}

	//计算runtime均值和标准差

	double mean_time = total_runtime / experiment_max;
	double variance_time = 0.0;

	for (int i = 0; i < experiment_max; i++)
	{
		variance_time += (runtimes[i] - mean_time) * (runtimes[i] - mean_time);
	}

	variance_time /= (experiment_max - 1);

	double stddev_time = sqrt(variance_time);

	double success_num = (double)experiment_max - (double)fail_num;
	double success_rate = success_num / experiment_max;

	fprintf(fp1, "%d,%f,%f,%.6f,%.6f,%.6f", trace_num, success_rate, (double)(count + evolution_max * fail_num) / (double)(k - 1), total_runtime, mean_time, stddev_time);
	for (int i = 0; i < experiment_max; i++)
	{
		fprintf(fp1, ",%.6f", runtimes[i]);
	}
	fprintf(fp1, "\n");

	printf("实验成功次数为：%d, 平均进化代数为：%f，实验运行时间：%f\n", (int)success_num, (double)(count + evolution_max * fail_num) / (double)(k - 1), total_runtime);  //输出平均代数

}


vector<int> evolution_set = { 200, 300, 400 };

int main() {

	string start = currentTimetoStr();
	srand((unsigned)time(NULL));

	int trace_number_min = 50, trace_number_max = 1500, trace_number_step = 50, pn = 20, pn_tournament = 16, evolution_max, experiment_max = 200;

	for (int k = 0; k < evolution_set.size(); k++) {
		evolution_max = evolution_set[k];
		
		string dir = "./GA-CPA-Simul/" + start + "/" + "pn" + to_string(pn) + "/evolution_max" + to_string(evolution_max) + "/";
		CreateFolder(dir);


		string res = "./GA-CPA-Simul/";
		res += start;
		res += "/pn" + to_string(pn) + "/evolution_max" + to_string(evolution_max) + "/Result_GA_" + to_string(evolution_max) + "_pn" + to_string(pn) + "_" + start + ".csv";

		FILE* fp1 = fopen(res.c_str(), "w");
		if (fp1 == NULL) {
			perror("fopen");
			return -1;
		}

		fprintf(fp1, "trace_num,success_rate,ave_iter_num,total_runtime,mean_time,stddev_time");
		for (int i = 1; i <= experiment_max; i++)
			fprintf(fp1, ",run%d", i);

		fprintf(fp1, "\n");


		for (int i = trace_number_min; i <= trace_number_max; i += trace_number_step) {

			const int n = i;
			unsigned char best_key[16];
			printf("Trace number is %d\n", n);

			string dir_n = dir + to_string(n) + "/";
			CreateFolder(dir_n);

			GA(best_key, n, pn, evolution_max, pn_tournament, experiment_max, fp1, dir_n);
		}
		fclose(fp1);
	}

	printf("--------------Experiments End---------------\n");

	return 0;
}
