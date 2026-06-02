#include "HybridTools_real_linux.h"
#include <sys/time.h>

void bGWO(unsigned char* key, int trace_num, string start, int pn, int evolution_max, int flag_num, FILE* fp1, string dir, int experiment_max) {
	int count = 0, count_2 = 0, count_fail = 0, fail_num = 0;
	int k, k_2 = 0;

	int start_trace_set[200 + 5];
	'''
	Choose random power traces data in your own dataset.
	'''

	vector<double> runtimes(experiment_max);

	struct timeval tpstart, tpend;
	double timeuse;
	double total_runtime = 0.0;


	node w[max_I], w_precise[128 + 1];
	PlainData Data[maxn];
	double trace[1][maxn];
	int HD_key[maxn];

	for (k = 1; k <= experiment_max; k++) {   

		int start_trace = start_trace_set[k];

		PrintData(k, dir, trace_num, Data, trace, start_trace);
		int i, flag = 0;

		gettimeofday(&tpstart, NULL);

		creat_node(w, pn, trace_num, Data, HD_key, trace);

		node alpha, beta, delta;

		alpha.corr = 0;
		beta.corr = 0;
		delta.corr = 0;



		memcpy(alpha.key, w[1].key, 16);
		memcpy(beta.key, w[2].key, 16);
		memcpy(delta.key, w[3].key, 16);

		
		for (i = 1; memcmp(alpha.key, Data[1].key, 16) && memcmp(beta.key, Data[1].key, 16) && memcmp(delta.key, Data[1].key, 16) && flag != flag_num; i++) {

			work_bGWO(trace_num, w, Data, HD_key, trace, pn, evolution_max, i);

			if (i >= 2 && memcmp(w[1].key, alpha.key, 16) == 0) {
				flag++;
			}
			else {
				flag = 0;
			}

			memcpy(alpha.key, w[1].key, 16);
			memcpy(beta.key, w[2].key, 16);
			memcpy(delta.key, w[3].key, 16);
		}
		memcpy(key, alpha.key, 16);

		int j = 1;

		if (flag == flag_num) {
			k_2 = k_2 + 1;
			for (j = 1; memcmp(alpha.key, Data[1].key, 16) && memcmp(beta.key, Data[1].key, 16) && memcmp(delta.key, Data[1].key, 16) && j <= 10; j++) {
				precise_search(alpha, w_precise, pn, trace_num, Data, HD_key, trace);

				memcpy(alpha.key, w_precise[1].key, 16);
				memcpy(beta.key, w_precise[2].key, 16);
				memcpy(delta.key, w_precise[3].key, 16);

			}

			if (j == 10 + 1) {
				fail_num++;
				count_fail += i - 1;
			}
			else {
				count += i - 1;
				count_2 += j - 1;
			}
		}

		gettimeofday(&tpend, NULL);

		timeuse = 1000000 * (tpend.tv_sec - tpstart.tv_sec) + tpend.tv_usec - tpstart.tv_usec;
		timeuse /= 1000000.0;
		runtimes[k - 1] = timeuse;
		total_runtime += timeuse;

		printf("第%d次实验的代数为：%d + %d，运行时间为：%f s\n", k, i - 1, j - 1, timeuse);

	}
	

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

	double ave_onebitflip = 0.0;
	if (k_2 > 0) {
		ave_onebitflip = (double)(count_2 + 10 * fail_num) / (double)(k_2);
	}

	fprintf(fp1, "%d,%f,%f,%f,%.6f,%.6f,%.6f", trace_num, success_rate, (double)(count + count_fail) / (double)(k - 1), ave_onebitflip, total_runtime, mean_time, stddev_time);
	for (int i = 0; i < experiment_max; i++)
	{
		fprintf(fp1, ",%.6f", runtimes[i]);
	}
	fprintf(fp1, "\n");
	printf("实验成功次数为：%d, k_2次数：%d，平均进化代数为：%f + %f，实验运行时间：%f s\n", (int)success_num, k_2, (double)(count + count_fail) / (double)(k - 1), ave_onebitflip, total_runtime);  

}

vector<int> pn_set = { 4, 5, 6, 7, 8 };

int main() {

	string start = currentTimetoStr();
	srand((unsigned)time(NULL));

	int trace_number_min = 50, trace_number_max = 1500, trace_number_step = 50, pn, evolution_max = 3000, flag_num = 100, experiment_max = 200;   // 50 1000 50 

	for (int k = 0; k < pn_set.size(); k++) {
		pn = pn_set[k];

		string dir = "bGWO-CPA-Real/" + start + "/" + "pn" + to_string(pn) + "/";  
		CreateFolder(dir);

		string res = "bGWO-CPA-Real/";
		res += start;
		res += "/pn" + to_string(pn) + "/Result_bGWO_flag" + to_string(flag_num) + "_pn" + to_string(pn) + "_" + start + ".csv";

		FILE* fp1 = fopen(res.c_str(), "w");
		if (fp1 == NULL) {
			perror("fopen");
			return -1;
		}

		fprintf(fp1, "trace_num,success_rate,ave_iter_num,ave_onebitflip_num,total_runtime,mean_time,stddev_time");
		for (int i = 1; i <= experiment_max; i++)
			fprintf(fp1, ",run%d", i);

		fprintf(fp1, "\n");

		for (int i = trace_number_min; i <= trace_number_max; i += trace_number_step) {

			const int n = i;
			unsigned char best_key[16];
			printf("Trace number is %d\n", n);

			string dir_n = dir + to_string(n) + "/";
			CreateFolder(dir_n);

			bGWO(best_key, n, start, pn, evolution_max, flag_num, fp1, dir_n, experiment_max);
		}

		fclose(fp1);

	}

	printf("--------------Experiments End---------------\n");

	return 0;
}
