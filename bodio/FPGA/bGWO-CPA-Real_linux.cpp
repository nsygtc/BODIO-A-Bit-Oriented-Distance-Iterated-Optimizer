//#include "uTools.h"
#include "HybridTools_real_linux.h"
#include <sys/time.h>

void bGWO(unsigned char* key, int trace_num, string start, int pn, int evolution_max, int flag_num, FILE* fp1, string dir, int experiment_max) {
	int count = 0, count_2 = 0, count_fail = 0, fail_num = 0;
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

	vector<double> runtimes(experiment_max);

	struct timeval tpstart, tpend;
	double timeuse;
	double total_runtime = 0.0;


	node w[max_I], w_precise[128 + 1];
	PlainData Data[maxn];
	double trace[1][maxn];//trace[x][y]在第y个波形里  第x个点的值
	int HD_key[maxn];

	for (k = 1; k <= experiment_max; k++) {  //进行k次实验
		//srand(time(0));  //把当前的时间作为随机数种子

		int start_trace = start_trace_set[k];

		PrintData(k, dir, trace_num, Data, trace, start_trace);
		int i, flag = 0;

		gettimeofday(&tpstart, NULL);

		creat_node(w, pn, trace_num, Data, HD_key, trace);
		/*for (int m = 1; m < pn; m++) {
			printf("%lf   ", w[m].corr);
		}*/
		node alpha, beta, delta;
		//best.corr = 0;
		alpha.corr = 0;
		beta.corr = 0;
		delta.corr = 0;

		//node correct;
		//memcpy(correct.key, Data[1].key, 16);
		//correct.corr = calc_corr_SboxPre(correct.key, trace_num, Data, HD_key, trace);

		memcpy(alpha.key, w[1].key, 16);
		memcpy(beta.key, w[2].key, 16);
		memcpy(delta.key, w[3].key, 16);

		//枚举
		for (i = 1; memcmp(alpha.key, Data[1].key, 16) && memcmp(beta.key, Data[1].key, 16) && memcmp(delta.key, Data[1].key, 16) && flag != flag_num; i++) {

			work_bGWO(trace_num, w, Data, HD_key, trace, pn, evolution_max, i);

			//printf("%d ", i);

			//枚举使用

			if (i >= 2 && memcmp(w[1].key, alpha.key, 16) == 0) {
				flag++;
			}
			else {
				flag = 0;
			}

			//best = min(best, w[1]);
			memcpy(alpha.key, w[1].key, 16);
			memcpy(beta.key, w[2].key, 16);
			memcpy(delta.key, w[3].key, 16);

			//printf("%lf  %lf  %lf  %f  %f ", correct.corr, w[1].corr, w[2].corr, w[3].corr, w[4].corr);
			//for (int j = 0; j < 16; j++) {
			//	printf("%02X", alpha.key[j]);  //%02X输出为16进制数
			//}
			//printf(" %d\n", dis(alpha.key, Data[1].key));//比较正确密钥和当前密钥x.key差几比特*/

		}
		memcpy(key, alpha.key, 16);


		//枚举

		int j = 1;

		if (flag == flag_num) {
			k_2 = k_2 + 1;
			for (j = 1; memcmp(alpha.key, Data[1].key, 16) && memcmp(beta.key, Data[1].key, 16) && memcmp(delta.key, Data[1].key, 16) && j <= 10; j++) {
				//printf("%d ", j);
				precise_search(alpha, w_precise, pn, trace_num, Data, HD_key, trace);

				//printf("%lf  %lf  %lf  %lf  ", w_precise[1].corr, w_precise[2].corr, w_precise[3].corr, w_precise[4].corr);
				//best = min(best, w_precise[1]);
				//memcpy(best.key, w_precise[1].key, 16);
				//best.corr = calc_corr_SboxPre(best.key, trace_num, Data, HW_key, trace);

				memcpy(alpha.key, w_precise[1].key, 16);
				memcpy(beta.key, w_precise[2].key, 16);
				memcpy(delta.key, w_precise[3].key, 16);
				//for (int j = 0; j < 16; j++) {
				//	printf("%02X", alpha.key[j]);  //%02X输出为16进制数
				//}
				//printf(" %d\n", dis(alpha.key, Data[1].key));//比较正确密钥和当前密钥x.key差几比特
			}

			if (j == 10 + 1) {
				//fprintf(fp, "fail\n");
				fail_num++;
				count_fail += i - 1;
			}
			else {
				count += i - 1;
				count_2 += j - 1;
				//fprintf(fp, "success\n");
			}
		}

		gettimeofday(&tpend, NULL);

		timeuse = 1000000 * (tpend.tv_sec - tpstart.tv_sec) + tpend.tv_usec - tpstart.tv_usec;
		timeuse /= 1000000.0;
		runtimes[k - 1] = timeuse;
		total_runtime += timeuse;

		printf("第%d次实验的代数为：%d + %d，运行时间为：%f s\n", k, i - 1, j - 1, timeuse);

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

	double ave_onebitflip = 0.0;
	if (k_2 > 0) {
		ave_onebitflip = (double)(count_2 + 10 * fail_num) / (double)(k_2);
	}

	fprintf(fp1, "%d,%f,%f,%f,%.6f,%.6f,%.6f", trace_num, success_rate, (double)(count + count_fail) / (double)(k - 1), ave_onebitflip, total_runtime, mean_time, stddev_time);
	//fprintf(fp1, "%d,%f,%f,%d,%d\n", trace_num, result / (double)experiment_max, (double)(count + evolution_max * fail_num) / (double)(k - 1), computation_fitness, computation_sort);
	for (int i = 0; i < experiment_max; i++)
	{
		fprintf(fp1, ",%.6f", runtimes[i]);
	}
	fprintf(fp1, "\n");

	printf("实验成功次数为：%d, k_2次数：%d，平均进化代数为：%f + %f，实验运行时间：%f s\n", (int)success_num, k_2, (double)(count + count_fail) / (double)(k - 1), ave_onebitflip, total_runtime);  //输出平均代数


}

//vector<int> pn_set = { 3, 4, 5, 6, 7, 8 }; 
vector<int> pn_set = { 4, 5, 6, 7, 8 };

int main() {

	string start = currentTimetoStr();
	srand((unsigned)time(NULL));

	int trace_number_min = 50, trace_number_max = 1500, trace_number_step = 50, pn, evolution_max = 3000, flag_num = 100, experiment_max = 200;   // 50 1000 50 

	for (int k = 0; k < pn_set.size(); k++) {
		pn = pn_set[k];
		//pn = 8;

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

			const int n = i;////n表示生成的波形数量，读取n次加密结果
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
