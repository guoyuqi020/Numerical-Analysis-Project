
#if defined(BUILD_LAPACK) || defined(BUILD_GOTO2) || defined(BUILD_ATLAS)
#define PACKAGE_SPECIFIED 1
#else
#define PACKAGE_SPECIFIED 0
#endif

static_assert(PACKAGE_SPECIFIED == 1, "No Package Specified!\n");

#include <stdio.h>
#ifdef BUILD_LAPACK
#include "lapacke.h"
#endif
extern "C"
{
#ifdef BUILD_GOTO2
#include "common.h"
#include "clapack.h"
#endif
#include "cblas.h"

#ifdef BUILD_ATLAS
#include "clapack.h"
#endif
}
#include <stdlib.h>
#include <stdint.h>
#include <vector>
#include <assert.h>
#include <algorithm>
#include <math.h>
#include <random>
#include <iostream>
#include <filesystem>
#include <boost/program_options.hpp>
#include <time.h>

/// @brief 生成随机数据
/// @param M 随机数据的个数
/// @param k 生成参数k
/// @param ratio 随机误差比率，以小数格式表示，例如0.05表示5%随机误差
/// @return std::pair<x, y>
std::pair<std::vector<double>, std::vector<double>> generateRandomData(uint32_t M, uint32_t k, double ratio)
{
	// assert(M >= 10000);
	// assert(k >= 100);
	assert(ratio >= 0);
	std::vector<double> ret_x, ret_y;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int64_t> dis(-INT64_MAX, INT64_MAX); // 一个正态“分布器”，高斯分布器是 std::normal_distribution

	for (uint32_t i = 1; i <= M; i++)
	{
		double standard_x = ((double)(2 * i - 1)) / (2 * M);
		double standard_y = cos(((double)((2 * i - 1) * k)) / (2 * M));
		ret_x.emplace_back(standard_x);
		ret_y.emplace_back(standard_y * (1 + ratio * ((double)dis(gen) / (double)INT64_MAX)));
	}
	return std::make_pair(ret_x, ret_y);
}

/// @brief 将rawData输出到文件
/// @param x
/// @param y
/// @param dirPath 输出到dirPath文件夹下
/// @return 输出文件名的index
uint32_t outputRawData(std::vector<double> &x, std::vector<double> &y, std::filesystem::path dirPath)
{
	assert(std::filesystem::exists(dirPath) && std::filesystem::is_directory(dirPath));
	uint32_t idx = 0;
	std::filesystem::path outputPath;
	outputPath = dirPath / ("rawdata" + std::to_string(idx) + "_x.txt");

	while (std::filesystem::exists(outputPath))
	{
		idx += 1;
		outputPath = dirPath / ("rawdata" + std::to_string(idx) + "_x.txt");
	}
	outputPath = dirPath / ("rawdata" + std::to_string(idx) + "_x" + ".txt");
	printf("Save Raw X data to %s\n", outputPath.c_str());
	FILE *fd = fopen(outputPath.c_str(), "w+");
	assert(fd);
	assert(x.size() == y.size());
	uint32_t N = x.size();
	for (uint32_t idx = 0; idx < N; idx++)
	{
		fprintf(fd, "%.8lf\n", x.at(idx));
	}
	fclose(fd);
	outputPath = dirPath / ("rawdata" + std::to_string(idx) + "_y" + ".txt");
	printf("Save Raw Y data to %s\n", outputPath.c_str());
	fd = fopen(outputPath.c_str(), "w+");
	for (uint32_t idx = 0; idx < N; idx++)
	{
		fprintf(fd, "%.8lf\n", y.at(idx));
	}
	fclose(fd);
	return idx;
}

/// @brief 将拟合结果输出到文件
/// @param x
/// @param dirPath 输出到dirPath文件夹下
/// @param idx 输出文件名的index
void outputResults(std::vector<double> &x, std::filesystem::path dirPath, uint32_t idx)
{
	assert(std::filesystem::exists(dirPath) && std::filesystem::is_directory(dirPath));
	std::filesystem::path outputPath;
	outputPath = dirPath / ("fitresults" + std::to_string(idx) + ".txt");
	printf("Save polyfit data to %s\n", outputPath.c_str());
	FILE *fd = fopen(outputPath.c_str(), "w+");
	assert(fd);
	uint32_t N = x.size();
	for (uint32_t idx = 0; idx < N; idx++)
	{
		fprintf(fd, "%.40lf\n", x.at(idx));
	}
	fclose(fd);
	return;
}

// 算法参考：最小二乘法拟合曲线，正则方程组解法，https://blog.csdn.net/weixin_52544906/article/details/121341257

/// @brief 对于M维度向量x，计算 (x_j)^i, i=0, 1, ..., N, j=1, 2, ..., M
/// @param x 观察数据x
/// @param N 从0次方计算到N次方
/// @param M 向量x的维数
/// @return 一个(N+1)*M维的，范德蒙矩阵的转置，第(i, j)元素存储(x_j)^i
std::vector<double> PowerNT(std::vector<double> &x, uint32_t N, uint32_t M)
{
	assert(x.size() == M);
	std::vector<double> results;
	results.resize((N + 1) * M, 1);
	for (uint32_t i = 1; i <= N; i++)
	{
		for (uint32_t j = 0; j < M; j++)
		{
			results.at(i * M + j) = results.at((i - 1) * M + j) * x.at(j);
		}
	}
	return results;
}

/// @brief 对于M维度向量x，计算 (x_i)^j, i=1, 2, ..., M, j=0, 1, ..., N
/// @param x 观察数据x
/// @param N 从0次方计算到N次方
/// @param M 向量x的维数
/// @return 一个M*(N+1)维的，范德蒙矩阵，第(i, j)元素存储(x_i)^j
std::vector<double> PowerN(std::vector<double> &x, uint32_t N, uint32_t M)
{
	assert(x.size() == M);
	std::vector<double> results;
	results.resize((N + 1) * M, 1);
	for (uint32_t i = 0; i < M; i++)
	{
		for (uint32_t j = 1; j <= N; j++)
		{
			results.at(i * (N + 1) + j) = results.at(i * (N + 1) + j - 1) * x.at(i);
		}
	}
	return results;
}

/// @brief 拟合多项式曲线a_0+a_1*x+a_2*x^2+...+a_n*x^n
/// @param x
/// @param y
/// @param M 观察数据的个数
/// @param N 待拟合的多项式的个数
/// @return 参数向量<a_0, a_1, ... a_n>
std::vector<double> fitPoints(std::vector<double> &x, std::vector<double> &y, uint32_t M, uint32_t N)
{

	time_t start_time = clock();

#if defined(BUILD_LAPACK)
	auto vonMatrixT = PowerN(x, N, M);
	LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', M, N + 1, 1, vonMatrixT.data(), N + 1, y.data(), 1);
#elif defined(BUILD_ATLAS) || defined(BUILD_GOTO2)
	auto vonMatrixT = PowerNT(x, N, M);
	clapack_dgels(CblasRowMajor, CblasTrans, N + 1, M, 1, vonMatrixT.data(), M, y.data(), M);
#endif
	y.resize(N + 1);
	time_t end_time = clock();
	printf("Time usage: %.4lfs\n", ((double)end_time - start_time) / CLOCKS_PER_SEC);
	return y;
}

int main(int argc, char **argv)
{
	using namespace boost::program_options;
	options_description desc("Options");
	uint32_t M,
		k, N;
	std::string rawDir, fitDir;
	desc.add_options()("M,M", value<uint32_t>(&M)->required())("k,k", value<uint32_t>(&k)->required())("N,N", value<uint32_t>(&N)->required())("rawdir,r", value<std::string>(&rawDir)->default_value("."))("fitdir,f", value<std::string>(&fitDir)->default_value("."));
	variables_map vm;
	store(parse_command_line(argc, argv, desc), vm);
	notify(vm);

	printf("Current Parameters: M=%u k=%u N=%u\n", M, k, N);

	auto rawData = generateRandomData(M, k, 0.005);

	auto filenameIdx = outputRawData(rawData.first, rawData.second, rawDir);
	auto res = fitPoints(rawData.first, rawData.second, M, N);
	assert(res.size() == N + 1);
	outputResults(res, fitDir, filenameIdx);
	return 0;
}