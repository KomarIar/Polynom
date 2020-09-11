#include <iostream>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <functional>
#include <thread>
#include <queue>
#include <mutex>
#include <chrono>
#include <condition_variable>

template <typename T>
class Queue
{
	std::queue<T> queue_;
	std::mutex mutex_;
	std::condition_variable cv_;

public:
	T pop()
	{
		std::unique_lock<std::mutex> mlock(mutex_);
		while (queue_.empty())
		{
			cv_.wait(mlock);
		}
		auto item = queue_.front();
		queue_.pop();
		mlock.unlock();
		return item;
	}

	void pop(T& item)
	{
		std::unique_lock<std::mutex> mlock(mutex_);
		while (queue_.empty())
		{
			cv_.wait(mlock);
		}
		item = queue_.front();
		queue_.pop();
		mlock.unlock();
	}

	void push(const T& item)
	{
		std::unique_lock<std::mutex> mlock(mutex_);
		queue_.push(item);
		cv_.notify_one();
	}

	int rough_size() { return queue_.size(); }
};

std::vector <int> num2vec(int num, int base, int len)
{
	std::vector <int> vec(len, 0);
	for (int i = len - 1; i >= 0; i--) {
		vec[i] = num % base;
		num /= base;
	}
	return vec;
}

int vec2num(std::vector <int> vec, int base)
{
	int num = 0;
	for (int i = vec.size() - 1; i >= 0; i--) {
		num *= base;
		num += vec[i];
	}
	return num;
}

class Func_represent
{
	int function_;
	int qvar_;
	std::vector <int> polynom_;

public:
	Func_represent(int qvar, int function, std::vector <int> polynom) :
		qvar_(qvar),
		function_(function),
		polynom_(polynom)
	{}

	void print(std::ofstream &output)
	{
		output << '(';
		std::vector <int> vector = num2vec(function_, 2, pow(2, qvar_));
		for (int i = 0; i < pow(2, qvar_); i++) {
			output << vector[i];
		}
		output << ')' << " = ";
		for (int i = 0; i < polynom_.size() - 1; i++) {
			if (i != 0) output << " + ";
			if (polynom_[i] == 0) {
				output << "1";
				continue;
			}
			std::vector <int> add = num2vec(polynom_[i], 3, qvar_);
			for (int j = 0; j < qvar_; j++) {
				if (add[j] == 2) output << "(-x" << j + 1 << ")";
				if (add[j] == 1) output << "(x" << j + 1 << ")";
			}
		}
		output << std::endl;
	}
};

class Polynom_gen
{
	std::vector <int> polynom_;
	int qvar_;
	int ec_max;

public:
	Polynom_gen(int qvar) :
		qvar_(qvar)
	{
		ec_max = pow(3, qvar_);
		polynom_.push_back(-1);
		polynom_.push_back(-1);
	}

	std::vector <int> next_polynom()
	{
		int i = 0;
		++polynom_[0];
		while (polynom_[i] == ec_max - i) {
			++i;
			++polynom_[i];
		}
		if (polynom_.back() != -1) {
			polynom_.push_back(-1);
		}
		while (i > 0) {
			--i;
			polynom_[i] = polynom_[i + 1] + 1;
		}
		return polynom_;
	}

	int get_cur_size()
	{
		return polynom_.size() - 1;
	}
};

class Problem_data
{
	int qvar_;                                   //количество переменных
	int func_left_;								 //осталось посчитать функций
	std::fstream calc_fs;						 //бинарник со статусами всех функций
	std::vector <Func_represent> Max_polynoms_;	 //вектор представлений функций с максимальной длиной полинома
	int shannon_function_;                       //текущая максимальная длина полинома

public:
	std::vector <int> ec_table;                  //таблица конъюнкций
	std::vector <char> bitmask;					 //маска битов для получения статуса функции из файла
	Problem_data(int qvar) :
		qvar_(qvar),
		shannon_function_(0),
		calc_fs("calculated.bin", std::fstream::binary | std::fstream::trunc | std::fstream::in | std::fstream::out)
	{
		int mask = pow(2, 7);
		for (int i = 0; i < 8; i++) {
			bitmask.push_back(mask);
			mask /= 2;
		}
		init_comp_func();
		make_ec_table();
	}

	~Problem_data()
	{
		calc_fs.close();
	}

	int get_qvar() { return qvar_; }
	int get_func_left() { return func_left_; }

	void write_func(std::vector <int> polynom)
	{
		int func = 0;
		int polynom_size = polynom.size() - 1;
		for (unsigned int i = 0; i < polynom_size; i++) {
			func ^= ec_table[polynom[i]];
		}
		int offset = func / 8;
		int bit = func % 8;
		char byte = 0;
		calc_fs.seekg(offset);
		calc_fs.get(byte);

		if ((byte & bitmask[bit]) != 0) return;
		if (polynom_size > shannon_function_) {
			shannon_function_ = polynom_size;
			std::cout << polynom_size << std::endl;
			Max_polynoms_.clear();
		}
		byte |= bitmask[bit];
		calc_fs.seekp(offset);
		calc_fs.put(byte);

		--func_left_;
		Max_polynoms_.push_back(Func_represent(qvar_, func, polynom));
	}

	void print_max_pol()
	{
		std::ofstream output("result.txt", std::ostream::trunc);
		output << "L(" << qvar_ << ") = " << shannon_function_ << std::endl;
		for (int i = 0; i < Max_polynoms_.size(); i++) {
			output << i << '\t';
			Max_polynoms_[i].print(output);
		}
		output.close();
	}

private:
	void init_comp_func()
	{
		func_left_ = pow(2, pow(2, qvar_));
		char byte = 0b10000000;
		calc_fs.put(byte);
		byte = 0b00000000;
		for (int i = 1; i < func_left_ / 8; i++) {
			calc_fs.put(byte);
		}
		--func_left_;
	}

	void make_ec_table()
	{
		ec_table.push_back(pow(2, pow(2, qvar_)) - 1);
		int ec_quantity = pow(3, qvar_);
		int set_quantity = pow(2, qvar_);
		for (int i = 1; i < ec_quantity; i++) {
			int result = 0;
			std::vector <int> conjunction = num2vec(i, 3, qvar_);
			for (int j = 0; j < set_quantity; j++) {
				std::vector <int> set = num2vec(j, 2, qvar_);
				int set_val = 1;
				for (int cur = 0; cur < qvar_; cur++) {
					if ((conjunction[cur] == 1 && set[cur] != 1) || (conjunction[cur] == 2 && set[cur] != 0)) { set_val = 0; break; }
				}
				result = result * 2 + set_val;
			}
			ec_table.push_back(result);
		}
	}
};

class Bruteforce
{
	Problem_data* data_;
	std::vector<std::thread> threads_;
	Polynom_gen generator_;
	Queue <std::vector <int>> pol_queue_;
	std::mutex file_mutex_;

public:
	Bruteforce(int qthread, Problem_data* data) :
		threads_(qthread),
		data_(data),
		generator_(data_->get_qvar())
	{}

	void check_polynoms()
	{
		while (data_->get_func_left() > 0) {
			std::vector <int> polynom = pol_queue_.pop();

			std::ifstream calc_fs ("calculated.bin", std::ifstream::binary);
			int func = 0;
			int polynom_size = polynom.size() - 1;
			for (unsigned int i = 0; i < polynom_size; i++) {
				func ^= data_->ec_table[polynom[i]];
				if (!func) continue;
			}
			int offset = func / 8;
			int bit = func % 8;
			char byte = 0;
			calc_fs.seekg(offset);
			calc_fs.get(byte);
			if ((byte & data_->bitmask[bit]) == 0) {
				file_mutex_.lock();
				data_->write_func(polynom);
				file_mutex_.unlock();
			}
			calc_fs.close();
		}
	}

	void gen_polynoms()
	{
		while (data_->get_func_left() > 0) {
			pol_queue_.push(generator_.next_polynom());
			while (pol_queue_.rough_size() > 1000000)
				std::this_thread::sleep_for(std::chrono::seconds(1));
		}
	}

	void brute_force()
	{
		for (int i = 0; i < threads_.size(); i++) {
			threads_[i] = std::thread(&Bruteforce::check_polynoms, this);
		}
		gen_polynoms();

		for (auto&& thread : threads_)
			thread.join();
	}
};

int main(int argc, char* argv[])
{
	std::istringstream arg(argv[1]);
	int proc;
	arg >> proc;
	Problem_data data(3);
	Bruteforce bf(proc, &data);
	bf.brute_force();
	data.print_max_pol();
}