#include <iostream>
#include <random>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <time.h>

using namespace std;

struct read
{
	int start;
	int last;
	int haplotype;
};

void read_pos_map(string filename, map<int, int>& pos_map);
void read_genotype(string filename, vector<vector<int> >& data);
void eliminate_homozygous(vector<vector<int> >& data, map<int, int>& pos_map);
int normal_dist_number(double mean, double dev);
int uniform_dist_number(int start, int last);
void gen_reads(vector<struct read>& reads, map<int, int> heter_pos_map);
void output_reads(string filename, vector<struct read>& reads, vector<int>& haplotype1, vector<int>& haplotype2);

#define NO_READS 500
#define ERR_RATE 100
#define MEAN 1000
#define VARIANCE 100
#define READ_LENGTH 10000
#define NO_SNP 34026
#define NO_PEOPLE 60
#define GENE_START_POS 14431347
#define GENE_LAST_POS 49519949

int main()
{
	srand(time(NULL));

	vector<vector<int> > raw_data;
	map<int, int> raw_pos_map;

	read_pos_map("genotypes_chr22_CEU_r21_nr_fwd_legend.txt", raw_pos_map);
	read_genotype("genotypes_chr22_CEU_r21_nr_fwd_phased", raw_data);
	eliminate_homozygous(raw_data, raw_pos_map);

	// reads
	for (int i=0; i<NO_PEOPLE; i++)
	{
		vector<vector<int> > heter_data;
		string heter_data_filename = "reads/phase"+to_string(i);
		read_genotype(heter_data_filename, heter_data);

		map<int, int> heter_pos_map;
		string heter_legend_filename = "reads/legend"+to_string(i);
		read_pos_map(heter_legend_filename, heter_pos_map);

		vector<struct read> reads;
		vector<int> haplotype1 = heter_data[i*2];
		vector<int> haplotype2 = heter_data[i*2+1];

		gen_reads(reads, heter_pos_map);

		string reads_filename = "reads/read"+to_string(i);
		output_reads(reads_filename, reads, haplotype1, haplotype2);
	}

//	vector<int> chrom1 = data[0];
//	vector<int> chrom2 = data[1];
//	gen_reads(chrom1, reads);
//	output_reads("reads", reads);
}

void read_pos_map(string filename, map<int, int>& pos_map)
{
	ifstream file;
	file.open(filename);
	string line;
	// skip first line
	getline(file, line);
	int idx = 0;
	while (getline(file, line))
	{
		stringstream sstream(line);
		string pos;
		sstream>>pos>>pos;
		pos_map[idx] = atoi(pos.c_str());
		idx++;
	}
	file.close();
}

void eliminate_homozygous(vector<vector<int> >& data, map<int, int>& pos_map)
{
	for (int i=0; i<NO_PEOPLE; i++)
	{
		vector<int> chrome1 = data[2*i];
		vector<int> chrome2 = data[2*i+1];
		map<int, int> heter_pos_map;
		for (int j=0; j<NO_SNP; j++)
		{
			if (chrome1[j] != chrome2[j])
				heter_pos_map[j] = pos_map[j];
		}
		string out_phase_filename = "reads/phase"+to_string(i);
		string out_legend_filename = "reads/legend"+to_string(i);
		cout<<out_phase_filename<<endl;
		ofstream ophase_file;
		ofstream olegend_file;
		ophase_file.open(out_phase_filename);
		olegend_file.open(out_legend_filename);
		for (int j=0; j<2*NO_PEOPLE; j++)
		{
			vector<int> chrome = data[j];
			map<int, int>::iterator itr = heter_pos_map.begin();
			while (itr != heter_pos_map.end())
			{
				int phase_idx = itr->first;
				int phase_val = chrome[phase_idx];
				ophase_file<<phase_val<<" ";
				++itr;
			}
			ophase_file<<endl;
		}
		ophase_file.close();
		// output legend
		map<int, int>::iterator itr = heter_pos_map.begin();
		int idx = 0;
		while (itr != heter_pos_map.end())
		{
			olegend_file<<idx<<"\t"<<itr->second<<endl;
			++itr;
			++idx;
		}
		olegend_file.close();
	}
}

int normal_dist_number(double mean, double dev)
{
	default_random_engine generator;
	generator.seed(rand());
	normal_distribution<double> distribution(mean, dev);
	return distribution(generator);
}

int uniform_dist_number(int start, int last)
{
	default_random_engine generator;
	generator.seed(rand());
	uniform_int_distribution<int> distribution(start, last);
	return distribution(generator);
}

void gen_reads(vector<struct read>& reads, map<int, int> heter_pos_map)
{
	int no_reads = 0;
	while (no_reads<NO_READS)
	{
//		cout<<"test1\n";
		int snp_start_pos = uniform_dist_number(GENE_START_POS, GENE_LAST_POS-READ_LENGTH+1);
		int snp_last_pos = snp_start_pos+READ_LENGTH-1;
		int haplo_idx = rand()%2;
		map<int, int>::iterator itr = heter_pos_map.begin();
		int start_idx = -1, last_idx = -1, prev_idx = -1;
		bool find_first = true;
//		cout<<snp_start_pos<<" "<<snp_last_pos<<endl;
		while (itr != heter_pos_map.end())
		{
			if (find_first)
			{
				if (itr->second > snp_start_pos)
				{
					start_idx = itr->first;
					prev_idx = start_idx;
					find_first = false;
				} else
				{
					++itr;
				}
			} else
			{
				if (itr->second > snp_last_pos)
				{
					last_idx = prev_idx;
					break;
				} else {
					prev_idx = itr->first;
					++itr;
				}
			}
		}
		if (start_idx != -1 && last_idx == -1)
			last_idx = prev_idx;
		struct read read;
		read.start = start_idx;
		read.last = last_idx;
		read.haplotype = haplo_idx+1;
		reads.push_back(read);
		no_reads++;
		if (no_reads == NO_READS)
			break;
		int gap = normal_dist_number(MEAN, VARIANCE);
		gap = gap>=0?gap:0;
		snp_start_pos = snp_last_pos+gap+1;
		if (snp_start_pos>GENE_LAST_POS-READ_LENGTH+1)
			continue;
		snp_last_pos = snp_start_pos+READ_LENGTH-1;

		find_first = true;
		while (itr != heter_pos_map.end())
		{
			if (find_first)
			{
				if (itr->second > snp_start_pos)
				{
					start_idx = itr->first;
					prev_idx = start_idx;
					find_first = false;
				} else
				{
					++itr;
				}
			} else
			{
				if (itr->second > snp_last_pos)
				{
					last_idx = prev_idx;
					break;
				} else {
					prev_idx = itr->first;
					++itr;
				}
			}
		}
		if (start_idx != -1 && last_idx == -1)
			last_idx = prev_idx;

		struct read read2;
		read2.start = start_idx;
		read2.last = last_idx;
		read2.haplotype = haplo_idx+1;
		reads.push_back(read);
		no_reads++;
	}
	cout<<"test2\n";
}

void read_genotype(string filename, vector<vector<int> >& data)
{
	ifstream file;
	file.open(filename);
	string line;
	while (getline(file, line))
	{
		vector<int> chrome;
		// parse 0-1 binary sequence to chrome
		stringstream sstream(line);
		string val;
		while (sstream.good())
		{
			sstream>>val;
			chrome.push_back(atoi(val.c_str()));
		}
		data.push_back(chrome);
	}
	file.close();
}

void output_reads(string filename, vector<struct read>& reads, vector<int>& haplotype1, vector<int>& haplotype2)
{
	ofstream file;
	file.open(filename);
	int length = haplotype1.size();
	for (int i=0; i<reads.size(); i++)
	{
		struct read read = reads[i];
		int start = read.start;
		int last = read.last;
		vector<int> haplotype = read.haplotype==1?haplotype1:haplotype2;
		for (int j=0; j<start; j++)
			file<<"- ";
		for (int j=start; j<=last; j++)
		{
			// simulate some error
			int err_flag = rand()%ERR_RATE;
			if (err_flag == 0)
				file<<1-haplotype[j]<<" ";
			else
				file<<haplotype[j]<<" ";
		}
		for (int j=last+1; j<length; j++)
			file<<"- ";
		file<<endl;
	}
	file.close();
}
