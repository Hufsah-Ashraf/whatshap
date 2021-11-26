#ifndef CALLER_H
#define CALLER_H
#include<iostream>
#include <deque>
#include <string>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>
#include <thread>
#include <utility>
class Caller{
public:
	Caller(std::string reference, int k);
	void add_read(int bam_alignment_pos,std::vector<std::vector<int>> bam_alignment_cigartuples, std::string bam_alignment_query);
	void finish();
	std::pair<int, std::map<int,int>> get_column(int pos);
	std::pair<std::pair<int,int>,std::map<int,int>> pop_column();
	std::vector<std::pair<std::pair<int,int>,std::map<int,int>>> process_complete_columns();
	void advance_to(int target_pos);
	std::deque<std::pair<int,int>> enumerate_reference_kmers(std::string reference, int k);
	std::pair<int,int> enumerate_kmers(std::string query_string, int k, std::vector<std::vector<int>> cigartuples);
	int i3_init;
	int i;
	int kmer;
	int pos;
	int k;
	int ref_pos;
	int target_pos_pre;
	int ref_kmer;
	int enum_kmer_i;
	int consecutive;
	int cigar_index;
	int enum_kmer_pos;
	int n_kmer;
	std::string query_bam;
	std::vector<std::vector<int>> cigar_bam;
	//std::deque<int> bam_records;
	std::deque<std::pair<int,int>> kmer_generators;
	std::deque<int> kmer_generators_finished;
	std::deque<std::pair<int,int>> current_kmers;
	std::deque<std::pair<int,int>> ref_kmer_generator;
	std::deque<std::map<int,int>> pileup_columns;
	std::deque<int> ref_kmers;
	std::map<int,int> pileup_column;
	//std::deque<std::deque<std::pair<int,int>>>::iterator i2;
	std::deque<std::pair<int,int>>::iterator i2a;
	std::deque<std::pair<int,int>>::iterator i1;
	std::deque<std::pair<int,int>>::iterator i3;
	std::pair<std::pair<int,int>,std::map<int,int>> result;
	std::pair<int,int> temp_pair1;
	//std::deque<std::pair<int,int>> temp_deque;
	std::pair<int,int> temp_pair2;
	std::pair<int,int> temp_pair3;
	
};
#endif
