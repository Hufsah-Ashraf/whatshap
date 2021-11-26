#include <cassert>
#include "caller.h"
using namespace std;

std::unordered_map<int,int> empty_dict;
std::unordered_map<int,int> none_dict= {{-1, -1}};
std::unordered_map<int,string> hashed;
std::deque<int> variantslist;
int window= 15;
std::deque<std::pair<int,int>> enum_refkmers;
std::deque<std::pair<int,int>> enum_kmers;
Caller::Caller(std::string &reference, int k){

	this->k=k;
	this->enumerate_reference_kmers(reference,k);
	this->ref_kmer_generator = enum_refkmers;
	this->i1= this->ref_kmer_generator.begin();
	//temp_pair1= *this->i1;
	kmer= (*this->i1).first;
	pos=(*this->i1).second;
	if (this->i1 != this->ref_kmer_generator.end()){
		this->i1++;
	}
	this->pileup_columns.push_back(empty_dict);
	this->ref_kmers.push_back(kmer);
	this->ref_pos = pos;
}
void Caller::all_variants(std::deque<int> & variant_list){
	variantslist =variant_list;
}
void Caller::add_read(int bam_alignment_pos, std::vector<std::vector<int>> &bam_alignment_cigartuples, std::string &bam_alignment_query){
	this->target_pos_pre= bam_alignment_pos;
	this->enumerate_kmers(bam_alignment_pos,bam_alignment_query, this->k, bam_alignment_cigartuples);
	this->kmer_generators.push_back(enum_kmers);
	this->kmer_generators_finished.push_back(false);
	int iterator_index = this->kmer_generators.size()-1;
	this->iterators.push_back(this->kmer_generators[iterator_index].begin());
	temp_pair2= *this->iterators[this->iterators.size()-1];
	kmer= temp_pair2.first;
	pos= temp_pair2.second;
	this->current_kmers.push_back(std::pair<int,int> (kmer,pos));
	if(this->iterators[this->iterators.size()-1]!= this->kmer_generators[iterator_index].end()-1){
			this->iterators[this->iterators.size()-1]++;
	}
	std::pair<int, int> pair_getcolumn = this->get_column(pos);
	ref_kmer= pair_getcolumn.first;
	signed index_getcolumn= pair_getcolumn.second;
	if (index_getcolumn>=0){
		this->pileup_columns[index_getcolumn][temp_pair2.first]+=1;
	}
	this->process_complete_columns();
}

void Caller::finish(){
	;
}
std::pair<int, int> Caller::get_column(int pos){

	int index = pos - this->ref_pos;
	if (index >=0){
		while (int (this->pileup_columns.size()) <= index){
			//temp_pair1= *this->i1;
			kmer= (*this->i1).first;
			pos= (*this->i1).second;
			if (this->i1 != this->ref_kmer_generator.end()){
				this->i1++;
			}
			this->ref_kmers.push_back(kmer);
			this->pileup_columns.push_back(empty_dict);
		}
		return std::pair<int,int> (this->ref_kmers[index], index);
		//cout << "exited if" << endl;
	} else {
		//cout << "entered else" << endl;
		//temp_pair1= *this->i1;
		kmer= (*this->i1).first;
		pos= (*this->i1).second;
		if (this->i1 != this->ref_kmer_generator.end()){
				this->i1++;
			}

		this->ref_kmers.push_back(kmer);

		this->pileup_columns.push_back(empty_dict);

		return std::pair<int,int> (kmer,-1);
	}

}
void Caller::pop_column(){
	int result_ref_pos;
	int result_ref_kmer;
	std::unordered_map<int,int> result_pileup_kmers;
	int result_kmer;
	int result_count;
	std::string hashed_ref_kmer;
	std::string hashed_kmer;
	if (int (this->pileup_columns.size()) > 0){
		//cout<<"entered if"<<endl;
		result_ref_pos= this->ref_pos;
		result_ref_kmer= this->ref_kmers.front();
		result_pileup_kmers= this->pileup_columns.front();
		this->ref_kmers.pop_front();
		this->pileup_columns.pop_front();
	}
	else {
		//temp_pair1= *this->i1;
		kmer= (*this->i1).first;
		pos=(*this->i1).second;
		if (this->i1 != this->ref_kmer_generator.end()){
				this->i1++;
			}
		assert (pos == this->ref_pos);
		result_ref_pos= this->ref_pos;
		result_ref_kmer=kmer;
		result_pileup_kmers= none_dict;
	}
	this->ref_pos += 1;
	int variantposition=variantslist.front();
	int varstart= variantposition-window;
	int varend= variantposition+window;
	if (result_ref_pos>=varstart and result_ref_pos<=varend){
		;
	}
	else if (int(variantslist.size())>0 and result_ref_pos>= (variantslist[1]-window) and result_ref_pos<= (variantslist[1]+window)){
			variantslist.pop_front();
	}
	else{
	if (result_pileup_kmers!=none_dict and int (result_pileup_kmers.size()) > 0){
		for(std::unordered_map<int,int> ::iterator it = result_pileup_kmers.begin(); it != result_pileup_kmers.end(); ++it){
			result_kmer=it->first;
			result_count=it->second;
			if (hashed.find(result_ref_kmer)==hashed.end()){
				hashed_ref_kmer= hash_to_dna(result_ref_kmer,this->k);
				hashed[result_ref_kmer]=hashed_ref_kmer;
			}
			else{
				hashed_ref_kmer= hashed[result_ref_kmer];
			}
			if (hashed.find(result_kmer)==hashed.end()){
				hashed_kmer=hash_to_dna(result_kmer, this->k);
				hashed[result_kmer]=hashed_kmer;
			}
			else{
				hashed_kmer=hashed[result_kmer];
		}
			//hashed_ref_kmer= hash_to_dna(result_ref_kmer,this->k);
			//hashed_kmer=hash_to_dna(result_kmer, this->k);
			cout<<result_ref_pos<<"\t"<<hashed_ref_kmer<<"\t"<<hashed_kmer<<"\t"<<result_count<<endl;
	}}}}
 void Caller::process_complete_columns(){
	/*
	Perform calling of columns that are complete, i.e. they cannot receive
	more reads because subsequent reads are further to the right'''
	# compute largest position for which k-mer pileup can safely be generated*/

	int target_pos= this->target_pos_pre + this->k - 1;
	std::vector<std::pair<std::pair<int,int>,std::unordered_map<int,int>>> complete_columns;
	this->advance_to(target_pos);
	//std::deque<std::pair<std::pair<int,int>,std::unordered_map<int,int>>> complete_columns;
	while (this->ref_pos < target_pos){
		this->pop_column();
		}

	}

	//return complete_columns;



void Caller::advance_to(int target_pos){
	/*
	Add all k-mer from all reads up to target_pos to pileup_columns.
	*/
	//print('advance_to', target_pos)
	//cout << "advance_to_started with target_pos " <<target_pos << endl;
	int i;
	//int size_kmer_generators = this->kmer_generators.size();



	for (i=0; i< int (this->kmer_generators.size()) ;i++){
		//cout << "entered for loop with i being "<< i << endl;


		std::deque<std::pair<int,int>> kmer_generator= this->kmer_generators[i];
		//this->i3= kmer_generator.begin();



		try{
			temp_pair3 = this->current_kmers[i];
			kmer = temp_pair3.first;
			pos= temp_pair3.second;
			//this->i3= this->iterators[i];
			//cout<<"within try and outside while kmer is "<<kmer<<endl;
			//cout<<"within try and outside while pos is "<<pos<<endl;
			while (pos < target_pos){
				//cout<<"entered while"<<endl;

				//temp_pair3= *this->i3;
				if(this->iterators[i]!= this->kmer_generators[i].end()){
					temp_pair3= *this->iterators[i];
					kmer=temp_pair3.first;
					pos=temp_pair3.second;
					//cout<<"within while kmer is "<<kmer<<endl;
					//cout<<"within while pos is "<<pos<<endl;
					//this->get_column(pos).second[kmer]+=1;
					//pileup_column[kmer] += 1;
					//std::pair<int, int> (ref_kmer, index_getcolumn) = this->get_column(pos);
					std::pair<int, int> pair_getcolumn = this->get_column(pos);
					ref_kmer= pair_getcolumn.first;
					signed index_getcolumn= pair_getcolumn.second;
					if (index_getcolumn>=0){
						this->pileup_columns[index_getcolumn][temp_pair3.first]+=1;
						//cout<<"advance_to "<< this->ref_pos<< " pileup addition at "<<temp_pair3.first <<" with count being "<<this->pileup_columns[index_getcolumn][temp_pair3.first]<<endl;
						//pileup_column=this->pileup_columns[index_getcolumn];
					}

				////cout<<"ref_kmer is "<< ref_kmer<<" and kmer is  "<< temp_pair3.first<<endl;
				////for(std::unordered_map<int,int> ::iterator it = pileup_column.begin(); it != pileup_column.end(); ++it){
				////	cout << it->first << " " << it->second << endl;

				////	}


					this->iterators[i]++;
				}
				else{

					throw(1);
				}


			}

			this->current_kmers[i] = std::pair<int,int> (kmer,pos);
			//cout<< kmer<<" "<<pos <<" is added to current_kmers of "<< i <<endl;

		}
		catch(int e){
			//this->current_kmers[i] = std::make_pair(kmer,pos);

			//cout<<"entered except "<<endl;
			this->kmer_generators_finished[i] = true;

		}


	}

	while ((this->kmer_generators.size() > 0) && this->kmer_generators_finished[0]){
						//std::cout<<"current_kmers popping "<< this->current_kmers.front().first<<" "<<this->current_kmers.front().second<<endl;
                        this->current_kmers.pop_front();
                        //this->bam_records.pop_back();
						//std::cout<<"kmer_generators popping "<< this->kmer_generators.front().first<<" "<<this->kmer_generators.front().second<<endl;
                        this->kmer_generators.pop_front();
						this->iterators.pop_front();
                        this->kmer_generators_finished.pop_front();
	}
}
void Caller::enumerate_reference_kmers(std::string &reference, int k){

	int A = 'A';
        int C = 'C';
        int G = 'G';
        int T = 'T';
        int c = 0;
        int h = 0;
        int mask = (1 << (2*k)) - 1;
        int i=0;

				enum_refkmers.clear();

	int length_reference = reference.length();

        for (i=0; i<length_reference; i++){
                c = reference[i];
                if (c == A){
                        h = ((h << 2) | 0) & mask;
				}
		else if (c == C){
                        h = ((h << 2) | 1) & mask;
				}
		else if (c == G){
                        h = ((h << 2) | 2) & mask;
				}
		else if (c == T){
                        h = ((h << 2) | 3) & mask;
				}
                else{
                        h = ((h << 2) | 0) & mask;
				}
                if (i >= k-1){
                        enum_refkmers.push_back(std::pair<int,int> (h,i+1));
                        //yield(h, i+1)
				}
		}
        //return enum_refkmers;
		//cout << "enumerate_reference_kmers_end" << endl;
}

void Caller::enumerate_kmers(int pos, std::string &query_string, int k, std::vector<std::vector<int>> &cigartuples){

        int h = 0;
        int n_kmer=0;
        int mask = (1 << (2*k)) - 1;
        int cigar_index = 0;
        int cigar_op = cigartuples[cigar_index][0];
        int cigar_length = cigartuples[cigar_index][1];
        //cdef int cigar_length = cigartuples[cigar_index][1]
        int BAM_CMATCH = 0;     // M
        int BAM_CINS = 1 ;      // I
        int BAM_CDEL = 2 ;      // D
        int BAM_CREF_SKIP = 3;  // N
        int BAM_CSOFT_CLIP = 4; // S
        int BAM_CHARD_CLIP = 5; // H
        int BAM_CPAD = 6 ;     // P
        int BAM_CEQUAL = 7 ;   //=
        int BAM_CDIFF = 8  ;    // X
        int BAM_CBACK = 9 ;     // B
        int A = 'A';
        int C = 'C';
        int G = 'G';
        int T = 'T';
        int i =0;
        int consecutive = 0;

				enum_kmers.clear();
		int length_query_string= query_string.length();
        while (i < length_query_string){
                //process cigar entries that don't consume a character from the read
                while (true){
                        if ((cigar_op == BAM_CDEL) or (cigar_op == BAM_CREF_SKIP)){
                                //print(pos, '-->', pos + cigar_length)
                                pos += cigar_length;
						}
						else if (cigar_op == BAM_CSOFT_CLIP){
                                //i += cigar_length
                                consecutive = 0;
						}
						else if ((cigar_length == 0) or (cigar_op == BAM_CHARD_CLIP)){

						}
                        else{
                                break;
						}
                        cigar_index += 1;
						cigar_op = cigartuples[cigar_index][0];
						cigar_length = cigartuples[cigar_index][1];
                        //cigar_op, cigar_length = cigartuples[cigar_index];
				}
                if (i >= length_query_string){
                        break;
				}

                // update hash

                if (query_string[i] == A){
                        h = ((h << 2) | 0) & mask;
                        consecutive += 1;
				}
				else if (query_string[i] == C){
                        h = ((h << 2) | 1) & mask;
                        consecutive += 1;
				}
				else if (query_string[i] == G){
                        h = ((h << 2) | 2) & mask;
                        consecutive += 1;
				}
				else if (query_string[i] == T){
                        h = ((h << 2) | 3) & mask;
                        consecutive += 1;
				}
                else {
                        n_kmer=1;
                        consecutive += 1;
                        //print( "N encountered at", i)
				}
                if (consecutive >= k){
                        if (n_kmer ==0){
								//cout <<"enum yielded "<< h << " and " <<pos+1<<endl;
								//cout<<" pos+1 "<<pos+1 <<endl;
                                enum_kmers.push_back(std::pair<int,int> (h,pos+1));
                                //yield (h, pos+1)
						}
						else if (n_kmer==1){
                                //print (h, pos+1, "this kmer had N")
                                n_kmer=0;
						}
				}

                //consume one character of read
                assert (cigar_length > 0);
                if ((cigar_op == BAM_CMATCH) or (cigar_op == BAM_CEQUAL) or (cigar_op == BAM_CDIFF)){
                        cigar_length -= 1;
                        pos += 1;
				}
				else if (cigar_op == BAM_CINS){
                        cigar_length -= 1;
				}
                else {
                        assert (false);
				}

                i += 1;
				}
        //return enum_kmers;

		}

std::string Caller::hash_to_dna(int h, int k){
        std::vector<char> l(k);
		int hashed;
		//cout<<"hash  started with h "<<h<<endl;
        int l_i;
        for (l_i =k-1;l_i>=0;l_i--){
				//cout<<"within for h is "<<h<<" and h&3 is "<< (h&3)<<endl;
                if ((h&3) == 0){
						//cout<<"entered if "<<endl;
                        l[l_i]='A'  ;
                                }
                else if ((h&3) == 1){
						//cout<<"entered firs else if"<<endl;
                        l[l_i]='C'  ;
                                }
                else if  ((h&3) == 2){
						//cout<<"entered seconond"<<endl;
                        l[l_i]='G';
                                }
                else if ((h&3) == 3){
						//cout<<"entered third"<<endl;
                        l[l_i]='T';
					}

                h = (h >> 2);
        }
		std::string hashed_string=" ";
		for (hashed=0; hashed<int (l.size()); hashed++){
			hashed_string+= l[hashed];
		}
		return hashed_string;
}
