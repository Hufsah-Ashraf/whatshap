#include <cassert>
#include "caller.h"
using namespace std;

std::map<int,int> empty_dict;
Caller::Caller(std::string reference, int k){
	this->k=k;
	//bam record for each active read
	//std::deque<int> this->bam_records;
	//k-mer generators for each active read
	//std::deque<int> this->kmer_generators;
	//std::deque<int> this->kmer_generators_finished;
	//.. corresponding queue with the latest (kmer, position) for each read
	//std::deque<int> this->current_kmers;
	// TODO provide implementation for enumerate_reference_kmers
	this->ref_kmer_generator = this->enumerate_reference_kmers(reference,k);
	//position that every other operation is relative to
	//std::deque<int> this->pileup_columns;
	//std::deque<int> this->ref_kmers;
	this->i1= this->ref_kmer_generator.begin();
	temp_pair1= *this->i1;
	kmer= temp_pair1.first;
	pos=temp_pair1.second;
	if (this->i1 != this->ref_kmer_generator.end()){
		this->i1++;
	}
	//cout<<"initial kmer"<<kmer<<endl;
	//cout<<"initial pos"<<pos<<endl;
	this->pileup_columns.push_back(empty_dict);
	this->ref_kmers.push_back(kmer);
	this->ref_pos = pos;
	this->i3_init=0;
	
}
void Caller::add_read(int bam_alignment_pos, std::vector<std::vector<int>> bam_alignment_cigartuples, std::string bam_alignment_query){
	this->enum_kmer_i=0;
	this->consecutive=0;
	this->cigar_index=0;
	this->n_kmer=0;
	//cout << "add_read_start" << endl;
	//this->bam_records.push_back(bam_alignment_pos);
	this->target_pos_pre= bam_alignment_pos;
	this->enum_kmer_pos= bam_alignment_pos;
	this ->cigar_bam = bam_alignment_cigartuples;
	this->query_bam = bam_alignment_query;
	// TODO provide implementation for enumerate_kmers
	this->kmer_generators.push_back(this->enumerate_kmers(bam_alignment_query, this->k, bam_alignment_cigartuples));
	//cout<<"one call to enumerate_kmers done"<<endl;
	this->kmer_generators_finished.push_back(false);
	if(this->i3_init==0){
		this->i3= this->kmer_generators.begin();
		this->i3_init=1;
	}
	
	//temp_deque= *this->i2;
	//for (std::deque<std::pair<int,int>>::iterator it = this->kmer_generators.begin(); it!=this->kmer_generators.end(); ++it){
		//std::pair<int,int> temp = *it;
		//std::cout << ' ' << temp.first << ' ' <<temp.second;
		//std::cout << '\n';
		//}
	//this->i2a= temp_deque.begin();
	temp_pair2= *this->i3;
	kmer=temp_pair2.first;
	pos=temp_pair2.second;
	//if (this->i2 != this->kmer_generators.end()){
		//		this->i2++;
			//}
	if (this->i3 != this->kmer_generators.end()){
			this->i3++;
			}
	this->current_kmers.push_back(std::make_pair(kmer,pos));
	
	std::pair<int,map<int,int>> (ref_kmer, pileup_column) = this->get_column(pos);
	pileup_column[kmer] += 1;
}
void Caller::finish(){
	;
}

std::pair<int, std::map<int,int>> Caller::get_column(int pos){	
	int index = pos - this->ref_pos;
	int size_pileup_columns = this->pileup_columns.size();
	if (index >=0){
		while (size_pileup_columns <= index){
			temp_pair1= *this->i1;
			kmer= temp_pair1.first;
			pos=temp_pair1.second;
			if (this->i1 != this->ref_kmer_generator.end()){
				this->i1++;
			}
			
			this->ref_kmers.push_back(kmer);
			this->pileup_columns.push_back(empty_dict);
			size_pileup_columns = this->pileup_columns.size();
			
		} 
	
		return std::pair<int, std::map<int,int>> (this->ref_kmers[index], this->pileup_columns[index]);

	} else {
		
		temp_pair1= *this->i1;
		kmer= temp_pair1.first;
		pos=temp_pair1.second;
		if (this->i1 != this->ref_kmer_generator.end()){
				this->i1++;
			}
		this->ref_kmers.push_back(kmer);
		this->pileup_columns.push_back(empty_dict);
		return std::pair<int, std::map<int,int>> (kmer,empty_dict);
	}
	
}
std::pair<std::pair<int,int>,std::map<int,int>> Caller::pop_column(){
	
	if (this->pileup_columns.size() > 0){
		result = std::make_pair(std::make_pair(this->ref_pos, this->ref_kmers.front()),this->pileup_columns.front());
		this->ref_kmers.pop_front();
		this->pileup_columns.pop_front();
	} else {
		temp_pair1= *this->i1;
		kmer= temp_pair1.first;
		pos=temp_pair1.second;
		if (this->i1 != this->ref_kmer_generator.end()){
				this->i1++;
			}
		assert (pos == this->ref_pos);
		result = std::pair<std::pair<int,int>,std::map<int,int>> (std::pair<int,int>(this->ref_pos, kmer), empty_dict);
	}
	this->ref_pos += 1;
	return result;
}
std::vector<std::pair<std::pair<int,int>,std::map<int,int>>> Caller::process_complete_columns(){
	/*
	Perform calling of columns that are complete, i.e. they cannot receive
	more reads because subsequent reads are further to the right'''
	# compute largest position for which k-mer pileup can safely be generated*/
	//int target_pos = *this->bam_records.begin() + this->k - 1;
	int target_pos= this->target_pos_pre + this->k - 1;
	this->advance_to(target_pos);
	std::vector<std::pair<std::pair<int,int>,std::map<int,int>>> complete_columns;
	while (this->ref_pos < target_pos){
		complete_columns.push_back(this->pop_column());
	}
	return complete_columns;
	
}
        
void Caller::advance_to(int target_pos){
	/*
	Add all k-mer from all reads up to target_pos to pileup_columns.
	*/
	//cout << "advance_to_started with target_pos " <<target_pos << endl;
	int size_kmer_generators = this->kmer_generators.size();
	//cout<<"length of self.kmer_generators is "<< size_kmer_generators<<endl;

	for (i=0;i<size_kmer_generators;i++){
		//cout << "entered for loop with i being "<< i << endl;
		//std::deque<std::pair<int,int>> kmer_generator;
		//std::deque<std::pair<int,int>>::iterator i3= kmer_generator.begin();
		//kmer_generator= this->kmer_generators[i];
		//this->i3= kmer_generators.begin();
		//if (this->i3 != kmer_generators.end()){
			//		this->i3++;
				//}
		
			//cout << "entered try" << endl;
			
		temp_pair3 = this->current_kmers[i];
		kmer = temp_pair3.first;
		pos= temp_pair3.second;
		cout<< kmer<<" "<<pos<<" " <<"came out of current_kmers of "<< i <<endl;
		cout << "within try and outside while kmer is "<<kmer << endl;
		cout << "within try and outside while pos is "<<pos << endl;		
		while (pos < target_pos){
			//cout<<"entered while"<<endl;
			this->kmer_generators.push_back(this->enumerate_kmers(this->query_bam, this->k, this->cigar_bam));
			int size_kmer_generators = this->kmer_generators.size();
			temp_pair3= *this->i3;
			kmer=temp_pair3.first;
			pos=temp_pair3.second;
			//cout<<"within while kmer is "<<kmer<<endl;
			//cout<<"within while pos is "<<pos<<endl;
			std::make_pair(ref_kmer, pileup_column) = this->get_column(pos);
			pileup_column[kmer] += 1;
			if (this->i3 != this->kmer_generators.end()-1){
				this->i3++;
			}
			else{
				this->kmer_generators_finished[i] = true; 
				break;
			}
		}
	
		if (this->kmer_generators_finished[i] != true){
			this->current_kmers[i] = std::make_pair(kmer,pos);
			cout<< kmer<<" "<<pos <<" is added to current_kmers of "<< i <<endl;
		}
		else{
			cout<<"entered except"<<endl;
			}
		
			
	}
	cout<<"self.current_kmers.front "<<this->current_kmers.front().first<<" ," <<this->current_kmers.front().second<< " and self.current_kmers.back "<<this->current_kmers.back().first<<" , "<<this->current_kmers.back().second<<endl;
	while ((size_kmer_generators > 0) && this->kmer_generators_finished[0]){
						cout<<"current_kmers popping "<<this->current_kmers.front().first<<","<<this->current_kmers.front().second<<endl;
                        this->current_kmers.pop_front();
                        //this->bam_records.pop_back();
						cout<<"kmer_generators popping "<< this->kmer_generators.front().first<<","<<this->kmer_generators.front().second<<endl;
                        this->kmer_generators.pop_front();
                        this->kmer_generators_finished.pop_front();
	}
	cout<<"at advance to end kmer is "<<kmer<<endl;
	cout<<"at advance to end pos is "<<pos<<endl;
	cout << "advance_to_end" << endl;
}
std::deque<std::pair<int,int>> Caller::enumerate_reference_kmers(std::string reference, int k){
       
	int A = 'A';
        int C = 'C';
        int G = 'G';
        int T = 'T';
        int c = 0;
        int h = 0;
        int mask = (1 << (2*k)) - 1;
        int i=0;
        std::deque<pair<int,int>> enum_refkmers;
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
                        enum_refkmers.push_back(std::make_pair(h,i+1));
                        //yield(h, i+1)
				}
		}
        return enum_refkmers;
		cout << "enumerate_reference_kmers_end" << endl;
}

std::pair<int,int> Caller::enumerate_kmers(std::string query_string, int k, std::vector<std::vector<int>> cigartuples){
		
        int h = 0;
        int mask = (1 << (2*k)) - 1;
        //int cigar_index = 0;
        int cigar_op = cigartuples[this->cigar_index][0];
        int cigar_length = cigartuples[this->cigar_index][1];
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
        //int i =0;
        //int consecutive = 0;
        std::pair<int,int> enum_kmers;
		int length_query_string= query_string.length();
        while (this->enum_kmer_i < length_query_string){
                //process cigar entries that don't consume a character from the read
				//cout<<"i "<<this->enum_kmer_i<<" consecutive "<<this->consecutive <<" cigar_index "<< this->cigar_index<<endl;
                while (true){
                        if ((cigar_op == BAM_CDEL) or (cigar_op == BAM_CREF_SKIP)){
                                //print(pos, '-->', pos + cigar_length)
                                this->enum_kmer_pos += cigar_length;
						}
						else if (cigar_op == BAM_CSOFT_CLIP){
                                //i += cigar_length
                                this->consecutive = 0;
						}
						else if ((cigar_length == 0) or (cigar_op == BAM_CHARD_CLIP))
                                ;
						
                        else
							break;
						
                        this->cigar_index += 1;
						cigar_op = cigartuples[this->cigar_index][0];
						cigar_length = cigartuples[this->cigar_index][1];
						cout<<"cigar_op "<<cigar_op<<" cigar_length "<<cigar_length<<endl;
                        //cigar_op, cigar_length = cigartuples[cigar_index];
				}
                if (this->enum_kmer_i >= length_query_string){
                        break;
				}

                // update hash
        
                if (query_string[this->enum_kmer_i] == A){
                        h = ((h << 2) | 0) & mask;
                        this->consecutive += 1;
				}
				else if (query_string[this->enum_kmer_i] == C){
                        h = ((h << 2) | 1) & mask;
                        this->consecutive += 1;
				}
				else if (query_string[this->enum_kmer_i] == G){
                        h = ((h << 2) | 2) & mask;
                        this->consecutive += 1;
				}
				else if (query_string[this->enum_kmer_i] == T){
                        h = ((h << 2) | 3) & mask;
                        this->consecutive += 1;
				}
                else {
                        this->n_kmer=1;
                        this->consecutive += 1;
                        //print( "N encountered at", i)
				}   
                if (this->consecutive >= k){
                        if (this->n_kmer ==0){
								
								//cout<<" pos+1 "<<pos+1 <<endl;
                                //enum_kmers.push_back(std::make_pair(h,pos+1));
                                enum_kmers= std::make_pair(h, this->enum_kmer_pos+1);
								//cout <<"enum yielded "<< h << " and " <<this->enum_kmer_pos+1<<endl;
								//return enum_kmers;
						}
						else if (this->n_kmer==1){
                                //print (h, pos+1, "this kmer had N")
                                this->n_kmer=0;
								//enum_kmers= std::make_pair(0,0);
								//return enum_kmers;
								
						}
						assert (cigar_length > 0);
						if ((cigar_op == BAM_CMATCH) or (cigar_op == BAM_CEQUAL) or (cigar_op == BAM_CDIFF)){
								cigar_length -= 1;
								this->enum_kmer_pos += 1;
						}
						else if (cigar_op == BAM_CINS){
								cigar_length -= 1;
						}
						else {
								assert (false);
						}

						this->enum_kmer_i += 1;
						return enum_kmers;
				}
                else{
					//consume one character of read
					assert (cigar_length > 0);
					if ((cigar_op == BAM_CMATCH) or (cigar_op == BAM_CEQUAL) or (cigar_op == BAM_CDIFF)){
							cigar_length -= 1;
							this->enum_kmer_pos += 1;
					}
					else if (cigar_op == BAM_CINS){
							cigar_length -= 1;
					}
					else {
							assert (false);
					}

					this->enum_kmer_i += 1;
					
				}
				
		}
		return enum_kmers;
			
}
                






