#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "PF_Tau.hpp"

using namespace std;
#define events_in_file (1)
#define tokens_per_event (50)


vector< vector<unsigned int> >  read_tracks_from_file(void){
	ifstream inFile;
	inFile.open("inputPatterns.txt");//open the input file

	stringstream strStream;
	strStream << inFile.rdbuf();//read the file
	string str = strStream.str();//str holds the content of the file

	string buf; // Have a buffer string
	stringstream ss(str); // Insert the string into a stream

	vector<string> tokens; // Create vector to hold our words

	while (ss >> buf)
		tokens.push_back(buf);

	vector< vector<unsigned int> > v;

	for (int i = 0; i < events_in_file; i++)
	{
		if (tokens[i * tokens_per_event] == "#Event")
		{
			cout << "Processing event #: " << tokens[i * tokens_per_event + 1] << endl;
		}
		//fix me
		vector<unsigned int> tracks;
		for (int j = 0; j < N_TRACKS; j++)
		{
			unsigned int idx = i * tokens_per_event + 2 + j;

			std::stringstream ss2;
			unsigned int x;

			ss2 << hex << tokens[idx].c_str() ;
			ss2 >> x;

			tracks.push_back(x);
		}
		v.push_back(tracks);
	}
	return v;
}

vector< vector<unsigned int> >  read_clusters_from_file(void){
	ifstream inFile;
	inFile.open("inputClusters.txt");//open the input file

	stringstream strStream;
	strStream << inFile.rdbuf();//read the file
	string str = strStream.str();//str holds the content of the file

	string buf; // Have a buffer string
	stringstream ss(str); // Insert the string into a stream

	vector<string> tokens; // Create vector to hold our words

	while (ss >> buf)
		tokens.push_back(buf);

	vector< vector<unsigned int> > v;

	for (int i = 0; i < events_in_file; i++)
	{
		if (tokens[i * tokens_per_event] == "#Event")
		{
			cout << "Processing event #: " << tokens[i * tokens_per_event + 1] << endl;
		}
		//fix me
		vector<unsigned int> tracks;
		for (int j = 0; j < N_CLUSTERS; j++)
		{
			unsigned int idx = i * tokens_per_event + 2 + j;

			std::stringstream ss2;
			unsigned int x;

			ss2 << hex << tokens[idx].c_str() ;
			ss2 >> x;

			tracks.push_back(x);
		}
		v.push_back(tracks);
	}
	return v;
}

// are tracks to be sorted by eta/phi?
void unpack_input_tracks(vector<unsigned int>  region_raw,
						track_t central_track [N_TRACKS]){
	//fix me for tracks

	for (unsigned int idx_out = 0; idx_out < N_TRACKS; idx_out++)
	{

		//unsigned int idx_in = 468 - ( phi + 1 ) * 26 + reg;

		central_track[idx_out].et       =   (region_raw[idx_out] & 0x7FF)     >> 0;
		central_track[idx_out].eta      =   (region_raw[idx_out] & 0x7FFFF)   >> 11;
		central_track[idx_out].eta_side =   (region_raw[idx_out] & 0xFFFFF)   >> 19;
		central_track[idx_out].phi      =   (region_raw[idx_out] & 0xFFFFFFF) >> 20;


	}

}

// are tracks to be sorted by eta/phi?
void unpack_input_clusters(vector<unsigned int>  region_raw,
						cluster_t central_cluster [N_CLUSTERS]){
	//fix me for tracks

	for (unsigned int idx_out = 0; idx_out < N_CLUSTERS; idx_out++)
	{

		//unsigned int idx_in = 468 - ( phi + 1 ) * 26 + reg;

		central_cluster[idx_out].et       =   (region_raw[idx_out] & 0x7FF)     >> 0;
		central_cluster[idx_out].eta      =   (region_raw[idx_out] & 0x7FFFF)   >> 11;
		central_cluster[idx_out].eta_side =   (region_raw[idx_out] & 0xFFFFF)   >> 19;
		central_cluster[idx_out].phi      =   (region_raw[idx_out] & 0xFFFFFFF) >> 20;

	}

}

int main()
{
	vector< vector< unsigned int> > tracks_raw = read_tracks_from_file();
	vector< vector< unsigned int> > clusters_raw = read_clusters_from_file();

	track_t central_tracks[N_TRACKS];
	cluster_t central_clusters[N_CLUSTERS];

	for (int i = 0; i < events_in_file; i++)
	{
		unpack_input_tracks(tracks_raw[i], central_tracks);
		unpack_input_clusters(clusters_raw[i], central_clusters);

	}

	algo_config_t algo_config;

	algo_config.three_prong_seed = 0x46;//corresponds to 7GeV

	algo_config.three_prong_delta_r = 0x6;

	algo_outputs_t algo_outputs;

	file_read_in(central_tracks, central_clusters, algo_config, algo_outputs);

	cout<<"Sum Pt "<<algo_outputs.sum_tracks<<std::endl;
}
