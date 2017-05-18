#include "file_read_in.hpp"

void file_read_in(
		track_t central_tracks[N_TRACKS],  // Number of tracks
		cluster_t central_clusters[N_CLUSTERS],  // Number of Clusters
		algo_config_t algo_config,           // algoritm configuration
		algo_outputs_t & algo_outputs        // algorithm outputs
		)
{
#pragma HLS ARRAY_PARTITION variable=central_tracks complete dim=1

	ap_uint<14> et_total_tmp = 0;

#pragma HLS PIPELINE II=6

	pfcharged_t charged_cands[N_TRACKS];

	////////////////////   PF Particle Matching   ////////////////////
	pf_match_alg( central_clusters, central_tracks, charged_cands );

	////////////////////Three Prong Tau Algorithm////////////////////
	track_t three_prong_tau_cand[3];
	tau_three_prong_alg( central_tracks, three_prong_tau_cand, algo_config);

	////////////////////    Algorithm Outputs    ////////////////////
	algo_outputs.three_prong_tau_et = three_prong_tau_cand[0].et + three_prong_tau_cand[1].et + three_prong_tau_cand[2].et;


}

/*
 * track index is from -20 to 20 in eta (excluding 0)
 * then from 0 to 71 in phi
 * -20	                        -1 | 1                  +20
 * |0   72  72*2      (20-eta)+phi |                       |
 * |1	73  72*2+1                 |                       |
 * |2	74  72*3+1                 |                       |
 */
void pf_match_alg(cluster_t central_clusters[N_CLUSTERS], track_t central_tracks[N_TRACKS] , pfcharged_t charged_cands[N_TRACKS]){
  
	for(int jdx = 0; jdx < N_TRACKS; jdx++)	//note, tracks are already sorted by PT
	  {
	    //#pragma HLS UNROLL
	    chargd_cands[jdx] = central_track[jdx];
	    int index = 0;
	    
	    index = find_the_index_crys( charged_cands.eta, charged_cands.eta_side, charged_cands.phi);
	    
	    // Take the Cluster and use H/E and E/H to determine if Hadron or electron/pi0
	    if( (central_clusters[index]).EoH > input_EoH_cut_ ){
   	       charged_cands[jdx].is_charged_hadron = 0;
	       charged_cands[jdx].is_electron = 1;
	     }
	     else {
	       charged_cands[jdx].is_charged_hadron = 1;
	       charged_cands[jdx].is_electron = 0;
	     }

  	     if(central_clusters[index].et() > charged_cands[jdx].et()){
	        central_clusters[index].et() = central_clusters[index].et() - charged_cands[jdx].et();
	     }
	     else{
	        central_clusters[index].et() = 0;
	     } 
         }


        for(int idx = 0; idx < N_CLUSTERS; idx++){
	  if(central_clusters[idx].EoH > input_EoH_cut_){
	    central_clusters[idx].is_photon = 1;
	    central_clusters[idx].is_neutral_hadron = 0;
	  }
	  else{
	    central_clusters[idx].is_photon = 0;
	    central_clusters[idx].is_neutral_hadron = 1;

	  }
         }
	
}



/*
 * Find the index given the crystal location input
 * First divide by 5 for 5 crystals per tower
 * Then match the detector geometry as seen above
 */

uint32_t find_the_index_crys( ap_uint<7> eta, ap_uint<1> eta_side, ap_uint<8> phi){
  // First go from cyrstal to tower
                eta = eta/5;
                phi = phi/5;
		//temporary caution while code is integrated 
		//the max number of towers is 20... put this in a define
		if(eta > 20)
		  return 0;

		if(eta == 0 && phi == 0)
		  return N_CLUSTERS-1;

		if(eta == 0 && phi == 1)
		  return N_CLUSTERS-2;

		if(eta_side < 1)
			index = (20-eta)*72+phi-2;
		else
			index = (20+eta)*72+phi;

		return index;
}

void tau_three_prong_alg(track_t central_tracks[N_TRACKS], track_t three_prong_tau_cand[3], algo_config_t algo_config){

	for (int idx = 0; idx < N_TRACKS; idx++)	//note, tracks are already sorted by PT
	{
#pragma HLS UNROLL

		int n_found_tracks = 0;
		track_t seedtrack = central_tracks[idx];

		if(seedtrack.et < algo_config.three_prong_seed) // 0x46 corresponds to a 7 GeV seed. Follow ales's config example in UCTSummary
			continue;
		n_found_tracks++;
		three_prong_tau_cand[0] = seedtrack;

		for (int jdx = 0; jdx < N_TRACKS; jdx++)
		{
#pragma HLS UNROLL
			track_t temptrack = central_tracks[jdx];
			if(idx==jdx)
				continue;

			if(Delta_R(seedtrack.eta, seedtrack.phi, temptrack.eta, temptrack.phi, algo_config.three_prong_delta_r)){// 0x6 corresponds to deltaR of 0.12
				if(n_found_tracks<3){
					n_found_tracks++;
					three_prong_tau_cand[2] = temptrack;
					break;
					}

				if(n_found_tracks<2){
					n_found_tracks++;
					three_prong_tau_cand[1] = temptrack;
					}


			}
		}
	}

}

ap_uint<1> Delta_R(ap_uint<8> eta_1, ap_uint<8> phi_1, ap_uint<8> eta_2, ap_uint<8> phi_2, ap_uint<8> maximum_delta_R){
	ap_uint<1> output;
	ap_uint<8> delta_eta = 0;
	ap_uint<8> delta_phi = 0;
	if(eta_2>eta_1)
		delta_eta = eta_2 - eta_1;
	else
		delta_eta = eta_1 - eta_2;

	if(phi_2>phi_1)
		delta_phi = phi_2 - phi_1;
	else
		delta_phi = phi_1 - phi_2;

	if(delta_phi + delta_eta < maximum_delta_R)
		output = 1;
	else
		output = 0;

	return output;
}
