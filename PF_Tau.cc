#include "file_read_in.hpp"

void file_read_in(
		track_t central_tracks[N_TRACKS],  // Number of tracks
		cluster_t central_clusters[N_CLUSTERS],  // Number of Clusters
		algo_config_t algo_config,           // algorithm configuration
		algo_outputs_t & algo_outputs        // algorithm outputs
		)
{
#pragma HLS ARRAY_PARTITION variable=central_tracks complete dim=1

	ap_uint<14> et_total_tmp = 0;

#pragma HLS PIPELINE II=6

	pf_charged_t charged_cands[N_TRACKS];

	////////////////////   PF Particle Matching   ////////////////////
	//changes central_clusters to neutral_clusters
	pf_match_alg( central_clusters, central_tracks, charged_cands, algo_config);

	////////////////////Three Prong Tau Algorithm////////////////////
	pftau_t tau_cands[12];

	// Tau_alg Takes in charged cands, central (neutral) clusters, the algorithm configuration 
	// Returns the tau cands
	tau_alg( charged_cands, central_clusters, algo_config, tau_cands);

	////////////////////    Algorithm Outputs    ////////////////////
	algo_outputs.taus = tau_cands;


}

/*
 * track index is from -20 to 20 in eta (excluding 0)
 * then from 0 to 71 in phi
 * -20	                        -1 | 1                  +20
 * |0   72  72*2      (20-eta)+phi |                       |
 * |1	73  72*2+1                 |                       |
 * |2	74  72*3+1                 |                       |
 */
void pf_match_alg(cluster_t central_clusters[N_CLUSTERS],
					track_t central_tracks[N_TRACKS] ,
					pf_charged_t pf_charged[N_TRACKS],
					algo_config_t algo_config){
  
	for(int jdx = 0; jdx < N_TRACKS; jdx++)	//note, tracks are already sorted by PT
	  {
	    //#pragma HLS UNROLL
	    pf_charged[jdx].et       = central_tracks[jdx].et;
	    pf_charged[jdx].eta      = central_tracks[jdx].eta;
	    pf_charged[jdx].phi      = central_tracks[jdx].phi;
	    pf_charged[jdx].eta_side = central_tracks[jdx].eta_side;
	    ap_uint<12> index = 0;
	    
	    index = find_the_index_crys( pf_charged[jdx].eta, pf_charged[jdx].eta_side, pf_charged[jdx].phi);
	    
	    // Take the Cluster and use H/E and E/H to determine if Hadron or electron/pi0
	    if( central_clusters[index].EoH > algo_config.input_EoH_cut ){
   	       pf_charged[jdx].is_charged_hadron = 0;
	       pf_charged[jdx].is_electron = 1;
	     }
	     else {
	       pf_charged[jdx].is_charged_hadron = 1;
	       pf_charged[jdx].is_electron = 0;
	     }

  	     if(central_clusters[index].et > pf_charged[jdx].et){
	        central_clusters[index].et = central_clusters[index].et - pf_charged[jdx].et;
	     }
	     else{
	        central_clusters[index].et() = 0;
	     } 
         }


        for(int idx = 0; idx < N_CLUSTERS; idx++){
	  if(central_clusters[idx].EoH > algo_config.input_EoH_cut){
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

ap_uint<12> find_the_index_crys( ap_uint<7> eta, ap_uint<1> eta_side, ap_uint<8> phi){
				// First go from cyrstal to tower
                eta = eta/5;
                phi = phi/5;
                //temporary caution while code is integrated
                //the max number of towers is 20... put this in a define
                ap_uint<12> index;
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


// Offset by -2 in eta and -2 in phi, special geometry to grab the grid
ap_uint<12> find_the_index_crys_offset( ap_uint<7> eta, ap_uint<1> eta_side, ap_uint<8> phi, ap_uint<3> eta_offset, ap_uint<3> phi_offset){
				// First go from cyrstal to tower
  ///IMPLEMENT ME this needs to take into account the edges, i.e. return 0 if on the edge
  if(eta > 2 ){
    eta = (eta + eta_offset - 2)/5;
  }
  else{
    eta = 0;
  }

  if(phi> 2){
    phi = (phi + phi_offset - 2)/5;
  }
  else{
    phi = 0;
  }
                //the max number of towers is 20... put this in a define
                ap_uint<12> index;
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
//add clusters to void
void tau_alg(pf_charged_t pf_charged[N_TRACKS], clusters_t neutral_clusters[N_CLUSTERS], pftau_t tau_cands[12], algo_config_t algo_config){

        ap_uint<4> n_taus = 0;

  	pf_charged_t electron_grid[12][5][5];

	for (unsigned int idx = 0; idx < N_TRACKS; idx++)	//note, tracks are already sorted by PT
	{
#pragma HLS UNROLL
	  
	  ap_uint<3> n_prongs_found = 0;
	  pf_charged_t seed_hadron = pf_charged[idx];
	  pf_charged_t prong_cands[3];
	  pf_charged_t electron_grid_temp[5][5];
	  uint32_t iso_sum_charged_hadron = 0;

	  if(seed_hadron.et < algo_config.one_prong_seed)
	    continue;
	  
	  prong_cands[0] = seed_hadron;
	  
	  // Now that we have 1 prong Loop through the remaining pf_charged objects to find more prongs
	  for (int jdx = 0; jdx < N_TRACKS; jdx++)
	    {
#pragma HLS UNROLL

	      if(idx <= jdx)
		continue;

	      ap_uint<8> seed_cand_dr = delta_r_pf_charged(pf_charged[jdx], seed_hadron);

	      //Build the remaining prongs, but only if the next track is a charged hadron
	      if(pf_charged[jdx].is_charged_hadron > 0){
		pf_charged_t pf_charged_hadron_signal_cand = pf_charged[jdx];
		n_prongs_found++;
		
		if(seed_cand_dr < algo_config.three_prong_delta_r ){
		  if(n_prongs_found<2){
		    prong_cands[1] = pf_charged_hadron_signal_cand;
		    n_prongs_found++;
		    continue;
		  }

		  if(n_prongs_found<3){
		    prong_cands[2] = pf_charged_hadron_signal_cand;
		    n_prongs_found++;
		  }
		} // check to make sure that the delta_r is less than 5
		else if(seed_cand_dr < algo_config.isolation_delta_r ){

		  // Sum Charged Hadron Isolation
		  iso_sum_charged_hadron += pf_charged_hadron_signal_cand.et;

		}// Less than isolation_delta_r

	      }// pf_charged[jdx] is NOT a charged_hadron
	      else if(pf_charged[jdx].is_electron > 0){
		if(seed_cand_dr < algo_config.isolation_delta_r){
		  pf_charged_t electron_cand;
		  ap_uint<8> index_eta = ieta_diff(electron_cand, seed_hadron);
		  ap_uint<8> index_phi = iphi_diff(electron_cand, seed_hadron);
		  //keep from breaking things, This check is to be removed when code is validated
		  if(index_eta > 4 || index_phi > 4){
		    index_eta = 0; 
		    index_phi = 0;
		  }
		  electron_grid_temp[index_eta][index_phi] = electron_cand;
		}
	      } 
	    }// Finished looking through all the charged pf candidates

	  //Create the 1 prong taus
	  if(n_prongs_found < 3 && n_taus < 12){
	    if(seed_hadron.et > algo_config.one_prong_seed){
	      pftau_t pf_tau_temp;
	      pf_tau_temp.et          = seed_hadron.et;
	      pf_tau_temp.eta         = seed_hadron.eta;
	      pf_tau_temp.phi         = seed_hadron.phi;
	      pf_tau_temp.iso_charged = iso_sum_charged_hadron + prong_cands[2]; // prong_cands[2] should be 0 if only 1 prong, check me
	      pf_tau_temp.tau_type    = 0;
	      electron_grid[n_taus]   = electron_grid_temp;
	      tau_cands[n_taus]       = pf_tau_temp;
	      n_taus++;
	    }
	  }

	  //Create the 3 prong taus
	  if(n_prongs_found == 3 && n_taus < 12)
	    if(seed_hadron.et > algo_config.three_prong_seed){
	      pftau_t pf_tau_temp;
	      pf_tau_temp.et          = prong_cands[0].et + prong_cands[1].et + prong_cands[2].et;
	      pf_tau_temp.eta         = weighted_avg_eta(prong_cands[0], prong_cands[1], prong_cands[2]);
	      pf_tau_temp.phi         = weighted_avg_phi(prong_cands[0], prong_cands[1], prong_cands[2]);
	      pf_tau_temp.iso_charged = iso_sum_charged_hadron;
	      pf_tau_temp.tau_type    = 10;
	      electron_grid[n_taus]   = electron_grid_temp;
	      tau_cands[n_taus]       = pf_tau_temp;
	      n_taus++;
	    }
	}
	//Process the strips in a separate module

	for(ap_uint<4> i = 0; i < n_taus ; i++){
	  if(tau_cands[i].tau_type == 10 || tau_cands[i].et == 0)
	    continue;
	  strip_alg(tau_cands[i], electron_grid[i], neutral_clusters, algo_config);
	}
	//FINISH ISOLATION CALCULATION


}

void strip_alg(pftau_t tau_cand, pf_charged_t electron_grid[5][5],  cluster_t neutral_clusters[N_CLUSTERS], algo_config_t algo_config){

    ap_uint<7> tau_eta      = tau_cand.eta;
    ap_uint<1> tau_eta_side = tau_cand.eta_side;
    ap_uint<8> tau_phi      = tau_cand.phi;

    // Create Grid of 5x5 neutral Clusters
    ap_uint<12> index_mm  = find_the_index_crys_offset(tau_eta, tau_eta_side, tau_phi,  0, 0);
    ap_uint<12> index_m   = find_the_index_crys_offset(tau_eta, tau_eta_side, tau_phi,  1, 0);
    ap_uint<12> index_cen = find_the_index_crys_offset(tau_eta, tau_eta_side, tau_phi,  2, 0);
    ap_uint<12> index_p   = find_the_index_crys_offset(tau_eta, tau_eta_side, tau_phi,  3, 0);
    ap_uint<12> index_pp  = find_the_index_crys_offset(tau_eta, tau_eta_side, tau_phi,  4, 0);

    strip_t temp_strip[5];    
    strip_t final_strip;

    ap_uint<12> index_m   = find_the_index_crys_offset(tau_eta, tau_eta_side, tau_phi,  1, 0);

    merge_strip_algo( clusters[index - 2], electron_grid[0][0], clusters[index - 1], electron_grid[0][1], temp_strip[0]);
    merge_strip_algo( clusters[index - 1], electron_grid[0][1], clusters[index - 0], electron_grid[0][2], temp_strip[0]);
    merge_strip_algo( clusters[index + 0], electron_grid[0][2], clusters[index + 1], electron_grid[0][3], temp_strip[0]);
    merge_strip_algo( clusters[index + 1], electron_grid[0][3], clusters[index + 2], electron_grid[0][4], temp_strip[0]);

    // Check if each tau_phi slice is strip like
    // This module looks at each strip starting from the top cluster and going to the bottom
    //  -> Geometry outlined above for index
    // grid for strip creation is below
    // strips can be at most two towers in phi, potentially 2 towers in eta if crystals pass neighbors requirements.
    //   -----------
    //   |0|0|0|0|0|
    //   |1|1|1|1|1|
    //   |2|2|2|2|2|
    //   |3|3|3|3|3|
    //   |4|4|4|4|4|
    //   -----------

    merge_strip_algo(neutral_clusters[index_mm  + tau_phi - 2], electron_grid[0][0], neutral_clusters[index_mm  + tau_phi - 1], electron_grid[0][1], temp_strip[0], algo_config);
    merge_strip_algo(neutral_clusters[index_mm  + tau_phi - 1], electron_grid[0][1], neutral_clusters[index_mm  + tau_phi - 0], electron_grid[0][2], temp_strip[0], algo_config);
    merge_strip_algo(neutral_clusters[index_mm  + tau_phi + 0], electron_grid[0][2], neutral_clusters[index_mm  + tau_phi + 1], electron_grid[0][3], temp_strip[0], algo_config);
    merge_strip_algo(neutral_clusters[index_mm  + tau_phi + 1], electron_grid[0][3], neutral_clusters[index_mm  + tau_phi + 2], electron_grid[0][4], temp_strip[0], algo_config);

    merge_strip_algo(neutral_clusters[index_m   + tau_phi - 2], electron_grid[1][0], neutral_clusters[index_m   + tau_phi - 1], electron_grid[1][1], temp_strip[1], algo_config);
    merge_strip_algo(neutral_clusters[index_m   + tau_phi - 1], electron_grid[1][1], neutral_clusters[index_m   + tau_phi - 0], electron_grid[1][2], temp_strip[1], algo_config);
    merge_strip_algo(neutral_clusters[index_m   + tau_phi + 0], electron_grid[1][2], neutral_clusters[index_m   + tau_phi + 1], electron_grid[1][3], temp_strip[1], algo_config);
    merge_strip_algo(neutral_clusters[index_m   + tau_phi + 1], electron_grid[1][3], neutral_clusters[index_m   + tau_phi + 2], electron_grid[1][4], temp_strip[1], algo_config);

    merge_strip_algo(neutral_clusters[index_cen + tau_phi - 2], electron_grid[2][0], neutral_clusters[index_cen + tau_phi - 1], electron_grid[2][1], temp_strip[2], algo_config);
    merge_strip_algo(neutral_clusters[index_cen + tau_phi - 1], electron_grid[2][1], neutral_clusters[index_cen + tau_phi - 0], electron_grid[2][2], temp_strip[2], algo_config);
    merge_strip_algo(neutral_clusters[index_cen + tau_phi + 0], electron_grid[2][2], neutral_clusters[index_cen + tau_phi + 1], electron_grid[2][3], temp_strip[2], algo_config);
    merge_strip_algo(neutral_clusters[index_cen + tau_phi + 1], electron_grid[2][3], neutral_clusters[index_cen + tau_phi + 2], electron_grid[2][4], temp_strip[2], algo_config);

    merge_strip_algo(neutral_clusters[index_p   + tau_phi - 2], electron_grid[3][0], neutral_clusters[index_p   + tau_phi - 1], electron_grid[3][1], temp_strip[3], algo_config);
    merge_strip_algo(neutral_clusters[index_p   + tau_phi - 1], electron_grid[3][1], neutral_clusters[index_p   + tau_phi - 0], electron_grid[3][2], temp_strip[3], algo_config);
    merge_strip_algo(neutral_clusters[index_p   + tau_phi + 0], electron_grid[3][2], neutral_clusters[index_p   + tau_phi + 1], electron_grid[3][3], temp_strip[3], algo_config);
    merge_strip_algo(neutral_clusters[index_p   + tau_phi + 1], electron_grid[3][3], neutral_clusters[index_p   + tau_phi + 2], electron_grid[3][4], temp_strip[3], algo_config);

    merge_strip_algo(neutral_clusters[index_pp  + tau_phi - 2], electron_grid[4][0], neutral_clusters[index_pp  + tau_phi - 1], electron_grid[4][1], temp_strip[4], algo_config);
    merge_strip_algo(neutral_clusters[index_pp  + tau_phi - 1], electron_grid[4][1], neutral_clusters[index_pp  + tau_phi - 0], electron_grid[4][2], temp_strip[4], algo_config);
    merge_strip_algo(neutral_clusters[index_pp  + tau_phi + 0], electron_grid[4][2], neutral_clusters[index_pp  + tau_phi + 1], electron_grid[4][3], temp_strip[4], algo_config);
    merge_strip_algo(neutral_clusters[index_pp  + tau_phi + 1], electron_grid[4][3], neutral_clusters[index_pp  + tau_phi + 2], electron_grid[4][4], temp_strip[4], algo_config);

    final_strip = temp_strip[0];
    for(ap_uint<3> j = 1; j < 5; j++){
      //first check if strip j is greater than final strip
      if(temp_strip[j] > final_strip){
	final_strip = temp_strip[j];
      }

      //merge strips if two are close by
      if(delta_r_strip( temp_strip[j], temp_strip[j-1]) < algo_config.max_neighbor_strip_dist){
	if(temp_strip[j].et + temp_strip[j-1].et > final_strip.et()){

	  ap_uint<12> et  = temp_strip[j].et + temp_strip[j-1].et;
	  final_strip.et  = et ;
	  final_strip.eta = weighted_avg_eta(temp_strip[j], temp_strip[j-1]);
	  final_strip.phi = weighted_avg_phi(temp_strip[j], temp_strip[j-1]);

	}
      }
    }

  //find 1 prong pi0's by combining the 1 prong and the pi0's
    if(tau_cand.tau_type == 0 && final_strip.et > algo_config.min_strip){
      //take care of eta_side
      ap_uint<12> temp_tau_et;
      ap_uint<8>  temp_tau_eta;
      ap_uint<8>  temp_tau_phi;
      temp_tau_et  = tau_cand.et + final_strip.et;
      temp_tau_eta = weighted_avg_eta(tau_cand,final_strip);
      temp_tau_phi = weighted_avg_phi(tau_cand,final_strip);

      tau_cand.et       = temp_tau_et;
      tau_cand.eta      = temp_tau_eta;
      tau_cand.phi      = temp_tau_phi;
      tau_cand.tau_type = 1;
    }
}


void tau_three_prong_alg(track_t central_tracks[N_TRACKS], track_t three_prong_tau_cand[3], algo_config_t algo_config){

	for (int idx = 0; idx < N_TRACKS; idx++)	//note, tracks are already sorted by PT
	{
#pragma HLS UNROLL

		int n_found_tracks = 0;
		track_t seedtrack = central_tracks[idx];

		if(seedtrack.et < algo_config.three_prong_seed) // 0x46 corresponds to a 7 GeV seed. Follow ales's config example in UCTSummary
		  continue;

		three_prong_tau_cand[0] = seedtrack;
		n_found_tracks++;

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

void merge_strip_algo(cluster_t cluster_1, pf_charged_t electron_1, cluster_t cluster_2, pf_charged_t electron_2, strip_t strip, algo_config_t algo_config){

  cluster_t temp_cluster_1.et = cluster_1.et + electron_1.et;
  temp_cluster_1.phi = weighted_avg_phi(cluster_1 + electron_1);
  temp_cluster_1.phi = weighted_avg_phi(cluster_1 + electron_1);

  cluster_t temp_cluster_2.et = cluster_2.et + electron_2.et;
  temp_cluster_2.phi = weighted_avg_phi(cluster_2 + electron_2);
  temp_cluster_2.phi = weighted_avg_phi(cluster_2 + electron_2);

  // Requirement that it is a photon is checking if E/H is high for the cluster
  if(temp_cluster_1.is_photon > 0 && temp_cluster_2.is_photon > 0
     && delta_r_cluster(temp_cluster_1, temp_cluster_2) < algo_config.eg_strip_merge
     && strip.et < temp_cluster_1.et + temp_cluster_2.et){
    
    strip.et  = temp_cluster_1.et + temp_cluster_2.et;
    strip.eta = weighted_avg_eta(temp_cluster_1, temp_cluster_2);
    strip.phi = weighted_avg_phi(temp_cluster_1, temp_cluster_2);
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

ap_uint<10> delta_r_strip(strip_t strip1, strip_t strip2){
  ap_uint<8> delta_eta = 0;
  ap_uint<8> delta_phi = 0;
  if(strip1.eta > strip2.eta)
    delta_eta = strip1.eta = strip2.eta;

  if(strip1.eta > strip2.eta)
    delta_eta = strip1.eta = strip2.eta;
  else
    delta_eta = strip2.eta = strip1.eta;

  if(strip1.phi > strip2.phi)
    delta_phi = strip1.phi = strip2.phi;
  else
    delta_phi = strip2.phi = strip1.phi;

  return(delta_eta + delta_phi);
}

ap_uint<10> delta_r_cluster(cluster_t cluster1, cluster_t cluster2){
  ap_uint<8> delta_eta = 0;
  ap_uint<8> delta_phi = 0;
  if(cluster1.eta > cluster2.eta)
    delta_eta = cluster1.eta = cluster2.eta;

  if(cluster1.eta > cluster2.eta)
    delta_eta = cluster1.eta = cluster2.eta;
  else
    delta_eta = cluster2.eta = cluster1.eta;

  if(cluster1.phi > cluster2.phi)
    delta_phi = cluster1.phi = cluster2.phi;
  else
    delta_phi = cluster2.phi = cluster1.phi;

  return(delta_eta + delta_phi);
}


//fix me - check if enough bits are allocated during operations
ap_uint<7> weighted_avg_eta(strip_t strip1, strip_t strip2){
  ap_uint<12> eta1 = strip1.eta*strip1.et;
  ap_uint<12> eta2 = strip2.eta*strip2.et;
  ap_uint<12> sum_et = strip1.et + strip2.et;
  ap_uint<7> output_eta = (eta1 + eta2)/sum_et;

  return output_eta;
}

ap_uint<7> weighted_avg_phi(strip_t strip1, strip_t strip2){
  ap_uint<12> phi1 = strip1.phi*strip1.et;
  ap_uint<12> phi2 = strip2.phi*strip2.et;
  ap_uint<12> sum_et = strip1.et + strip2.et;
  ap_uint<7> output_phi = (phi1 + phi2)/sum_et;

  return output_phi;
}

ap_uint<7> weighted_avg_eta(cluster_t cluster1, cluster_t cluster2){
  ap_uint<12> eta1 = cluster1.eta*cluster1.et;
  ap_uint<12> eta2 = cluster2.eta*cluster2.et;
  ap_uint<12> sum_et = cluster1.et + cluster2.et;
  ap_uint<7> output_eta = (eta1 + eta2)/sum_et;

  return output_eta;
}

ap_uint<8> weighted_avg_phi(cluster_t cluster1, cluster_t cluster2){
  ap_uint<12> phi1 = cluster1.phi*cluster1.et;
  ap_uint<12> phi2 = cluster2.phi*cluster2.et;
  ap_uint<12> sum_et = cluster1.et + cluster2.et;
  ap_uint<7> output_phi = (phi1 + phi2)/sum_et;

  return output_phi;
}

ap_uint<7> weighted_avg_eta(pf_charged_t strip1, strip_t strip2){
  ap_uint<12> eta1 = strip1.eta*strip1.et;
  ap_uint<12> eta2 = strip2.eta*strip2.et;
  ap_uint<12> sum_et = strip1.et + strip2.et;
  ap_uint<7> output_eta = (eta1 + eta2)/sum_et;

  return output_eta;
}

ap_uint<8> weighted_avg_phi(pf_charged_t strip1, strip_t strip2){
  ap_uint<12> phi1 = strip1.phi*strip1.et;
  ap_uint<12> phi2 = strip2.phi*strip2.et;
  ap_uint<12> sum_et = strip1.et + strip2.et;
  ap_uint<7> output_phi = (phi1 + phi2)/sum_et;

  return output_phi;
}

ap_uint<7> weighted_avg_eta(cluster_t cluster1, strip_t strip2){
  ap_uint<12> eta1 = cluster1.eta*cluster1.et;
  ap_uint<12> eta2 = strip2.eta*strip2.et;
  ap_uint<12> sum_et = cluster1.et + strip2.et;
  ap_uint<7> output_eta = (eta1 + eta2)/sum_et;

  return output_eta;
}

ap_uint<8> weighted_avg_phi(cluster_t cluster1, strip_t strip2){
  ap_uint<12> phi1 = cluster1.phi*cluster1.et;
  ap_uint<12> phi2 = strip2.phi*strip2.et;
  ap_uint<12> sum_et = cluster1.et + strip2.et;
  ap_uint<7> output_phi = (phi1 + phi2)/sum_et;

  return output_phi;
}

ap_uint<7> weighted_avg_eta(pf_charged_t pf1, pf_charged_t pf2, pf_charged_t pf3){
  ap_uint<12> eta1 = pf1.eta*pf1.et;
  ap_uint<12> eta2 = pf2.eta*pf2.et;
  ap_uint<12> eta3 = pf3.eta*pf3.et;
  ap_uint<12> sum_et = pf1.et + pf2.et + pf3.et;
  ap_uint<7> output_eta = (eta1 + eta2 + eta3)/sum_et;

  return output_eta;
}

ap_uint<7> weighted_avg_phi(pf_charged_t pf1, pf_charged_t pf2, pf_charged_t pf3){
  ap_uint<12> phi1 = pf1.phi*pf1.et;
  ap_uint<12> phi2 = pf2.phi*pf2.et;
  ap_uint<12> phi3 = pf3.phi*pf3.et;
  ap_uint<12> sum_et = pf1.et + pf2.et + pf3.et;
  ap_uint<7> output_phi = (phi1 + phi2 + phi3)/sum_et;

  return output_phi;
}

ap_uint<8> delta_r_pf_charged(pf_charged_t pf_1, pf_charged_t pf_2){
	ap_uint<8> delta_eta = 0;
	ap_uint<8> delta_phi = 0;
	ap_uint<8> eta_1 = pf_1.eta;
	ap_uint<8> eta_2 = pf_2.eta;
	ap_uint<8> phi_1 = pf_1.phi;
	ap_uint<8> phi_2 = pf_2.phi;
	
	if(eta_2>eta_1)
		delta_eta = eta_2 - eta_1;
	else
		delta_eta = eta_1 - eta_2;

	if(phi_2>phi_1)
		delta_phi = phi_2 - phi_1;
	else
		delta_phi = phi_1 - phi_2;

	return (delta_eta + delta_phi);
}


ap_uint<8> ieta_diff(pf_charged_t cand_1, pf_charged_t cand_2){
  //fix me for eta side
  if(cand_1.eta > cand_2.eta){
    return(cand_1.eta - cand_2.eta)
  }
  else
    return(cand_2.eta - cand_1.eta)

}

ap_uint<8> iphi_diff(pf_charged_t cand_1, pf_charged_t cand_2){
  //fix me for phi side
  if(cand_1.phi > cand_2.phi){
    return(cand_1.phi - cand_2.phi)
  }
  else
    return(cand_2.phi - cand_1.phi)

}
