 #include "PF_Tau.hpp"

void file_read_in(
		  track_t central_tracks[N_TRACKS],  // Number of tracks
		 cluster_t central_clusters[N_CLUSTERS],  // Number of Clusters
		  algo_config_t algo_config,           // algorithm configuration
		  algo_outputs_t & algo_outputs        // algorithm outputs
		  )
{
#pragma HLS ARRAY_PARTITION variable=central_tracks complete dim=1
  
  ap_uint<14> et_total_tmp = 0;
  
#pragma HLS PIPELINE II=8
  
  pf_charged_t charged_cands[N_TRACKS];
  
  ////////////////////   PF Particle Matching   ////////////////////
  //changes central_clusters to neutral_clusters
  pf_match_alg( central_clusters, central_tracks, charged_cands, algo_config);
  
  ////////////////////Three Prong Tau Algorithm////////////////////
  pftau_t tau_cands[N_TAUS];
  
  // Tau_alg Takes in charged cands, central (neutral) clusters, the algorithm configuration 
  // Returns the tau cands
  tau_alg( charged_cands, central_clusters, algo_config, tau_cands);
  
  ////////////////////    Algorithm Outputs    ////////////////////
  //algo_outputs.taus = tau_cands;
}

 /*
  * track index is from -20 to 20 in eta (excluding 0)
  * then from 0 to 71 in phi
  * -20	                        -1 | 1                  +20
  * |0   72  72*2      (20-eta)+phi |                       |
  * |1	73  72*2+1                 |                       |
  * |2	74  72*3+1                 |                       |
  *
  * Take as input the clusters, tracks and output pf_charged objects
  * This function should be replaced by cluster track linker!
  */
 void pf_match_alg(cluster_t central_clusters[N_CLUSTERS],
		   track_t central_tracks[N_TRACKS] ,
		   pf_charged_t pf_charged[N_TRACKS],
		   algo_config_t algo_config){
 #pragma HLS PIPELINE II=6 
 #pragma HLS ARRAY_PARTITION variable=central_clusters complete dim=0
 #pragma HLS ARRAY_PARTITION variable=central_tracks complete dim=0
 #pragma HLS ARRAY_PARTITION variable=pf_charged complete dim=0

   //Initialize pf_charged objects (represents all )
   for(int jdx = 0; jdx < N_TRACKS; jdx++)	//note, tracks are already sorted by PT
     {
#pragma HLS UNROLL
       pf_charged[jdx].et       = central_tracks[jdx].et;
       pf_charged[jdx].eta      = central_tracks[jdx].eta;
       pf_charged[jdx].phi      = central_tracks[jdx].phi;
       pf_charged[jdx].eta_side = central_tracks[jdx].eta_side;
     }
   
   for(ap_uint<9> index = 0; index < N_CLUSTERS; index++){
#pragma HLS UNROLL
     for(int jdx = 0; jdx < N_TRACKS; jdx++)	//note, tracks are already sorted by PT
       {
#pragma HLS UNROLL
	 check_pf_cand(pf_charged[jdx], central_clusters[index], algo_config);
	 
       }
   }
   
   for(int idx = 0; idx < N_CLUSTERS; idx++){
#pragma HLS UNROLL
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

 // This should also be replaced by the cluster track linker 
 void check_pf_cand(pf_charged_t &pf_charged, cluster_t &central_cluster, algo_config_t algo_config){
   if( delta_r_c_p(central_cluster, pf_charged) < 3){

     // Take the Cluster and use H/E and E/H to determine if Hadron or electron/pi0
     if( central_cluster.EoH > algo_config.input_EoH_cut ){
       pf_charged.is_charged_hadron = 0;
       pf_charged.is_electron = 1;
     }
     else {
       pf_charged.is_charged_hadron = 1;
       pf_charged.is_electron = 0;
     }
   }

   if(central_cluster.et > pf_charged.et){
     central_cluster.et = central_cluster.et - pf_charged.et;
   }
   else{
     central_cluster.et = 0;
   }

 }


/*
 *The Tau algorithm takes as input track-linked clusters, Neutral 
 * clusters (i.e. clusters without tracks linked) and output N_TAUS
 *
 */

void tau_alg(pf_charged_t pf_charged[N_TRACKS], cluster_t neutral_clusters[N_CLUSTERS], algo_config_t algo_config, pftau_t tau_cands_O[N_TAUS]){
  //Unroll variables to speed up algorithm
#pragma HLS ARRAY_PARTITION variable=neutral_clusters complete dim=0
#pragma HLS ARRAY_PARTITION variable=pf_charged complete dim=0
#pragma HLS ARRAY_PARTITION variable=tau_cands_O complete dim=0
#pragma HLS PIPELINE II=8
  //initialize output tau_cands
  pftau_t tau_cands[N_TAUS];
  
#pragma HLS ARRAY_PARTITION variable=tau_cands complete dim=0
  //found number of taus
   ap_uint<4> n_taus = 0;
   //electron_grid for isolation sum
   pf_charged_t electron_grid[N_TAUS][5][5];

#pragma HLS ARRAY_PARTITION variable=electron_grid complete dim=0

   //To build tau candidates we consider each track as 
   // a potential tau hadron seed. 
   for (unsigned int idx = 0; idx < N_TRACKS; idx++)	//note, tracks are already sorted by PT
     {
#pragma HLS UNROLL
       
       ap_uint<3> n_prongs_found = 0;
       pf_charged_t seed_hadron = pf_charged[idx];
       pf_charged_t prong_cands[3];
#pragma HLS ARRAY_PARTITION variable=prong_cands complete dim=0
       uint32_t iso_sum_charged_hadron = 0;
       // Check to see if the tau candidate seed passes the minimum criteria
       //  In particular, we have a very large number of low pt tracks! Only
       // tracks which pass minimum pt and charged hadron ID requirements 
       // are chose. (Note: Since this is firmware this does not actually
       // affect the number of times the rest of the algorithm is run)
       if(seed_hadron.et > algo_config.one_prong_seed && seed_hadron.is_charged_hadron > 0){
	 
	 // Minimal requirements for building a tau have been met, now we begin construction
	 // Note: Any additional minimal requirements for building tau should go above.)\
	 prong_cands[0] = seed_hadron;
	 
	 // Now that we have 1 prong Loop through the remaining pf_charged objects to 
	 // find three prong cands. If the track does not match 3 prong tau criteria
	 // then it can become part of the Isolation cone.
	 for (int jdx = 0; jdx < N_TRACKS; jdx++)
	   {
#pragma HLS UNROLL

	     if(jdx < idx){
	       // Calculate the delta R between current charged hadron and seed charged hadron
	       ap_uint<8> seed_cand_dr = delta_r_pf_charged(pf_charged[jdx], seed_hadron);

	       //Build the remaining prongs, but only if the track is a charged hadron
	       if(pf_charged[jdx].is_charged_hadron > 0){
		 
		 // Check DeltaR and charged hadron criteria to see if this object can be a 3 prong tau
		 find_tau_prongs( n_prongs_found,  prong_cands,  pf_charged[jdx],  seed_hadron, seed_cand_dr, iso_sum_charged_hadron, algo_config);
		 
		   }// pf_charged[jdx] is NOT a charged_hadron
	       else if(pf_charged[jdx].is_electron > 0){
		 // Electrons for pi0 Strip finding, Not used in 3 prong tau algorithm!
		 if(seed_cand_dr < algo_config.isolation_delta_r){
		   build_electron_grid(pf_charged[jdx], seed_cand_dr, seed_hadron, electron_grid,  n_taus);
		 }
	       } 
	     }// Finished looking through all the charged pf candidates
	   }
	 //If there are less than 3 prongs, we have not found N_Taus yet, and the seed 
	 //prong meets minimum transverse energy requirements then Create the 1 prong taus
	 if(n_prongs_found < 3 && n_taus < N_TAUS){
	   if(seed_hadron.et > algo_config.one_prong_seed){
	     
	     tau_cands[n_taus].et          = seed_hadron.et;
	     tau_cands[n_taus].eta         = seed_hadron.eta;
	     tau_cands[n_taus].phi         = seed_hadron.phi;
	     //TODO: prong_cands[2] should be 0 if only 1 prong, check me 
	     tau_cands[n_taus].iso_charged = iso_sum_charged_hadron + prong_cands[2].et; 
	     tau_cands[n_taus].tau_type    = 0;
	     n_taus++;
	   }
	     }// Sum the found prongs to Create the 3 prong taus
	     else if(n_prongs_found == 3 && n_taus < N_TAUS)
	       if(seed_hadron.et > algo_config.three_prong_seed){
		 tau_cands[n_taus].et          = prong_cands[0].et + prong_cands[1].et + prong_cands[2].et;
		 tau_cands[n_taus].eta         = weighted_avg_eta_p_p_p(prong_cands[0], prong_cands[1], prong_cands[2]);
		 tau_cands[n_taus].phi         = weighted_avg_phi_p_p_p(prong_cands[0], prong_cands[1], prong_cands[2]);
		 tau_cands[n_taus].iso_charged = iso_sum_charged_hadron;
		 tau_cands[n_taus].tau_type    = 10;
		 n_taus++;
	       }
	   }
	   //Process the strips in a separate module
	 }
   // This is the strip processing. Ignore for 3 prong tau finding!!
   for(ap_uint<4> i = 0; i < N_TAUS ; i++){
#pragma HLS UNROLL
     if(tau_cands[i].tau_type == 0){
       strip_alg(tau_cands[i], electron_grid[i], neutral_clusters, algo_config);}
     else{
       tau_cands[i].tau_type = 0;}
   }
   
   for(ap_uint<4> i = 0; i < N_TAUS ; i++){
#pragma HLS UNROLL
     tau_cands_O[i] = tau_cands[i];

   }
   //TODO: FINISH ISOLATION CALCULATION
    
}

// Simple function to check 
//1.) Check to see if candidate prong matches deltaR requirements to be a signal prong
//2.) Check to see if candidate prong can be put in place [2] or [2]
//3.) Check to see if the seed cand can be put in the isolation cone
//    If yes then increment isolation cone
 void find_tau_prongs(  ap_uint<3> &n_prongs_found, pf_charged_t prong_cands[3], pf_charged_t pf_charged_hadron_signal_cand, pf_charged_t seed_hadron, ap_uint<8> seed_cand_dr,ap_uint<12> iso_sum_charged_hadron, algo_config_t algo_config){
 #pragma HLS ARRAY_PARTITION variable=prong_cands complete dim=0
	   n_prongs_found++;
	   
	   //If charged hadron matches deltaR requirements then it can be one of the prongs for a 3 prong tau
	   if(seed_cand_dr < algo_config.three_prong_delta_r ){
	       //if this is the third prong found then place it at [2]
	     if(n_prongs_found==2){
	       prong_cands[2] = pf_charged_hadron_signal_cand;
	       n_prongs_found++;
	     }
	       //if this is the second prong found then place it at [2]
	     if(n_prongs_found==1){
	       prong_cands[1] = pf_charged_hadron_signal_cand;
	       n_prongs_found++;
	       
	     }
	   } // check to make sure that the delta_r is less than 5 -> Then the charged hadron goes in the isolation cone.
	   else if(seed_cand_dr < algo_config.isolation_delta_r ){
	     
	     // Sum Charged Hadron Isolation
	     iso_sum_charged_hadron += pf_charged_hadron_signal_cand.et;
	     
	   }// Less than isolation_delta_r
	   
 }
 // This needs the ieta and iphi from deltar function as an input
 void build_electron_grid(pf_charged_t electron_cand, ap_uint<8> seed_cand_dr, pf_charged_t seed_hadron, pf_charged_t electron_grid[N_TAUS][5][5], ap_uint<4> n_taus){
 #pragma HLS ARRAY_PARTITION variable=electron_grid complete dim=0

   ap_uint<7> index_eta = seed_cand_dr>>2;
   ap_uint<8> index_phi = seed_cand_dr>>2;
	 //keep from breaking things, This check is to be removed when code is validated
	 if(index_eta > 4 || index_phi > 4){
	   index_eta = 0;
	   index_phi = 0;
	 }
	 electron_grid[n_taus][index_eta][index_phi] = electron_cand;
 }

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
 void strip_alg(pftau_t &tau_cand, pf_charged_t electron_grid[5][5],  cluster_t neutral_clusters[N_CLUSTERS], algo_config_t algo_config){
 #pragma HLS ARRAY_PARTITION variable=electron_grid complete dim=0
 #pragma HLS ARRAY_PARTITION variable=neutral_clusters complete dim=0
 //#pragma HLS PIPELINE II=6

   ap_uint<7> tau_eta      = tau_cand.eta;
   ap_uint<1> tau_eta_side = tau_cand.eta_side;
   ap_uint<8> tau_phi      = tau_cand.phi;

     ap_uint<12> index[5];
 #pragma HLS ARRAY_PARTITION variable=index complete dim=0

     cluster_t neutral_cluster_grid[5][5];
 #pragma HLS ARRAY_PARTITION variable=neutral_cluster_grid complete dim=0

     for(ap_uint<3> i = 0; i<5; i++){
 #pragma HLS UNROLL
	   for(ap_uint<3> j = 0; j<5; j++){
 #pragma HLS UNROLL
	     ap_uint<7> grid_eta = tau_eta + i>>2;
	     ap_uint<7> grid_phi = tau_phi + j>>2;
		   neutral_cluster_grid[i][j] = find_matched_cluster(neutral_clusters, grid_eta, grid_phi);
	   }
     }

     strip_t temp_strip[5];    
 #pragma HLS ARRAY_PARTITION variable=temp_strip complete dim=0
     strip_t final_strip;

     cluster_t cluster;
     pf_charged_t charged;

     ap_uint<8> index1; ap_uint<8> index2;

     for(ap_uint<3> i = 0; i<5; i++){
 #pragma HLS UNROLL
       for(ap_uint<3> j = 0; j<4; j++){
 #pragma HLS UNROLL
	   ap_uint<3> jp = j+1;
	   merge_strip_algo(neutral_cluster_grid[i][j], electron_grid[i][j], neutral_cluster_grid[i][jp], electron_grid[i][jp], temp_strip[i], algo_config);
       }
     }

     final_strip = temp_strip[0];
     for(ap_uint<3> j = 1; j < 5; j++){
 #pragma HLS UNROLL
       //first check if strip j is greater than final strip
       if(temp_strip[j].et > final_strip.et){
	      final_strip = temp_strip[j];
       }

       //merge strips if two are close by
       if(delta_r_strip( temp_strip[j], temp_strip[j-1]) < algo_config.max_neighbor_strip_dist){
	     if(temp_strip[j].et + temp_strip[j-1].et > final_strip.et){

		 ap_uint<12> et  = temp_strip[j].et + temp_strip[j-1].et;
		 final_strip.et  = et ;
		 final_strip.eta = weighted_avg_eta_s_s(temp_strip[j], temp_strip[j-1]);
		 final_strip.phi = weighted_avg_phi_s_s(temp_strip[j], temp_strip[j-1]);

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
       temp_tau_eta = weighted_avg_eta_t_s(tau_cand,final_strip);
       temp_tau_phi = weighted_avg_phi_t_s(tau_cand,final_strip);

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

		 if(seedtrack.et < algo_config.three_prong_seed)
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

 void merge_strip_algo(cluster_t cluster_1, pf_charged_t electron_1, cluster_t cluster_2, pf_charged_t electron_2, strip_t &strip, algo_config_t algo_config){
   cluster_t temp_cluster_1;
   temp_cluster_1.et = cluster_1.et + electron_1.et;
   temp_cluster_1.phi = weighted_avg_phi_c_p(cluster_1, electron_1);
   temp_cluster_1.eta = weighted_avg_eta_c_p(cluster_1, electron_1);

   cluster_t temp_cluster_2;
   temp_cluster_2.et = cluster_2.et + electron_2.et;
   temp_cluster_2.phi = weighted_avg_phi_c_p(cluster_2, electron_2);
   temp_cluster_2.eta = weighted_avg_eta_c_p(cluster_2, electron_2);

   // Requirement that it is a photon is checking if E/H is high for the cluster
   if(temp_cluster_1.is_photon > 0 && temp_cluster_2.is_photon > 0
      && delta_r_cluster(temp_cluster_1, temp_cluster_2) < algo_config.eg_strip_merge
      && strip.et < temp_cluster_1.et + temp_cluster_2.et){

     strip.et  = temp_cluster_1.et + temp_cluster_2.et;
     strip.eta = weighted_avg_eta_c_c(temp_cluster_1, temp_cluster_2);
     strip.phi = weighted_avg_phi_c_c(temp_cluster_1, temp_cluster_2);
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
 ap_uint<7> weighted_avg_eta_t_s(pftau_t strip1, strip_t strip2){

   ap_uint<7> output_eta = (strip1.eta + strip2.eta)>>1;
   return output_eta;
 }

 ap_uint<8> weighted_avg_phi_t_s(pftau_t strip1, strip_t strip2){
   //ap_uint<12> phi1 = strip1.phi*strip1.et;
   //ap_uint<12> phi2 = strip2.phi*strip2.et;
   //ap_uint<12> sum_et = strip1.et + strip2.et;
   //ap_uint<7> output_phi = (phi1 + phi2)/sum_et;

   ap_uint<7> output_phi = (strip1.phi+strip2.phi)>>1;
   return output_phi;
 }

 //fix me - check if enough bits are allocated during operations
 ap_uint<7> weighted_avg_eta_s_s(strip_t strip1, strip_t strip2){
   //ap_uint<12> eta1 = strip1.eta*strip1.et;
   //ap_uint<12> eta2 = strip2.eta*strip2.et;
   //ap_uint<12> sum_et = strip1.et + strip2.et;
   //ap_uint<7> output_eta = (eta1 + eta2)/sum_et;
   ap_uint<7> output_eta = (strip1.eta + strip2.eta)>>1;
   return output_eta;
 }

 ap_uint<8> weighted_avg_phi_s_s(strip_t strip1, strip_t strip2){
   //ap_uint<12> phi1 = strip1.phi*strip1.et;
   //ap_uint<12> phi2 = strip2.phi*strip2.et;
   //ap_uint<12> sum_et = strip1.et + strip2.et;
   //ap_uint<7> output_phi = (phi1 + phi2)/sum_et;
   ap_uint<7> output_phi  = (strip1.phi+ strip2.phi)>>1;
   return output_phi;
 }

 ap_uint<7> weighted_avg_eta_c_p(cluster_t cluster1, pf_charged_t cluster2){
   //ap_uint<12> eta1 = cluster1.eta*cluster1.et;
   //ap_uint<12> eta2 = cluster2.eta*cluster2.et;
   //ap_uint<12> sum_et = cluster1.et + cluster2.et;
   //ap_uint<7> output_eta = (eta1 + eta2)/sum_et;

   ap_uint<7> output_eta = (cluster1.eta + cluster2.eta)>>1;
   return output_eta;
 }

 ap_uint<8> weighted_avg_phi_c_p(cluster_t cluster1, pf_charged_t cluster2){
   //ap_uint<12> phi1 = cluster1.phi*cluster1.et;
   //ap_uint<12> phi2 = cluster2.phi*cluster2.et;
   //ap_uint<12> sum_et = cluster1.et + cluster2.et;
   ap_uint<8> output_phi = (cluster1.phi + cluster2.phi)>>1;
   return output_phi;
 }

 ap_uint<7> weighted_avg_eta_c_c(cluster_t cluster1, cluster_t cluster2){
   //ap_uint<12> eta1 = cluster1.eta*cluster1.et;
   //ap_uint<12> eta2 = cluster2.eta*cluster2.et;
   // ap_uint<12> sum_et = cluster1.et + cluster2.et;
   // ap_uint<7> output_eta = (eta1 + eta2)/sum_et;

   //return output_eta;

   ap_uint<7> output_eta = (cluster1.eta + cluster2.eta)>>1;
   return output_eta;
 }

 ap_uint<8> weighted_avg_phi_c_c(cluster_t cluster1, cluster_t cluster2){
   //ap_uint<12> phi1 = cluster1.phi*cluster1.et;
   //ap_uint<12> phi2 = cluster2.phi*cluster2.et;
   //ap_uint<12> sum_et = cluster1.et + cluster2.et;
   //ap_uint<7> output_phi = (phi1 + phi2)/sum_et;

   //return output_phi;

   ap_uint<8> output_phi = (cluster1.phi + cluster2.phi)>>1;
   return output_phi;
 }

 ap_uint<7> weighted_avg_eta_p_s(pf_charged_t strip1, strip_t strip2){
   //ap_uint<12> eta1 = strip1.eta*strip1.et;
   //ap_uint<12> eta2 = strip2.eta*strip2.et;
   //ap_uint<12> sum_et = strip1.et + strip2.et;
   //ap_uint<7> output_eta = (eta1 + eta2)/sum_et;

   ap_uint<7> output_eta = (strip1.eta + strip2.eta)>>1;
   return output_eta;
 }

 ap_uint<8> weighted_avg_phi_p_s(pf_charged_t strip1, strip_t strip2){
   //ap_uint<12> phi1 = strip1.phi*strip1.et;
   //ap_uint<12> phi2 = strip2.phi*strip2.et;
   //ap_uint<12> sum_et = strip1.et + strip2.et;
   //ap_uint<7> output_phi = (phi1 + phi2)/sum_et;
   ap_uint<7> output_phi = (strip1.phi + strip2.phi)>>1;
   return output_phi;
 }

 ap_uint<7> weighted_avg_eta_c_s(cluster_t cluster1, strip_t strip2){
   //ap_uint<12> eta1 = cluster1.eta*cluster1.et;
   //ap_uint<12> eta2 = strip2.eta*strip2.et;
   //ap_uint<12> sum_et = cluster1.et + strip2.et;
   //ap_uint<7> output_eta = (eta1 + eta2)/sum_et;
   ap_uint<7> output_eta = (cluster1.eta + strip2.eta)>>1;
   return output_eta;
 }

 ap_uint<8> weighted_avg_phi_c_s(cluster_t cluster1, strip_t strip2){
   //ap_uint<12> phi1 = cluster1.phi*cluster1.et;
   //ap_uint<12> phi2 = strip2.phi*strip2.et;
   //ap_uint<12> sum_et = cluster1.et + strip2.et;
   //ap_uint<7> output_phi = (phi1 + phi2)/sum_et;
   ap_uint<7> output_phi = (cluster1.phi+strip2.phi)>>1;
   return output_phi;
 }

 ap_uint<7> weighted_avg_eta_p_p_p(pf_charged_t pf1, pf_charged_t pf2, pf_charged_t pf3){
   //ap_uint<12> eta1 = pf1.eta*pf1.et;
   //ap_uint<12> eta2 = pf2.eta*pf2.et;
   //ap_uint<12> eta3 = pf3.eta*pf3.et;
   //ap_uint<12> sum_et = pf1.et + pf2.et + pf3.et;
   //ap_uint<7> output_eta = (eta1 + eta2 + eta3)/sum_et;

   ap_uint<7> output_eta = (pf1.eta + ((pf2.eta + pf3.eta)>>1 ))>>1;
   return output_eta;
 }

 ap_uint<8> weighted_avg_phi_p_p_p(pf_charged_t pf1, pf_charged_t pf2, pf_charged_t pf3){
   //ap_uint<12> phi1 = pf1.phi*pf1.et;
   //ap_uint<12> phi2 = pf2.phi*pf2.et;
   //ap_uint<12> phi3 = pf3.phi*pf3.et;
   //ap_uint<12> sum_et = pf1.et + pf2.et + pf3.et;
   ap_uint<8> output_phi = (pf1.phi + ((pf2.phi + pf3.phi)>>1))>>1;

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

 ap_uint<8> delta_r_c_p(cluster_t pf_1, pf_charged_t pf_2){
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


 ap_uint<7> ieta_diff(pf_charged_t cand_1, pf_charged_t cand_2){
   //fix me for eta side
   if(cand_1.eta > cand_2.eta){
     return(cand_1.eta - cand_2.eta);
   }

    return(cand_2.eta - cand_1.eta);

 }

 ap_uint<8> iphi_diff(pf_charged_t cand_1, pf_charged_t cand_2){
   //fix me for phi side
   if(cand_1.phi > cand_2.phi){
     return(cand_1.phi - cand_2.phi);
   }

     return(cand_2.phi - cand_1.phi);

 }
 /*
  * Find the index given the crystal location input
  * First divide by 5 for 5 crystals per tower
  * Then match the detector geometry as seen above
  */

 ap_uint<12> find_the_index_crys( ap_uint<7> eta, ap_uint<1> eta_side, ap_uint<8> phi){
		 // First go from cyrstal to tower
		 //eta = eta/5;
		 //phi = phi/5;
	 return 5;
		 ap_uint<5> eta_calc = eta/5;
		 ap_uint<5> phi_calc = phi/5;
	  //temporary caution while code is integrated
	  //the max number of towers is 20... put this in a define
		 ap_uint<12> index;
		 /*
		if(eta > 20)
		    return 0;

		if(eta == 0 && phi == 0)
		    return N_CLUSTERS-1;

		if(eta == 0 && phi == 1)
		    return N_CLUSTERS-2;
 */
		if(eta_side < 1){
		    index = (20-eta_calc)*72+phi_calc-2;
		}
		else{
		    index = (20+eta_calc)*72+phi_calc;
		}

		return index;
 }


 cluster_t find_matched_cluster(cluster_t neutral_clusters[N_CLUSTERS], ap_uint<8> eta_1, ap_uint<8> phi_1){
 #pragma HLS ARRAY_PARTITION variable=neutral_clusters complete dim=0
	 cluster_t cluster_to_return;
	 cluster_to_return = neutral_clusters[0];
   for(ap_uint<9> i =0; i<N_CLUSTERS; i++){
 #pragma HLS UNROLL    
     if(Delta_R(neutral_clusters[i].eta, neutral_clusters[i].phi, eta_1, phi_1, 4))
	 cluster_to_return = neutral_clusters[i];
   }
   //fix me
   return cluster_to_return;

 }

 // Offset by -2 in eta and -2 in phi, special geometry to grab the grid
 // IMPLEMENTE ME We are now making the cluster grid 2^7 in phi (128)
 ap_uint<12> find_index_crys_offset( ap_uint<7> eta, ap_uint<1> eta_side, ap_uint<8> phi, ap_uint<3> eta_offset, ap_uint<3> phi_offset){
				 // First go from cyrstal to tower
   ///IMPLEMENT ME this needs to take into account the edges, i.e. return 0 if on the edge
   //#pragma HLS PIPELINE II=1
   if(eta > 2 ){
     eta = (eta + eta_offset - 2)>>2;
   }
   else{
     eta = 0;
   }

   if(phi> 2){
     phi = (phi + phi_offset - 2)>>2;
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
     index = (20-eta)<<7+phi-2;
   else
     index = (20+eta)<<7+phi;
   
   return index;
 }
