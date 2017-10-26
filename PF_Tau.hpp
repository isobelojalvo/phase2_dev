#ifndef pf_tau_H_
#define pf_tau_H_

#include <stdint.h>
#include <ap_int.h>

//define number of tracks
#define N_TRACKS (15)
#define N_CLUSTERS (100)
#define N_TAUS (4)

typedef struct
{
	ap_uint<11> et;
	ap_uint<1> eta_side;
	ap_uint<7> eta;
	ap_uint<8> phi;

} track_t;

typedef struct
{
	ap_uint<11> et;
	ap_uint<1> eta_side;
	ap_uint<7> eta;
	ap_uint<8> phi;
	ap_uint<8> EoH;
	ap_uint<8> HoE;
    ap_uint<1> is_photon;
    ap_uint<1> is_neutral_hadron;

} cluster_t;


typedef struct
{
	ap_uint<11> et;
	ap_uint<1> eta_side;
	ap_uint<7> eta;
	ap_uint<8> phi;
        ap_uint<1> is_photon;
        ap_uint<1> is_neutral_hadron;
} pf_neutral_t;

typedef struct
{
		ap_uint<11> et;
		ap_uint<1> eta_side;
		ap_uint<7> eta;
		ap_uint<8> phi;
        ap_uint<1> is_electron;
        ap_uint<1> is_charged_hadron;
        ap_uint<1> is_muon;
} pf_charged_t;

typedef struct
{
		ap_uint<11> et;
		ap_uint<1> eta_side;
		ap_uint<7> eta;
		ap_uint<8> phi;
        ap_uint<4> tau_type;
        ap_uint<12> iso_charged;

} pftau_t;

typedef struct
{
  	ap_uint<11> et;
	ap_uint<1> eta_side;
	ap_uint<7> eta;
	ap_uint<8> phi;

} strip_t;

typedef struct{
	ap_uint<11> three_prong_seed;
	ap_uint<11> three_prong_delta_r;
	ap_uint<11> isolation_delta_r;
	ap_uint<11> one_prong_seed;
	ap_uint<1>  dummy;
	ap_uint<5>  input_EoH_cut;
	ap_uint<5>  max_neighbor_strip_dist;
	ap_uint<11> min_strip;
	ap_uint<5>  eg_strip_merge;
} algo_config_t;

typedef struct{
	ap_uint<14> sum_tracks;
	ap_uint<14> three_prong_tau_et;
	pftau_t taus[N_TAUS];
} algo_outputs_t;

void file_read_in(track_t tracks[N_TRACKS],
		  cluster_t clusters[N_CLUSTERS],
		  algo_config_t algo_config,
		  algo_outputs_t & algo_outputs
);


ap_uint<1> Delta_R(ap_uint<8> eta_1,
			ap_uint<8> phi_1,
			ap_uint<8> eta_2,
			ap_uint<8> phi_2,
			ap_uint<8> maximum_delta_R);

ap_uint<12> find_the_index_crys( ap_uint<7> eta, ap_uint<1> eta_side, ap_uint<8> phi);

ap_uint<12> find_the_index_crys( ap_uint<7> eta, ap_uint<1> eta_side, ap_uint<8> phi);

ap_uint<8> delta_r_pf_charged(pf_charged_t pf_1, pf_charged_t pf_2);

ap_uint<7> weighted_avg_eta_p_p_p(pf_charged_t pf1, pf_charged_t pf2, pf_charged_t pf3);

ap_uint<8> weighted_avg_phi_p_p_p(pf_charged_t pf1, pf_charged_t pf2, pf_charged_t pf3);

ap_uint<7> weighted_avg_eta_s_s(strip_t strip1, strip_t strip2);

ap_uint<8> weighted_avg_phi_s_s(strip_t strip1, strip_t strip2);

ap_uint<7> weighted_avg_eta_c_c(cluster_t cluster1, cluster_t cluster2);

ap_uint<8> weighted_avg_phi_c_c(cluster_t cluster1, cluster_t cluster2);

ap_uint<7> weighted_avg_eta_p_s(pf_charged_t strip1, strip_t strip2);

ap_uint<8> weighted_avg_phi_p_s(pf_charged_t strip1, strip_t strip2);

ap_uint<7> weighted_avg_eta_c_s(cluster_t cluster1, strip_t strip2);

ap_uint<8> weighted_avg_phi_c_s(cluster_t cluster1, strip_t strip2);

ap_uint<8> weighted_avg_phi_c_p(cluster_t cluster1, pf_charged_t cluster2);

ap_uint<7> weighted_avg_eta_c_p(cluster_t cluster1, pf_charged_t cluster2);

ap_uint<8> weighted_avg_phi_t_s(pftau_t strip1, strip_t strip2);

ap_uint<7> weighted_avg_eta_t_s(pftau_t strip1, strip_t strip2);

ap_uint<7> ieta_diff(pf_charged_t cand_1, pf_charged_t cand_2);

ap_uint<8> iphi_diff(pf_charged_t cand_1, pf_charged_t cand_2);

ap_uint<10> delta_r_strip(strip_t strip1, strip_t strip2);

ap_uint<10> delta_r_cluster(cluster_t cluster1, cluster_t cluster2);

ap_uint<8> delta_r_c_p(cluster_t pf_1, pf_charged_t pf_2);

void merge_strip_algo(cluster_t cluster_1, pf_charged_t electron_1, cluster_t cluster_2, pf_charged_t electron_2, strip_t &strip, algo_config_t algo_config);

void strip_alg(pftau_t &tau_cand, pf_charged_t electron_grid[5][5], cluster_t neutral_clusters[N_CLUSTERS], algo_config_t algo_config);

void tau_three_prong_alg(track_t central_tracks[N_TRACKS], track_t three_prong_tau_cand[3], algo_config_t algo_config);

void pf_match_alg(cluster_t central_clusters[N_CLUSTERS], track_t central_tracks[N_TRACKS] , pf_charged_t charged_cands[N_TRACKS],  algo_config_t algo_config);

void check_pf_cand(pf_charged_t &pf_charged, cluster_t &central_cluster, algo_config_t algo_config);
// Tau_alg Takes in charged cands, central (neutral) clusters, the algorithm configuration 
// Returns the tau cands

void tau_alg(pf_charged_t pf_charged[N_TRACKS], cluster_t neutral_clusters[N_CLUSTERS], algo_config_t algo_config, pftau_t tau_cands[N_TAUS]);


void build_electron_grid(pf_charged_t electron_cand, ap_uint<8> seed_cand_dr, pf_charged_t seed_hadron, pf_charged_t electron_grid[N_TAUS][5][5], ap_uint<4> n_taus);
void find_tau_prongs(ap_uint<3> &n_prongs_found, pf_charged_t prong_cands[3], pf_charged_t pf_charged_hadron_signal_cand, pf_charged_t seed_hadron, ap_uint<8> seed_cand_dr,ap_uint<12> iso_sum_charged_hadron, algo_config_t algo_config);

#endif
