#ifndef file_read_in_H_
#define file_read_in_H_

#include <stdint.h>
#include <ap_int.h>

//define number of tracks
#define N_TRACKS (50)
#define N_CLUSTERS (2880)

typedef struct
{
	ap_uint<11> et;
	ap_uint<1> eta_side;
	ap_uint<7> eta;
	ap_uint<8> phi;
	ap_uint<1> charged_hadron;
} track_t;

typedef struct
{
	ap_uint<11> et;
	ap_uint<1> eta_side;
	ap_uint<7> eta;
	ap_uint<8> phi;
	ap_uint<8> EoH;
	ap_uint<8> HoE;
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
	ap_uint<1> dummy;

} algo_config_t;

typedef struct{
	ap_uint<14> sum_tracks;
	ap_uint<14> three_prong_tau_et;
} algo_outputs_t;

void file_read_in(track_t tracks[N_TRACKS],
		algo_config_t algo_config,
		algo_outputs_t & algo_outputs
);


ap_uint<1> Delta_R(ap_uint<8> eta_1,
			ap_uint<8> phi_1,
			ap_uint<8> eta_2,
			ap_uint<8> phi_2,
			ap_uint<8> maximum_delta_R);

void tau_three_prong_alg(track_t central_tracks[N_TRACKS], track_t three_prong_tau_cand[3], algo_config_t algo_config);

void pf_match_alg(cluster_t central_clusters[N_CLUSTERS], track_t central_tracks[N_TRACKS] , pfcharged_t charged_cands[N_TRACKS]);

void strip_alg(pftau_t tau_cand, pfcharged_t electron_grid[5][5], pfneutral_t neutral_clusters[N_CLUSTERS]);

#endif
