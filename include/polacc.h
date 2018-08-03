#ifndef _POLACC_H_
#define _POLACC_H_

typedef struct t_polacc_data {
	int deg;		//!< degree of polynomial 
	double cf[2];	//!< coefficients of mapping \f$ p(x) \f$
	t_op op;
	void* data;		//!< op internal data, e.g. wilson \f$ \kappa \f$
	double sigma1;	//!< \f$ \sigma_1 \f$ for Chebyshev iteration (see Saad)
	double c, e;	//!< see Y.Saad page 223
} t_polacc_data;


void build_mapping(double m, double l, t_polacc_data* data);

void build_mapping_ov(double lambda1, double xb, t_polacc_data* data);

void pol_acc(void* data, t_gauge_vector* gv, t_cds_vector* out, t_cds_vector* in);

void pol_acc_ov(void* data, t_gauge_vector* gv, t_cds_vector* out, t_cds_vector* in);

void orig_eval_hw(t_gauge_vector* gv, t_op_ev* oev);

#endif
