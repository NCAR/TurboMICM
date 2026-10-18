#ifndef CHEM_ADAPTIVE_H
#define CHEM_ADAPTIVE_H
#include "chem.h"

/* Experimental autonomous ODE integrators. No algebraic constraints. */
typedef enum { CHEM_ADAPT_ROS, CHEM_ADAPT_EXTRAP, CHEM_ADAPT_SDIRK } chem_adaptive_method;
typedef struct {
    chem_adaptive_method method;
    chem_ros_method ros;
    size_t max_steps, jacobian_age, extrap_order;
    double safety, growth, h_start, h_max, reuse_window;
    int controller, error_norm, clip_stages, warm_start;
    size_t newton_iterations;
    double newton_tolerance;
    int adaptive_order;
} chem_adaptive_parameters;
typedef struct {
    const chem_model *model;
    size_t cells;
    double *jacobian, *vectors;
    double next_h;
    uint64_t jacobian_reuses, factor_reuses, positivity_rejections;
    uint64_t nonfinite_rejections, clipped_values;
    uint64_t nonlinear_iterations, nonlinear_rejections;
    uint64_t extrap_orders[3];
    double clipped_mass;
} chem_adaptive_workspace;
chem_adaptive_parameters chem_adaptive_defaults(void);
chem_status chem_adaptive_create(chem_adaptive_workspace *, const chem_state *);
void chem_adaptive_destroy(chem_adaptive_workspace *);
chem_result chem_adaptive_solve(chem_state *, double,
    const chem_adaptive_parameters *, chem_adaptive_workspace *);
#endif
