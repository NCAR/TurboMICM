#ifndef CK_TUTORIALS_H
#define CK_TUTORIALS_H
#include "chem.h"
#include <stdio.h>
const ck_mechanism *tutorial_hand_rates(int user);
const ck_mechanism *tutorial_hand_robertson(void);
int run_tutorial(const char *name, FILE *output, int configuration);
int test_configuration_equivalence(int which);
int test_tutorial(const char *name);
#endif
