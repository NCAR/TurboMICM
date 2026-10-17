/*
 * tmc: MICM's State::PrintHeader and State::PrintState, which MICM's tutorial
 * programs use to print a state as a table. The header is "  time" (width 6),
 * then ",  grid" when the state has more than one cell, then each variable
 * name right-aligned in a column; the column width is 11, or the longest name
 * plus 2 when a name has 10 or more characters. Each state row is the time
 * (%g, width 6), the cell index (width 6) when there is more than one cell,
 * and each variable's value in scientific notation with two decimals, all
 * separated by commas and ending in a newline; one row per cell.
 */
#ifndef TMC_TUTORIAL_EXT_H
#define TMC_TUTORIAL_EXT_H

#include "tmc.h"

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* The header or the state rows formatted into buffer: at most size - 1
 * characters and a terminating NUL are written (nothing when size is 0), and
 * the length of the whole text is returned, as snprintf does. */
size_t tmc_state_format_header(const tmc_state *state, char *buffer, size_t size);
size_t tmc_state_format_state(const tmc_state *state, double time, char *buffer, size_t size);

/* The same text written to a stream (MICM writes to std::cout); returns the
 * number of characters written, or a negative value on an output error. */
int tmc_state_print_header(const tmc_state *state, FILE *stream);
int tmc_state_print_state(const tmc_state *state, double time, FILE *stream);

#ifdef __cplusplus
}
#endif

#endif /* TMC_TUTORIAL_EXT_H */
