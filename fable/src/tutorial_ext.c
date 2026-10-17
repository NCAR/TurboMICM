/* MICM's State::PrintHeader and PrintState (solver/state.inl), see
 * tmc_tutorial_ext.h. */
#include "tmc_tutorial_ext.h"

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

/* snprintf-style accumulation into a buffer that may be too small or NULL */
typedef struct writer
{
  char *buffer;
  size_t size;
  size_t length;
} writer;

static void write_format(writer *w, const char *format, ...)
{
  va_list args;
  va_start(args, format);
  char *at = w->buffer && w->length < w->size ? w->buffer + w->length : NULL;
  size_t room = at ? w->size - w->length : 0;
  int n = vsnprintf(at, room, format, args);
  va_end(args);
  if (n > 0)
    w->length += (size_t)n;
}

/* MICM's column width: 11, or the longest name plus 2 from 10 characters up */
static int column_width(const tmc_state *state)
{
  size_t largest = 0;
  const int n = tmc_state_num_variables(state);
  for (int i = 0; i < n; ++i)
  {
    const size_t length = strlen(tmc_state_variable_name(state, i));
    if (length > largest)
      largest = length;
  }
  return largest < 10 ? 11 : (int)largest + 2;
}

size_t tmc_state_format_header(const tmc_state *state, char *buffer, size_t size)
{
  writer w = { buffer, size, 0 };
  if (buffer && size > 0)
    buffer[0] = '\0';
  const int width = column_width(state);
  write_format(&w, "%6s", "time");
  if (tmc_state_num_cells(state) > 1)
    write_format(&w, ",%6s", "grid");
  const int n = tmc_state_num_variables(state);
  for (int i = 0; i < n; ++i)
    write_format(&w, ",%*s", width, tmc_state_variable_name(state, i));
  write_format(&w, "\n");
  return w.length;
}

size_t tmc_state_format_state(const tmc_state *state, double time, char *buffer, size_t size)
{
  writer w = { buffer, size, 0 };
  if (buffer && size > 0)
    buffer[0] = '\0';
  const int width = column_width(state);
  const int n_cells = tmc_state_num_cells(state), n = tmc_state_num_variables(state);
  for (int cell = 0; cell < n_cells; ++cell)
  {
    write_format(&w, "%6g,", time);
    if (n_cells > 1)
      write_format(&w, "%6d,", cell);
    for (int i = 0; i < n; ++i)
      write_format(&w, "%s%*.2e", i > 0 ? "," : "", width, tmc_state_get(state, cell, i));
    write_format(&w, "\n");
  }
  return w.length;
}

static int print_text(FILE *stream, size_t (*format)(const tmc_state *, double, char *, size_t), const tmc_state *state, double time)
{
  const size_t length = format(state, time, NULL, 0);
  char *text = malloc(length + 1);
  if (!text)
    return -1;
  format(state, time, text, length + 1);
  const int written = fputs(text, stream) < 0 ? -1 : (int)length;
  free(text);
  return written;
}

static size_t format_header_adapter(const tmc_state *state, double time, char *buffer, size_t size)
{
  (void)time;
  return tmc_state_format_header(state, buffer, size);
}

int tmc_state_print_header(const tmc_state *state, FILE *stream)
{
  return print_text(stream, format_header_adapter, state, 0.0);
}

int tmc_state_print_state(const tmc_state *state, double time, FILE *stream)
{
  return print_text(stream, tmc_state_format_state, state, time);
}
