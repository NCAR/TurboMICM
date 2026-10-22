/* Shared standard-C CSV and argument helpers for the benchmark drivers. */
static char *read_line(FILE *file)
{
    size_t length=0, capacity=4096;
    char *s=malloc(capacity);
    if (!s) return NULL;
    int c;
    while ((c=fgetc(file))!=EOF && c!='\n') {
        if (length+1==capacity) {
            if (capacity>SIZE_MAX/2) { free(s); return NULL; }
            capacity*=2;
            char *next=realloc(s,capacity);
            if (!next) { free(s); return NULL; }
            s=next;
        }
        s[length++]=(char)c;
    }
    if (!length && c==EOF) { free(s); return NULL; }
    if (length && s[length-1]=='\r') --length;
    s[length]='\0'; return s;
}

/* CSV fields are compacted in place; double quotes and escaped quotes are accepted. */
static size_t split_csv(char *s, char **fields, size_t limit)
{
    size_t count=0;
    char *read=s,*write=s;
    for (;;) {
        if (count==limit) return CHEM_NONE;
        fields[count++]=write;
        int quoted=*read=='"';
        if (quoted) ++read;
        while (*read) {
            if (quoted && *read=='"') {
                if (read[1]=='"') { *write++='"'; read+=2; continue; }
                ++read; quoted=0; break;
            }
            if (!quoted && *read==',') break;
            *write++=*read++;
        }
        if (quoted || (*read && *read!=',')) return CHEM_NONE;
        int more=*read==',';
        if (more) ++read;
        *write++='\0';
        if (!more) break;
    }
    return count;
}

static int read_cells(const char *path, chem_state *states, size_t group_size, size_t cells)
{
    FILE *file=fopen(path,"r");
    if (!file) { fprintf(stderr,"Cannot open input: %s\n",path); return 0; }
    const chem_model *m=states[0].model;
    size_t limit=m->n_species+m->n_custom+2;
    char **fields=calloc(limit,sizeof(char *)); size_t *indices=calloc(limit,sizeof(size_t));
    int *kind=calloc(limit,sizeof(int)); char *line=read_line(file);
    int success=0;
    if (!fields || !indices || !kind || !line) goto done;
    size_t columns=split_csv(line,fields,limit);
    if (columns==CHEM_NONE || !columns) goto done;
    for (size_t i=0; i<columns; ++i) {
        for (size_t j=0; j<i; ++j) if (!strcmp(fields[i],fields[j])) goto done;
        if (!strcmp(fields[i],"temperature")) kind[i]=0;
        else if (!strcmp(fields[i],"pressure")) kind[i]=1;
        else if ((indices[i]=chem_species_index(m,fields[i]))!=CHEM_NONE) kind[i]=2;
        else if ((indices[i]=chem_custom_index(m,fields[i]))!=CHEM_NONE) kind[i]=3;
        else { fprintf(stderr,"Unknown input column: %s\n",fields[i]); goto done; }
    }
    free(line); line=NULL;
    for (size_t row=0; row<cells;) {
        line=read_line(file);
        if (!line) { fprintf(stderr,"Input has fewer than %zu cells\n",cells); goto done; }
        if (!*line) { free(line); line=NULL; continue; }
        if (split_csv(line,fields,limit)!=columns) goto done;
        chem_state *s=&states[row/group_size]; size_t c=row%group_size;
        for (size_t i=0; i<columns; ++i) {
            char *end; errno=0; double value=strtod(fields[i],&end);
            while (isspace((unsigned char)*end)) ++end;
            if (fields[i]==end || *end || errno==ERANGE || !isfinite(value)) goto done;
            if (kind[i]==0) s->conditions[c].temperature=value;
            else if (kind[i]==1) s->conditions[c].pressure=value;
            else if (kind[i]==2) s->y[indices[i]*s->cells+c]=value;
            else {
                for (size_t j=indices[i]; j<m->n_custom; ++j)
                    if (!strcmp(m->custom_names[j],m->custom_names[indices[i]])) s->custom[j*s->cells+c]=value;
            }
        }
        if (s->conditions[c].temperature<=0 || s->conditions[c].pressure<=0) goto done;
        chem_ideal_air_density(&s->conditions[c]);
        ++row; free(line); line=NULL;
    }
    success=1;
done:
    free(line); free(fields); free(indices); free(kind); fclose(file);
    if (!success) fprintf(stderr,"Invalid CSV input\n");
    return success;
}

static int write_final(chem_state *states, size_t groups)
{
    FILE *file=fopen("ts1_final.csv","w");
    if (!file) return 0;
    const chem_model *m=states[0].model;
    fprintf(file,"temperature,pressure");
    for (size_t i=0; i<m->n_species; ++i) fprintf(file,",%s",m->mechanism->species[m->var_to_spec[i]].name);
    fputc('\n',file);
    for (size_t g=0; g<groups; ++g) {
        chem_state *s=&states[g];
        for (size_t c=0; c<s->cells; ++c) {
            fprintf(file,"%.17g,%.17g",s->conditions[c].temperature,s->conditions[c].pressure);
            for (size_t i=0; i<m->n_species; ++i) fprintf(file,",%.17g",s->y[i*s->cells+c]);
            fputc('\n',file);
        }
    }
    int success=!ferror(file);
    if (fclose(file)) success=0;
    return success;
}

static int integer(const char *s, size_t *value, int zero_ok)
{
    if (*s=='-' || !*s) return 0;
    char *end; errno=0; uintmax_t v=strtoumax(s,&end,10);
    if (*end || errno || v>SIZE_MAX || (!zero_ok && !v)) return 0;
    *value=(size_t)v; return 1;
}

