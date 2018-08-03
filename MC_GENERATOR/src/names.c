#include <names.h>

static char lattice_name[MAX_BUF];
static int  lattice_name_length = 0;

/* ==================================================================== */
void lattice_set_name( double beta, uchar T, uchar L, uchar D ){
  lattice_name[0] = 0;
  sprintf( &lattice_name[strlen(lattice_name)], "b%.4f",  beta );
  sprintf( &lattice_name[strlen(lattice_name)], "_D%d",   D );
  sprintf( &lattice_name[strlen(lattice_name)], "_T%d",  T );
  sprintf( &lattice_name[strlen(lattice_name)], "_L%d",  L );
  lattice_name_length = strlen(lattice_name);
};
/* ==================================================================== */
char *lattice_get_name( char *ext ){
  char *ret;
  if( !lattice_name_length ) return NULL;
  sprintf( &lattice_name[lattice_name_length], ".%s", ext );
  ret = strdup( lattice_name );
  lattice_name[lattice_name_length] = 0;
  return ret;
};
/* ==================================================================== */
char *lattice_parse_name( char *name, char *ext,
			  double *beta, int *size0, int *size1, int *dim ){
  char scratch[MAX_BUF];
  char *p;
  int n;
  bzero( &scratch[0], MAX_BUF * sizeof(char));
  if( !name || !*name ) return NULL;
  if( strlen(name) > MAX_BUF/2 ) panic("Filename too long");
  if( ext && !*ext ) ext = NULL;
  strcat( scratch, name );
  n = strlen(name) - 1;
  if( ext ) n -= strlen(ext);
  if( n <= 0 ) return NULL;
  p = &scratch[n];
  if( ext ){
    if( !p || *p++ != '.' || strcmp( p, ext ) ) return NULL;
    *(--p) = 0;
  };
  while( p > scratch && *p != '/' ) p--;
  if( *p == '/' ) p++;
  //----------------------------------
  if( sscanf( p, "b%lf_D%d_T%d_L%d.%s", beta, dim, size0, size1, &scratch[MAX_BUF/2] ) == 5 ){
    p = strdup( &scratch[MAX_BUF/2] );
    return p;
  };
  if( sscanf( p, "b%lf_D%d_T%d_L%d", beta, dim, size0, size1 ) == 4 ){
    p = strdup("");
    return p;
  };
  return NULL;
};
