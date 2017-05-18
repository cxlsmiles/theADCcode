
#ifndef READ_SCF_DATA_H
#define READ_SCF_DATA_H

/* $Id: read_scf_data.h,v 1.2 2003/12/01 14:55:46 joergb Exp $
 *
 * $Log: read_scf_data.h,v $
 * Revision 1.2  2003/12/01 14:55:46  joergb
 * avoid double parsing using "#ifndef"
 *
 * Revision 1.1  2003/10/05 08:49:46  joergb
 * splitted into globals.h, adc_macros.h, get_input.h,
 *
 */



/* scf_t:
 * structure for the SCF information
 */
struct scf_t {
  int    nSym;           /* number of symmetries */
  int    nBas;           /* number of orbitals ? */
  int    nAtoms;         /* number of atoms */
  int    *loa;           /* list of active orbitals */
  double e_hf;           /* total HF energy */
  double *epsi;          /* list of orbital energies */
  double *occ;           /* list of occupation numbers */ 
  int    *sym;           /* list of symmetries */
};

/* void read_scf_data(char *backend, int debug, FILE *data_fp, struct scf_t *scf); */
/* //void read_scf_data(struct inp_t *inp, struct scf_t *scf); */

#endif
