/*
 * 
 * Some notes:
 *  - Functions return a simple Scheme type or a list if multiple values are expected.
 *  - Functions return #f if the feature is not present in the backend.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <libguile.h>
#include <guile/gh.h>

#include "phis.h"

static long phis_cap;

/* Prototypes for the Scheme wrappers */     
SCM scm_phis_init(SCM backend, SCM flags);
SCM scm_phis_load_Vpqrs();
SCM scm_phis_get_info();
SCM scm_phis_get_epsi();
SCM scm_phis_get_sym();
SCM scm_phis_get_occ();
SCM scm_phis_Vpqrs(SCM p, SCM q, SCM r, SCM s);
SCM scm_phis_get_ao();
SCM scm_phis_get_scfvec();
SCM scm_phis_get_overlap();
SCM scm_phis_get_geometry();


void     
scm_phis_mod_init()
{
        /* Define Scheme interfaces to C functions. */
        gh_new_procedure("c-init", scm_phis_init, 2, 0, 0);
        gh_new_procedure("load-vpqrs", scm_phis_load_Vpqrs, 0 ,0, 0);
        gh_new_procedure("get-info", scm_phis_get_info, 0, 0, 0);
        gh_new_procedure("get-epsi", scm_phis_get_epsi, 0, 0, 0);
        gh_new_procedure("get-sym", scm_phis_get_sym, 0, 0, 0);
        gh_new_procedure("get-occ", scm_phis_get_occ, 0, 0, 0);
        gh_new_procedure("vpqrs", scm_phis_Vpqrs, 4, 0, 0);
        gh_new_procedure("get-ao", scm_phis_get_ao, 0, 0, 0);
        gh_new_procedure("c-get-scfvec", scm_phis_get_scfvec, 0, 0, 0);
        gh_new_procedure("get-overlap", scm_phis_get_overlap, 0, 0, 0);
        gh_new_procedure("get-geometry", scm_phis_get_geometry, 0, 0, 0);

}

/*************************/

SCM 
scm_phis_init(SCM backend, SCM flags)
{
        int b = gh_scm2long(backend);
        int f = gh_scm2long(flags);
        int c_flags = (f<<16) + b;
        phis_cap = phis_init(&c_flags); 

        return gh_ulong2scm(phis_cap);
}


SCM
scm_phis_load_Vpqrs()
{
        phis_load_Vpqrs();
        return SCM_UNDEFINED;
}


SCM 
scm_phis_get_info()
{
        int nSym, nBas, nCenters;

        if (!(phis_cap & HAVE_INFO))
                return SCM_BOOL_F;

        phis_get_info(&nSym, &nBas, &nCenters);
        return gh_list(
                gh_long2scm(nSym),
                gh_long2scm(nBas),
                gh_long2scm(nCenters),
                SCM_UNDEFINED);
}

SCM
scm_phis_get_epsi()
{
        SCM epsi;
        int n = 0;
        double Ehf;
        double *e;

        if (!(phis_cap & HAVE_EPSI))
                return SCM_BOOL_F;

        phis_get_epsi(NULL,NULL,&n);
        n = -n;
        e = (double *)malloc_we(n*sizeof(double));
        phis_get_epsi(&Ehf,e, &n);
        if (n<0) {
                perror("phis-get-epsi");
                exit(-1);
        }

        epsi = gh_doubles2scm(e, n);

        free(e);
        return gh_list(
                gh_double2scm(Ehf),
                epsi,
                SCM_UNDEFINED);
}

SCM
scm_phis_get_sym()
{
        SCM sym;
        int n = 0;
        int *s;

        if (!(phis_cap & HAVE_SYM))
                return SCM_BOOL_F;

        phis_get_sym(NULL,&n);
        n = -n;
        s = (int *)malloc_we(n*sizeof(int));
        phis_get_sym(s, &n);
        if (n<0) {
                perror("phis-get-sym");
                exit(-1);
        }

        sym = gh_ints2scm(s, n);

        free(s);
        return sym;
}

SCM
scm_phis_get_occ()
{
        SCM occ;
        int n = 0;
        double *o;

        if (!(phis_cap & HAVE_EPSI))
                return SCM_BOOL_F;

        phis_get_occ(NULL,&n);
        n = -n;
        o = (double *)malloc_we(n*sizeof(double));
        phis_get_occ(o, &n);
        if (n<0) {
                perror("phis-get-occ");
                exit(-1);
        }

        occ = gh_doubles2scm(o, n);

        free(o);
        return occ;
}

SCM
scm_phis_Vpqrs(SCM p, SCM q, SCM r, SCM s) 
{
        int i = gh_scm2ulong(p);
        int j = gh_scm2ulong(q);
        int k = gh_scm2ulong(r);
        int l = gh_scm2ulong(s);

        return gh_double2scm(Vpqrs(&i,&j,&k,&l));
}

SCM
scm_phis_get_ao()
{
        int i, nao = 0, max_ncc = 0;

        int *ncc, *c_poly, *c_center;
        double *c_cc, *c_alpha;

        SCM poly, cc, alpha, center;

        if (!(phis_cap & HAVE_AO))
                return SCM_BOOL_F;

        /* ask the backend for the dimensions */
        phis_get_ao(NULL,NULL,NULL,NULL,NULL,&nao,&max_ncc);
        nao = abs(nao); max_ncc = abs(max_ncc);

        /* allocate some memory for the arrays */
        if ((ncc = (int*) malloc(nao*sizeof(int))) == NULL ||
            (c_poly = (int*) malloc(nao*3*sizeof(int))) == NULL ||
            (c_cc = (double*) malloc(nao*max_ncc*sizeof(double))) == NULL ||
            (c_alpha = (double*) malloc(nao*max_ncc*sizeof(double))) == NULL ||
            (c_center = (int*) malloc(nao*sizeof(int))) == NULL) {
                perror("malloc");
                exit(-1);
        }

        /* ask the backend for the AOs */
        phis_get_ao(c_poly, ncc, c_cc, c_alpha, c_center, &nao, &max_ncc);
        
        /* error check */
        if(nao<=0 || max_ncc<=0) {
                fprintf(stderr,"phis-get-ao: Problem getting AOs.\n"
                        "             backend returned %i %i\n",
                        nao, max_ncc);
                exit(-1);
        }

        /* convert everything to Scheme */
        poly = gh_make_vector(gh_long2scm(nao), SCM_BOOL_F);
        cc = gh_make_vector(gh_long2scm(nao), SCM_BOOL_F);
        alpha = gh_make_vector(gh_long2scm(nao), SCM_BOOL_F);
        for (i=0; i<nao; i++) {
                gh_vector_set_x (poly, 
                                 gh_long2scm(i),
                                 gh_ints2scm((c_poly + i*3), 3));
                gh_vector_set_x (cc,
                                 gh_long2scm(i),
                                 gh_doubles2scm((c_cc + i*max_ncc), ncc[i]));
                gh_vector_set_x (alpha,
                                 gh_long2scm(i),
                                 gh_doubles2scm((c_alpha + i*max_ncc), ncc[i]));
        }
        center = gh_ints2scm(c_center, nao);

        /* collecting garbage by hand :( */
        free(ncc);
        free(c_poly);
        free(c_cc);
        free(c_alpha);
        free(c_center);

        return gh_list(poly, cc, alpha, center, SCM_UNDEFINED);
}

SCM
scm_phis_get_scfvec()
{
        int n = 0, len = 0;

        double *c_scfvec;
        SCM scfvec;

        if (!(phis_cap & HAVE_SCFVEC))
                return SCM_BOOL(0);

        phis_get_scfvec(NULL, &n, &len);

        /* allocate memory */
        n = abs(n); len = abs(len);
        if ((c_scfvec = (double*) malloc(n*len * sizeof(double))) == NULL) {
                perror("malloc");
                exit(-1);
        }

        /* get vectors */
        phis_get_scfvec(c_scfvec, &n, &len);

        /* error check */
        if (n<=0 || len<=0) {
                fprintf(stderr,"phis-get-scfvec: Problem getting SCF vectors.\n"
                        "             backend returned %i %i\n",
                        n, len);
                exit(-1);
        }

        /* transform the matrix to Scheme somehow. */
        /* also: use a nice array wrapper in C */

        /* This makes a vector of vectors.
           Not very practical, better use a single vector and an affine mapping
           as in SRFI-25.
        scfvec = gh_make_vector(gh_long2scm(n), SCM_BOOL_F);
        for (i=0; i<n; i++)
                gh_vector_set_x (scfvec, 
                                 gh_long2scm(i),
                                 gh_doubles2scm((c_scfvec + i*len), len));
        */

        /* Return a single vector; is transformed into an array in the Scheme part.
           FIXME: should implement arrays natively in guile. */
        scfvec = gh_doubles2scm(c_scfvec, n*len);
                         
        free(c_scfvec);
        return scm_values(
                scm_list_3(gh_long2scm(n), gh_long2scm(len), scfvec));
}

SCM
scm_phis_get_overlap()
{
        int n = 0;

        double *c_overlap;
        SCM overlap;

        if (!(phis_cap & HAVE_OVERLAP))
                return SCM_BOOL(0);

        phis_get_overlap(NULL, &n);

        /* allocate memory */
        n = abs(n);
        if ((c_overlap = (double*) malloc(n*n * sizeof(double))) == NULL) {
                perror("malloc");
                exit(-1);
        }

        /* get vectors */
        phis_get_overlap(c_overlap, &n);

        /* error check */
        if (n<=0 ) {
                fprintf(stderr, "phis-get-overlap: Problem getting overlap matrix.\n"
                                "                  backend returned %i\n", n);
                exit(-1);
        }

        /* FIXME: see get-scfvec */
        overlap = gh_doubles2scm(c_overlap, n*n);

        free(c_overlap);
        return overlap;
}

SCM
scm_phis_get_geometry() {

        int i, n = 0;
        double *geo, *Z_nuc;
        double E_nuc;

        SCM Z, geometry;

        phis_get_geometry(&n, NULL, NULL, NULL);

        n = abs(n);
        if((geo = (double *) malloc(3*n*sizeof(double))) == NULL ||
           (Z_nuc = (double *) malloc(n*sizeof(double))) == NULL) {
                perror("malloc");
                exit(-1);
        }

        phis_get_geometry(&n, geo, Z_nuc, &E_nuc);

        geometry = gh_make_vector(gh_long2scm(n), SCM_BOOL_F);
        for (i=0; i<n; i++)
                gh_vector_set_x (geometry, 
                                 gh_long2scm(i),
                                 gh_doubles2scm((geo + i*3), 3));
        Z = gh_doubles2scm(Z_nuc, n);
                

        free(geo);
        free(Z_nuc);

        return gh_list(gh_double2scm(E_nuc), Z, geometry, SCM_UNDEFINED);
}
