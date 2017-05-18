;;; The exported bindings are mostly C functions defined in phis-scm.c.
;;; Some bindings need some additional scheme code. In this case, the C wrapper
;;; has the prefix c-.
;;; I have omitted a "phis" prefix from the bindings, but the application may 
;;; add one on import.

(define-module (phis)
; FIXME: these are needed for the array stuff only.
;  :use-module (srfi srfi-25) ; array library: not part of the guile distribution
;  :use-module (srfi srfi-11) ; multiple value handling
  :export (multab
           init load-vpqrs
           get-info get-epsi get-sym get-occ
           vpqrs 
;--- FIXME: deal with arrays later
;           get-ao get-scfvec get-overlap
           get-geometry))

(load-extension "phis-scm" "scm_phis_mod_init")

;--- I prefer to use the closed-form for the multiplication table below
;    instead of a wrapper around the table in phis/init.c.
(define (multab i j)
  (+ (logxor (- i 1) (- j 1) 1)))


(define (init backend flags)
  ;; The backend should be given as a single symbol (cf. backend-alist below),
  ;; flags as a list of symbols (cf. flags-alist).
  ;; FIXME: flags and backends should be implemented as a wrapper around
  ;;        the C definitions to avoid inconsistencies.
  (let* ((backend-alist '((MOLCAS 0)
                          (GUK 1)
                          (STUB 2)))
         (flags-alist '((SHOW_INACTIVE #x0001)
                        (SYM_BLOCKED   #x0002)))
         ;; Should I check for duplicates here?
         (encode-flags (lambda (lst)
                         (apply + (map 
                                   (lambda (x) 
                                     (cadr (assoc x flags-alist)))
                                   lst)))))
    (c-init
     (cadr (assoc backend backend-alist))
     (encode-flags flags))))


;--- FIXME: There is no proper support for arrays (SRFI-25),
;           so I omit the functions that are supposed to return arrays.
;(define (get-scfvec)
;  (let-values (((n len v) (c-get-scfvec)))
;              (apply array (shape 0 n 0 len)
;                     (vector->list v))))
