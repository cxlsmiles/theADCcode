;;; Scheme version of the MP2 example.
;;; Written basically to learn more about Scheme. This code needs a lot of clean-up
;;; and is untested. Don't look too closely ;)
;;; FIXME: Document the general layout here.
;;; 
;;; NOTE: this needs LD_LIBRARY_PATH and GUILE_LOAD_PATH to contain <phisdir>/scm.

(use-modules ((phis) :renamer (symbol-prefix-proc 'phis:)))

;--- utility functions ---
(define (++ i) (+ i 1))
(define (-- i) (- i 1))

(define (vector-for-each/index f . vecs)
  (let ((len (vector-length (car vecs))))
    (do ((i 0 (++ i)))
	((= i len))
      (apply f (cons i (map (lambda (v) (vector-ref v i))
			    vecs))))))

;--- end utility functions ---


(define (sort-orbitals sym occ nsym)
  "Returns a vector of length nSym. 
   Each element is a list of two lists, first the occupied then the virtual orbitals."
  (let ((sorted (make-vector nsym '(() ()))))
    (vector-for-each/index
     (lambda (i s o)
       (vector-set! sorted (-- s)
                    (case o
                      ((2.0) (list 
                              (cons i (car (vector-ref sorted (-- s))))
                              (cadr (vector-ref sorted (-- s)))))
                      ((0.0) (list
                              (car (vector-ref sorted (-- s)))
                              (cons i (cadr (vector-ref sorted (-- s))))))
                      (else (error "Only closed shell allowed.")))))
     sym occ)
    sorted))

; Needs a tail-recursive version or nsym=8 causes a stack overflow.
(define (generate-sym-blocks nsym)
  "Generates a list of 4-tupels of irreps the product of which transforms
as the totally symmetric representation."

  (define (4-tupel p q r)
    "Returns a totally symmetric 4-tupel when three irreps are given."
    (list p q r (phis:multab (phis:multab p q) r)))

  (let loop ((psym (iota nsym))
             (qsym (iota nsym))
             (rsym (iota nsym))
             (seed '()))
    (cond
     ((null? psym) seed)
     ((null? qsym) (loop (cdr psym) (iota nsym) (iota nsym)
                         seed))
     ((null? rsym) (loop psym (cdr qsym) (iota nsym)
                         seed))
     (else (loop psym qsym (cdr rsym)
                 (cons (4-tupel (car psym) (car qsym) (car rsym))
                       seed))))))


(define (ijab-loop/sym func sym)
  "Applies func to all orbital tupels (ij|ab) of the given symmetry
and sums up the results"
  ;;--- FIXME: is it really clever to call the PHIS wrapper for every symmetry block?
  (let* ((orbitals-sorted (sort-orbitals (phis:get-sym) (phis:get-occ) (car (phis:get-info))))
         (all-iorb (car (vector-ref orbitals-sorted (car sym))))
         (all-jorb (car (vector-ref orbitals-sorted (cadr sym))))
         (all-aorb (cadr (vector-ref orbitals-sorted (caddr sym))))
         (all-borb (cadr (vector-ref orbitals-sorted (caddr sym)))))
    (let loop ((iorb all-iorb)
               (jorb all-jorb)
               (aorb all-aorb)
               (borb all-borb)
               (sum 0))
      (cond
       ((null? iorb) sum)
       ((null? jorb) (loop (cdr iorb) all-jorb all-aorb all-borb sum))
       ((null? aorb) (loop iorb (cdr jorb) all-aorb all-borb sum))
       ((null? borb) (loop iorb jorb (cdr aorb) all-borb sum))
       (else (loop iorb jorb aorb (cdr borb) 
                   (+ (func (car iorb) (car jorb) (car aorb) (car borb)) sum)))))))


(define (calculate-mp2)

  ;; FIXME: an evil and fragile hack!!
  (case (string->symbol (cadr (command-line)))
    ((MOLCAS) (phis:init 'MOLCAS '(SYM_BLOCKED)))
    ((GUK) (phis:init 'GUK '()))
    (else (error "Sorry, no help!")))
  (phis:load-vpqrs)

  (let ((epsi (cadr (phis:get-epsi))))

    (define (mp2-contribution i j a b)
      "Calculate the second-order MP-correction for the given set of orbital indices."
      (let ((iajb (phis:vpqrs i a j b))
            (ibja (phis:vpqrs i b j a)))
        (/ 
         (- (* 2 iajb ibja) (* iajb ibja))
         (+ (vector-ref epsi i)
            (vector-ref epsi j)
            (- (vector-ref epsi a))
            (- (vector-ref epsi b))))))

    (apply + (map 
              (lambda (sym-block) (ijab-loop/sym mp2-contribution sym-block))
              (generate-sym-blocks (car (phis:get-info)))))))
