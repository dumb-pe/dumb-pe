

(define-module (dumb-pe field-units)
  #:export (api-from-sg
            r-from-f))

(define (api-from-sg sg)
  (- (/ 141.5 sg) 131.5))

(define (r-from-f f)
  (+ f 459.67))
