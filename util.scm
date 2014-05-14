(define-module (dumb-pe util)
  #:export (test-bound-throw))


(define-syntax test-bound-throw
  (syntax-rules ()
    ((test-bound-throw var lower upper)
     (if (or (< var lower) (> var upper))
         (throw 'value-out-of-bound
                (format #f "~a = ~a out of (~a, ~a)" 'var var lower upper)
                var lower upper)))))
