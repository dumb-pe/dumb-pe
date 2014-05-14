;;; Commentary:

;;; These correlations are based on the paper "Reservoir-Fluid
;;; Property Correlations - State of the Art" by W.D. McCain Jr. (SPE
;;; 18571)

;;; Code:

(define-module (dumb-pe bo-corr)
  #:use-module ((dumb-pe field-units))
  #:use-module ((dumb-pe util))
  #:export (stock-tank-gor
            bubble-point-pressure
            solution-gor
            saturated-oil-fvf
            saturated-oil-fvf2))


(define (stock-tank-gor stock-tank-oil-sg
                        separator-gas-sg
                        separator-pressure
                        separator-temp)
  "Estimation of stock-tank GOR.
Arguments:
    stock-tank-oil-sg: stock tank oil specific gravity (dimensionless)
    separator-gas-sg: separator gas specific gravity (dimensionless)
    separator-pressure: separator pressure (psia)
    separator-temp: separator temperature (F)
Return:
    stock tank GOR value (scf/STB)

Often the reported values of producing GOR do not include stock-tank
vent gas. In this case, the use of initial producing GOR for solution
GOR results in values that are low by 10% or more. Addition of the
estimate of stock-tank GOR from this correlation to the separater GOR
results in an estimate of solution GOR acurate to within 3%. This
correlation should not be used if the separator temperature is >140F."
  (test-bound-throw separator-temp 0 140)
  (expt 10 (+ 0.3818
              (* -5.506 (log10 stock-tank-oil-sg))
              (* 2.902 (log10 separator-gas-sg))
              (* 1.327 (log10 separator-pressure))
              (* -0.7355 (log10 separator-temp)))))


(define (bubble-point-pressure initial-solution-gor
                               gas-sg
                               stock-tank-oil-api
                               reservoir-temp)
  "Estimation of reservoir oil bubble point pressure to an accuracy of 15%.
Arguments:
    initial-solution-gor: initial solution gas oil ratio (scf/STB)
    gas-sg: (separator) gas specific gravity (dimensionless)
    stock-tank-oil-api: stock-tank oil API gravity (API)
    reservoir-temp: reservoir temperature (F)
Return:
    bubble point pressure (psia)

The initial solution gas oil ratio should include stock-tank vent gas.
"
  (test-bound-throw reservoir-temp 0 325)
  (* 18.2 (- (* (expt (/ initial-solution-gor gas-sg) 0.83)
                (expt 10 (- (* 0.00091 reservoir-temp)
                            (* 0.0125 stock-tank-oil-api)))) 1.4)))


(define (solution-gor gas-sg
                      stock-tank-oil-api
                      reservoir-temp
                      pressure)
  "Estimation of solution gas oil ratio for pressures below the bubble
point pressure to an accuracy of 15%.
Arguments:
    gas-sg: (separator) gas specific gravity (dimensionless)
    stock-tank-oil-api: stock-tank oil API gravity (API)
    reservoir-temp: reservoir temperature (F)
    pressure: pressure (psia)
Return:
    solution gas oil ratio (scf/STB)"
  (test-bound-throw reservoir-temp 0 325)
  (* gas-sg (expt (/ (+ (/ pressure 18.2) 1.4)
                               (expt 10 (- (* 0.00091 reservoir-temp)
                                           (* 0.0125 stock-tank-oil-api))))
                            (/ 1 0.83))))


(define (solution-gor2 gas-sg
                       stock-tank-oil-api
                       reservoir-temp
                       initial-solution-gor
                       field-derived-pb
                       pressure)
  "An improved estimation of solution GOR if a field-derived
bubble-point pressure has been obtained from pressure measurement.
Arguments:
    gas-sg: (separator) gas specific gravity (dimensionless)
    stock-tank-oil-api: stock-tank oil API gravity (API)
    reservoir-temp: reservoir temperature (F)
    initial-solution-gor: initial solution GOR (scf/STB)
    field-derived-pb: field derived bubble point pressure (psia)
    pressure: pressure (below bubble point) (psia)
Return:
    solution GOR (scf/STB)

This function works very well for pressures near the bubble point. It
is less accurate at low pressures.
"
  (test-bound-throw pressure 0 field-derived-pb)
  (let ((calculated-pb (bubble-point-pressure initial-solution-gor
                                              gas-sg
                                              stock-tank-oil-api
                                              reservoir-temp)))
    (solution-gor gas-sg stock-tank-oil-api reservoir-temp
                  (+ pressure (- calculated-pb field-derived-pb)))))


(define (saturated-oil-fvf gas-sg
                           stock-tank-oil-sg
                           reservoir-temp
                           pressure)
  "The oil FVF for use at pressures equal to or below bubblepoint.
Arguments:
    gas-sg: (separator) gas specific gravity (dimensionless)
    stock-tank-oil-api: stock-tank oil API gravity (API)
    reservoir-temp: reservoir temperature (F)
    pressure: pressure (psia)
Return:
    saturated oil formation volume factor (RB/STB)"
  (let ((rs (solution-gor gas-sg
                          (api-from-sg stock-tank-oil-sg)
                          reservoir-temp
                          pressure)))
    (saturated-oil-fvf2 gas-sg
                        stock-tank-oil-sg
                        reservoir-temp
                        rs)))


(define (saturated-oil-fvf2 gas-sg
                            stock-tank-oil-sg
                            reservoir-temp
                            solution-gor)
  "The oil FVF for use at pressures equal to or below bubblepoint
if accurate solution GOR is known.
Argments:
    gas-sg: (separator) gas specific gravity (dimensionless)
    stock-tank-oil-api: stock-tank oil API gravity (API)
    reservoir-temp: reservoir temperature (F)
    solution-gor: solution GOR (scf/STB)
Return:
    saturated oil formation volume factor (RB/STB)

It is less accurate to use correlations for solution-gor, which is
implemented in saturated-oil-fvf."
  (test-bound-throw reservoir-temp 0 325)
  (+ 0.9759 (* 0.00012 (expt (+ (* solution-gor
                                   (expt (/ gas-sg stock-tank-oil-sg)
                                         0.5))
                                (* 1.25 reservoir-temp))
                             1.2))))

