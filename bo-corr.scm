;;; Commentary:

;;; These correlations are based on the paper "Reservoir-Fluid
;;; Property Correlations - State of the Art" by W.D. McCain Jr. (SPE
;;; 18571)

;;; Code:

(define-module (dumb-pe bo-corr)
  #:use-module ((dumb-pe field-units))
  #:export (stock-tank-gor
            bubble-point-pressure
            solution-gor
            saturated-oil-fvf))


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
  (if (> separator-temp 140)
      (throw 'value-out-of-range "separator-temp<140" separator-temp)
      (expt 10 (+ 0.3818
                  (* -5.506 (log10 stock-tank-oil-sg))
                  (* 2.902 (log10 separator-gas-sg))
                  (* 1.327 (log10 separator-pressure))
                  (* -0.7355 (log10 separator-temp))))))


(define (bubble-point-pressure initial-solution-gor
                               separator-gas-sg
                               stock-tank-oil-api
                               reservoir-temp)
  "Estimation of reservoir oil bubble point pressure to an accuracy of 15%.
Arguments:
    initial-solution-gor: initial solution gas oil ratio (scf/STB)
    separator-gas-sg: separator gas specific gravity (dimensionless)
    stock-tank-oil-api: stock-tank oil API gravity (API)
    reservoir-temp: reservoir temperature (F)
Return:
    bubble point pressure (psia)

The initial solution gas oil ratio should include stock-tank vent gas.
"
  (if (> reservoir-temp 350)
      (throw 'value-out-of-range "reservoir-temp<350" reservoir-temp)
      (* 18.2 (- (* (expt (/ initial-solution-gor separator-gas-sg) 0.83)
                    (expt 10 (- (* 0.00091 reservoir-temp)
                                (* 0.0125 stock-tank-oil-api)))) 1.4))))


(define (solution-gor separator-gas-sg
                      stock-tank-oil-api
                      reservoir-temp
                      pressure)
  "Estimation of solution gas oil ratio for pressures below the bubble point pressure.
Arguments:
    separator-gas-sg: separator gas specific gravity (dimensionless)
    stock-tank-oil-api: stock-tank oil API gravity (API)
    reservoir-temp: reservoir temperature (R)
    pressure: pressure (psia)
Return:
    solution gas oil ratio (scf/STB)"
  (* separator-gas-sg (expt (/ (+ (/ pressure 18.2) 1.4)
                               (expt 10 (- (* 0.00091 reservoir-temp)
                                           (* 0.0125 stock-tank-oil-api))))
                            (/ 1 0.83))))


(define (saturated-oil-fvf separator-gas-sg
                           stock-tank-oil-sg
                           reservoir-temp
                           pressure)
  "The oil FVF for use at pressures equal to or below bubblepoint.
Arguments:
    separator-gas-sg: separator gas specific gravity (dimensionless)
    stock-tank-oil-api: stock-tank oil API gravity (API)
    reservoir-temp: reservoir temperature (R)
    pressure: pressure (psia)
Return:
    saturated oil formation volume factor (RB/STB)"
  (let ((rs (solution-gor separator-gas-sg
                          (api-from-sg stock-tank-oil-sg)
                          reservoir-temp
                          pressure)))
    (+ 0.9759 (* 0.00012 (expt (+ (* rs (expt (/ separator-gas-sg stock-tank-oil-sg)
                                              0.5))
                                  (* 1.25 reservoir-temp))
                               1.2))))) 
