;;; Commentary:

;;; These correlations are based on the paper "Reservoir-Fluid
;;; Property Correlations - State of the Art" by W.D. McCain Jr. (SPE
;;; 18571)

;;; Code:


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


(define (bubble-point-pressure solution-gor
                               separator-gas-sg
                               stock-tank-oil-api
                               reservoir-temp)
  (if ))


