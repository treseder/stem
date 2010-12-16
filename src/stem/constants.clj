(ns stem.constants)

(def *stem-version* 2.0)

;; if this value is true, then a caught exception will exit the
;; system.  If in dev, we don't want it to exit since it destroys the
;; REPL session.  This should be set to true before generating the
;; uber jar
(def *in-production* true)

(def *internal-node-name* "internal")
(def *root-name* "stem-root")

(def *num-saved-trees-default* 10)
(def *burnin-default* 100)
(def *beta-default* 0.0005)
(def *bound-total-iter-default* 200000)
(def *mle-filename-default* "mle.tre")
(def *bootstrap-filename* "boot.tre")
(def *search-filename-default* "search.tre")

