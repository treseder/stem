(ns stem.constants)

(def *stem-version* 2.0)

;; if this value is true, then a caught exception will exit the
;; system.  If in dev, we don't want it to exit since it destroys the
;; REPL session.  This should be set to true before generating the
;; uber jar
(def *in-production* false)

(def *internal-node-name* "internal")
(def *root-name* "stem-root")

