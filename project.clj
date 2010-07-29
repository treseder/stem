(defproject stem "0.1.0"
  :description "STEM is a program for inferring maximum likelihood species trees from a collection of estimated gene trees under the coalescent model."
  :url ""
;  :dependencies [[org.clojure/clojure "1.1.0"]
;                 [org.clojure/clojure-contrib "1.1.0-master-SNAPSHOT"]
  :dependencies [[org.clojure/clojure "1.2.0-beta1"]
                 [org.clojure/clojure-contrib "1.2.0-beta1"]
                 [org.clojars.kjw/snakeyaml "1.5"]]
  :dev-dependencies [[swank-clojure "1.3.0-SNAPSHOT"]
                                        ; [vijual "0.1.0-SNAPSHOT"]
                     ]
  :main stem.core)
