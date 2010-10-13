(defproject stem "0.1.0"
  :description "STEM is a program for inferring maximum likelihood species trees from a collection of estimated gene trees under the coalescent model."
  :url ""
  :dependencies [[org.clojure/clojure "1.2.0"]
                 [org.clojure/clojure-contrib "1.2.0"]
                 [org.clojars.kjw/snakeyaml "1.5"]
                 [org.clojars.overtone/vijual "0.2.1"]]
  :dev-dependencies [[swank-clojure "1.3.0-SNAPSHOT"]]
;  :warn-on-reflection true
  :main stem.core)
