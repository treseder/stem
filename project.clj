(defproject stem "0.1.0"
  :description "STEM is a program for inferring maximum likelihood species trees from a collection of estimated gene trees under the coalescent model."
  :url ""
  :dependencies [[org.clojure/clojure "1.1.0"]
                 [org.clojure/clojure-contrib "1.1.0-master-SNAPSHOT"]
                 [org.clojars.kjw/snakeyaml "1.5"]]
  :dev-dependencies [[leiningen/lein-swank "1.0.0-SNAPSHOT"]
                     [vijual "0.1.0-SNAPSHOT"]]
  :main stem.core)
