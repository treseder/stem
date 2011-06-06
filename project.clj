(defproject stem "0.1.0"
  :description "STEM is a program for inferring maximum likelihood species trees from a collection of estimated gene trees under the coalescent model."
  :url ""
  :dependencies [[org.clojure/clojure "1.2.0"]
                 [org.clojure/clojure-contrib "1.2.0"]
                 [org.clojars.kjw/snakeyaml "1.5"]]
  :dev-dependencies [[swank-clojure "1.4.0-SNAPSHOT"]
                     [clojure-source "1.2.0"]]
  :jvm-opts ["-agentlib:jdwp=transport=dt_socket,server=y,suspend=n"]
;  :warn-on-reflection true
  :main stem.core)
