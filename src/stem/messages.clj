(ns stem.messages
  (:require [clojure.contrib.str-utils2 :as s]
            [clojure.contrib.seq-utils :as s-utils]
            [stem.util :as util]))

(def e-strs {:yaml "An error occurred parsing the settings file."})

(defmacro println-if [test form]
  `(if ~test
    (println ~form)
    (flush)))

(defn yaml-message [prop-map]
  (println "Successfully parsed the settings file")
  (flush))

(defn theta-message [theta]
  (println "Using theta =" theta) (println) (flush))

(defn lin-set-message [s]
  (println "There are" (count s) "lineages:")
  (println (apply str (interpose ", " s)))
  (println) (flush))

(defn spec-set-message [s]
  (println "There are" (count s) "species:")
  (println (apply str (interpose ", " s)))
  (println) (flush))

(defn spec-newick-message [s num-trees]
  (println "The Maximum Likelihood tree (written to file 'mle.tre') is:")
  (println)
  (println s)
  (println-if (> num-trees 1) (str "There were also " num-trees " trees that have the same likelihood estimate as the tree above.\nThese trees are also output to 'mle.tre'."))
  (flush))

(defn mle-message [mle]
  (println (str "The likelihood estimate is: " mle))
  (flush))
