(ns stem.core
  (:use [clojure.contrib pprint] [clojure set])
  (:require [clojure.contrib.str-utils2 :as s]
            [clojure.contrib.seq-utils :as s-utils]
            [stem.newick :as newick]
            [stem.gene-tree :as g-tree]
            [stem.util :as util])
  (:import [java.io BufferedReader FileReader]
           [org.yaml.snakeyaml Yaml])
  (:gen-class))

(defn parse-yaml-config-file
  "Reads in the yaml config file and returns a map with the following keys:
  :files = a map - keys are filenames, values are the constant
  :species = a map - keys are specie names, values are comma separated lineage names
  :properties = a map with user configurable properties"
  [file-name]
  (let [yaml (Yaml.)
        yaml-map (.load yaml (slurp file-name))]
    (zipmap (map keyword (.keySet yaml-map))
            (map #(into {} %) (.values yaml-map)))))

(defn get-lineages-to-spec-map
  "From the generic property map, builds the lineages to species map"
  [map]
  (let [s-map (:species map)
        specs (keys s-map)]
    (reduce
     (fn [m [k v]]
       (let [lineages (s/split v #",")]
         (merge m (zipmap lineages (repeat k)))))
     {} s-map)))

(defn lins->indexes [lins index-map]
  (map index-map lins))

(defn fill-matrix-for-lineage
  [i js matrix c-time]
  (doseq [j js]
    (if (> i j) ;; matrix is symmetric - only fill half
      (util/aset! matrix j i c-time)
      (util/aset! matrix i j c-time))))

(defn rec-fill-time-matrix
  [node matrix lin-index parent-desc-set parent-c-time]
  (if (newick/is-leaf? node)
    (let [[node-map _ _] node
          n-name (node-map :name)]
      (fill-matrix-for-lineage (lin-index n-name)
                               (lins->indexes (difference parent-desc-set #{n-name}) lin-index)
                               matrix
                               parent-c-time))
    (let [[m left right] node
          c-time (m :c-time)
          desc (m :desc)]
      (rec-fill-time-matrix left matrix lin-index desc c-time)
      (rec-fill-time-matrix right matrix lin-index desc c-time))))

(defn fill-time-matrix-for-tree
  [tree matrix lin-set lin-index]
  (rec-fill-time-matrix tree matrix lin-set lin-index (atom {})))





(defn -main [& args]
  (let [prop-map (parse-yaml-config-file "settings.yaml")
        line-to-spec-map (get-lineages-to-spec-map prop-map)
        line-set (set (keys line-to-spec-map))
        spec-set (set (vals line-to-spec-map))
        gene-trees (g-tree/get-gene-trees (:files prop-map))]
    (pprint line-set)
    (pprint spec-set)))



