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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;functions for building the coalescent matrix ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn get-sorted-lins-by-index [l-index]
  (sort-by #(l-index %) (keys l-index)))

(defn create-lin-index-map [l-set]
  (zipmap l-set (range 0 10000)))

(defn lins->indexes
  "Takes a seq of lineage names and turns them into the indexes of the matrix"
  [lins index-map]
  (map index-map lins))

(defn fill-matrix-for-lineage
  [i js matrix c-time]
  (doseq [j js]
    (if (> i j) ;; matrix is symmetric - only fill half
      (util/aset! matrix i j c-time)
      (util/aset! matrix j i c-time))))

(defn rec-fill-time-matrix
  [node matrix lin-index parent-desc-set parent-c-time]
  (let [[{name :name  desc-set :desc c-time :c-time} l-node r-node] node
        ;; if leaf consider self the only desc
        checked-set (if (empty? desc-set) #{name} desc-set)
        diff-set (difference parent-desc-set checked-set)
        js (lins->indexes diff-set lin-index)]
    (doseq [i checked-set] (fill-matrix-for-lineage (lin-index i) js matrix parent-c-time))
    (if-not (newick/is-leaf? node)
      (do (rec-fill-time-matrix l-node matrix lin-index desc-set c-time)
          (rec-fill-time-matrix r-node matrix lin-index desc-set c-time)))))

(defn fill-time-matrix-for-tree
  [tree matrix l-set]
  (rec-fill-time-matrix tree
                        matrix
                        (create-lin-index-map l-set)
                        ((first tree) :desc)
                        0))

(defn -main [& args]
  (let [prop-map (parse-yaml-config-file "settings.yaml")
        line-to-spec-map (get-lineages-to-spec-map prop-map)
        line-set (set (keys line-to-spec-map))
        spec-set (set (vals line-to-spec-map))
        gene-trees (g-tree/get-gene-trees (:files prop-map))]
))
