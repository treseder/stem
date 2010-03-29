(ns stem.core
  (:use [clojure.contrib pprint])
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


(defn -main [& args]
  (let [prop-map (parse-yaml-config-file "settings.yaml")
        lin-to-spec-map (get-lineages-to-spec-map prop-map)
        gene-trees (g-tree/get-gene-trees (prop-map :files))]
    (pprint (count gene-trees))))



