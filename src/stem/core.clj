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
  :files
  :species
  :properties

  Each value of the map is another map of key value pairs"
  [file-name]
  (let [yaml (Yaml.)
        yaml-map (.load yaml (slurp file-name))]
    (zipmap (map keyword (.keySet yaml-map))
            (map #(into {} %) (.values yaml-map)))))

(defn -main [& args]
  (let [prop-map (parse-yaml-config-file "settings.yaml")
        gene-trees (g-tree/parse-gene-tree-file "genetrees.tre")]
    (pprint (count gene-trees))))
