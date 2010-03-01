(ns stem.core
  (:use [clojure.contrib duck-streams pprint])
  (:require [clojure.contrib.str-utils2 :as s]
            [clojure.contrib.seq-utils :as s-utils]
            [stem.newick :as newick]
            [stem.gene-tree :as g-tree]
            [stem.util :as util])
  (:import [java.io BufferedReader FileReader]
           [org.yaml.snakeyaml Yaml])
  (:gen-class))

(defn map-from-line [line]
  (if-not (empty? line)
    (let [[k v] (.split line ":")]
      {(keyword (.trim k)) (.trim v)})))

(defn lines-into-map [seq]
  (apply merge (map map-from-line seq)))


(defn parse-prop-file
  "Parses a property file with the following format:
  prop1:value1
  prop2:value2
  ...

  Returns a map of the properties, with keys as the property names."
  [f-name]
  (with-open [rdr (reader (util/get-file f-name))]
    (lines-into-map (line-seq rdr))))


(defn parse-yaml-config-file
  "Reads in the yaml config file and returns a map with the following keys:
  :files
  :species
  :properties

  Each value of the map is another map of key value pairs"
  []
  (let [yaml (Yaml.)
        yaml-map (.load yaml (slurp "test.yaml"))]
    (zipmap (map keyword (.keySet yaml-map))
            (map #(into {} %) (.values yaml-map)))))

(defn -main [& args]
  (let [prop-map (parse-prop-file "test.prop")
        gene-trees (g-tree/parse-gene-tree-file "genetrees.tre")]
    (pprint prop-map)
    (pprint (count gene-trees))))


