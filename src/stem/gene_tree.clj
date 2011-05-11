(ns stem.gene-tree
  (:require [stem.util :as util]
            [stem.newick :as newick]
            [clojure.string :as str]))

(defrecord GeneTree [vec-tree rate])

(defn create-gene-tree
  "The struct that holds the tree structure of the gene-tree
  and its rate (double)"
  [vec-tree rate]
  (GeneTree. vec-tree rate))

(defn parse-rate-and-tree [s]
  "s could be of the form
  [rate](newick-tree)
       or
  (newick-tree)"
  (re-find #"(\[(\d*\.\d*)\])?(.*)" s))


(defn create-gene-tree-from-parts
  [rate-str newick-str theta]
  (let [rate (util/to-double rate-str)]
    (create-gene-tree
     (newick/newick->tree newick-str rate theta)
     rate)))

(defn parse-gene-tree-str
  ([s theta & [r]]
     ;; if the user has input the rate as zero in the settings file,
     ;; then the rates are defined in the gene-tree file
     (let [[_ _ rate n-str] (parse-rate-and-tree s)
           true-rate (if (nil? r) rate (if (zero? r) rate r))]
       (create-gene-tree-from-parts true-rate n-str theta))))

(defn parse-gene-tree-file
  "Returns a sequence of gene-tree structs found in file-name"
  [file-name theta & [rate]]
  (let [file-str (util/remove-whitespace (util/read-file file-name))
        newick-strs (str/split file-str #";")]
    (map #(parse-gene-tree-str % theta rate) newick-strs)))

(defn get-gene-trees [file-map, theta]
  "file-map is a map with keys of file names, and values of rates.
  The user has the option of omitting this portion of the settings file
  and putting all gene trees in a file called genetrees.tre.  This
  function checks for both"
  (if file-map
    (reduce (fn [s [file-name rate]] (concat s (parse-gene-tree-file file-name theta rate))) [] file-map)
    ;assumes all genetrees are in one file called genetrees.tre where
    ;each newick tree is preceded by the rate
    (parse-gene-tree-file "genetrees.tre" theta)))
